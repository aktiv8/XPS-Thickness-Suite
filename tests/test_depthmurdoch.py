"""
Murdoch and McCulloch regularised depth-profile tests (depthmurdoch.py). Synthetic data only.

    python -m unittest discover -s tests -v
"""
import math
import os
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

import depthmurdoch as DM  # noqa: E402
import depthprofile as DP  # noqa: E402

ANG = np.array([0.0, 30.0, 45.0, 55.0, 62.0, 68.0])
LAM = np.array([2.0, 2.4, 2.8])
NAMES = ["C", "Ox", "Sub"]
DZ, DEPTH = 0.1, 6.0


def _truth(thick=(0.6, 1.5)):
    z = DM.depths(int(round(DEPTH / DZ)), DZ)
    return z, DM.layer_guess([0, 1, 2], list(thick), 3, z, floor=1e-6)


def _data(thick=(0.6, 1.5)):
    z, c = _truth(thick)
    return DM.apparent_from_profile(c, DM.attenuation_kernel(z, LAM, ANG))


class MurdochTests(unittest.TestCase):
    def test_ealformula_matches_the_spreadsheet_expression(self):
        for e in (400.0, 1000.0, 1200.0):
            self.assertAlmostEqual(DM.eal_murdoch_nm(e), 143.0 / e ** 2 + 0.054 * math.sqrt(e))
        with self.assertRaises(ValueError):
            DM.eal_murdoch_nm(0)

    def test_apparent_concentrations_close_to_one(self):
        x = DM.apparent_at_pct([[1, 2], [3, 6]], [1.0, 3.0])
        self.assertTrue(np.allclose(x.sum(axis=0), 1.0))

    def test_gradient_matches_finite_differences(self):
        z, c0 = _truth()
        K = DM.attenuation_kernel(z, LAM, ANG)
        prob = DM._Problem(K, _data((0.7, 1.2)), 1.0, 1e-4, DZ)
        rng = np.random.default_rng(3)
        v = rng.normal(size=len(z) * 3)
        f0, g = prob.value_and_grad(v)
        for i in rng.integers(0, len(v), 6):
            e = np.zeros_like(v)
            e[i] = 1e-6
            num = (prob.value_and_grad(v + e)[0] - prob.value_and_grad(v - e)[0]) / 2e-6
            self.assertAlmostEqual(g[i], num, delta=1e-5 * max(1.0, abs(num)))

    def test_profile_is_closed_and_chi2_vanishes_at_the_truth(self):
        z, c = _truth()
        x = _data()
        f = DM.fit(NAMES, ANG, x, LAM, DZ, DEPTH, 1e-8, 1.0, c)
        self.assertTrue(np.allclose(f.c_pct.sum(axis=1), 100.0))
        self.assertLess(f.chi2, 1e-6)
        self.assertLess(f.peak_position_nm(0), 0.6)
        self.assertGreater(f.peak_position_nm(2), 2.1)

    def test_recovers_the_layer_structure_from_a_perturbed_guess(self):
        x = _data()
        z = DM.depths(int(round(DEPTH / DZ)), DZ)
        c0 = DM.layer_guess([0, 1, 2], [0.4, 2.0], 3, z)
        f = DM.fit(NAMES, ANG, x, LAM, DZ, DEPTH, 1e-6, 1.0, c0)
        self.assertLess(f.chi2, 2e-3)
        self.assertGreater(f.c_pct[0, 0], 60)                         # carbon on top
        self.assertGreater(f.c_pct[-1, 2], 90)                        # substrate at depth
        self.assertTrue(0.3 < f.peak_position_nm(1) < 2.2)

    def test_large_alpha_flattens_and_small_alpha_fits_better(self):
        x = _data()
        z = DM.depths(int(round(DEPTH / DZ)), DZ)
        c0 = DM.layer_guess([0, 1, 2], [0.4, 2.0], 3, z)
        smooth = DM.fit(NAMES, ANG, x, LAM, DZ, DEPTH, 1e2, 1.0, c0)
        sharp = DM.fit(NAMES, ANG, x, LAM, DZ, DEPTH, 1e-8, 1.0, c0)
        self.assertLess(smooth.r_term, sharp.r_term)
        self.assertGreater(smooth.chi2, sharp.chi2)

    def test_gradient_regulariser_keeps_the_substrate(self):
        # a profile that is the substrate everywhere is not penalised, unlike a squared-concentration penalty
        z = DM.depths(60, DZ)
        K = DM.attenuation_kernel(z, LAM, ANG)
        c = np.zeros((60, 3))
        c[:, 2] = 1.0
        prob = DM._Problem(K, DM.apparent_from_profile(np.maximum(c, 1e-6), K), 1.0, 1.0, DZ)
        v = np.log(np.maximum(c, 1e-6)).ravel()
        self.assertLess(prob.parts(v)[7], 1e-6)

    def test_s_curve_and_suggested_alpha(self):
        x = _data()
        z = DM.depths(int(round(DEPTH / DZ)), DZ)
        c0 = DM.layer_guess([0, 1, 2], [0.5, 1.8], 3, z)
        curve = DM.s_curve(NAMES, ANG, x, LAM, [-2, -4, -6, -8], DZ, DEPTH, 1.0, c0)
        self.assertEqual([c[0] for c in curve], [-8, -6, -4, -2])
        self.assertLessEqual(curve[0][1], curve[-1][1] + 1e-9)       # chi2 rises with alpha
        self.assertIn(DM.suggest_alpha(curve), (-8.0, -6.0, -4.0, -2.0))
        self.assertEqual(DM.suggest_alpha([(-2, 1.0, 0), (-4, 0.1, 0), (-6, 0.09, 0)]), -4.0)

    def test_ensemble_reports_mean_and_spread(self):
        x = _data()
        z, mean, sd, fits = DM.ensemble(NAMES, ANG, x, LAM, [0, 1, 2], [0.6, 1.5], 0.5, 4, 1, DZ, DEPTH, 1e-6)
        self.assertEqual(len(fits), 4)
        self.assertEqual(mean.shape, (len(z), 3))
        self.assertTrue(np.all(sd >= 0))
        self.assertTrue(np.allclose(mean.sum(axis=1), 100.0))

    def test_csv_reader(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "x.csv")
            open(p, "w").write("angle,C 1s,O 1s\n0,20,80\n30,30,70\n60,50,50\n")
            names, ang, x = DM.load_at_pct_csv(p)
            self.assertEqual(names, ["C 1s", "O 1s"])
            self.assertTrue(np.allclose(ang, [0, 30, 60]))
            self.assertTrue(np.allclose(x[0], [0.2, 0.3, 0.5]))
            open(p, "w").write("angle,C\n0,1\n")
            with self.assertRaises(ValueError):
                DM.load_at_pct_csv(p)

    def test_inputs_from_species(self):
        ang = np.array([0.0, 30.0, 60.0, 75.0])
        a = DP.Species("A", "r", ang, np.array([100.0, 90, 70, 60]), sigma=np.full(4, 1.0), lam_nm=2.0, rsf=2.0)
        b = DP.Species("B", "r", ang, np.array([50.0, 60, 80, 90]), sigma=np.full(4, 1.0), lam_nm=2.0, rsf=1.0)
        names, angles, x, sig = DM.inputs_from_species([a, b], 70.0)
        self.assertEqual(list(angles), [0.0, 30.0, 60.0])
        self.assertTrue(np.allclose(x.sum(axis=0), 1.0))
        self.assertAlmostEqual(x[0, 0], 0.5)                          # 100/2 against 50/1
        self.assertEqual(sig.shape, x.shape)
        with self.assertRaises(ValueError):
            DM.inputs_from_species([a], 70.0)


if __name__ == "__main__":
    unittest.main()
