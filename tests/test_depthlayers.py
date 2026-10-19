"""
Layer-stack (Gries signal-ratio) tests: the forward model, the least-squares fit, the thickness-range scan and the first
guess of which layer holds which species. Synthetic data only.

    python -m unittest discover -s tests -v
"""
import math
import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

import depthlayers as DL  # noqa: E402
import depthprofile as DP  # noqa: E402
import synthetic as S  # noqa: E402

ANGLES = np.array(S.DEFAULT_TILTS)
SUB = DL.SUBSTRATE_LAYER


def _species(name, layer, d, lam=3.0, noise=None, seed=0):
    intensity = 1000.0 * DL.layer_intensity(ANGLES, d, layer, [lam] * len(d))
    sigma = None
    if noise:
        rng = np.random.default_rng(seed)
        counts = rng.poisson(intensity * noise) / noise
        sigma, intensity = np.sqrt(np.maximum(counts, 1.0) / noise), counts
    return DP.Species(name, "r", ANGLES.copy(), intensity, sigma=sigma, lam_nm=lam, ke_ev=1000.0)


class ForwardModel(unittest.TestCase):
    def test_the_substrate_is_attenuated_by_every_layer(self):
        th = np.array([0.0, 60.0])
        got = DL.layer_intensity(th, [1.0, 2.0], SUB, [3.0, 4.0])
        self.assertTrue(np.allclose(got, np.exp(-(1.0 / 3.0 + 2.0 / 4.0) / np.cos(np.radians(th)))))

    def test_a_layer_species_is_attenuated_only_by_the_layers_above_and_grows_with_its_own_thickness(self):
        th = np.array([0.0, 45.0])
        c = np.cos(np.radians(th))
        got = DL.layer_intensity(th, [1.0, 2.0], 1, [3.0, 4.0])
        self.assertTrue(np.allclose(got, np.exp(-1.0 / 3.0 / c) * (1 - np.exp(-2.0 / 4.0 / c))))
        top = DL.layer_intensity(th, [1.0, 2.0], 0, [3.0, 4.0])
        self.assertTrue(np.allclose(top, 1 - np.exp(-1.0 / 3.0 / c)))

    def test_an_infinitely_thick_top_layer_gives_unit_intensity_and_no_substrate_signal(self):
        th = np.array([0.0, 60.0])
        self.assertTrue(np.allclose(DL.layer_intensity(th, [200.0], 0, [3.0]), 1.0))
        self.assertTrue(np.allclose(DL.layer_intensity(th, [200.0], SUB, [3.0]), 0.0, atol=1e-12))

    def test_model_curves_use_each_species_own_lambda_unless_a_matrix_is_given(self):
        a = _species("A", 0, [1.0], lam=2.0)
        b = _species("B", SUB, [1.0], lam=5.0)
        curves = DL.model_curves([a, b], [0, SUB], [1.0])
        self.assertTrue(np.allclose(curves[0], DL.layer_intensity(ANGLES, [1.0], 0, [2.0])))
        self.assertTrue(np.allclose(curves[1], DL.layer_intensity(ANGLES, [1.0], SUB, [5.0])))
        alt = DL.model_curves([a, b], [0, SUB], [1.0], lam_layers=[[9.0], [9.0]])
        self.assertTrue(np.allclose(alt[1], DL.layer_intensity(ANGLES, [1.0], SUB, [9.0])))


class Fit(unittest.TestCase):
    D = [0.8, 1.5]

    def _three(self, **kw):
        return [_species("C", 0, self.D, **kw), _species("O", 1, self.D, **kw), _species("Al", SUB, self.D, **kw)]

    def test_noise_free_data_gives_back_the_thicknesses_exactly(self):
        fit = DL.fit_layers(self._three(), [0, 1, SUB], 2, max_angle=80)
        self.assertTrue(np.allclose(fit.d_nm, self.D, atol=1e-4))
        self.assertAlmostEqual(fit.chi2, 0.0, places=6)
        self.assertEqual(fit.warnings, [])

    def test_chi_squared_is_zero_at_the_truth_and_grows_away_from_it(self):
        sps = self._three()
        at = DL.gries_chi2(self.D, sps, [0, 1, SUB], 80)[0]
        off = DL.gries_chi2([0.8, 3.0], sps, [0, 1, SUB], 80)[0]
        self.assertAlmostEqual(at, 0.0, places=8)
        self.assertGreater(off, 100.0)

    def test_the_fit_with_noise_stays_within_its_own_range_and_the_range_brackets_the_truth(self):
        sps = self._three(noise=0.05, seed=3)
        fit = DL.fit_layers(sps, [0, 1, SUB], 2, max_angle=80)
        DL.scan_ranges(fit, sps, 2, max_angle=80, delta=9.0)           # about 3 sigma, so the truth is nearly always inside
        for k, truth in enumerate(self.D):
            lo, hi = fit.ranges[k]
            self.assertLessEqual(lo, fit.d_nm[k] + 1e-9)
            self.assertGreaterEqual(hi, fit.d_nm[k] - 1e-9)
            self.assertTrue(lo - 0.05 <= truth <= hi + 0.05, f"layer {k + 1}: {truth} outside {lo}-{hi}")

    def test_information_fades_with_depth_so_the_deeper_layer_has_the_wider_range(self):
        sps = self._three(noise=0.05, seed=3)
        fit = DL.fit_layers(sps, [0, 1, SUB], 2, max_angle=80)
        DL.scan_ranges(fit, sps, 2, max_angle=80)
        width = [r[1] - r[0] for r in fit.ranges]
        self.assertGreater(width[1], width[0])

    def test_an_unassigned_layer_and_a_run_to_the_limit_are_reported(self):
        sps = [_species("C", 0, [0.8]), _species("Al", SUB, [0.8])]
        fit = DL.fit_layers(sps, [0, SUB], 2, max_angle=80)
        self.assertTrue(any("Layer(s) 2 hold no species" in w for w in fit.warnings))
        thick = [_species("C", 0, [30.0]), _species("Al", SUB, [30.0])]
        hit = DL.fit_layers(thick, [0, SUB], 1, d_max=5.0, max_angle=80)
        self.assertTrue(any("upper limit" in w for w in hit.warnings))

    def test_nothing_to_fit_is_an_error(self):
        sps = [_species("Al", SUB, [0.8])]
        with self.assertRaisesRegex(ValueError, "at least one species to a layer"):
            DL.fit_layers(sps, [SUB], 1)
        short = _species("C", 0, [0.8])
        short.intensity[1:] = -1.0
        with self.assertRaisesRegex(ValueError, "two usable angles"):
            DL.fit_layers([short], [0], 1)

    def test_only_the_shape_matters_a_common_scale_changes_nothing(self):
        sps = self._three()
        for sp in sps:
            sp.intensity = sp.intensity * 7.3
        fit = DL.fit_layers(sps, [0, 1, SUB], 2, max_angle=80)
        self.assertTrue(np.allclose(fit.d_nm, self.D, atol=1e-4))


class Helpers(unittest.TestCase):
    def test_the_simplex_finds_the_minimum_of_a_bounded_quadratic(self):
        x, v = DL.nelder_mead(lambda p: (p[0] - 1.2) ** 2 + 3 * (p[1] - 0.4) ** 2, [4.0, 4.0], [0.5, 0.5], [0, 0], [10, 10])
        self.assertTrue(np.allclose(x, [1.2, 0.4], atol=1e-4))
        x, _ = DL.nelder_mead(lambda p: (p[0] + 5.0) ** 2, [1.0], [0.5], [0.0], [10.0])      # the minimum is outside: stop at the bound
        self.assertAlmostEqual(float(x[0]), 0.0, places=6)

    def test_the_default_assignment_orders_species_by_depth_and_puts_the_deepest_in_the_substrate(self):
        d = [0.8, 1.5]
        sps = [_species("Al", SUB, d), _species("C", 0, d), _species("O", 1, d)]
        self.assertEqual(DL.default_assignment(sps, 2, 80), [SUB, 0, 1])

    def test_measured_curves_are_normalised_to_the_first_usable_angle(self):
        sp = _species("C", 0, [0.8])
        th, y, sig, ok = DL.measured_normalised(sp, 70)
        self.assertEqual(y[0], 1.0)
        self.assertEqual(list(th), [0.0, 15.0, 30.0, 45.0, 60.0])
        self.assertIsNone(DL.measured_normalised(_species("C", 0, [0.8]), -1.0))


class LayeredSynthetic(unittest.TestCase):
    def test_the_synthetic_layer_stack_round_trips_through_regions(self):
        blocks = S.layered_series([0.8, 1.5], [dict(name="C 1s", be=285.0, layer=0), dict(name="O 1s", be=532.0, layer=1),
                                               dict(name="Al 2p", be=73.0, layer=SUB)])
        by_region = {}
        for b in blocks:
            by_region.setdefault(b.block_id, []).append((b.exp_values[0], b))
        sps = [DP.species_from_region(by_region[r], transmission=False) for r in ("C 1s Scan", "O 1s Scan", "Al 2p Scan")]
        for sp in sps:
            sp.lam_nm = 3.0
        fit = DL.fit_layers(sps, [0, 1, SUB], 2, max_angle=80)
        self.assertTrue(np.allclose(fit.d_nm, [0.8, 1.5], atol=0.1))


if __name__ == "__main__":
    unittest.main()
