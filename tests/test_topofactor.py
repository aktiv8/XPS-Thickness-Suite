"""
Tests for topofactor.py (Shard, Wang, Spencer, Surf. Interface Anal. 41, 541 (2009)). Pure maths, no display.

    python -m unittest discover -s tests -v
"""
import math
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

import thickogram as T  # noqa: E402
import topofactor as TF  # noqa: E402


def _eq6_by_hand(c, t, b):
    c0, c1, c2, c3, c4, k = c
    return c0 + (c1 + c2 * math.log(b)) * t ** k + c3 * math.exp(-c4 * b)


class PaperConstants(unittest.TestCase):
    def test_table_1_is_verbatim(self):
        self.assertEqual(TF.TABLE_1["cylinder"], (0.46, 0.25, 0.09, 0.16, 1.1, 0.22))
        self.assertEqual(TF.TABLE_1["sphere"], (0.40, 0.20, 0.09, 0.14, 1.5, 0.33))
        self.assertEqual(TF.TABLE_1["sphere_shadowed"], (0.46, 0.19, 0.08, 0.13, 1.4, 0.30))

    def test_thin_limits_and_constants(self):
        self.assertAlmostEqual(TF.SHAPES["cylinder"].thin_limit, 2 / math.pi)       # 0.64
        self.assertEqual(TF.SHAPES["sphere"].thin_limit, 0.50)
        self.assertEqual(TF.SHAPES["sphere_shadowed"].thin_limit, 0.57)
        self.assertEqual(TF.SHAPES["cylinder"].constant, 0.79)
        self.assertEqual(TF.SHAPES["sphere"].constant, 0.67)
        self.assertEqual(TF.SHAPES["sphere_shadowed"].constant, 0.71)
        self.assertEqual(TF.SHAPES["flat"].constant, 1.0)

    def test_ranges(self):
        self.assertEqual((TF.T_THIN, TF.T_MAX), (0.1, 3.0))
        self.assertEqual(TF.B_USUAL, (0.7, 1.4))


class Equation6(unittest.TestCase):
    def test_matches_a_hand_calculation(self):
        for key in ("cylinder", "sphere", "sphere_shadowed"):
            for t, b in ((0.2, 1.0), (1.5, 1.0), (2.9, 0.8), (1.0, 1.4)):
                self.assertAlmostEqual(TF.eq6(key, t, b), _eq6_by_hand(TF.TABLE_1[key], t, b), places=12)

    def test_values_at_t_1p5_b_1_agree_with_the_papers_quoted_averages(self):
        """The paper quotes 0.79 (cylinder), 0.71 (shadowed sphere), 0.67 (sphere) at t = 1.5, B = 1."""
        cyl, sph, shad = (TF.eq6(k, 1.5, 1.0) for k in ("cylinder", "sphere", "sphere_shadowed"))
        self.assertAlmostEqual(cyl, 0.7866, delta=5e-4)
        self.assertAlmostEqual(sph, 0.6598, delta=5e-4)
        self.assertAlmostEqual(shad, 0.7067, delta=5e-4)
        self.assertAlmostEqual(cyl, 0.79, delta=0.79 * 0.01)
        self.assertAlmostEqual(shad, 0.71, delta=0.71 * 0.01)
        self.assertAlmostEqual(sph, 0.67, delta=0.67 * 0.02)      # the paper's own t-vs-T inconsistency, see module doc

    def test_increases_with_thickness_and_stays_below_one(self):
        for key in ("cylinder", "sphere", "sphere_shadowed"):
            for b in (0.7, 1.0, 1.4):
                vals = [TF.eq6(key, t, b) for t in (0.1, 0.3, 0.6, 1.0, 2.0, 3.0)]
                self.assertEqual(vals, sorted(vals), f"{key} B={b}")
                self.assertTrue(all(0.3 < v < 1.0 for v in vals), f"{key} B={b}: {vals}")

    def test_eq6_approaches_the_thin_limits_near_t_0p1(self):
        for key in ("cylinder", "sphere", "sphere_shadowed"):
            self.assertAlmostEqual(TF.eq6(key, 0.1, 1.0), TF.SHAPES[key].thin_limit, delta=0.06, msg=key)


class Topofactor(unittest.TestCase):
    def test_flat_is_exactly_one_everywhere(self):
        for t in (0.0, 0.05, 1.0, 5.0):
            r = TF.topofactor("flat", t, 1.0, "eq6", theta_normal_deg=40.0)
            self.assertEqual((r.tf, r.regime, r.warnings), (1.0, "planar", []))

    def test_regimes_at_the_boundaries(self):
        self.assertEqual(TF.topofactor("cylinder", 0.0999, 1.0).regime, "thin limit")
        self.assertEqual(TF.topofactor("cylinder", 0.1, 1.0).regime, "Eq. 6")
        self.assertEqual(TF.topofactor("cylinder", 3.0, 1.0).regime, "Eq. 6")
        self.assertEqual(TF.topofactor("cylinder", 3.0001, 1.0).regime, "Eq. 6 (held at t = 3)")
        self.assertAlmostEqual(TF.topofactor("cylinder", 0.05, 1.0).tf, 2 / math.pi)
        self.assertAlmostEqual(TF.topofactor("sphere", 0.05, 1.0).tf, 0.50)
        self.assertAlmostEqual(TF.topofactor("sphere_shadowed", 0.05, 1.0).tf, 0.57)

    def test_above_t_3_the_factor_is_held_and_warned(self):
        r = TF.topofactor("sphere", 8.0, 1.0)
        self.assertAlmostEqual(r.tf, TF.eq6("sphere", 3.0, 1.0), places=12)
        self.assertTrue(any("above 3" in w for w in r.warnings))

    def test_constant_mode(self):
        for key in ("cylinder", "sphere", "sphere_shadowed"):
            for t in (0.01, 1.0, 9.0):
                r = TF.topofactor(key, t, 1.0, "constant")
                self.assertEqual((r.tf, r.regime), (TF.SHAPES[key].constant, "constant"))

    def test_warnings(self):
        self.assertEqual(TF.topofactor("cylinder", 1.0, 1.0, "eq6", 0.0).warnings, [])
        self.assertTrue(any("theta = 0" in w for w in TF.topofactor("cylinder", 1.0, 1.0, "eq6", 30.0).warnings))
        self.assertTrue(any("outside the usual range" in w for w in TF.topofactor("cylinder", 1.0, 2.0).warnings))
        # B is irrelevant to the constant values, so no B warning there
        self.assertFalse(any("usual range" in w for w in TF.topofactor("cylinder", 1.0, 2.0, "constant").warnings))

    def test_the_nanoscale_note_is_always_present_for_curved_samples(self):
        for key in ("cylinder", "sphere", "sphere_shadowed"):
            for mode in TF.MODES:
                self.assertIn(TF.NANOSCALE_NOTE, TF.topofactor(key, 1.0, 1.0, mode).notes)

    def test_invalid_input(self):
        for args in (("torus", 1.0, 1.0), ("cylinder", -0.1, 1.0), ("cylinder", float("nan"), 1.0),
                     ("cylinder", 1.0, 0.0), ("cylinder", 1.0, -1.0)):
            with self.assertRaises(ValueError, msg=str(args)):
                TF.topofactor(*args)
        with self.assertRaises(ValueError):
            TF.topofactor("cylinder", 1.0, 1.0, mode="nope")


class BConversion(unittest.TestCase):
    def test_b_paper_is_the_energy_ratio_to_the_power_0p75(self):
        self.assertAlmostEqual(TF.b_paper(1.0), 1.0)
        self.assertAlmostEqual(TF.b_paper(2.0), 2 ** 0.75)
        for bad in (0.0, -1.0, float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                TF.b_paper(bad)

    def test_reproduces_the_papers_experimental_b_values(self):
        """Paper (Al Ka): B(C 1s : O 1s) = 1.19 and B(C 1s : Si 2p) = 0.90."""
        b_co = TF.b_paper(T.energy_ratio_b(285.0, 532.0, 1486.6))
        b_csi = TF.b_paper(T.energy_ratio_b(285.0, 99.5, 1486.6))
        self.assertAlmostEqual(b_co, 1.19, delta=0.01)
        self.assertAlmostEqual(b_csi, 0.90, delta=0.01)


class Apply(unittest.TestCase):
    @staticmethod
    def _thickogram(theta=0.0, io=20000.0, lam=1.6):
        # Fe 2p overlayer on Si 2p demo numbers (Al Ka); B = (KEo/KEs) = 0.562 -> paper's B = 0.65
        return T.thickness(io=io, is_=6000.0, so=10.8, ss=0.82, be_o=707.0, be_s=100.0, source_ev=1486.6,
                           lam_nm=lam, theta_normal_deg=theta)

    def test_corrected_thickness_is_tf_times_the_planar_one(self):
        res = self._thickogram()
        out = TF.apply_topofactor(res, "cylinder", "eq6")
        t_over_l = res.t_nm / 1.6
        want = TF.topofactor("cylinder", t_over_l, res.b ** 0.75).tf
        self.assertAlmostEqual(out.tf, want, places=12)
        self.assertAlmostEqual(out.t_nm, want * res.t_nm, places=12)
        self.assertAlmostEqual(out.planar_t_nm, res.t_nm)
        self.assertAlmostEqual(out.t_over_l, t_over_l)
        self.assertAlmostEqual(out.b_paper, res.b ** 0.75, places=12)

    def test_flat_leaves_the_result_unchanged(self):
        res = self._thickogram()
        out = TF.apply_topofactor(res, "flat")
        self.assertEqual((out.tf, out.t_nm, out.regime), (1.0, res.t_nm, "planar"))

    def test_constant_mode_applies_the_papers_average(self):
        res = self._thickogram()
        self.assertAlmostEqual(TF.apply_topofactor(res, "sphere", "constant").t_nm, 0.67 * res.t_nm, places=12)

    def test_corrected_thickness_is_always_smaller_than_the_planar_one(self):
        res = self._thickogram()
        for key in ("cylinder", "sphere", "sphere_shadowed"):
            for mode in TF.MODES:
                self.assertLess(TF.apply_topofactor(res, key, mode).t_nm, res.t_nm)

    def test_off_normal_angle_warns(self):
        out = TF.apply_topofactor(self._thickogram(theta=30.0), "cylinder")
        self.assertTrue(any("theta = 0" in w for w in out.warnings))

    def test_thin_and_thick_regimes_through_the_whole_pipeline(self):
        thin = TF.apply_topofactor(self._thickogram(io=300.0), "cylinder")
        self.assertEqual(thin.regime, "thin limit")
        thick = TF.apply_topofactor(self._thickogram(lam=100.0), "cylinder")        # t/L large -> still defined
        self.assertIn(thick.regime, ("Eq. 6", "Eq. 6 (held at t = 3)", "thin limit"))


if __name__ == "__main__":
    unittest.main()
