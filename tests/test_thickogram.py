"""
Tests for thickogram.py (Cumpson's method). Pure maths + a matplotlib smoke test; no display needed.

    python -m unittest discover -s tests -v

GOLDEN values were produced by the ORIGINAL, unmodified solver from the user's THickogram_v1_3.py (kept in
legacy/), so these tests prove the rewrite reproduces it.
"""
import math
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

import numpy as np  # noqa: E402

import app as A  # noqa: E402
import thickogram as T  # noqa: E402

try:
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.figure import Figure
    HAVE_MPL = True
except ImportError:
    HAVE_MPL = False

# (A, B, x) from the original solve_fAB_0
GOLDEN = [
    (0.05, 0.5, 0.049760045408000046),
    (0.5, 0.9, 0.41611521926899997),
    (1.0, 1.0, 0.6931471805590003),
    (2.0, 1.8, 0.8197637662969993),
    (5.4851, 1.8, 1.3000024150409997),
    (10.0, 3.0, 1.1726230294710005),
    (0.2, 0.4, 0.1995501140590001),
    (20.0, 2.2, 1.762602929446002),
    (0.8, 1.35, 0.5304752263649996),
]


class Solver(unittest.TestCase):
    def test_matches_the_original_solver(self):
        for a, b, x_orig in GOLDEN:
            self.assertAlmostEqual(T.solve_x(a, b), x_orig, delta=1e-9, msg=f"A={a} B={b}")

    def test_round_trip_over_a_grid(self):
        for b in (0.4, 0.6, 0.9, 1.0, 1.4, 2.0, 3.0):
            for x in (0.02, 0.1, 0.5, 1.3, 3.0, 5.9):
                a = T.forward_a(x, b)
                self.assertAlmostEqual(T.solve_x(a, b), x, delta=1e-9 * max(1, x), msg=f"x={x} B={b}")

    def test_extreme_intensity_ratios(self):
        for a in (1e-10, 1e-4, 1e4, 1e8):
            x = T.solve_x(a, 1.2)
            # judge the root by how well the forward model reproduces A (relative), not by the steep f near x -> 0
            self.assertTrue(math.isclose(T.forward_a(x, 1.2), a, rel_tol=1e-9), f"A={a}: x={x}")

    def test_equal_energies_reduce_to_ln_1_plus_a(self):
        for a in (0.05, 0.5, 2.0, 10.0, 300.0):
            self.assertAlmostEqual(T.solve_x(a, 1.0), math.log1p(a), places=10)

    def test_solution_increases_with_a(self):
        xs = [T.solve_x(a, 1.5) for a in (0.1, 0.3, 1, 3, 10)]
        self.assertEqual(xs, sorted(xs))

    def test_invalid_inputs_raise_instead_of_hanging(self):
        for a, b in ((0, 1), (-1, 1), (1, 0), (1, -2), (float("nan"), 1), (1, float("inf")), (float("inf"), 1)):
            with self.assertRaises(ValueError, msg=f"A={a} B={b}"):
                T.solve_x(a, b)


class CrossCheckWithTheCalculator(unittest.TestCase):
    def test_equal_energy_case_is_the_strohmeier_equation(self):
        """B = 1 (same KE): Thickogram t = x L cos(theta_normal) must equal the Calculator's d with theta from the surface."""
        lam = 2.8
        for theta_normal in (0.0, 30.0, 60.0):
            for a in (0.3, 1.0, 4.0):
                x = T.solve_x(a, 1.0)
                t_thick = x * lam * math.cos(math.radians(theta_normal))
                res = A.strohmeier_thickness(
                    T.surface_from_normal(theta_normal), 1.0, lam, 1.0, 1.0, [A.OverlayerInput("x", 1.0, lam, 1.0, a)])
                self.assertAlmostEqual(t_thick, res["d_nm"], places=9, msg=f"A={a} theta={theta_normal}")

    def test_angle_conversion(self):
        for th in (0.0, 15.5, 45.0, 89.0):
            self.assertAlmostEqual(T.surface_from_normal(T.normal_from_surface(th)), th)
        self.assertEqual(T.normal_from_surface(90.0), 0.0)


class Thickness(unittest.TestCase):
    # Fe 2p overlayer on Si 2p substrate, as in the original UI hints (Al Ka)
    KW = dict(io=1000.0, is_=3000.0, so=10.0, ss=0.9, be_o=707.0, be_s=100.0, source_ev=1486.6, lam_nm=1.2,
              theta_normal_deg=0.0)

    def test_pipeline_matches_the_original_formulae(self):
        r = T.thickness(**self.KW)
        ke_o, ke_s = 1486.6 - 707.0, 1486.6 - 100.0
        b = ke_o / ke_s
        a = 1000.0 / (3000.0 * (10.0 / 0.9))
        self.assertAlmostEqual(r.b, b, places=12)
        self.assertAlmostEqual(r.a, a, places=12)
        self.assertAlmostEqual(r.ro, 10.0 / 0.9, places=12)
        self.assertAlmostEqual(r.x, T.solve_x(a, b), places=12)
        self.assertAlmostEqual(r.t_nm, r.x * 1.2, places=12)
        self.assertAlmostEqual(r.ke_o, ke_o)
        self.assertAlmostEqual(r.ke_s, ke_s)

    def test_cos_theta_scaling(self):
        t0 = T.thickness(**self.KW).t_nm
        kw = dict(self.KW, theta_normal_deg=60.0)
        self.assertAlmostEqual(T.thickness(**kw).t_nm, t0 * 0.5, places=12)

    def test_warnings_for_out_of_range_b_and_x(self):
        r = T.thickness(**dict(self.KW, be_o=1400.0))             # KE_o = 86.6 eV, KE_s = 1386.6 eV -> B far below 0.4
        self.assertTrue(any("outside the printed nomogram scale" in w for w in r.warnings))
        r2 = T.thickness(**dict(self.KW, io=1e10))                 # huge A -> x beyond the marked 6
        self.assertTrue(any("x =" in w for w in r2.warnings))
        # in range: B = 0.86 and x ~ 1.5 -> no warnings
        ok = T.thickness(**dict(self.KW, io=1000.0, is_=300.0, so=1.0, ss=1.0, be_o=300.0, be_s=100.0))
        self.assertGreater(ok.x, 0.1)
        self.assertEqual(ok.warnings, [])
        # a very thin layer: x below the green curve's first marked point is flagged too
        self.assertTrue(any("x =" in w for w in T.thickness(**dict(self.KW, be_o=300.0, be_s=100.0)).warnings))

    def test_validation(self):
        bad = [dict(lam_nm=0.0), dict(lam_nm=-1.0), dict(theta_normal_deg=90.0), dict(theta_normal_deg=-1.0),
               dict(io=0.0), dict(is_=-5.0), dict(so=0.0), dict(ss=0.0), dict(be_o=1486.6), dict(be_s=2000.0),
               dict(source_ev=500.0)]
        for change in bad:
            with self.assertRaises(ValueError, msg=str(change)):
                T.thickness(**dict(self.KW, **change))


class Nomogram(unittest.TestCase):
    def test_constants_are_identical_to_the_original_script(self):
        np.testing.assert_allclose(T.AXP, np.arange(0.1, 6.1, 0.1))
        self.assertEqual(T.AXD, [0.1, 0.2, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        self.assertEqual(T.AYP, [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6,
                                 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30])
        self.assertEqual(T.AYD, [0.01, 0.1, 1.0, 10.0])
        np.testing.assert_allclose(T.BX, np.arange(0.4, 3.2, 0.2))
        self.assertEqual((T.B_INTERCEPT, T.B_STEP), (4.2, 0.6))

    def test_original_helper_func_equals_xb_times_b_term(self):
        for b in (0.4, 1.0, 2.2, 3.0):
            xb = T.x_of_b(b)
            self.assertAlmostEqual(T.func(xb), xb * (b ** 0.75 - 0.5), places=12)

    def test_red_line_crosses_the_green_curve_exactly_at_the_solution(self):
        """The redline fix: from (0, ln(A/2)) to X = 2*x_B, the line must meet y = ln sinh(X/4) at X = 2x."""
        for a, b, x in GOLDEN:
            (x0, y0), (x1, y1) = T.solution_line(a, b)
            self.assertEqual(x0, 0.0)
            self.assertAlmostEqual(y0, math.log(a / 2), places=12)
            self.assertAlmostEqual(x1, 2 * T.x_of_b(b), places=12)
            y_on_line = y0 + (y1 - y0) * (2 * x) / x1
            self.assertAlmostEqual(y_on_line, math.log(math.sinh(x / 2)), places=6, msg=f"A={a} B={b}")

    @unittest.skipUnless(HAVE_MPL, "matplotlib not installed")
    def test_draw_without_and_with_a_result(self):
        fig = Figure(figsize=(6, 4))
        ax = fig.add_subplot(111)
        T.draw_nomogram(ax)
        n_plain = len(ax.lines)
        self.assertFalse(any(l.get_color() == "red" for l in ax.lines))
        r = T.thickness(io=1000.0, is_=3000.0, so=10.0, ss=0.9, be_o=707.0, be_s=100.0, source_ev=1486.6,
                        lam_nm=1.2, theta_normal_deg=0.0)
        fig2 = Figure(figsize=(6, 4))
        ax2 = fig2.add_subplot(111)
        T.draw_nomogram(ax2, r)
        red = [l for l in ax2.lines if l.get_color() == "red"]
        self.assertEqual(len(red), 2)                       # marker + solution line
        self.assertEqual(len(ax2.lines), n_plain + 2)
        fig2.canvas.draw()                                   # renders without error


if __name__ == "__main__":
    unittest.main()
