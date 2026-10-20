"""
DAXPES maths tests (daxpes.py): the thickness of a film from the peak areas of one core level measured with a soft and a
hard X-ray source (Murdoch 2026), its coverage independence, the coverage, error handling and uncertainty propagation.
Synthetic numbers only; the attenuation-length pairs are those of Murdoch's Table 1.

    python -m unittest discover -s tests -v
"""
import math
import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

import daxpes as DX  # noqa: E402

# (soft Al Ka, hard Ag La) attenuation lengths in nm from Murdoch's Table 1
PAIRS = {"Ti 2p": (1.87, 3.77), "Ga 2p": (0.84, 2.73), "Sn 3d": (1.69, 3.93), "Si 2p": (2.94, 5.29)}


class Forward(unittest.TestCase):
    def test_a_thick_film_gives_one_and_a_thin_one_is_proportional_to_the_thickness(self):
        self.assertAlmostEqual(DX.film_factor(1000.0, 2.0, 0.0), 1.0, places=12)
        self.assertAlmostEqual(DX.film_factor(1e-6, 2.0, 0.0), 0.5e-6, places=12)

    def test_the_forward_ratios_carry_the_coverage_as_a_common_factor(self):
        rx, rh = DX.forward_ratios(3.0, 0.6, 1.87, 3.77, 30.0)
        fx, fh = DX.forward_ratios(3.0, 1.0, 1.87, 3.77, 30.0)
        self.assertAlmostEqual(rx, 0.6 * fx, places=12)
        self.assertAlmostEqual(rh / rx, fh / fx, places=12)


class Thickness(unittest.TestCase):
    def test_equation_6_recovers_the_thickness_for_every_table_1_pair_and_angle(self):
        for name, (lx, lh) in PAIRS.items():
            for theta in (0.0, 30.0, 45.0, 60.0):
                c = math.cos(math.radians(theta))
                for frac in (0.3, 0.7, 1.0, 1.5, 2.0):
                    d = frac * lx * c
                    rx, rh = DX.forward_ratios(d, 1.0, lx, lh, theta)
                    got = DX.daxpes_thickness(rh / rx, lx, lh, theta)
                    self.assertAlmostEqual(got, d, delta=1e-9 * max(d, 1.0), msg=f"{name} {theta} deg d={d:.3g}")

    def test_the_coverage_cancels_in_the_ratio_so_d_is_the_same_for_any_coverage(self):
        lx, lh = PAIRS["Ti 2p"]
        base = None
        for cov in (1.0, 0.8, 0.5, 0.1):
            i_x, i_h = DX.forward_ratios(2.5, cov, lx, lh, 0.0)
            res = DX.analyse(i_x * 1000.0, i_h * 5000.0, lx, lh, 0.0, iinf_soft=1000.0, iinf_hard=5000.0)
            base = base or res.d_nm
            self.assertAlmostEqual(res.d_nm, base, places=9)
            self.assertAlmostEqual(res.d_nm, 2.5, places=8)
            self.assertAlmostEqual(res.coverage, cov, places=9)

    def test_only_the_reference_ratio_is_needed_for_the_thickness(self):
        lx, lh = PAIRS["Ga 2p"]
        i_x, i_h = DX.forward_ratios(1.2, 0.4, lx, lh, 0.0)
        res = DX.analyse(i_x * 800.0, i_h * 2400.0, lx, lh, 0.0, rho=3.0)
        self.assertAlmostEqual(res.d_nm, 1.2, places=8)
        self.assertIsNone(res.coverage)
        self.assertIsNone(res.t_soft_nm)

    def test_single_source_thickness_is_eq_3_inverted_and_underestimates_a_patchy_film(self):
        lx, lh = PAIRS["Si 2p"]
        d, cov = 4.0, 0.5
        rx, rh = DX.forward_ratios(d, cov, lx, lh, 0.0)
        self.assertAlmostEqual(DX.single_source_thickness(DX.film_factor(d, lx, 0.0), lx), d, places=10)
        self.assertLess(DX.single_source_thickness(rx, lx), d)
        self.assertLess(DX.single_source_thickness(rh, lh), d)

    def test_the_paper_coverage_approaches_the_exact_one_for_thin_films(self):
        lx, lh = PAIRS["Ti 2p"]
        for cov in (0.3, 0.6, 0.9):
            i_x, i_h = DX.forward_ratios(0.2, cov, lx, lh, 0.0)
            res = DX.analyse(i_x, i_h, lx, lh, 0.0, iinf_soft=1.0, iinf_hard=1.0)
            self.assertAlmostEqual(res.coverage, cov, places=9)
            self.assertAlmostEqual(res.coverage_paper, cov, delta=0.03)

    def test_a_thick_film_warns_about_the_elastic_limit_and_a_saturated_soft_source(self):
        lx, lh = PAIRS["Ti 2p"]
        i_x, i_h = DX.forward_ratios(11.0, 1.0, lx, lh, 0.0)
        res = DX.analyse(i_x, i_h, lx, lh, 0.0, iinf_soft=1.0, iinf_hard=1.0)
        self.assertAlmostEqual(res.limit_nm, 3 * lh, places=9)
        self.assertTrue(any("elastic limit" in w for w in res.warnings))
        self.assertTrue(any("saturated" in w for w in res.warnings))


class LambdaBand(unittest.TestCase):
    """How far the answer moves if the attenuation lengths are off by +-10 %."""
    LX, LH = 1.87, 3.77

    def _res(self, d=3.0, cov=0.6, **kw):
        rx, rh = DX.forward_ratios(d, cov, self.LX, self.LH, 0.0)
        return DX.analyse(rx * 1000.0, rh * 1000.0, self.LX, self.LH, 0.0, iinf_soft=1000.0, iinf_hard=1000.0, **kw)

    def test_changing_both_lengths_together_scales_the_thickness_and_leaves_the_coverage(self):
        res = self._res()
        lo, hi, c_lo, c_hi = res.lam_band["common"]
        self.assertAlmostEqual(lo, 0.9 * 3.0, places=9)
        self.assertAlmostEqual(hi, 1.1 * 3.0, places=9)
        self.assertAlmostEqual(c_lo, 0.6, places=9)
        self.assertAlmostEqual(c_hi, 0.6, places=9)

    def test_changing_the_hard_length_alone_moves_both_the_thickness_and_the_coverage(self):
        res = self._res()
        lo, hi, c_lo, c_hi = res.lam_band["ratio"]
        self.assertLess(lo, res.d_nm)
        self.assertGreater(hi, res.d_nm)
        self.assertGreater(c_lo, c_hi)                          # a thicker film means fewer covered atoms for the same signal
        self.assertGreater(hi - lo, 0.2 * res.d_nm)             # the ratio matters more than the common scale

    def test_a_variant_with_no_solution_is_reported_as_none_and_band_zero_switches_it_off(self):
        low, _ = DX.ratio_limits(self.LX, self.LH)
        res = DX.analyse(1.0, 1.0 * 1.05 * low, self.LX, self.LH, 0.0, rho=1.0)
        self.assertIsNone(res.lam_band["ratio"])                # lowering lambda_hard lifts the thin-film limit above q
        self.assertIsNotNone(res.lam_band["common"])
        self.assertEqual(DX.analyse(1.0, 1.05 * low, self.LX, self.LH, 0.0, rho=1.0, band=0.0).lam_band, {})

    def test_without_reference_areas_there_is_no_coverage_in_the_band(self):
        rx, rh = DX.forward_ratios(2.0, 1.0, self.LX, self.LH, 0.0)
        res = DX.analyse(rx, rh, self.LX, self.LH, 0.0, rho=1.0)
        self.assertIsNone(res.lam_band["common"][2])


class Errors(unittest.TestCase):
    def test_the_hard_source_must_have_the_longer_attenuation_length(self):
        with self.assertRaises(ValueError):
            DX.ratio_limits(3.0, 2.0)

    def test_ratios_outside_the_resolvable_range_are_refused_with_a_reason(self):
        lx, lh = PAIRS["Ti 2p"]
        low, high = DX.ratio_limits(lx, lh)
        with self.assertRaisesRegex(ValueError, "thin-film limit"):
            DX.daxpes_thickness(low * 0.99, lx, lh)
        with self.assertRaisesRegex(ValueError, "above 1"):
            DX.daxpes_thickness(1.01, lx, lh)
        with self.assertRaises(ValueError):
            DX.daxpes_thickness(float("nan"), lx, lh)

    def test_a_reference_larger_than_the_sample_and_bad_inputs_are_refused(self):
        with self.assertRaises(ValueError):
            DX.single_source_thickness(1.0, 2.0)
        with self.assertRaises(ValueError):
            DX.single_source_thickness(-0.1, 2.0)
        with self.assertRaises(ValueError):
            DX.analyse(0.0, 1.0, 2.0, 4.0, rho=1.0)
        with self.assertRaises(ValueError):
            DX.analyse(1.0, 1.0, 2.0, 4.0)                       # neither the references nor their ratio

    def test_a_coverage_above_100_percent_is_flagged(self):
        lx, lh = PAIRS["Ti 2p"]
        i_x, i_h = DX.forward_ratios(2.0, 1.0, lx, lh, 0.0)
        res = DX.analyse(i_x * 1.3, i_h * 1.3, lx, lh, 0.0, iinf_soft=1.0, iinf_hard=1.0)
        self.assertGreater(res.coverage, 1.05)
        self.assertTrue(any("above 100 %" in w for w in res.warnings))


class Uncertainty(unittest.TestCase):
    def test_the_error_matches_a_monte_carlo_estimate(self):
        lx, lh = PAIRS["Ti 2p"]
        rng = np.random.default_rng(5)
        rx, rh = DX.forward_ratios(2.0, 0.7, lx, lh, 0.0)
        i_x, i_h, inf_x, inf_h = rx * 5000.0, rh * 3000.0, 5000.0, 3000.0
        rel = 0.01
        res = DX.analyse(i_x, i_h, lx, lh, 0.0, iinf_soft=inf_x, iinf_hard=inf_h, sig_soft=rel * i_x, sig_hard=rel * i_h,
                         sig_inf_soft=rel * inf_x, sig_inf_hard=rel * inf_h)
        ds, covs = [], []
        for _ in range(3000):
            n = rng.normal(1.0, rel, 4)
            r = DX.analyse(i_x * n[0], i_h * n[1], lx, lh, 0.0, iinf_soft=inf_x * n[2], iinf_hard=inf_h * n[3])
            ds.append(r.d_nm)
            covs.append(r.coverage)
        self.assertAlmostEqual(res.d_se_nm / float(np.std(ds)), 1.0, delta=0.1)
        self.assertAlmostEqual(res.coverage_se / float(np.std(covs)), 1.0, delta=0.1)

    def test_no_sigma_gives_no_error(self):
        lx, lh = PAIRS["Ti 2p"]
        i_x, i_h = DX.forward_ratios(2.0, 1.0, lx, lh, 0.0)
        res = DX.analyse(i_x, i_h, lx, lh, 0.0, iinf_soft=1.0, iinf_hard=1.0)
        self.assertIsNone(res.d_se_nm)
        self.assertIsNone(res.coverage_se)


class Helpers(unittest.TestCase):
    def test_the_sources_and_the_tpp2m_start_value(self):
        self.assertAlmostEqual(DX.SOURCES["Al Kα"], 1486.7)
        self.assertGreater(DX.SOURCES["Ag Lα"], 2 * DX.SOURCES["Mg Kα"])

        def fake(ke, rho, mass, nv, gap):
            return 0.001 * ke + rho / 10.0 + nv / 100.0 + gap / 100.0

        want = 0.9 * (0.001 * 1000.0 + 4.26 / 10.0 + 16 / 100.0 + 3.2 / 100.0)
        self.assertAlmostEqual(DX.tpp2m_start_nm(fake, "TiO2", 1000.0, 0.9), want, places=12)


class SeahCurve(unittest.TestCase):
    def test_eq_7_atomic_size_of_sio2_is_close_to_the_025_nm_average(self):
        self.assertAlmostEqual(DX.seah_a_nm(60.08, 2.196, 3), 0.2475, delta=0.001)       # Seah: a = 0.25 +- 0.03 nm

    def test_eq_11_by_hand_for_an_element(self):
        # Si, Z = 14, a = 0.272 nm, E = 1388 eV, W = 0
        z, a, e = 14.0, DX.seah_a_nm(28.09, 2.33, 1), 1388.0
        want = (5.8 + 0.0041 * z ** 1.7 + 0.088 * e ** 0.93) * a ** 1.82 / z ** 0.38
        self.assertAlmostEqual(DX.seah_al_nm(e, z, a), want, places=12)
        self.assertAlmostEqual(DX.seah_start_nm("Si", e), want, places=9)

    def test_w_lengthens_the_attenuation_length_and_the_hard_source_is_longer(self):
        a = DX.seah_a_nm(79.87, 4.26, 3)
        self.assertGreater(DX.seah_al_nm(1000.0, 12.7, a, 0.2), DX.seah_al_nm(1000.0, 12.7, a, 0.0))
        self.assertGreater(DX.seah_start_nm("TiO2", 2525.0), DX.seah_start_nm("TiO2", 1027.0))

    def test_it_agrees_with_murdochs_table_1_within_the_scatter_of_the_choices_made(self):
        table = {"TiO2": ((1027.0, 1.87), (2525.0, 3.77)), "Ga2O3": ((369.0, 0.84), (1867.0, 2.73)),
                 "SnO2": ((1000.0, 1.69), (2497.0, 3.93)), "Si": ((1388.0, 2.94), (2885.0, 5.29))}
        for material, pairs in table.items():
            for ke, want in pairs:
                got = DX.seah_start_nm(material, ke)
                self.assertAlmostEqual(got / want, 1.0, delta=0.2, msg=f"{material} KE {ke:.0f}: {got:.2f} vs {want}")

    def test_sio2_at_the_al_ka_energy_is_near_seahs_figure_7(self):
        self.assertAlmostEqual(DX.seah_start_nm("SiO2", 1387.0), 3.3, delta=0.3)

    def test_bad_inputs_are_refused(self):
        with self.assertRaises(ValueError):
            DX.seah_al_nm(0.0, 14.0, 0.25)
        with self.assertRaises(ValueError):
            DX.seah_al_nm(1000.0, 14.0, 0.25, 1.0)


class AttenuationStart(unittest.TestCase):
    def test_tpp2m_times_the_default_factor_is_close_to_murdochs_table_1(self):
        import app
        table = {"TiO2": ((1486.7 - 459.0, 1.87), (2984.3 - 459.0, 3.77)),
                 "Ga2O3": ((1486.7 - 1118.0, 0.84), (2984.3 - 1118.0, 2.73)),
                 "SnO2": ((1486.7 - 487.0, 1.69), (2984.3 - 487.0, 3.93)),
                 "Si": ((1486.7 - 99.0, 2.94), (2984.3 - 99.0, 5.29))}
        for material, pairs in table.items():
            for ke, want in pairs:
                got = DX.tpp2m_start_nm(app.tpp2m_imfp_nm, material, ke, 0.9)
                self.assertAlmostEqual(got / want, 1.0, delta=0.2, msg=f"{material} KE {ke:.0f}: {got:.2f} vs {want}")


if __name__ == "__main__":
    unittest.main()
