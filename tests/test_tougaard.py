"""
Tougaard background tests (tougaard.py): the cross-section, the kernel matrix, the recovery of B1 from spectra that satisfy
Eq. 1 exactly, errors, and the Monte Carlo uncertainty. Synthetic data only.

    python -m unittest discover -s tests -v
"""
import math
import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

import tougaard as TG  # noqa: E402

OXIDE = TG.CROSS_SECTIONS["Metal oxides (C = -1643, D = 0)"]
SILICON = TG.CROSS_SECTIONS["Silicon (C = 325, D = 96)"]


def _peak(x, centre, fwhm):
    sigma = fwhm / 2.3548
    return np.exp(-0.5 * ((x - centre) / sigma) ** 2) / (sigma * math.sqrt(2 * math.pi))


def _spectrum(b1, cd=OXIDE, centre=458.8, area=1000.0, baseline=5.0, step=0.1, lo=-12.0, hi=75.0):
    """A spectrum that satisfies Eq. 1 exactly: J = (1 - B1 M)^-1 F, plus a flat zero level."""
    x = np.arange(centre + lo, centre + hi + step / 2, step)
    m = TG.kernel_matrix(x, *cd)
    f = area * _peak(x, centre, 1.2)
    return x, np.linalg.solve(np.eye(len(x)) - b1 * m, f) + baseline, f, m


class Kernel(unittest.TestCase):
    def test_the_cross_section_is_zero_at_zero_and_peaks_near_the_plasmon_for_silicon(self):
        self.assertEqual(float(TG.cross_section(0.0, *OXIDE)), 0.0)
        t = np.linspace(1, 60, 600)
        k = TG.cross_section(t, *SILICON)
        self.assertAlmostEqual(float(t[int(np.argmax(k))]), 18.0, delta=2.0)             # sqrt(325) = 18 eV
        self.assertAlmostEqual(float(TG.cross_section(30.0, *OXIDE)), 30.0 / (1643.0 + 900.0) ** 2, places=12)

    def test_the_kernel_matrix_is_the_trapezoid_of_the_convolution(self):
        x = np.arange(0.0, 40.0, 0.5)
        y = np.exp(-0.5 * ((x - 10.0) / 2.0) ** 2)
        m = TG.kernel_matrix(x, *OXIDE)
        for i in (5, 30, 70):
            t = x[i] - x[:i + 1]
            want = float(np.trapezoid(y[:i + 1] * TG.cross_section(t, *OXIDE), x[:i + 1]))
            self.assertAlmostEqual(float((m @ y)[i]), want, places=12)
        self.assertTrue(np.allclose(np.triu(m), 0.0))                                 # only lower binding energies contribute

    def test_smoothing_keeps_a_straight_line_away_from_the_ends_and_the_length(self):
        y = np.linspace(0, 10, 50)
        s = TG.smooth(y)
        self.assertEqual(len(s), 50)
        self.assertTrue(np.allclose(s[3:-3], y[3:-3]))
        self.assertTrue(np.allclose(TG.smooth(y, 1), y))


class Recovery(unittest.TestCase):
    def test_b1_is_recovered_exactly_for_both_cross_sections(self):
        for cd in (OXIDE, SILICON):
            for b1 in (60.0, 300.0, 1000.0):
                x, j, f, m = _spectrum(b1, cd)
                fit = TG.fit_b1(x, j, (452.0, 466.0), *cd, smooth_n=1)
                self.assertAlmostEqual(fit.b1, b1, delta=1e-6 * b1)
                self.assertAlmostEqual(fit.area, 1000.0, delta=0.5)
                self.assertAlmostEqual(fit.centroid_be, 458.8, delta=0.01)

    def test_it_does_not_depend_on_the_zero_level_the_peak_area_or_the_energy_step(self):
        ref = TG.fit_b1(*_spectrum(500.0)[:2], (452.0, 466.0), *OXIDE, smooth_n=1).b1
        for kw in (dict(baseline=0.0), dict(baseline=40.0), dict(area=123.0), dict(step=0.2)):
            x, j, _, _ = _spectrum(500.0, **kw)
            self.assertAlmostEqual(TG.fit_b1(x, j, (452.0, 466.0), *OXIDE, smooth_n=1).b1, ref, delta=1e-3 * ref)

    def test_delta_e_between_15_and_30_ev_gives_the_same_b1_for_an_exact_spectrum(self):
        x, j, _, _ = _spectrum(700.0)
        for de in (15.0, 20.0, 30.0):
            self.assertAlmostEqual(TG.fit_b1(x, j, (452.0, 466.0), *OXIDE, delta_e=de, smooth_n=1).b1, 700.0, delta=0.01)
        self.assertAlmostEqual(TG.fit_b1(x, j, (452.0, 466.0), *OXIDE, eval_be=490.0, smooth_n=1).b1, 700.0, delta=0.01)

    def test_five_point_smoothing_changes_b1_only_slightly_for_a_wide_peak(self):
        x, j, _, _ = _spectrum(700.0)
        self.assertAlmostEqual(TG.fit_b1(x, j, (452.0, 466.0), *OXIDE).b1, 700.0, delta=700.0 * 0.01)

    def test_a_thicker_film_has_a_larger_b1_in_the_model_of_murdochs_eq_5(self):
        b0, lam = 1500.0, 1.87
        vals = [b0 * (1 - math.exp(-d / (2 * lam))) for d in (1.0, 3.0, 6.0)]
        fits = [TG.fit_b1(*_spectrum(v)[:2], (452.0, 466.0), *OXIDE, smooth_n=1).b1 for v in vals]
        self.assertTrue(fits[0] < fits[1] < fits[2])


class Band(unittest.TestCase):
    def test_the_band_fit_gives_the_same_b1_as_the_point_for_an_exact_spectrum(self):
        for b1 in (100.0, 700.0):
            x, j, _, _ = _spectrum(b1)
            fit = TG.fit_b1(x, j, (452.0, 466.0), *OXIDE, smooth_n=1, band=TG.DEFAULT_BAND)
            self.assertAlmostEqual(fit.b1, b1, delta=1e-3 * b1)
            lo, hi = fit.band_be
            self.assertAlmostEqual(lo, fit.centroid_be + 15.0, places=6)
            self.assertAlmostEqual(hi, fit.centroid_be + 45.0, places=6)
            self.assertAlmostEqual(fit.eval_be, 0.5 * (lo + hi), places=9)

    def test_the_band_fit_is_much_less_noisy_than_the_single_point(self):
        x, j, _, _ = _spectrum(700.0, area=20000.0, baseline=3.0)
        counts = j * 40.0
        window = (452.0, 466.0)
        point = TG.b1_sigma(x, counts, 1 / 40.0, window, *OXIDE, draws=60)
        band = TG.b1_sigma(x, counts, 1 / 40.0, window, *OXIDE, draws=60, band=TG.DEFAULT_BAND)
        self.assertLess(band, 0.6 * point)

    def test_a_band_that_runs_off_the_spectrum_is_refused(self):
        x, j, _, _ = _spectrum(500.0, hi=30.0)
        with self.assertRaisesRegex(ValueError, "wider range"):
            TG.fit_b1(x, j, (452.0, 466.0), *OXIDE, band=(15.0, 45.0))


class Sensitivity(unittest.TestCase):
    def test_the_table_lists_each_cross_section_class_delta_e_and_the_band(self):
        xs, js, _, _ = _spectrum(1200.0)
        xh, jh, _, _ = _spectrum(800.0)
        rows = TG.sensitivity_table(xs, js, xh, jh, (452.0, 466.0))
        labels = [r[0] for r in rows]
        self.assertEqual(len(rows), 8)                                              # 2 classes x (3 delta E + band)
        self.assertIn("Metal oxides, DeltaE 30 eV", labels)
        self.assertIn("Silicon, band 15-45 eV", labels)
        oxide = {r[0]: r for r in rows if r[0].startswith("Metal oxides")}
        for label, bs, bh in oxide.values():                                       # exact oxide spectra: every ΔE agrees
            self.assertAlmostEqual(bs, 1200.0, delta=5.0)
            self.assertAlmostEqual(bh, 800.0, delta=5.0)


class Errors(unittest.TestCase):
    def test_a_spectrum_that_stops_too_soon_is_refused_with_the_energy_it_needs(self):
        x, j, _, _ = _spectrum(500.0, hi=20.0)
        with self.assertRaisesRegex(ValueError, "wider range"):
            TG.fit_b1(x, j, (452.0, 466.0), *OXIDE)

    def test_no_tail_or_a_bad_window_is_refused(self):
        x = np.arange(440.0, 530.0, 0.1)
        flat = 5.0 + 1000.0 * _peak(x, 458.8, 1.2)                                  # a peak with no inelastic tail at all
        with self.assertRaisesRegex(ValueError, "inelastic tail"):
            TG.fit_b1(x, flat, (452.0, 466.0), *OXIDE, smooth_n=1)
        with self.assertRaises(ValueError):
            TG.fit_b1(x, flat, (300.0, 466.0), *OXIDE)
        with self.assertRaises(ValueError):
            TG.fit_b1(x[:10], flat[:10], (441.0, 442.0), *OXIDE)


class Noise(unittest.TestCase):
    def test_the_monte_carlo_sigma_matches_the_scatter_of_independent_noisy_fits(self):
        x, j, _, _ = _spectrum(700.0, area=20000.0, baseline=3.0)
        counts = j * 40.0                                                            # 40 counts per unit rate
        window = (452.0, 466.0)
        sigma = TG.b1_sigma(x, counts, 1 / 40.0, window, *OXIDE, draws=80)
        rng = np.random.default_rng(11)
        spread = np.std([TG.fit_b1(x, rng.poisson(counts) / 40.0, window, *OXIDE).b1 for _ in range(80)], ddof=1)
        self.assertGreater(sigma, 0.0)
        self.assertAlmostEqual(sigma / spread, 1.0, delta=0.35)


if __name__ == "__main__":
    unittest.main()
