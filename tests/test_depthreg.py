"""
Regularised depth-profile tests (depthreg.py): the kernel, Tikhonov limits, automatic alpha, non-negativity, bulk
constraints and the covariance bands. Synthetic data only.

    python -m unittest discover -s tests -v
"""
import math
import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

import depthprofile as DP  # noqa: E402
import depthreg as DR  # noqa: E402

ANG = np.array([0.0, 20.0, 40.0, 55.0, 65.0, 72.0, 78.0])
LAM, NZ, ZMAX = 3.0, 60, 12.0
Z, DZ = DR.depth_grid(ZMAX, NZ)


def _species(c_true, noise=0.01, seed=1, name="s", rsf=1.0):
    A = DR.kernel(ANG, LAM, NZ, DZ)
    clean = A @ c_true
    rng = np.random.default_rng(seed)
    data = clean * (1 + rng.normal(0, noise, len(clean))) if noise else clean.copy()
    return DP.Species(name, "r", ANG.copy(), data * rsf, sigma=max(noise, 1e-3) * clean * rsf, lam_nm=LAM, rsf=rsf), clean


def _rms(a, b):
    return float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))


class Pieces(unittest.TestCase):
    def test_the_kernel_is_the_slab_integral_of_the_attenuation(self):
        A = DR.kernel([0.0, 60.0], 3.0, 4, 0.5)
        for m, th in enumerate((0.0, 60.0)):
            lc = 3.0 * math.cos(math.radians(th))
            for n in range(4):
                zs = np.linspace(n * 0.5, (n + 1) * 0.5, 20001)
                self.assertAlmostEqual(A[m, n], float(np.trapezoid(np.exp(-zs / lc), zs)) if hasattr(np, "trapezoid")
                                       else float(np.trapz(np.exp(-zs / lc), zs)), places=6)

    def test_all_slabs_to_infinity_give_lambda_cos_theta(self):
        A = DR.kernel([0.0, 60.0], 3.0, 400, 0.1)
        self.assertAlmostEqual(A[0].sum(), 3.0, places=4)       # 40 nm is not quite infinite
        self.assertAlmostEqual(A[1].sum(), 1.5, places=4)

    def test_second_difference_ignores_a_straight_line_and_sees_a_kink(self):
        L = DR.second_difference(10)
        self.assertEqual(L.shape, (8, 10))
        line = 2.0 + 0.5 * np.arange(10)
        self.assertTrue(np.allclose(L @ line, 0.0))
        kink = np.abs(np.arange(10) - 5.0)
        self.assertGreater(np.abs(L @ kink).max(), 1.5)

    def test_nnls_equals_least_squares_when_it_is_positive_and_clamps_when_it_is_not(self):
        rng = np.random.default_rng(2)
        A = rng.normal(size=(30, 6))
        x_true = np.array([1.0, 0.0, 2.0, 0.5, 0.0, 3.0])
        self.assertTrue(np.allclose(DR.nnls(A, A @ x_true), x_true, atol=1e-8))
        b = A @ x_true + rng.normal(0, 3.0, 30)
        x = DR.nnls(A, b)
        self.assertTrue((x >= 0).all())
        ls = np.linalg.lstsq(A, b, rcond=None)[0]
        self.assertLessEqual(float(np.sum((A @ ls - b) ** 2)), float(np.sum((A @ x - b) ** 2)) + 1e-9)
        # the KKT conditions: the gradient is non-negative where x is zero and zero where x is positive
        g = A.T @ (A @ x - b)
        self.assertTrue((g[x == 0] >= -1e-8).all())
        self.assertTrue(np.allclose(g[x > 0], 0.0, atol=1e-6))


class Profiles(unittest.TestCase):
    def test_a_layer_and_an_exponential_profile_are_recovered(self):
        for c_true, limit in ((np.where(Z < 2.0, 1.0, 0.0), 0.25), (np.exp(-Z / 2.0), 0.06)):
            sp, _ = _species(c_true)
            f = DR.fit_profile(sp, ZMAX, 30, max_angle=80)
            self.assertLess(_rms(f.c, np.interp(f.z_nm, Z, c_true)), limit)
            self.assertLess(f.chi2_red, 3.0)

    def test_the_automatic_alpha_fits_within_the_noise_and_gcv_is_less_smooth(self):
        sp, _ = _species(np.exp(-Z / 2.0))
        disc = DR.fit_profile(sp, ZMAX, 30, max_angle=80)
        gcv = DR.fit_profile(sp, ZMAX, 30, max_angle=80, method="gcv")
        self.assertLess(disc.chi2_red, 3.0)
        self.assertLessEqual(gcv.alpha, disc.alpha)

    def test_no_regularisation_oscillates_and_a_huge_one_gives_a_straight_line(self):
        sp, _ = _species(np.exp(-Z / 2.0))
        wild = DR.fit_profile(sp, ZMAX, 30, alpha=1e-12, max_angle=80)
        self.assertGreater(np.ptp(wild.c), 20.0)
        flat = DR.fit_profile(sp, ZMAX, 30, alpha=1e12, max_angle=80)
        self.assertLess(float(np.abs(np.diff(flat.c, 2)).max()), 1e-6)

    def test_negative_parts_are_flagged_and_the_non_negative_option_removes_them(self):
        sp, _ = _species(np.where(Z < 2.0, 1.0, 0.0))
        free = DR.fit_profile(sp, ZMAX, 30, max_angle=80)
        self.assertTrue(free.negative)
        self.assertTrue(any("negative concentrations" in n for n in free.notes))
        pos = DR.fit_profile(sp, ZMAX, 30, nonneg=True, max_angle=80)
        self.assertGreaterEqual(float(pos.c.min()), 0.0)

    def test_holding_the_bulk_constant_or_zero_below_a_depth(self):
        c_true = np.where(Z < 2.0, 1.0, 0.0)
        sp, _ = _species(c_true)
        zero = DR.fit_profile(sp, ZMAX, 30, bulk=DR.BULK_ZERO, z_c_nm=4.0, max_angle=80)
        deep = zero.z_nm > 4.0
        self.assertTrue(np.all(zero.c[deep] == 0.0))
        self.assertLess(_rms(zero.c, np.interp(zero.z_nm, Z, c_true)), 0.18)           # the prior sharpens the result
        free = DR.fit_profile(sp, ZMAX, 30, max_angle=80)
        self.assertLess(_rms(zero.c, np.interp(zero.z_nm, Z, c_true)), _rms(free.c, np.interp(free.z_nm, Z, c_true)))
        const = DR.fit_profile(_species(np.exp(-Z / 50.0) * 0 + 1.0)[0], ZMAX, 30, bulk=DR.BULK_CONSTANT, z_c_nm=5.0,
                               max_angle=80)
        self.assertTrue(np.allclose(const.c[const.z_nm > 5.0], const.c[const.z_nm > 5.0][0]))

    def test_the_noise_bands_match_the_scatter_of_repeated_noisy_fits(self):
        c_true = np.exp(-Z / 2.0)
        A = DR.kernel(ANG, LAM, NZ, DZ)
        clean = A @ c_true
        rng = np.random.default_rng(5)
        first = None
        draws = []
        for _ in range(400):
            sp = DP.Species("s", "r", ANG.copy(), clean * (1 + rng.normal(0, 0.01, len(clean))), sigma=0.01 * clean,
                            lam_nm=LAM)
            f = DR.fit_profile(sp, ZMAX, 30, alpha=2e3, max_angle=80)
            first = first or f
            draws.append(f.c)
        sample_sd = np.std(draws, axis=0)
        self.assertTrue(np.all(first.c_se > 0))
        ratio = sample_sd[:15] / first.c_se[:15]
        self.assertTrue(np.all((ratio > 0.8) & (ratio < 1.2)), ratio)

    def test_sensitivity_factors_scale_the_concentration_not_the_shape(self):
        c_true = np.exp(-Z / 2.0)
        a, _ = _species(c_true, rsf=1.0)
        b, _ = _species(c_true, rsf=4.0)
        fa, fb = DR.fit_profile(a, ZMAX, 30, alpha=2e3, max_angle=80), DR.fit_profile(b, ZMAX, 30, alpha=2e3, max_angle=80)
        self.assertTrue(np.allclose(fa.c, fb.c, rtol=1e-6))

    def test_too_few_angles_is_an_error_and_the_atomic_fraction_adds_up_to_one(self):
        sp, _ = _species(np.exp(-Z / 2.0))
        sp.intensity[2:] = -1.0
        with self.assertRaisesRegex(ValueError, "at least three usable angles"):
            DR.fit_profile(sp, ZMAX, 30)
        a = DR.fit_profile(_species(np.exp(-Z / 2.0))[0], ZMAX, 30, max_angle=80)
        b = DR.fit_profile(_species(np.where(Z < 4.0, 0.3, 1.0))[0], ZMAX, 30, max_angle=80)
        frac = DR.atomic_fraction([a, b])
        self.assertTrue(np.allclose(frac.sum(axis=0), 1.0))


if __name__ == "__main__":
    unittest.main()
