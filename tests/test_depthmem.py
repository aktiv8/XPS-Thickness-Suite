"""
Maximum-entropy depth-profile tests (depthmem.py): the model, the entropy, the gradient, the fit on the synthetic structures of
Scorciapino et al. (2009), the composition-dependent attenuation lengths and the G-1 equation of Gries (1996). Synthetic data only.

    python -m unittest discover -s tests -v
"""
import math
import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

import depthmem as MM  # noqa: E402

Z = MM.depths(150, 0.1)                                   # the paper's 150 layers of 1 angstrom


def _structure(name):
    st = MM.PAPER_STRUCTURES[name]
    th = np.array(st["thickness"])
    return list(st["species"]), th, st["compositions"], st["lam"]


def _layer_errors(fit, th, comp, margin=0.05):
    """Largest absolute composition error within each true layer (away from its edges)."""
    edges = np.concatenate([[0.0], np.cumsum(th)])
    out = []
    for li in range(len(th)):
        m = (Z >= edges[li] + margin) & (Z < edges[li + 1] - margin)
        if m.any():
            true = np.array(comp[li], dtype=float) / np.sum(comp[li])
            out.append(float(np.max(np.abs(fit.n[m].mean(axis=0) - true))))
    return out


class Model(unittest.TestCase):
    def test_lambda_matrix_scales_with_the_electron_energy_to_the_three_quarters(self):
        lm = MM.lambda_matrix([3.0, 2.0], [1000.0, 500.0])
        self.assertAlmostEqual(lm[0, 0], 3.0)
        self.assertAlmostEqual(lm[1, 1], 2.0)
        self.assertAlmostEqual(lm[0, 1], 2.0 * (1000.0 / 500.0) ** 0.75)
        self.assertAlmostEqual(lm[1, 0], 3.0 * (500.0 / 1000.0) ** 0.75)
        self.assertTrue(np.allclose(MM.lambda_matrix([3.0, 2.0], None), [[3.0, 2.0], [3.0, 2.0]]))

    def test_the_layer_length_is_the_own_one_in_a_pure_layer_and_the_thickness_weighted_mean_in_a_bilayer(self):
        lm = MM.lambda_matrix([4.0, 2.0])
        pure = np.tile([1.0, 0.0], (6, 1))
        self.assertTrue(np.allclose(MM.layer_lambda(pure, lm), 4.0))
        bilayer = np.vstack([np.tile([1.0, 0.0], (2, 1)), np.tile([0.0, 1.0], (2, 1))])      # 2 layers of A, then 2 of B
        lam = MM.layer_lambda(bilayer, lm)
        self.assertTrue(np.allclose(lam[:, :2], 4.0))
        self.assertAlmostEqual(lam[0, 2], (4 + 4 + 2) / 3)
        self.assertAlmostEqual(lam[0, 3], (4 + 4 + 2 + 2) / 4)
        self.assertAlmostEqual(lam[1, 3], (4 + 4 + 2 + 2) / 4)                                # same materials, same mean

    def test_the_apparent_concentration_matches_a_hand_calculation(self):
        # A in layer 0, B in layer 1 and below (z = 0, t, 2t, ...), lambda A = B = 2 nm, t = 1 nm, theta = 0
        z = MM.depths(40, 1.0)
        c = np.zeros((40, 2))
        c[0, 0], c[1:, 1] = 1.0, 1.0
        k = MM.kernel(z, np.full((2, 40), 2.0), [0.0])
        a_sum = sum(math.exp(-i / 2.0) for i in range(40))
        ra = 1.0 / a_sum
        rb = (a_sum - 1.0) / a_sum
        x = MM.apparent(c, k)
        self.assertAlmostEqual(x[0, 0], ra / (ra + rb), places=12)
        self.assertAlmostEqual(x[0, 0] + x[1, 0], 1.0, places=12)

    def test_the_entropy_is_zero_at_the_prior_and_negative_elsewhere(self):
        m = MM.layered_model([1.0], [[1, 0], [0.3, 0.7]], Z, floor=0.02)
        self.assertAlmostEqual(MM.entropy(m, m), 0.0, places=12)
        rng = np.random.default_rng(0)
        for _ in range(5):
            c = np.maximum(m + 0.1 * rng.random(m.shape), 1e-6)
            c = c / c.sum(axis=1, keepdims=True)
            self.assertLess(MM.entropy(c, m), 0.0)

    def test_a_layered_model_has_the_compositions_in_the_right_layers(self):
        n0 = MM.layered_model([0.5, 0.5], [[2, 2, 0], [0, 1, 0], [0, 0, 5]], Z)
        self.assertTrue(np.allclose(n0[0], [0.5, 0.5, 0]))
        self.assertTrue(np.allclose(n0[7], [0, 1, 0]))
        self.assertTrue(np.allclose(n0[-1], [0, 0, 1]))
        with self.assertRaises(ValueError):
            MM.layered_model([1.0], [[1, 0]], Z)


class Gradient(unittest.TestCase):
    def test_the_analytic_gradient_agrees_with_finite_differences(self):
        rng = np.random.default_rng(1)
        s, n = 3, 30
        z = MM.depths(n, 0.2)
        ang = np.array([0.0, 20.0, 40.0, 55.0, 60.0])
        lm = MM.lambda_matrix([3.0, 2.0, 1.5], [1400.0, 900.0, 1200.0])
        c0 = MM.layered_model([1.0], [[1, 0, 0], [0, 0.5, 0.5]], z, floor=0.02)
        k = MM.kernel(z, MM.layer_lambda(c0, lm), ang)
        xm = MM.apparent(MM.layered_model([1.5], [[0.7, 0.3, 0], [0, 0.2, 0.8]], z, floor=0.02), k)
        for fixed in (None, c0[-5:]):
            prob = MM.MemProblem(k, xm, 0.05 * xm, 3.0, np.maximum(c0, 1e-4), fixed)
            v = rng.normal(size=prob.n_free * s)
            _, g = prob.value_and_grad(v)
            gn = np.empty_like(g)
            for i in range(len(v)):
                vp, vm = v.copy(), v.copy()
                vp[i] += 1e-6
                vm[i] -= 1e-6
                gn[i] = (prob.value(vp) - prob.value(vm)) / 2e-6
            self.assertLess(float(np.max(np.abs(g - gn))), 1e-5 * max(float(np.max(np.abs(g))), 1.0))


class Fit(unittest.TestCase):
    def test_the_result_sums_to_one_in_every_layer_and_is_never_negative(self):
        names, th, comp, lam = _structure("4_2+2")
        x = MM.acd_from_layers(th, comp, lam, MM.PAPER_ANGLES)
        fit = MM.fit_mem(names, MM.PAPER_ANGLES, x, lam, MM.layered_model(th * 1.2, comp, Z, floor=0.01), sigma=0.02,
                         iterate_imfp=False)
        self.assertTrue(np.allclose(fit.n.sum(axis=1), 1.0))
        self.assertGreaterEqual(float(fit.n.min()), 0.0)
        self.assertEqual(len(fit.cycles), 5)
        self.assertEqual([c["alpha"] for c in fit.cycles], list(MM.DEFAULT_ALPHAS))

    def test_a_huge_alpha_returns_the_starting_model(self):
        names, th, comp, lam = _structure("4_2+2")
        x = MM.acd_from_layers(th, comp, lam, MM.PAPER_ANGLES)
        start = MM.layered_model(th * 1.3, ((0.5, 0.5, 0, 0), (0, 0, 0.5, 0.5)), Z, floor=0.01)
        fit = MM.fit_mem(names, MM.PAPER_ANGLES, x, lam, start, sigma=0.02, iterate_imfp=False, alphas=(1e7,))
        self.assertLess(float(np.max(np.abs(fit.n - start))), 5e-3)

    def test_a_small_alpha_fits_the_data_from_a_wrong_start(self):
        names, th, comp, lam = _structure("4_2+2")
        x = MM.acd_from_layers(th, comp, lam, MM.PAPER_ANGLES)
        start = MM.layered_model(th * 1.3, ((0.5, 0.5, 0, 0), (0, 0, 0.5, 0.5)), Z, floor=0.01)
        wrong = MM.fit_mem(names, MM.PAPER_ANGLES, x, lam, start, sigma=0.02, iterate_imfp=False, alphas=(1e7,))
        fit = MM.fit_mem(names, MM.PAPER_ANGLES, x, lam, start, sigma=0.02, iterate_imfp=False,
                         alphas=(100, 10, 1, 0.1, 0.01))
        self.assertLess(fit.chi2, 0.01 * wrong.chi2)

    def test_the_stopping_rule_takes_the_largest_alpha_that_fits_within_the_noise(self):
        names, th, comp, lam = _structure("4_2+2")
        x = MM.add_acd_error(MM.acd_from_layers(th, comp, lam, MM.PAPER_ANGLES), 0.10, seed=3)
        start = MM.layered_model(th * 1.1, comp, Z, floor=0.01)
        fit = MM.fit_mem(names, MM.PAPER_ANGLES, x, lam, start, sigma=0.06, iterate_imfp=False,
                         alphas=(100, 30, 10, 3, 1, 0.3, 0.1), target_chi2_per_point=1.0)
        self.assertLessEqual(fit.chi2 / fit.n_points, 1.0)
        self.assertLess(len(fit.cycles), 7)
        self.assertTrue(all(c["chi2_per_point"] > 1.0 for c in fit.cycles[:-1]))

    def test_a_known_bulk_below_the_free_depth_is_held_at_the_starting_composition(self):
        names, th, comp, lam = _structure("8")
        x = MM.acd_from_layers(th, comp, lam, MM.PAPER_ANGLES)
        start = MM.layered_model(th, comp, Z, floor=0.01)
        fit = MM.fit_mem(names, MM.PAPER_ANGLES, x, lam, start, sigma=0.03, iterate_imfp=False, free_depth_nm=8.0)
        deep = fit.n[Z >= 8.0]
        self.assertTrue(np.allclose(deep, fit.prior[Z >= 8.0]))

    def test_bad_inputs_are_refused(self):
        names, th, comp, lam = _structure("4_2+2")
        x = MM.acd_from_layers(th, comp, lam, MM.PAPER_ANGLES)
        start = MM.layered_model(th, comp, Z)
        with self.assertRaises(ValueError):
            MM.fit_mem(names, MM.PAPER_ANGLES, x[:3], lam, start)
        with self.assertRaises(ValueError):
            MM.fit_mem(names, MM.PAPER_ANGLES, x, lam, start, alphas=(0.0,))
        with self.assertRaises(ValueError):
            MM.fit_mem(names, MM.PAPER_ANGLES, x, lam, start, lam_mat=np.ones((2, 2)))


class PaperStructures(unittest.TestCase):
    """Numerical experiments of Scorciapino et al.: start models within 10 % of the truth, as the paper's simulator gives."""

    def test_the_4_2_plus_2_structure_is_recovered_with_and_without_error(self):
        names, th, comp, lam = _structure("4_2+2")
        x0 = MM.acd_from_layers(th, comp, lam, MM.PAPER_ANGLES)
        for tag, x, sg in (("clean", x0, 0.02), ("error", MM.add_acd_error(x0, 0.10, seed=3), 0.06)):
            for scale in (0.9, 1.0, 1.1):
                fit = MM.fit_mem(names, MM.PAPER_ANGLES, x, lam, MM.layered_model(th * scale, comp, Z, floor=0.01),
                                 sigma=sg, iterate_imfp=False, free_depth_nm=8.0, alphas=(100, 30, 10, 3, 1, 0.3, 0.1),
                                 target_chi2_per_point=None if tag == "clean" else 1.0)
                over = fit.crossing_depth_nm([0, 1])
                self.assertAlmostEqual(over / th[0], 1.0, delta=0.12, msg=f"{tag} x{scale}: overlayer {over:.2f} nm")
                self.assertLess(_layer_errors(fit, th, comp)[0], 0.06, msg=f"{tag} x{scale}")

    def test_the_8_component_structure_is_recovered_within_the_papers_limits(self):
        names, th, comp, lam = _structure("8")
        x0 = MM.acd_from_layers(th, comp, lam, MM.PAPER_ANGLES)
        for tag, x, sg in (("clean", x0, 0.02), ("error", MM.add_acd_error(x0, 0.10, seed=3), 0.06)):
            for scale in (0.9, 1.1):
                fit = MM.fit_mem(names, MM.PAPER_ANGLES, x, lam, MM.layered_model(th * scale, comp, Z, floor=0.01),
                                 sigma=sg, iterate_imfp=False, free_depth_nm=8.0, alphas=(100, 30, 10, 3, 1, 0.3, 0.1),
                                 target_chi2_per_point=None if tag == "clean" else 1.0)
                over = fit.crossing_depth_nm([1, 7])                                    # the B:H overlayer
                self.assertAlmostEqual(over / th[0], 1.0, delta=0.12, msg=f"{tag} x{scale}: overlayer {over:.2f} nm")
                self.assertLess(max(_layer_errors(fit, th, comp)), 0.16, msg=f"{tag} x{scale}")

    def test_the_7_b_structure_is_recovered_from_the_exact_start_and_roughly_from_a_close_one(self):
        names, th, comp, lam = _structure("7_b")
        x0 = MM.acd_from_layers(th, comp, lam, MM.PAPER_ANGLES)
        exact = MM.fit_mem(names, MM.PAPER_ANGLES, x0, lam, MM.layered_model(th, comp, Z, floor=0.01), sigma=0.02,
                           iterate_imfp=False, free_depth_nm=8.0, alphas=(100, 30, 10, 3, 1, 0.3, 0.1))
        self.assertLess(max(_layer_errors(exact, th, comp)), 0.06)
        noisy = MM.add_acd_error(x0, 0.10, seed=3)
        close = MM.fit_mem(names, MM.PAPER_ANGLES, noisy, lam, MM.layered_model(th * 0.9, comp, Z, floor=0.01), sigma=0.06,
                           iterate_imfp=False, free_depth_nm=8.0, alphas=(100, 30, 10, 3, 1, 0.3, 0.1),
                           target_chi2_per_point=1.0)
        errs = _layer_errors(close, th, comp)
        self.assertLess(max(errs[:2]), 0.15)                      # the 0.6 nm pure-A layer is near the resolution limit

    def test_a_noise_ensemble_over_the_starting_thicknesses_has_a_positive_spread(self):
        names, th, comp, lam = _structure("4_2+2")
        x = MM.acd_from_layers(th, comp, lam, MM.PAPER_ANGLES)
        z, mean, sd, fits = MM.ensemble_mem(names, MM.PAPER_ANGLES, x, lam, th, comp, spread=0.3, n=4, seed=2,
                                            sigma=0.02, iterate_imfp=False, alphas=(10, 1, 0.1))
        self.assertEqual(len(fits), 4)
        self.assertEqual(mean.shape, sd.shape)
        self.assertGreater(float(sd.max()), 0.0)
        self.assertTrue(np.allclose(mean.sum(axis=1), 1.0))


class ImfpIteration(unittest.TestCase):
    def test_lengths_that_follow_the_composition_are_recovered_only_by_the_iteration(self):
        names, th, comp, lam = _structure("4_2+2")
        x_fixed = MM.acd_from_layers(th, comp, lam, MM.PAPER_ANGLES, iterate_imfp=False)
        x = MM.acd_from_layers(th, comp, lam, MM.PAPER_ANGLES, iterate_imfp=True)
        self.assertGreater(float(np.max(np.abs(x - x_fixed))), 0.05)                     # the dependence is not small
        start = MM.layered_model(th * 1.1, comp, Z, floor=0.01)
        kw = dict(sigma=0.02, alphas=(100, 30, 10, 3, 1, 0.3, 0.1))
        fixed = MM.fit_mem(names, MM.PAPER_ANGLES, x, lam, start, iterate_imfp=False, **kw)
        itr = MM.fit_mem(names, MM.PAPER_ANGLES, x, lam, start, iterate_imfp=True, **kw)
        self.assertLess(abs(itr.crossing_depth_nm([0, 1]) - th[0]), abs(fixed.crossing_depth_nm([0, 1]) - th[0]))
        self.assertLess(abs(itr.crossing_depth_nm([0, 1]) - th[0]) / th[0], 0.03)
        self.assertLess(itr.chi2, fixed.chi2)
        self.assertTrue(all(c["imfp_iterations"] <= MM.MAX_IMFP_ITERATIONS for c in itr.cycles))
        self.assertTrue(any(c["imfp_iterations"] > 1 for c in itr.cycles))               # it did iterate, and it settled
        self.assertGreater(float(np.ptp(itr.lam_ji[0])), 0.5)                            # the length varies with depth


class G1(unittest.TestCase):
    def test_water_reproduces_table_3_of_gries(self):
        # E (eV) -> lambda G1 (nm) printed in Table 3
        table = {200: 1.42, 300: 1.81, 400: 2.18, 500: 2.53, 600: 2.88, 800: 3.54, 1000: 4.17, 1500: 5.66, 2000: 7.08}
        water = MM.Material("H2O", 1.0)
        for e, want in table.items():
            self.assertAlmostEqual(MM.g1_imfp_nm(e, water), want, delta=0.02, msg=f"{e} eV")

    def test_the_equation_by_hand_for_silicon_dioxide(self):
        z_star = (math.sqrt(14) + 2 * math.sqrt(8)) / 3
        v_a = (28.086 + 2 * 15.999) / (2.196 * 3)
        want = 0.0019 * (v_a / z_star) * 1387.0 / (math.log10(1387.0) - 1.30)
        self.assertAlmostEqual(MM.g1_imfp_nm(1387.0, MM.Material("SiO2", 2.196)), want, places=10)

    def test_the_material_class_follows_the_periodic_table_and_the_formula(self):
        cls = lambda f, cat="auto": MM.g1_class(MM.parse_formula(f), cat)             # noqa: E731
        self.assertEqual(cls("Cu"), "3d")
        self.assertEqual(cls("Ag"), "4d")
        self.assertEqual(cls("W"), "5d")
        self.assertEqual(cls("Si"), "main")
        self.assertEqual(cls("SiO2"), "inorganic")
        self.assertEqual(cls("CH2"), "organic")
        self.assertEqual(cls("CH2", "inorganic"), "inorganic")
        self.assertEqual(cls("Al2O3", "organic"), "organic")

    def test_density_scales_the_length_inversely_and_energy_raises_it(self):
        a, b = MM.Material("Al2O3", 3.97), MM.Material("Al2O3", 1.985)
        self.assertAlmostEqual(MM.g1_imfp_nm(1000.0, b) / MM.g1_imfp_nm(1000.0, a), 2.0, places=10)
        self.assertGreater(MM.g1_imfp_nm(1500.0, a), MM.g1_imfp_nm(1000.0, a))

    def test_formulas_are_parsed_and_bad_ones_refused(self):
        self.assertEqual(MM.parse_formula("Al2O3"), [("Al", 2.0), ("O", 3.0)])
        self.assertEqual(MM.parse_formula("SiO1.5"), [("Si", 1.0), ("O", 1.5)])
        for bad in ("", "xyz", "Qq2", "Si0", "SiO2)"):
            with self.assertRaises(ValueError):
                MM.parse_formula(bad)
        with self.assertRaises(ValueError):
            MM.g1_imfp_nm(20.0, MM.Material("Si", 2.33))
        with self.assertRaises(ValueError):
            MM.g1_imfp_nm(1000.0, MM.Material("Si", 0.0))

    def test_the_matrix_has_the_energy_of_the_row_and_the_material_of_the_column(self):
        mats = [MM.Material("Si", 2.33), MM.Material("SiO2", 2.196)]
        lm = MM.g1_lambda_matrix([1400.0, 900.0], mats)
        self.assertAlmostEqual(lm[0, 1], MM.g1_imfp_nm(1400.0, mats[1]))
        self.assertAlmostEqual(lm[1, 0], MM.g1_imfp_nm(900.0, mats[0]))
        with self.assertRaises(ValueError):
            MM.g1_lambda_matrix([1400.0, None], mats)

    def test_the_range_warning_names_the_energies_outside_the_fitted_range(self):
        self.assertIsNone(MM.g1_range_warning([300.0, 1500.0]))
        self.assertIn("2985", MM.g1_range_warning([300.0, 2985.0]))

    def test_a_fit_can_take_the_g1_matrix(self):
        names, th, comp, lam = _structure("4_2+2")
        mats = [MM.Material("C6H8O2", 1.19), MM.Material("SiO2", 2.196), MM.Material("Al2O3", 3.97), MM.Material("Si", 2.33)]
        ke = [1200.0, 1100.0, 1400.0, 1380.0]
        lm = MM.g1_lambda_matrix(ke, mats)
        x = MM.acd_from_layers(th, comp, np.diag(lm), MM.PAPER_ANGLES, ke_ev=ke, iterate_imfp=False)
        fit = MM.fit_mem(names, MM.PAPER_ANGLES, x, np.diag(lm), MM.layered_model(th * 1.1, comp, Z, floor=0.01),
                         sigma=0.02, iterate_imfp=False, lam_mat=lm)
        self.assertTrue(np.allclose(fit.lam_nm, np.diag(lm)))
        self.assertAlmostEqual(fit.crossing_depth_nm([0, 1]) / th[0], 1.0, delta=0.1)


if __name__ == "__main__":
    unittest.main()
