"""
Depth-profile maths tests (depthprofile.py): species from fits and regions, the Surrey Trend Plot statistics, and the
robust Shirley background they depend on. Synthetic data only.

    python -m unittest discover -s tests -v
"""
import math
import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

import arxps as AR  # noqa: E402
import depthprofile as DP  # noqa: E402
import spectrum_sim as sim  # noqa: E402
import synthetic as S  # noqa: E402

ANGLES = np.array(S.DEFAULT_TILTS)


def _species(name, z_nm, lam_nm, ke=1000.0, i0=1000.0):
    sec = 1.0 / np.cos(np.radians(ANGLES))
    return DP.Species(name, "r", ANGLES.copy(), i0 * np.exp(-z_nm * sec / lam_nm), lam_nm=lam_nm, ke_ev=ke)


def _series(blocks, region=None):
    return [(b.exp_values[0], b) for b in blocks if region is None or b.block_id == region]


class TrendFit(unittest.TestCase):
    def test_a_delta_layer_gives_slope_minus_z_over_lambda(self):
        sp = _species("a", 2.0, 3.0)
        slope, se, icpt, rms, n = DP.trend_fit(sp.angles, sp.intensity, 80)
        self.assertAlmostEqual(slope, -2.0 / 3.0, places=10)
        self.assertAlmostEqual(icpt, math.log(1000.0), places=10)
        self.assertAlmostEqual(rms, 0.0, places=10)
        self.assertEqual(n, 6)
        row = DP.trend_table([sp], 80)[0]
        self.assertAlmostEqual(row.depth_nm, 2.0, places=9)

    def test_noise_gives_a_scatter_and_a_standard_error(self):
        rng = np.random.default_rng(3)
        sp = _species("a", 2.0, 3.0)
        sp.intensity = sp.intensity * np.exp(rng.normal(0, 0.03, len(ANGLES)))
        slope, se, _, rms, _ = DP.trend_fit(sp.angles, sp.intensity, 80)
        self.assertGreater(rms, 0.005)
        self.assertGreater(se, 0)
        self.assertAlmostEqual(slope, -2.0 / 3.0, delta=4 * se + 0.05)

    def test_angles_above_the_maximum_and_non_positive_intensities_are_left_out(self):
        sp = _species("a", 2.0, 3.0)
        sp.intensity[-1] *= 5                              # a bad grazing point
        slope, _, _, _, n = DP.trend_fit(sp.angles, sp.intensity, 70)
        self.assertEqual(n, 5)
        self.assertAlmostEqual(slope, -2.0 / 3.0, places=10)
        sp.intensity[0] = -3.0
        self.assertEqual(DP.trend_fit(sp.angles, sp.intensity, 70)[4], 4)

    def test_two_points_give_a_slope_but_no_scatter_and_one_gives_nothing(self):
        slope, se, icpt, rms, n = DP.trend_fit([0.0, 60.0], [100.0, 50.0], 70)
        self.assertEqual((n, se, rms), (2, None, None))
        self.assertAlmostEqual(slope, math.log(0.5) / (2.0 - 1.0), places=10)
        self.assertEqual(DP.trend_fit([0.0], [100.0], 70), (None, None, None, None, 1))

    def test_normalised_curves_start_at_one(self):
        sp = _species("a", 2.0, 3.0)
        self.assertAlmostEqual(DP.normalised(sp)[0], 1.0)


class CoLocation(unittest.TestCase):
    def _noisy(self, name, z, lam, ke, seed):
        sp = _species(name, z, lam, ke)
        rng = np.random.default_rng(seed)
        sp.intensity = sp.intensity * np.exp(rng.normal(0, 0.02, len(ANGLES)))
        return sp

    def test_species_at_the_same_depth_are_parallel_and_a_deeper_one_is_not(self):
        a, b = self._noisy("A", 1.0, 3.0, 1200, 1), self._noisy("B", 1.0, 3.0, 1100, 2)
        c = self._noisy("C", 4.0, 3.0, 1150, 3)
        rows = {r.name: r for r in DP.trend_table([a, b, c], 80)}
        self.assertEqual(rows["A"].similar, ["B"])
        self.assertEqual(rows["B"].similar, ["A"])
        self.assertEqual(rows["C"].similar, [])

    def test_very_different_kinetic_energies_are_not_called_similar_and_get_a_note(self):
        a, b = self._noisy("A", 1.0, 3.0, 1400, 1), self._noisy("B", 1.0, 3.0, 400, 2)
        rows = {r.name: r for r in DP.trend_table([a, b], 80)}
        self.assertEqual(rows["A"].similar, [])
        self.assertIn("kinetic energy differs", rows["A"].note)

    def test_a_single_point_fit_never_claims_a_match(self):
        a = _species("A", 1.0, 3.0)
        b = _species("B", 1.0, 3.0)
        a.intensity[1:] = -1.0
        rows = DP.trend_table([a, b], 80)
        self.assertEqual(rows[0].similar, [])


class Stratification(unittest.TestCase):
    SEC = 1.0 / np.cos(np.radians(ANGLES))

    def _delta(self, name, z, lam, n, rsf=1.0, sigma=None):
        sp = DP.Species(name, "r", ANGLES.copy(), n * rsf * np.exp(-z * self.SEC / lam), lam_nm=lam, rsf=rsf, ke_ev=1000.0)
        if sigma:
            sp.sigma = sp.intensity * sigma
        return sp

    def test_two_angles_recover_the_depth_and_the_amount_of_thin_layers(self):
        sps = [self._delta("A", 0.5, 3.0, 10.0), self._delta("B", 2.0, 3.0, 40.0, rsf=0.5),
               self._delta("C", 5.0, 3.5, 25.0)]
        rows = DP.stratification(sps, 0.0, 60.0)
        for r, z, n in zip(rows, (0.5, 2.0, 5.0), (10.0, 40.0, 25.0)):
            self.assertAlmostEqual(r.z_nm, z, places=9)
            self.assertAlmostEqual(r.amount, n, places=9)
        self.assertAlmostEqual(sum(r.amount_pct for r in rows), 100.0)
        self.assertAlmostEqual(rows[1].amount_pct, 100 * 40 / 75, places=9)

    def test_fitting_every_angle_gives_the_same_exact_answer(self):
        sps = [self._delta("A", 0.5, 3.0, 10.0), self._delta("B", 2.0, 3.0, 40.0)]
        for r, z in zip(DP.stratification(sps, use_all=True, max_angle=80), (0.5, 2.0)):
            self.assertAlmostEqual(r.z_nm, z, places=9)

    def test_the_default_angles_are_the_smallest_and_the_one_nearest_67_5_within_the_limit(self):
        sp = self._delta("A", 1.0, 3.0, 10.0)
        at_60 = math.exp(-1.0 * (2.0 - 1.0) / 3.0)
        at_75 = math.exp(-1.0 * (1 / math.cos(math.radians(75)) - 1) / 3.0)
        self.assertAlmostEqual(DP.stratification([sp], max_angle=80)[0].ratio, at_60, places=9)    # 60 and 75 tie: the lower
        self.assertAlmostEqual(DP.stratification([sp], max_angle=70)[0].ratio, at_60, places=9)
        self.assertAlmostEqual(DP.stratification([sp], theta_b=75.0, max_angle=80)[0].ratio, at_75, places=9)

    def test_a_species_that_rises_with_angle_is_at_the_surface_and_says_so(self):
        sp = DP.Species("S", "r", ANGLES.copy(), np.exp(+0.4 * self.SEC), lam_nm=3.0)
        r = DP.stratification([sp], 0.0, 60.0)[0]
        self.assertLess(r.z_nm, 0)
        self.assertIn("rises with angle", r.note)

    def test_a_species_in_two_layers_is_flagged_as_curved(self):
        two = DP.Species("T", "r", ANGLES.copy(), np.exp(-0.5 * self.SEC / 3) + np.exp(-6.0 * self.SEC / 3), lam_nm=3.0)
        r = DP.stratification([two], 0.0, 60.0)[0]
        self.assertGreater(abs(r.curved_t), 3)
        self.assertIn("more than one layer", r.note)
        self.assertTrue(0.5 < r.z_nm < 6.0)                                  # reported at one average depth, in between
        single = DP.stratification([self._delta("S", 2.0, 3.0, 1.0)], 0.0, 60.0)[0]
        self.assertEqual(single.note, "")

    def test_a_species_without_two_usable_angles_is_reported_not_crashed(self):
        sp = self._delta("A", 1.0, 3.0, 10.0)
        sp.intensity[:] = -1.0
        sp.intensity[0] = 5.0
        r = DP.stratification([sp])[0]
        self.assertIsNone(r.z_nm)
        self.assertIn("two usable angles", r.note)
        same = DP.stratification([self._delta("A", 1.0, 3.0, 10.0)], 0.0, 0.0)[0]
        self.assertIsNotNone(same.z_nm)                                      # the nearest OTHER angle is used

    def test_the_depth_uncertainty_from_counting_statistics_covers_the_truth(self):
        rng = np.random.default_rng(11)
        hits, trials = 0, 300
        for _ in range(trials):
            clean = self._delta("A", 2.0, 3.0, 5000.0)
            counts = rng.poisson(clean.intensity)
            sp = DP.Species("A", "r", ANGLES.copy(), counts.astype(float), sigma=np.sqrt(counts), lam_nm=3.0)
            r = DP.stratification([sp], 0.0, 60.0)[0]
            hits += abs(r.z_nm - 2.0) <= 2.0 * r.z_se_nm
        self.assertGreater(hits / trials, 0.90)

    def test_dividing_by_cos_theta_removes_a_geometry_factor_common_to_all_species(self):
        cos = np.cos(np.radians(ANGLES))
        sps = [self._delta("A", 0.5, 3.0, 10.0), self._delta("B", 2.0, 3.0, 40.0)]
        for sp in sps:
            sp.intensity = sp.intensity * cos                                # the analyser-defined-area projection
        wrong = DP.stratification(sps, 0.0, 60.0)
        right = DP.stratification(DP.apply_geometry(sps, 1.0), 0.0, 60.0)
        self.assertGreater(wrong[0].z_nm - 0.5, 0.5)
        for r, z in zip(right, (0.5, 2.0)):
            self.assertAlmostEqual(r.z_nm, z, places=9)
        self.assertIs(DP.apply_geometry(sps, 0.0)[0], sps[0])               # power 0 changes nothing


class SpeciesFromFits(unittest.TestCase):
    def test_states_become_species_with_the_exact_stored_areas(self):
        ser = _series(S.two_state_series(d_nm=1.5, fitted=True))
        sp = DP.species_from_fits(ser, transmission=False)
        self.assertEqual([s.name for s in sp], ["Si El", "Si Ox"])
        x = 1.5 / (3.89 * np.cos(np.radians(ANGLES)))
        self.assertTrue(np.allclose(sp[0].intensity, 60000.0 * np.exp(-x), rtol=1e-6))
        self.assertTrue(np.allclose(sp[1].intensity, 0.5547 * 60000.0 * (1 - np.exp(-x)), rtol=1e-6))
        self.assertAlmostEqual(sp[0].ke_ev, 1486.68 - (98.9 + 0.2), delta=0.2)        # doublet centroid, area-weighted over the drift
        self.assertEqual(sp[0].source, "fit")

    def test_the_substrate_state_follows_a_straight_line_in_sec_theta_and_the_overlayer_does_not(self):
        sp = DP.species_from_fits(_series(S.two_state_series(d_nm=1.5, fitted=True)), transmission=False)
        slope, _, _, rms, _ = DP.trend_fit(sp[0].angles, sp[0].intensity, 80)
        self.assertAlmostEqual(slope, -1.5 / 3.89, places=5)
        self.assertAlmostEqual(rms, 0.0, places=6)
        self.assertGreater(DP.trend_fit(sp[1].angles, sp[1].intensity, 80)[3], 0.05)

    def test_transmission_correction_divides_by_the_transmission_at_each_component(self):
        ser = _series(S.two_state_series(fitted=True))
        raw = DP.species_from_fits(ser, transmission=False)[0].intensity
        corr = DP.species_from_fits(ser, transmission=True)[0].intensity
        self.assertTrue(np.allclose(corr, raw / (S.TRANSMISSION / AR.TRANSMISSION_SCALE)))

    def test_angle_offset_and_a_block_without_a_fit(self):
        ser = _series(S.two_state_series(fitted=True))
        sp = DP.species_from_fits(ser, angle_offset=5.0)
        self.assertEqual(list(sp[0].angles), [a + 5.0 for a in S.DEFAULT_TILTS])
        with self.assertRaisesRegex(ValueError, "no CasaXPS fit"):
            DP.species_from_fits(_series(S.two_state_series()))

    def test_the_counting_uncertainty_is_step_times_root_counts_over_time(self):
        blk = S.two_state_series(noise=True)[0]
        step = float(np.mean(np.diff(blk.be_ev)))
        want = step * math.sqrt(float(blk.counts.sum())) / (blk.scans * blk.dwell_s)
        self.assertAlmostEqual(DP.area_sigma(blk), want, places=9)
        self.assertAlmostEqual(DP.area_sigma(blk, transmission_factor=2.0), want / 2.0, places=9)


class SpeciesFromRegions(unittest.TestCase):
    def test_a_whole_region_area_is_close_to_the_sum_of_its_states(self):
        blocks = S.two_state_series(d_nm=1.5)
        sp = DP.species_from_region(_series(blocks), transmission=False)
        x = 1.5 / (3.89 * np.cos(np.radians(ANGLES)))
        want = 60000.0 * np.exp(-x) + 0.5547 * 60000.0 * (1 - np.exp(-x))
        self.assertTrue(np.allclose(sp.intensity, want, rtol=0.06))
        self.assertEqual((sp.source, sp.name), ("region", "Si 2p"))

    def test_the_envelope_leaves_out_a_long_loss_tail(self):
        be = np.arange(90.0, 115.0, 0.1)
        y = 50000 * sim.pseudo_voigt(be, 99.0, 1.0, 0.0) + np.where(be > 105, 3000.0, 300.0)
        w = DP.auto_envelope(be, y)
        self.assertLess(w.hi, 104.0)
        self.assertGreater(w.lo, 90.0)


class RobustShirley(unittest.TestCase):
    def test_it_equals_the_classic_iteration_wherever_that_converges(self):
        be = np.arange(90.0, 110.0, 0.1)
        peak = 1000.0 * sim.pseudo_voigt(be, 100.0, 0.9, 0.0)
        y = peak + sim.shirley_background(peak, step=120.0, base=30.0)
        bg, ok = AR.shirley_info(y, 0.1)
        classic, done = AR._shirley_iterate(y, 0.1, float(y[:3].mean()), float(y[-3:].mean()), 200, 1e-8)
        self.assertTrue(ok and done)
        self.assertTrue(np.allclose(bg, classic))

    def test_when_the_iteration_diverges_the_result_is_still_sensible(self):
        be = np.arange(90.0, 115.0, 0.1)
        y = 50000 * sim.pseudo_voigt(be, 99.0, 1.0, 0.0) + np.where(be > 100, 9000.0, 300.0)
        bg, ok = AR.shirley_info(y, 0.1)
        self.assertTrue(np.all(np.isfinite(bg)))
        self.assertGreaterEqual(float(np.sum(y - bg)), 0.2 * float(np.sum(y - np.linspace(y[:3].mean(), y[-3:].mean(), len(y)))))


if __name__ == "__main__":
    unittest.main()
