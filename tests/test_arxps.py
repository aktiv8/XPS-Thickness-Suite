"""
ARXPS maths tests: background and areas, alignment, auto-detected windows, the fits, and that a thickness put into
synthetic data is recovered. Synthetic data only (see synthetic.py), so these always run.

    python -m unittest discover -s tests -v
"""
import math
import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

import app as A  # noqa: E402
import arxps as AR  # noqa: E402
import spectrum_sim as sim  # noqa: E402
import synthetic as S  # noqa: E402
import thickogram as TG  # noqa: E402

LAM, RINF = 3.89, 0.5547


def _series(blocks, region=None):
    return [(b.exp_values[0], b) for b in blocks if region is None or b.block_id == region]


def _auto(ser, sub=99.3, ov=103.4):
    return AR.auto_two_state_windows([(b.be_ev, AR.intensity(b)) for _, b in ser], sub, ov)


class Background(unittest.TestCase):
    def test_shirley_recovers_a_peak_on_a_known_shirley_background(self):
        be = np.arange(90.0, 110.0, 0.1)
        peak = 1000.0 * sim.pseudo_voigt(be, 100.0, 0.9, 0.0)
        y = peak + sim.shirley_background(peak, step=120.0, base=30.0)
        ex = AR.single_peak_area(be, y, AR.Window(94.0, 106.0), align=False)
        self.assertAlmostEqual(ex.area_sub, 1000.0, delta=15.0)      # 1.5%: the window cuts a sliver of the tails

    def test_shirley_is_flat_for_a_flat_spectrum(self):
        bg = AR.shirley(np.full(50, 7.0), 0.1)
        self.assertTrue(np.allclose(bg, 7.0))

    def test_intensity_per_second_makes_different_scan_counts_comparable(self):
        b = S.two_state_series(noise=False)[0]
        longer = S.two_state_series(noise=False)[0]
        longer.counts, longer.scans = b.counts * 1.5, 6
        self.assertTrue(np.allclose(AR.intensity(b), AR.intensity(longer)))
        self.assertFalse(np.allclose(AR.intensity(b, per_second=False), AR.intensity(longer, per_second=False)))


class Alignment(unittest.TestCase):
    def test_shift_of_a_shoulder_with_a_taller_neighbour_is_found(self):
        be = np.arange(10.0, 30.0, 0.1)

        def spec(shift):
            return (800 * sim.pseudo_voigt(be, 20.0 + shift, 1.2, 0.0)         # tall overlayer-like peak
                    + 300 * sim.pseudo_voigt(be, 18.7 + shift, 0.9, 0.0))      # the substrate is only a shoulder
        for shift in (-0.6, 0.0, 0.37, 0.8):
            got = AR.xcorr_shift(be, spec(0.0), be, spec(shift), 14.0, 25.0)
            self.assertAlmostEqual(got, shift, delta=0.02)

    def test_windows_follow_a_charging_drift_and_report_the_shifts(self):
        ser = _series(S.two_state_series(drift_ev=0.3))
        res = AR.analyse_two_states(ser, _auto(ser), LAM, RINF)
        got = [r.shift_ev for r in res.rows]
        want = [0.3 * t / 75.0 for t in S.DEFAULT_TILTS]
        for g, w in zip(got, want):
            self.assertAlmostEqual(g, w, delta=0.02)
        # the drift must not leak into the thickness
        self.assertAlmostEqual(res.d_nm, 1.5, delta=0.03)


class AutoWindows(unittest.TestCase):
    def test_finds_substrate_and_overlayer_states_and_a_split_between_them(self):
        w = _auto(_series(S.two_state_series()))
        self.assertAlmostEqual(w.sub_pos, 99.0, delta=0.15)
        self.assertAlmostEqual(w.ov_pos, 103.0, delta=0.2)
        self.assertTrue(w.sub_pos < w.split < w.ov_pos)
        self.assertLess(w.env_lo, w.sub_pos)
        self.assertGreater(w.env_hi, w.ov_pos)

    def test_a_thick_overlayer_that_is_the_taller_peak_is_not_mistaken_for_the_substrate(self):
        ser = _series(S.two_state_series(d_nm=4.0))
        w = _auto(ser)
        self.assertLess(w.sub_pos, w.ov_pos)
        self.assertAlmostEqual(w.sub_pos, 99.0, delta=0.2)
        res = AR.analyse_two_states(ser, w, LAM, RINF)
        self.assertAlmostEqual(res.d_nm, 4.0, delta=0.08)

    def test_needs_two_angles(self):
        b = S.two_state_series(tilts=(0.0,))[0]
        with self.assertRaises(ValueError):
            AR.auto_two_state_windows([(b.be_ev, AR.intensity(b))], 99.3, 103.4)


class Recovery(unittest.TestCase):
    def test_two_states_recovers_the_thickness_from_noise_free_data(self):
        for d in (0.8, 1.5, 2.6):
            ser = _series(S.two_state_series(d_nm=d))
            res = AR.analyse_two_states(ser, _auto(ser), LAM, RINF, max_angle=80)
            self.assertAlmostEqual(res.d_nm, d, delta=0.02 * d + 0.01, msg=f"d = {d}")
            self.assertLess(res.d_se_nm, 0.01)
            self.assertAlmostEqual(res.line[0], 0.0, delta=0.01)              # passes through the origin
            self.assertEqual(res.warnings, [])

    def test_two_states_with_noise_and_lorentzian_tails_is_close(self):
        ser = _series(S.two_state_series(d_nm=1.5, noise=True, lorentz=0.2, seed=3))
        res = AR.analyse_two_states(ser, _auto(ser), LAM, RINF)
        self.assertAlmostEqual(res.d_nm, 1.5, delta=0.15)

    def test_per_angle_thickness_equals_the_calculator_strohmeier_function(self):
        n_sub, n_ov, lam_sub, lam_ov = 0.0830, 0.0366, 3.09, 3.89
        rinf = n_ov * lam_ov / (n_sub * lam_sub)
        ser = _series(S.two_state_series(d_nm=1.5, lam_nm=lam_ov, rinf=rinf))
        res = AR.analyse_two_states(ser, _auto(ser), lam_ov, rinf, max_angle=80)
        for r in res.rows:
            want = A.strohmeier_thickness(90.0 - r.angle_deg, n_sub, lam_sub, 1.0, r.i_sub,
                                          [A.OverlayerInput("ox", n_ov, lam_ov, 1.0, r.i_ov)])["d_nm"]
            self.assertAlmostEqual(r.d_nm, want, places=9, msg=f"{r.tilt} deg")

    def test_two_peaks_recovers_the_thickness_and_matches_the_thickogram_per_angle(self):
        blocks = S.two_peak_series(d_nm=1.5)
        sub, ov = _series(blocks, "Si2p Scan"), _series(blocks, "O1s Scan")
        sw = AR.auto_peak_window([(b.be_ev, AR.intensity(b)) for _, b in sub], 99.3)
        ow = AR.auto_peak_window([(b.be_ev, AR.intensity(b)) for _, b in ov], 532.0)
        res = AR.analyse_two_peaks(sub, ov, sw, ow, 2.93, 0.82, 532.0, 99.3, 1486.68, 3.0, max_angle=80)
        self.assertAlmostEqual(res.d_nm, 1.5, delta=0.015)
        self.assertAlmostEqual(res.b, ((1486.68 - 532.0) / (1486.68 - 99.3)), places=9)
        for r in res.rows:
            want = TG.thickness(io=r.i_ov, is_=r.i_sub, so=2.93, ss=0.82, be_o=532.0, be_s=99.3, source_ev=1486.68,
                                lam_nm=3.0, theta_normal_deg=r.angle_deg).t_nm
            self.assertAlmostEqual(r.d_nm, want, places=9)

    def test_angle_offset_shifts_the_emission_angle(self):
        ser = _series(S.two_state_series(tilts=(0.0, 20.0, 40.0)))
        res = AR.analyse_two_states(ser, _auto(ser), LAM, RINF, angle_offset=10.0)
        self.assertEqual([r.angle_deg for r in res.rows], [10.0, 30.0, 50.0])


class Fits(unittest.TestCase):
    def test_through_origin_slope_and_error(self):
        u = [1.0, 1.4, 2.0, 3.0]
        k, se = AR.fit_through_origin(u, [0.4 * v for v in u])
        self.assertAlmostEqual(k, 0.4, places=12)
        self.assertAlmostEqual(se, 0.0, places=12)
        self.assertIsNone(AR.fit_through_origin([1.0], [0.4])[1])

    def test_free_intercept_line(self):
        u = np.array([1.0, 1.4, 2.0, 3.0])
        c, k, se_c, se_k = AR.fit_line(u, 0.2 + 0.3 * u)
        self.assertAlmostEqual(c, 0.2, places=10)
        self.assertAlmostEqual(k, 0.3, places=10)
        self.assertIsNone(AR.fit_line([1.0, 2.0], [0.1, 0.2]))

    def test_rinf_free_fit_recovers_both_parameters(self):
        u = np.array([1.0, 1.155, 1.414, 2.0, 2.9])
        r = 0.55 * np.expm1(0.42 * u)
        fit = AR.fit_rinf_free(u, r)
        self.assertAlmostEqual(fit["k"], 0.42, delta=0.002)
        self.assertAlmostEqual(fit["rinf"], 0.55, delta=0.005)
        self.assertFalse(fit["edge"])

    def test_rinf_free_fit_says_when_the_ratio_is_just_linear_in_u(self):
        u = np.array([1.0, 1.155, 1.414, 2.0])
        self.assertTrue(AR.fit_rinf_free(u, 0.3 * u)["edge"])

    def test_rinf_free_fit_says_when_the_exponential_is_nearly_linear_over_the_angles(self):
        u = np.array([1.0, 1.155, 1.414, 2.0, 2.9])
        fit = AR.fit_rinf_free(u, 4.0 * np.expm1(0.08 * u))              # k u_max = 0.23: R-inf and k are interchangeable
        self.assertTrue(fit["edge"])

    def _rows(self, xs, tilts, rinf=0.55):
        return [AR.AngleRow(t, t, 100.0, 100.0 * rinf * math.expm1(x), rinf * math.expm1(x))
                for t, x in zip(tilts, xs)]

    def test_a_series_that_misses_the_origin_raises_the_intercept_warning(self):
        tilts = [0, 20, 40, 55, 65]
        u = [1 / math.cos(math.radians(t)) for t in tilts]
        rows = self._rows([0.25 + 0.2 * v for v in u], tilts)
        res = AR.fit_series(rows, 3.0, lambda r: math.log1p(r / 0.55), 70.0, 1, rinf=0.55)
        self.assertTrue(any("free-intercept" in w for w in res.warnings))

    def test_a_clean_series_has_no_warnings(self):
        tilts = [0, 20, 40, 55, 65]
        u = [1 / math.cos(math.radians(t)) for t in tilts]
        res = AR.fit_series(self._rows([0.3 * v for v in u], tilts), 3.0, lambda r: math.log1p(r / 0.55), 70.0, 1,
                            rinf=0.55)
        self.assertEqual(res.warnings, [])
        self.assertAlmostEqual(res.d_nm, 0.9, places=9)

    def test_angles_above_the_maximum_are_shown_but_not_fitted(self):
        tilts = [0, 30, 60, 75]
        u = [1 / math.cos(math.radians(t)) for t in tilts]
        xs = [0.3 * v for v in u]
        xs[-1] *= 0.5                                      # a bad grazing point
        res = AR.fit_series(self._rows(xs, tilts), 3.0, lambda r: math.log1p(r / 0.55), 70.0, 1, rinf=0.55)
        self.assertEqual([r.used for r in res.rows], [True, True, True, False])
        self.assertIsNotNone(res.rows[-1].x)
        self.assertAlmostEqual(res.d_nm, 0.9, places=9)
        self.assertTrue(any("75 deg" in w for w in res.warnings))

    def test_zero_area_angles_are_skipped_with_a_warning_and_no_usable_angle_raises(self):
        rows = [AR.AngleRow(0.0, 0.0, 100.0, 0.0, 0.0), AR.AngleRow(30.0, 30.0, 100.0, 20.0, 0.2),
                AR.AngleRow(50.0, 50.0, 100.0, 30.0, 0.3)]
        res = AR.fit_series(rows, 3.0, lambda r: math.log1p(r / 0.55), 70.0, 1, rinf=0.55)
        self.assertEqual([r.used for r in res.rows], [False, True, True])
        self.assertTrue(any("0 deg skipped" in w for w in res.warnings))
        with self.assertRaises(ValueError):
            AR.fit_series([AR.AngleRow(0.0, 0.0, 100.0, 0.0, 0.0)], 3.0, lambda r: r, 70.0, 1)

    def test_few_angles_warn(self):
        tilts = [0, 40]
        u = [1 / math.cos(math.radians(t)) for t in tilts]
        res = AR.fit_series(self._rows([0.3 * v for v in u], tilts), 3.0, lambda r: math.log1p(r / 0.55), 70.0, 1,
                            rinf=0.55)
        self.assertTrue(any("3 or more" in w for w in res.warnings))

    def test_unresolved_states_raise_the_overlap_warning(self):
        ser = _series(S.two_state_series(sub_be=98.9, ov_be=100.3, fwhm_sub=0.9, fwhm_ov=1.1, drift_ev=0.0))
        res = AR.analyse_two_states(ser, AR.TwoStateWindows(96.4, 103.0, 99.6), LAM, RINF)
        self.assertTrue(any("less than 2 eV apart" in w for w in res.warnings))

    def test_resolved_states_do_not(self):
        ser = _series(S.two_state_series())
        res = AR.analyse_two_states(ser, _auto(ser), LAM, RINF)
        self.assertFalse(any("apart" in w for w in res.warnings))


class FittedAreas(unittest.TestCase):
    """Substrate and overlayer areas from CasaXPS-style fits stored in the VMS (synthetic, see synthetic.py)."""

    def _roles(self, ser):
        return AR.auto_roles(AR.read_fits(ser))

    def test_the_trend_rule_finds_the_substrate_and_the_overlayer(self):
        ser = _series(S.two_state_series(fitted=True))
        fits = AR.read_fits(ser)
        self.assertEqual(list(AR.state_names(fits).values()), ["Si El", "Si Ox"])
        self.assertEqual(AR.auto_roles(fits), {"i0": AR.SUBSTRATE, "i1": AR.OVERLAYER})
        shares = AR.state_shares(fits)
        self.assertGreater(shares["i0"][0], shares["i0"][1])             # falls with angle
        self.assertLess(shares["i1"][0], shares["i1"][1])                # rises

    def test_a_thick_overlayer_that_dominates_is_still_the_overlayer(self):
        ser = _series(S.two_state_series(d_nm=4.0, fitted=True))
        self.assertEqual(self._roles(ser), {"i0": AR.SUBSTRATE, "i1": AR.OVERLAYER})

    def test_a_state_whose_share_hardly_moves_is_not_assigned(self):
        fits = AR.read_fits(_series(S.two_state_series(fitted=True)))
        for f in fits:                                                   # make both states identical at every angle
            for st in f.states.values():
                for c in st.comps:
                    c.area = 100.0
        self.assertEqual(set(AR.auto_roles(fits).values()), {AR.IGNORE})

    def test_use_rule_for_two_peak_mode(self):
        fits = AR.read_fits(_series(S.two_state_series(fitted=True)))
        self.assertEqual(AR.auto_use(fits, AR.SUBSTRATE), {"i0": AR.SUBSTRATE, "i1": AR.IGNORE})
        self.assertEqual(AR.auto_use(fits, AR.OVERLAYER), {"i0": AR.IGNORE, "i1": AR.OVERLAYER})
        one = AR.read_fits(_series(S.two_peak_series(fitted=True), "O1s Scan"))
        self.assertEqual(AR.auto_use(one, AR.OVERLAYER), {"i0": AR.OVERLAYER})              # a lone state is always used

    def test_two_states_recovers_the_thickness_exactly_and_matches_the_window_method(self):
        blocks = S.two_state_series(d_nm=1.5, fitted=True)
        ser = _series(blocks)
        res = AR.analyse_fitted_states(ser, self._roles(ser), LAM, RINF, max_angle=80, transmission=False)
        self.assertAlmostEqual(res.d_nm, 1.5, places=5)
        for r in res.rows:
            self.assertAlmostEqual(r.d_nm, 1.5, places=5)
        self.assertEqual(res.warnings, [])
        windows = AR.analyse_two_states(ser, _auto(ser), LAM, RINF, max_angle=80)
        self.assertAlmostEqual(res.d_nm, windows.d_nm, delta=0.03)

    def test_the_energy_drift_is_read_from_the_state_positions(self):
        ser = _series(S.two_state_series(fitted=True, drift_ev=0.3))
        res = AR.analyse_fitted_states(ser, self._roles(ser), LAM, RINF)
        for r, tilt in zip(res.rows, S.DEFAULT_TILTS):
            self.assertAlmostEqual(r.shift_ev, 0.3 * tilt / 75.0, places=6)

    def test_the_charge_correction_changes_nothing_about_the_thickness(self):
        a = _series(S.two_state_series(fitted=True))
        b = _series(S.two_state_series(fitted=True, calib=0.7))
        ra = AR.analyse_fitted_states(a, self._roles(a), LAM, RINF)
        rb = AR.analyse_fitted_states(b, self._roles(b), LAM, RINF)
        self.assertAlmostEqual(ra.d_nm, rb.d_nm, places=9)

    def test_assigning_both_doublet_partners_is_what_makes_the_area(self):
        ser = _series(S.two_state_series(fitted=True))
        fits = AR.read_fits(ser)
        total = fits[0].states["i0"].area
        self.assertAlmostEqual(total, 60000.0 * math.exp(-1.5 / LAM), delta=1e-3 * total)     # 2p3/2 + 2p1/2

    def test_both_roles_are_needed_and_a_block_without_a_fit_is_refused(self):
        ser = _series(S.two_state_series(fitted=True))
        with self.assertRaisesRegex(ValueError, "at least one overlayer state"):
            AR.analyse_fitted_states(ser, {"i0": AR.SUBSTRATE, "i1": AR.IGNORE}, LAM, RINF)
        with self.assertRaisesRegex(ValueError, "at least one substrate state"):
            AR.analyse_fitted_states(ser, {"i0": AR.IGNORE, "i1": AR.OVERLAYER}, LAM, RINF)
        with self.assertRaisesRegex(ValueError, "no CasaXPS fit"):
            AR.analyse_fitted_states(_series(S.two_state_series()), {"i0": AR.SUBSTRATE}, LAM, RINF)

    def test_several_fitted_azimuths(self):
        from vamas import VamasFile
        blocks = S.two_state_series(azimuths=(0.0, 90.0), d_by_azimuth={0.0: 1.2, 90.0: 1.8}, fitted=True)
        series = VamasFile("", blocks).series("Si2p Scan", "exp1", "exp2")
        roles = self._roles(series[0.0])
        got = {az: AR.analyse_fitted_states(rows, roles, LAM, RINF, max_angle=80).d_nm for az, rows in series.items()}
        self.assertAlmostEqual(got[0.0], 1.2, places=4)
        self.assertAlmostEqual(got[90.0], 1.8, places=4)

    def test_two_peaks_recovers_the_thickness_and_matches_the_thickogram(self):
        blocks = S.two_peak_series(d_nm=1.5, fitted=True)
        sub, ov = _series(blocks, "Si2p Scan"), _series(blocks, "O1s Scan")
        fs, fo = AR.read_fits(sub), AR.read_fits(ov)
        res = AR.analyse_fitted_peaks(sub, ov, AR.auto_use(fs, AR.SUBSTRATE), AR.auto_use(fo, AR.OVERLAYER), 2.93, 0.82,
                                      532.0, 99.3, 1486.68, 3.0, max_angle=80)
        self.assertAlmostEqual(res.d_nm, 1.5, places=5)
        for r in res.rows:
            want = TG.thickness(io=r.i_ov, is_=r.i_sub, so=2.93, ss=0.82, be_o=532.0, be_s=99.3, source_ev=1486.68,
                                lam_nm=3.0, theta_normal_deg=r.angle_deg).t_nm
            self.assertAlmostEqual(r.d_nm, want, places=9)

    def test_the_transmission_function_is_applied_per_component_in_two_peak_mode(self):
        import re
        blocks = S.two_peak_series(d_nm=1.5, fitted=True)
        for b in blocks:
            if b.block_id == "O1s Scan":                      # a detector twice as sensitive at the O 1s energy
                b.transmission = b.transmission * 2.0
                b.counts = b.counts * 2.0
                b.comments = [re.sub(r"Area (\S+)", lambda m: f"Area {float(m.group(1)) * 2.0:.8g}", c) for c in b.comments]
        sub, ov = _series(blocks, "Si2p Scan"), _series(blocks, "O1s Scan")
        fs, fo = AR.read_fits(sub), AR.read_fits(ov)
        args = (AR.auto_use(fs, AR.SUBSTRATE), AR.auto_use(fo, AR.OVERLAYER), 2.93, 0.82, 532.0, 99.3, 1486.68, 3.0)
        corrected = AR.analyse_fitted_peaks(sub, ov, *args, max_angle=80, transmission=True)
        raw = AR.analyse_fitted_peaks(sub, ov, *args, max_angle=80, transmission=False)
        self.assertAlmostEqual(corrected.d_nm, 1.5, places=5)
        self.assertGreater(abs(raw.d_nm - 1.5), 0.1)

    def test_the_fit_region_can_be_chosen_when_a_block_holds_several(self):
        import synthetic
        b = S.two_state_series(fitted=True)
        extra = synthetic.casa_fit_lines("As 3d", 1440.0, 1450.0, 0.5,
                                         [dict(name="As", area=123.0, ke=1445.0, fwhm=1.0, index=-1, group="As", rsf=0.5)])
        for blk in b:
            blk.comments += extra
        ser = _series(b)
        self.assertEqual(AR.fit_regions(ser), ["Si 2p", "As 3d"])
        self.assertEqual(list(AR.state_names(AR.read_fits(ser, "As 3d")).values()), ["As"])
        self.assertEqual(list(AR.state_names(AR.read_fits(ser, "Si 2p")).values()), ["Si El", "Si Ox"])


class Azimuths(unittest.TestCase):
    def test_each_azimuth_is_fitted_separately_and_summarised(self):
        blocks = S.two_state_series(azimuths=(0.0, 90.0), d_by_azimuth={0.0: 1.2, 90.0: 1.8})
        from vamas import VamasFile
        f = VamasFile("", blocks)
        series = f.series("Si2p Scan", "exp1", "exp2")
        w = _auto(series[0.0])
        results = {az: AR.analyse_two_states(rows, w, LAM, RINF, max_angle=80) for az, rows in series.items()}
        self.assertAlmostEqual(results[0.0].d_nm, 1.2, delta=0.03)
        self.assertAlmostEqual(results[90.0].d_nm, 1.8, delta=0.04)
        mean, spread = AR.summarise_azimuths(results)
        self.assertAlmostEqual(mean, 1.5, delta=0.04)
        self.assertAlmostEqual(spread, 0.6 / math.sqrt(2), delta=0.05)

    def test_single_azimuth_has_no_spread(self):
        ser = _series(S.two_state_series())
        res = {None: AR.analyse_two_states(ser, _auto(ser), LAM, RINF)}
        self.assertIsNone(AR.summarise_azimuths(res)[1])


class Coverage(unittest.TestCase):
    """Partial coverage f (Sodhi 2017, Eq. 2) in the two-states mode."""
    RINF, LAM, D = 0.6, 3.0, 2.0
    TILTS = [0.0, 30.0, 45.0, 60.0, 70.0]

    def _rows(self, f, d=None):
        rows = []
        for t in self.TILTS:
            x = (d or self.D) / (self.LAM * math.cos(math.radians(t)))
            e = math.exp(-x)
            r = self.RINF * f * (1 - e) / ((1 - f) + f * e)
            rows.append(AR.AngleRow(t, t, 100.0, 100.0 * r, r))
        return rows

    def test_f_equal_to_one_is_the_ordinary_logarithm(self):
        for r in (0.05, 0.4, 3.0):
            self.assertEqual(AR.coverage_x(r, 0.55, 1.0), math.log1p(r / 0.55))

    def test_the_forward_model_and_its_inverse_agree(self):
        for f in (0.3, 0.7, 0.95):
            for x in (0.1, 0.8, 2.5):
                e = math.exp(-x)
                r = 0.55 * f * (1 - e) / ((1 - f) + f * e)
                self.assertAlmostEqual(AR.coverage_x(r, 0.55, f), x, places=10)

    def test_a_ratio_above_the_limit_for_that_coverage_is_refused(self):
        limit = 0.55 * 0.5 / 0.5
        with self.assertRaises(ValueError):
            AR.coverage_x(limit * 1.01, 0.55, 0.5)
        with self.assertRaises(ValueError):
            AR.coverage_x(0.3, 0.55, 0.0)

    def test_a_given_coverage_recovers_the_thickness_a_complete_layer_would_underestimate(self):
        rows = self._rows(0.7)
        solve, f, info = AR._coverage_solver(rows, self.RINF, 0.7, 70.0)
        res = AR.fit_series(rows, self.LAM, solve, 70.0, 1, rinf=self.RINF, coverage=f)
        self.assertAlmostEqual(res.d_nm, self.D, places=9)
        self.assertIsNone(res.rinf_free)                      # that fit assumes a complete layer
        full = AR.fit_series(self._rows(0.7), self.LAM, lambda r: math.log1p(r / self.RINF), 70.0, 1, rinf=self.RINF)
        self.assertLess(full.d_nm, 0.7 * self.D)

    def test_fitting_the_coverage_finds_it_and_the_thickness(self):
        for f in (0.5, 0.8):
            rows = self._rows(f)
            solve, fu, info = AR._coverage_solver(rows, self.RINF, None, 70.0)
            self.assertAlmostEqual(fu, f, delta=0.006)
            self.assertFalse(info["edge"])
            self.assertLess(info["cv"], 0.01)
            self.assertGreater(info["cv1"], 0.05)
            res = AR.fit_series(rows, self.LAM, solve, 70.0, 1, rinf=self.RINF, coverage=fu)
            self.assertAlmostEqual(res.d_nm, self.D, delta=0.02)

    def test_a_complete_layer_is_not_mistaken_for_patchy_and_the_fit_says_it_is_undetermined(self):
        rows = self._rows(1.0)
        solve, fu, info = AR._coverage_solver(rows, self.RINF, None, 70.0)
        self.assertEqual(fu, 1.0)
        self.assertTrue(info["edge"])
        res = AR.fit_series(rows, self.LAM, solve, 70.0, 1, rinf=self.RINF, coverage=fu)
        AR._note_coverage(res, info)
        self.assertTrue(any("not well determined" in w for w in res.warnings))

    def test_fitting_needs_three_angles(self):
        with self.assertRaises(ValueError):
            AR._coverage_solver(self._rows(0.7)[:2], self.RINF, None, 70.0)

    def test_the_window_method_accepts_a_coverage_and_the_default_is_unchanged(self):
        ser = _series(S.two_state_series(d_nm=1.5))
        base = AR.analyse_two_states(ser, _auto(ser), LAM, RINF, max_angle=80)
        same = AR.analyse_two_states(ser, _auto(ser), LAM, RINF, max_angle=80, coverage=1.0)
        self.assertEqual(base.d_nm, same.d_nm)
        self.assertEqual(base.coverage, 1.0)
        patchy = AR.analyse_two_states(ser, _auto(ser), LAM, RINF, max_angle=80, coverage=0.9)
        self.assertGreater(patchy.d_nm, base.d_nm)             # less cover for the same ratio means thicker patches


if __name__ == "__main__":
    unittest.main()
