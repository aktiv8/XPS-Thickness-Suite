"""
CasaXPS fit parsing tests (casa.py, ported from the author's eXPoSe project). The grammar lines are copied from real
CasaXPS files; the test against the user's own fitted tilt series runs only when that file exists on this machine.

    python -m unittest discover -s tests -v
"""
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

import arxps as AR  # noqa: E402
import casa  # noqa: E402
import synthetic as S  # noqa: E402
import vamas as VM  # noqa: E402

REAL = (r"C:\Users\david\Documents\Working Directory\20250626\Dave Morgan\X-Ray014 400um - FG  ON"
        r"\Stage Tilt ARXPS Si Wafer\Fitted Si tilt.vms")

SI_REGION = "CASA region (*Si 2p*) (*Shirley*) 1379.1765 1392.2142 0.817 2 0 0 0 -450 0 0 (*Si 2p*) 28.0858 0 0.817"
SI_EL = ("CASA comp (*Si 2p El*) (*LA(50)*) Area 73449.445 1e-20 748475.55 -1 1 MFWHM 1.0834043 0.24 6 -1 1 "
         "Position 1387.8626 1379.1765 1392.2142 -1 1 RSF 0.817 MASS 28.0858 INDEX 0 (*Si 2p*) CONST (**) UNCORR")
SI_OX = ("CASA comp (*Si 2p Ox*) (*LA(50)*) Area 24509.88 1e-20 364219.69 -1 1 MFWHM 1.7027027 0.19920319 4.9800797 -1 1 "
         "Position 1384.063 1379.1765 1392.2142 -1 1 RSF 0.817 MASS 28.0858 INDEX 1 (*Si 2p*) CONST")
HV = 1486.68


def _block(lines):
    class B:
        comments = lines
    return B()


class Grammar(unittest.TestCase):
    def test_the_real_casaxps_lines(self):
        fit = casa.parse(["Casa Info Follows", "0", "0", "1", SI_REGION, "2", SI_EL, SI_OX, "X-ray spot-size: 400 um"])
        self.assertEqual([r.name for r in fit.regions], ["Si 2p"])
        reg = fit.regions[0]
        self.assertEqual((reg.background, reg.start_ke, reg.end_ke, reg.rsf), ("Shirley", 1379.1765, 1392.2142, 0.817))
        el, ox = fit.components
        self.assertEqual((el.name, el.shape, el.area, el.fwhm, el.pos_ke, el.rsf, el.index, el.group, el.region),
                         ("Si 2p El", "LA(50)", 73449.445, 1.0834043, 1387.8626, 0.817, 0, "Si 2p", "Si 2p"))
        self.assertEqual((ox.name, ox.area, ox.pos_ke, ox.index), ("Si 2p Ox", 24509.88, 1384.063, 1))

    def test_no_fit_gives_none(self):
        self.assertIsNone(casa.parse(["Casa Info Follows", "0", "0", "1", "X-ray spot-size: 400 um by 400 um"]))
        self.assertIsNone(casa.parse([]))
        self.assertIsNone(casa.parse(None))

    def test_library_line_shapes_constraints_and_the_uncorrected_rsf_tail(self):
        line = ("CASA comp (*Ru 3p Metal*) (*%&<>LS-RUMET3P*) Area 1345.4905 1e-020 10000000 -1 1 MFWHM 1 1 1 -1 1 "
                "Position 1025.3055 981.53 1046.53 -1 1 RSF 2.04 MASS 101.069 INDEX 2 (*Ru 3p*) CONST (**) "
                "UNCORRECTEDRSF 2.04")
        c = casa.parse(["CASA region (*Ru 3p*) (*Shirley*) 1008.4492 1036.0906 2.04 2 0 0 29493.218 -3650 0 0 "
                        "(*Ru 3p*) 101.069 0 2.04", "1", line]).components[0]
        self.assertEqual((c.name, c.shape, c.area, c.pos_ke, c.rsf, c.index), ("Ru 3p Metal", "%&<>LS-RUMET3P",
                                                                                  1345.4905, 1025.3055, 2.04, 2))
        linked = ("CASA comp (*Ru 3d met*) (*LF(1.01,1.25,1000,50)*) Area 1202.2565 0 10000000 4 0.6666667 MFWHM "
                  "1.0218928 0.2 5 4 1.6385406 Position 1202.4337 0 0 4 -4.2175 RSF 4.27 MASS 101.069 INDEX 5 "
                  "(*Ru 3d*) CONST (**) UNCORRECTEDRSF 4.27")
        c = casa.parse([SI_REGION, "1", linked]).components[0]
        self.assertEqual((c.area, c.pos_ke, c.index), (1202.2565, 1202.4337, 5))

    def test_a_comp_line_without_area_or_position_is_ignored(self):
        self.assertEqual(len(casa.parse([SI_REGION, "2", "CASA comp (*x*) (*GL(30)*) MFWHM 1", SI_EL]).components), 1)

    def test_region_attribution_when_regions_are_listed_before_their_comps(self):
        s2p = "CASA region (*S 2p*) (*Shirley*) 1310 1330 1.67 1 0 0 0 0 0 0 (*S 2p*) 32.06"
        p2s = "CASA region (*P 2s*) (*Shirley*) 1240 1260 1.2 1 0 0 0 0 0 0 (*P 2s*) 30.97"
        in_s = "CASA comp (*S a*) (*GL(30)*) Area 10 0 1e7 -1 1 MFWHM 1 0 5 -1 1 Position 1320 0 0 -1 1 RSF 1.67 MASS 32 INDEX -1"
        in_p = "CASA comp (*P a*) (*GL(30)*) Area 20 0 1e7 -1 1 MFWHM 1 0 5 -1 1 Position 1250 0 0 -1 1 RSF 1.2 MASS 31 INDEX -1"
        fit = casa.parse([s2p, p2s, "2", in_s, in_p])           # both comps follow the LAST region line
        self.assertEqual({c.name: c.region for c in fit.components}, {"S a": "S 2p", "P a": "P 2s"})
        self.assertEqual([c.name for c in fit.region_components("S 2p")], ["S a"])


class States(unittest.TestCase):
    def test_a_state_in_a_group_of_its_own_is_named_by_its_component(self):
        fit = casa.parse([SI_REGION, "2", SI_EL, SI_OX])                 # INDEX 0 and 1, both tagged (*Si 2p*)
        st = casa.states_of(fit, HV)
        self.assertEqual([s.name for s in st.values()], ["Si 2p El", "Si 2p Ox"])     # the tag is only the region label
        self.assertAlmostEqual(st["i0"].area, 73449.445)
        self.assertAlmostEqual(st["i0"].be, HV - 1387.8626)

    def test_doublet_partners_sharing_an_index_add_up_and_use_the_group_tag(self):
        a = ("CASA comp (*Ti 2p3/2*) (*GL(30)*) Area 300 0 1e7 -1 1 MFWHM 1 0 5 -1 1 Position 1028 0 0 -1 1 RSF 2 MASS 47.9 "
             "INDEX 0 (*Metal*) CONST")
        b = ("CASA comp (*Ti 2p1/2*) (*GL(30)*) Area 150 0 1e7 0 0.5 MFWHM 1 0 5 0 1 Position 1022 0 0 0 -6 RSF 0 MASS 47.9 "
             "INDEX 0 (*Metal*) CONST")
        c = ("CASA comp (*TiO2*) (*GL(30)*) Area 90 0 1e7 -1 1 MFWHM 1 0 5 -1 1 Position 1024 0 0 -1 1 RSF 2 MASS 47.9 "
             "INDEX -1")
        st = casa.states_of(casa.parse(["CASA region (*Ti 2p*) (*Shirley*) 1010 1040 2 1 0 0 0 0 0 0 (*Ti 2p*) 47.9", "3",
                                        a, b, c]), HV)
        self.assertEqual([s.name for s in st.values()], ["Metal", "TiO2"])
        self.assertAlmostEqual(st["i0"].area, 450.0)
        self.assertAlmostEqual(st["i0"].be, HV - (300 * 1028 + 150 * 1022) / 450)          # area-weighted
        self.assertEqual(len(st["i0"].comps), 2)

    def test_the_component_rsf_is_never_used_so_a_zero_rsf_partner_still_counts(self):
        sb = casa.states_of(casa.parse(S.two_state_series(fitted=True)[0].comments), HV)
        self.assertEqual([c.name for c in sb["i0"].comps], ["Si El 2p3/2", "Si El 2p1/2"])
        self.assertGreater(sb["i0"].comps[1].area, 0)


class Calibration(unittest.TestCase):
    CALIB = "Calib M = 281.7289 A = 282 BE ADD"

    def test_the_charge_shift_moves_the_displayed_binding_energy_not_the_area(self):
        fit = casa.parse([self.CALIB, SI_REGION, "1", SI_EL])
        self.assertAlmostEqual(fit.calib_shift, 282 - 281.7289)
        c = fit.components[0]
        self.assertAlmostEqual(casa.component_be(c, HV, fit), HV - 1387.8626)         # no flags: stored = displayed frame
        self.assertAlmostEqual(casa.states_of(fit, HV)["i0"].comps[0].ke_raw, 1387.8626 + 282 - 281.7289)

    def test_the_flags_say_the_fit_was_already_moved_with_the_spectrum(self):
        fit = casa.parse(["Calib M = 281.7289 A = 282 Regions Comps BE ADD", SI_REGION, "1", SI_EL])
        self.assertEqual(fit.comp_shift, 0.0)
        self.assertAlmostEqual(casa.component_be(fit.components[0], HV, fit), HV - 1387.8626 + (282 - 281.7289))

    def test_a_calibrated_synthetic_file_reports_the_shifted_binding_energies(self):
        b = S.two_state_series(fitted=True, calib=0.7, drift_ev=0.0)[0]
        st = casa.states_of(casa.fits_of(b), b.hv_ev)
        self.assertAlmostEqual(st["i0"].be, 98.9 + 0.2 + 0.7, places=6)                   # centroid of the doublet + shift


class Files(unittest.TestCase):
    def test_fits_survive_writing_and_reading_a_vamas_file(self):
        fd, path = tempfile.mkstemp(suffix=".vms")
        os.close(fd)
        try:
            VM.write(path, S.two_state_series(fitted=True))
            f = VM.read(path)
        finally:
            os.remove(path)
        series = f.series("Si2p Scan", "exp1")[None]
        fits = AR.read_fits(series)
        self.assertEqual(list(AR.state_names(fits).values()), ["Si El", "Si Ox"])
        self.assertEqual(len(fits), 6)

    @unittest.skipUnless(os.path.isfile(REAL), "the user's fitted CasaXPS tilt series is not on this machine")
    def test_the_real_fitted_si_tilt_series(self):
        f = VM.read(REAL)
        series = f.series("Si2p Scan", "exp1")[None]
        fits = AR.read_fits(series)
        names = AR.state_names(fits)
        self.assertEqual(list(names.values()), ["Si 2p El", "Si 2p Ox"])
        el, ox = (fits[0].states[k] for k in names)
        self.assertAlmostEqual(el.area, 73449.45, places=1)                    # Si.csv, cycle 0
        self.assertAlmostEqual(ox.area, 24509.88, places=1)
        self.assertAlmostEqual(el.be, 98.82, places=2)
        self.assertAlmostEqual(ox.be, 102.62, places=2)
        self.assertAlmostEqual(fits[-1].states[list(names)[0]].area, 8143.61, places=1)       # cycle 75
        self.assertEqual(AR.fit_regions(f.series("C1s Scan", "exp1")[None]), [])             # only Si 2p was fitted
        self.assertEqual(AR.auto_roles(fits), {"i0": AR.SUBSTRATE, "i1": AR.OVERLAYER})
        res = AR.analyse_fitted_states(series, AR.auto_roles(fits), 3.89, 0.5547, fits=fits)
        self.assertAlmostEqual(res.d_nm, 1.72, delta=0.01)


if __name__ == "__main__":
    unittest.main()
