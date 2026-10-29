"""
Batch thickness tests (batch.py, batch_tab.py): synthetic fitted VAMAS files, plus the user's own Al oxide file when it
exists on this machine.

    python -m unittest discover -s tests -v
"""
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

import numpy as np  # noqa: E402

import app as A  # noqa: E402
import batch as B  # noqa: E402
import synthetic as S  # noqa: E402
import vamas as VM  # noqa: E402

REAL = r"D:\Temp\AlOxide\bae al samples october 2026 (Calibrated to Al Metal).vms"
HV = 1486.68


def _block(sample, region, comps):
    be = np.linspace(80.0, 68.0, 61)[::-1]
    lines = S.casa_fit_lines(region, HV - 80.0, HV - 68.0, 0.537, [
        dict(name=n, area=a, ke=HV - e, fwhm=1.0, index=-1, group=region, rsf=0.537) for n, a, e in comps])
    return VM.Block(block_id=region, sample_id=sample, species="Al", transition="2p", hv_ev=HV, be_ev=np.sort(be),
                    counts=np.ones(61), dwell_s=0.05, scans=6, comments=lines)


SAMPLES = [("A", 10000.0, 45000.0), ("B", 9000.0, 20000.0), ("C", 8000.0, 8000.0)]


def _blocks():
    out = [_block(n, "Al2p Scan", [("Al 2p Metal", m, 72.5), ("Al 2p Oxide", o, 75.0)]) for n, m, o in SAMPLES]
    out.append(VM.Block(block_id="Al2p Scan", sample_id="no fit", hv_ev=HV, be_ev=np.linspace(68, 80, 61),
                        counts=np.ones(61), dwell_s=0.05, scans=6))
    return out


def _file(testcase):
    d = tempfile.TemporaryDirectory()
    testcase.addCleanup(d.cleanup)
    path = os.path.join(d.name, "fitted.vms")
    VM.write(path, _blocks())
    return path, VM.read(path)


class Maths(unittest.TestCase):
    def setUp(self):
        _, self.vf = _file(self)
        self.samples, self.skipped = B.collect(self.vf, "Al2p Scan")
        self.states = B.fit_states(self.samples)

    def test_collect_keeps_fitted_blocks_and_counts_the_rest(self):
        self.assertEqual([s.sample_id for s in self.samples], ["A", "B", "C"])
        self.assertEqual(self.skipped, 1)
        self.assertEqual(B.fitted_regions(self.vf), ["Al2p Scan"])

    def test_roles_are_guessed_from_binding_energy(self):
        self.assertEqual([s.name for s in self.states], ["Al 2p Metal", "Al 2p Oxide"])
        roles = B.guess_roles(self.states)
        self.assertEqual(sorted(roles.values()), [B.OVERLAYER, B.SUBSTRATE])
        self.assertEqual(roles[self.states[0].key], B.SUBSTRATE)

    def test_thickness_equals_the_direct_strohmeier_call(self):
        roles = B.guess_roles(self.states)
        sub, ov = A.SUBSTRATES["Al (metal)"], A.OVERLAYERS["Al2O3"]
        mats = {self.states[1].key: ov}
        rows = B.compute(self.samples, roles, sub, mats, 90.0, A.volume_density)
        for row, (_, m, o) in zip(rows, SAMPLES):
            ref = A.strohmeier_thickness(
                90.0, A.volume_density(sub.density_g_cm3, sub.molar_mass_g_mol, sub.atoms_per_formula), sub.imfp_nm,
                sub.rsf, m, [A.OverlayerInput("Al2O3", A.volume_density(ov.density_g_cm3, ov.molar_mass_g_mol,
                                                                         ov.atoms_per_formula), ov.imfp_nm, ov.rsf, o)])
            self.assertAlmostEqual(row.d_nm, ref["d_nm"], places=9)
        self.assertGreater(rows[0].d_nm, rows[1].d_nm)       # more oxide per metal -> thicker

    def test_a_sample_missing_a_state_gives_an_error_row(self):
        roles = B.guess_roles(self.states)
        samples = list(self.samples)
        samples.append(B.Sample(9, "broken", {self.states[1].key: ("Al 2p Oxide", 5.0, 75.0)}))
        rows = B.compute(samples, roles, A.SUBSTRATES["Al (metal)"], {self.states[1].key: A.OVERLAYERS["Al2O3"]},
                         90.0, A.volume_density)
        self.assertIsNone(rows[-1].d_nm)
        self.assertIn("no", rows[-1].error)
        self.assertIsNotNone(rows[0].d_nm)

    def test_roles_must_be_complete(self):
        with self.assertRaises(ValueError):
            B.compute(self.samples, {s.key: B.OVERLAYER for s in self.states}, A.SUBSTRATES["Al (metal)"], {}, 90.0,
                      A.volume_density)


@unittest.skipUnless(os.path.exists(REAL), "the user's Al oxide file is not on this machine")
class RealFile(unittest.TestCase):
    def test_al_oxide_set(self):
        vf = VM.read(REAL)
        samples, _ = B.collect(vf, "Al2p Scan")
        self.assertEqual(len(samples), 14)
        first = samples[0]
        areas = sorted(round(v[1], 2) for v in first.states.values())
        self.assertEqual(areas, [10031.64, 45231.28])
        states = B.fit_states(samples)
        rows = B.compute(samples, B.guess_roles(states), A.SUBSTRATES["Al (metal)"],
                         {states[1].key: A.OVERLAYERS["Al2O3"]}, 90.0, A.volume_density)
        self.assertTrue(all(r.d_nm and r.d_nm > 0 for r in rows))


SI_REAL = (r"C:\Users\david\Documents\Working Directory\20250626\Dave Morgan\X-Ray014 400um - FG  ON"
           r"\Stage Tilt ARXPS Si Wafer\Fitted Si tilt.vms")


class SiliconPresets(unittest.TestCase):
    """SiO2 on Si works like Al2O3 on Al: the preset guess must not land on SiC or the blank entry."""

    def test_material_guess_for_si_and_al(self):
        try:
            from batch_tab import BatchTab
        except Exception:
            self.skipTest("needs numpy + matplotlib")
        for el, ox, sub, ovl in (("Si 2p El", "Si 2p Ox", "Si (metal/elemental)", "SiO2"),
                                 ("Al 2p Metal", "Al 2p Oxide", "Al (metal)", "Al2O3")):
            self.assertEqual(BatchTab._guess_material(el, A.SUBSTRATES), sub)
            self.assertEqual(BatchTab._guess_material(ox, A.OVERLAYERS), ovl)

    def test_synthetic_si_series_recovers_the_oxide_thickness(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "si.vms")
            VM.write(path, S.two_state_series(d_nm=1.5, fitted=True))
            vf = VM.read(path)
        samples, _ = B.collect(vf, "Si2p Scan")
        states = B.fit_states(samples)
        self.assertEqual(len(states), 2)
        rows = B.compute(samples[:1], B.guess_roles(states), A.SUBSTRATES["Si (metal/elemental)"],
                         {states[1].key: A.OVERLAYERS["SiO2"]}, 90.0 - 0.0, A.volume_density)
        self.assertTrue(0.5 < rows[0].d_nm < 4.0)

    @unittest.skipUnless(os.path.exists(SI_REAL), "the user's Si tilt file is not on this machine")
    def test_real_si_wafer_zero_tilt(self):
        vf = VM.read(SI_REAL)
        samples, _ = B.collect(vf, "Si2p Scan")
        states = B.fit_states(samples)
        self.assertEqual([s.name for s in states], ["Si 2p El", "Si 2p Ox"])
        rows = B.compute(samples[:1], B.guess_roles(states), A.SUBSTRATES["Si (metal/elemental)"],
                         {states[1].key: A.OVERLAYERS["SiO2"]}, 90.0, A.volume_density)
        self.assertAlmostEqual(rows[0].d_nm, 1.83, delta=0.05)      # the ARXPS tab gave 1.72 nm


try:
    import tkinter as tk
    _probe = tk.Tk()
    _probe.destroy()
    import batch_tab  # noqa: F401
    HAVE_GUI = True
except Exception:
    HAVE_GUI = False


@unittest.skipUnless(HAVE_GUI, "needs a display and numpy + matplotlib")
class TabSmoke(unittest.TestCase):
    def test_load_and_calculate(self):
        path, _ = _file(self)
        app = A.App()
        app.withdraw()
        self.addCleanup(app.destroy)
        tab = app.batch
        self.assertTrue(tab.load_file(path))
        self.assertEqual(len(tab.samples), 3)
        self.assertEqual(tab.sub_preset.get(), "Al (metal)")
        self.assertEqual(sorted(v.get() for v in tab.role_vars.values()), ["Overlayer", "Substrate"])
        tab.calculate()
        self.assertEqual(len(tab.rows), 3)
        self.assertEqual(len(tab.table.get_children()), 3)
        self.assertTrue(all(r.d_nm > 0 for r in tab.rows))


if __name__ == "__main__":
    unittest.main()
