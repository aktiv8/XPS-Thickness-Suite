"""
GUI smoke tests: build the real window and drive the Thickogram tab. They need a display and skip themselves
(rather than fail) when Tk cannot open one, e.g. on a headless CI box.

    python -m unittest discover -s tests -v
"""
import os
import sys
import unittest
from unittest import mock

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

import tkinter as tk  # noqa: E402

try:
    _probe = tk.Tk()
    _probe.destroy()
    HAVE_DISPLAY = True
except tk.TclError:
    HAVE_DISPLAY = False

import app as A  # noqa: E402
import thickogram as T  # noqa: E402

try:
    import thickogram_tab as TT
except ImportError:
    TT = None


def _make_app(testcase):
    app = A.App()
    app.withdraw()
    testcase.addCleanup(app.destroy)
    return app


@unittest.skipUnless(HAVE_DISPLAY and TT is not None, "needs a display and numpy + matplotlib")
class ThickogramTabSmoke(unittest.TestCase):
    def test_window_has_the_three_tabs_and_the_suite_title(self):
        app = _make_app(self)
        self.assertEqual(app.title(), "XPS Thickness Suite")
        nb = app.nametowidget(app.calc_tab.winfo_parent())
        self.assertEqual([nb.tab(i, "text") for i in range(nb.index("end"))],
                         ["Calculator & simulation", "Thickogram (different core levels)", "Method / References"])

    def test_example_reproduces_the_library_result(self):
        app = _make_app(self)
        tab = app.thickogram
        tab.load_example()
        tab.calculate()
        v = {k: float(tab.vars[k].get()) for k in TT.EXAMPLE}
        want = T.thickness(io=v["Io"], is_=v["Is"], so=v["So"], ss=v["Ss"], be_o=v["Eo"], be_s=v["Es"],
                           source_ev=v["Esource"], lam_nm=v["Lambda"], theta_normal_deg=v["Theta"])
        self.assertAlmostEqual(tab.result.t_nm, want.t_nm, places=12)
        self.assertEqual(tab.readout_var.get(), f"{want.t_nm:.2f} nm")
        self.assertEqual(tab.badge.cget("text"), "Within the printed nomogram scales")
        red = [ln for ln in tab.fig.axes[0].lines if ln.get_color() == "red"]
        self.assertEqual(len(red), 2)

    def test_blank_nomogram_is_drawn_before_any_calculation(self):
        app = _make_app(self)
        self.assertIsNone(app.thickogram.result)
        self.assertGreater(len(app.thickogram.fig.axes[0].lines), 100)

    def test_presets_prefill_binding_energy_and_lambda(self):
        app = _make_app(self)
        tab = app.thickogram
        tab.ov_preset.set("SiO2")
        tab._apply_preset("ov")
        self.assertEqual(float(tab.vars["Eo"].get()), A.OVERLAYERS["SiO2"].be_ev)
        self.assertEqual(float(tab.vars["Lambda"].get()), A.OVERLAYERS["SiO2"].imfp_nm)
        tab.sub_preset.set("Al (metal)")
        tab._apply_preset("sub")
        self.assertEqual(float(tab.vars["Es"].get()), A.SUBSTRATES["Al (metal)"].be_ev)

    def test_auger_preset_clears_the_be_and_says_why(self):
        app = _make_app(self)
        tab = app.thickogram
        tab.vars["Eo"].set("123")
        tab.ov_preset.set("CuO (LMM)")
        tab._apply_preset("ov")
        self.assertEqual(tab.vars["Eo"].get(), "")
        self.assertIn("Auger", tab.badge.cget("text"))

    def test_copy_from_calculator_converts_the_angle_and_copies_energies(self):
        app = _make_app(self)
        app.theta_var.set("60")                         # from the surface -> 30 from the normal
        s = app.substrate_row
        s.choice_var.set("Al (metal)")
        s._apply_preset(A.SUBSTRATES)
        for i, row in enumerate(app.overlayer_rows):
            row.choice_var.set("Al2O3" if i == 0 else "None")
            row._apply_preset(A.OVERLAYERS)
        tab = app.thickogram
        tab.copy_from_calculator()
        self.assertEqual(float(tab.vars["Theta"].get()), 30.0)
        self.assertEqual(float(tab.vars["Es"].get()), A.SUBSTRATES["Al (metal)"].be_ev)
        self.assertEqual(float(tab.vars["Eo"].get()), A.OVERLAYERS["Al2O3"].be_ev)
        self.assertEqual(float(tab.vars["Lambda"].get()), A.OVERLAYERS["Al2O3"].imfp_nm)
        self.assertIn("= 60", tab.theta_note.cget("text"))      # the live read-out converts back to the surface angle
        self.assertEqual(tab.ov_preset.get(), "Al2O3")

    def test_bad_input_shows_a_message_instead_of_crashing(self):
        app = _make_app(self)
        tab = app.thickogram
        tab.load_example()
        tab.vars["Io"].set("abc")
        with mock.patch.object(TT.messagebox, "showerror") as showerror:
            tab.calculate()
        self.assertEqual(showerror.call_count, 1)
        self.assertIn("Overlayer peak area Io", showerror.call_args[0][1])
        self.assertIsNone(tab.result)
        tab.vars["Io"].set("-5")
        with mock.patch.object(TT.messagebox, "showerror") as showerror:
            tab.calculate()
        self.assertEqual(showerror.call_count, 1)

    def test_clear_resets_the_form(self):
        app = _make_app(self)
        tab = app.thickogram
        tab.load_example()
        tab.clear()
        self.assertEqual(tab.vars["Io"].get(), "")
        self.assertEqual(tab.vars["Esource"].get(), "1486.6")
        self.assertEqual(tab.vars["Theta"].get(), "0")


@unittest.skipUnless(HAVE_DISPLAY, "needs a display")
class WithoutPlottingPackages(unittest.TestCase):
    def test_calculator_still_works_and_the_tab_shows_an_install_hint(self):
        with mock.patch.object(A, "HAVE_THICKOGRAM", False), mock.patch.object(A, "HAVE_PLOT", False):
            app = _make_app(self)
            self.assertFalse(hasattr(app, "thickogram"))
            texts = [w.cget("text") for w in app.thick_tab.winfo_children() if w.winfo_class() == "TLabel"]
            self.assertTrue(any("pip install numpy matplotlib" in t for t in texts))
            s = app.substrate_row
            s.choice_var.set("Al (metal)")
            s._apply_preset(A.SUBSTRATES)
            for i, row in enumerate(app.overlayer_rows):
                row.choice_var.set("Al2O3" if i == 0 else "None")
                row._apply_preset(A.OVERLAYERS)
                row.area_var.set("65" if i == 0 else "")
            s.area_var.set("35")
            app.theta_var.set("90")
            app._on_calculate()
            self.assertAlmostEqual(app.last_calc["res"]["d_nm"], 3.577, places=2)


if __name__ == "__main__":
    unittest.main()
