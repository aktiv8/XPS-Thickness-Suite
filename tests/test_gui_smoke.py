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


@unittest.skipUnless(HAVE_DISPLAY and TT is not None, "needs a display and numpy + matplotlib")
class TopofactorInTheTab(unittest.TestCase):
    """The Topofactor controls of the Thickogram tab (Shard, Wang & Spencer 2009)."""

    @staticmethod
    def _label(key):
        return TT.TFAC.SHAPES[key].label

    def _calc(self, shape="flat", mode="eq6"):
        app = _make_app(self)
        tab = app.thickogram
        tab.load_example()
        tab.shape_var.set(self._label(shape))
        tab.mode_var.set(tab._mode_label(mode))
        tab.calculate()
        return tab

    def _figure_text(self, tab):
        return " ".join(t.get_text() for t in tab.fig.axes[0].texts)

    def test_flat_is_the_default_and_leaves_the_planar_result(self):
        tab = self._calc()
        self.assertIsNone(tab.corrected)
        self.assertEqual(tab.readout_var.get(), f"{tab.result.t_nm:.2f} nm")
        self.assertNotIn("Topofactor", self._figure_text(tab))

    def test_cylinder_shows_the_corrected_thickness_everywhere(self):
        tab = self._calc("cylinder")
        want = TT.TFAC.apply_topofactor(tab.result, "cylinder", "eq6")
        self.assertAlmostEqual(tab.corrected.t_nm, want.t_nm, places=12)
        self.assertLess(tab.corrected.t_nm, tab.result.t_nm)
        self.assertEqual(tab.readout_var.get(), f"{want.t_nm:.2f} nm")
        self.assertIn("Topofactor", tab.readout_sub_var.get())
        self.assertIn("CORRECTED THICKNESS", tab.details.get("1.0", "end"))
        self.assertIn(f"T = {want.t_nm:.2f} nm", self._figure_text(tab))
        self.assertIn(f"planar t = {tab.result.t_nm:.2f} nm", self._figure_text(tab))
        # the demo's B (paper) is 0.65, outside the usual 0.7-1.4, so the badge says so rather than "applied"
        self.assertIn("outside the usual range", tab.badge.cget("text"))

    def test_an_in_range_b_gives_a_clean_applied_badge(self):
        app = _make_app(self)
        tab = app.thickogram
        for key, value in {"Eo": "285", "Es": "99.5", "Io": "800", "Is": "3000", "So": "1.0", "Ss": "3.0",
                           "Lambda": "3.0", "Theta": "0", "Esource": "1486.6"}.items():
            tab.vars[key].set(value)                                          # C 1s overlayer on Si 2p: B(paper) = 0.90
        tab.shape_var.set(self._label("cylinder"))
        tab.calculate()
        self.assertAlmostEqual(tab.corrected.b_paper, 0.90, delta=0.01)
        self.assertTrue(tab.badge.cget("text").startswith("Topofactor applied"), tab.badge.cget("text"))

    def test_changing_the_shape_afterwards_reapplies_without_recalculating(self):
        tab = self._calc("cylinder")
        planar = tab.result
        tab.shape_var.set(self._label("sphere"))
        tab._morphology_changed()
        self.assertIs(tab.result, planar)                                   # the planar solution is untouched
        self.assertAlmostEqual(tab.corrected.t_nm, TT.TFAC.apply_topofactor(planar, "sphere").t_nm, places=12)
        tab.shape_var.set(self._label("flat"))
        tab._morphology_changed()
        self.assertIsNone(tab.corrected)
        self.assertEqual(tab.readout_var.get(), f"{planar.t_nm:.2f} nm")

    def test_constant_mode_uses_the_papers_average(self):
        tab = self._calc("sphere_shadowed", "constant")
        self.assertAlmostEqual(tab.corrected.tf, 0.71)
        self.assertEqual(tab.corrected.regime, "constant")

    def test_off_normal_angle_raises_a_warning_badge(self):
        app = _make_app(self)
        tab = app.thickogram
        tab.load_example()
        tab.vars["Theta"].set("30")
        tab.shape_var.set(self._label("cylinder"))
        tab.calculate()
        self.assertIn("theta = 0", tab.badge.cget("text"))

    def test_clear_returns_to_a_flat_sample(self):
        tab = self._calc("cylinder", "constant")
        tab.clear()
        self.assertEqual(tab._shape_key(), "flat")
        self.assertEqual(tab._mode_key(), "eq6")


@unittest.skipUnless(HAVE_DISPLAY, "needs a display")
class HelpViewerSmoke(unittest.TestCase):
    """The in-app documentation (docs/*.md rendered in a Tk Text with typeset equations)."""

    @staticmethod
    def _viewer(testcase, **kw):
        import help_viewer
        root = tk.Tk()
        root.withdraw()
        testcase.addCleanup(root.destroy)
        return help_viewer, help_viewer.build_viewer(root, **kw)

    def test_help_tab_builds_the_viewer_when_it_is_opened(self):
        app = _make_app(self)
        self.assertIsNone(app._help_view)                       # not built at start-up (keeps launch quick)
        app.nb.select(app.help_tab)
        app.update()
        self.assertIsNotNone(app._help_view)
        self.assertTrue(app._help_view.text.get("1.0", "end").strip())

    def test_every_page_renders_in_strict_mode(self):
        hv, view = self._viewer(self, strict=True)
        for _, page in hv.PAGES:
            view.show(page)
            self.assertEqual(view.page, page)
            self.assertGreater(len(view.text.get("1.0", "end")), 500, page)

    def test_equations_are_typeset_as_images_when_matplotlib_is_present(self):
        hv, view = self._viewer(self, strict=True)
        if not view.math.available:
            self.skipTest("matplotlib not installed")
        view.show("methods.md")
        self.assertGreater(len(view.text.image_names()), 30)
        self.assertTrue(all(v is not None for v in view.math.cache.values()), "an equation fell back to plain text")

    def test_links_navigate_between_pages_and_to_anchors(self):
        hv, view = self._viewer(self, strict=True)
        view.open_link("methods.md#spectrum-simulation")
        self.assertEqual(view.page, "methods.md")
        self.assertIn("spectrum-simulation", view._anchors)
        view.open_link("#imfps-from-tpp-2m")
        self.assertEqual(view.page, "methods.md")
        view.open_link("references.md")
        self.assertEqual(view.page, "references.md")
        with mock.patch.object(hv.webbrowser, "open") as opened:
            view.open_link("https://doi.org/10.1002/sia.3044")
        opened.assert_called_once_with("https://doi.org/10.1002/sia.3044")

    def test_missing_docs_folder_gives_a_pointer_instead_of_a_crash(self):
        import help_viewer
        with mock.patch.object(help_viewer, "docs_dir", return_value=os.path.join(os.path.dirname(__file__), "nope")):
            _, view = self._viewer(self)
        text = view.text.get("1.0", "end")
        self.assertIn("not found", text)
        self.assertIn("github.com/aktiv8/XPS-Thickness-Suite", text)

    def test_without_matplotlib_the_equations_are_readable_text(self):
        import help_viewer
        math = help_viewer.MathRenderer()
        math.available = False
        _, view = self._viewer(self, math=math)
        view.show("methods.md")
        text = view.text.get("1.0", "end")
        self.assertIn("λ_ox sinθ ln", text)                     # the Strohmeier equation, as Unicode text
        self.assertNotIn("\\frac", text)
        self.assertEqual(len(view.text.image_names()), 2)       # only the two figures remain images


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
