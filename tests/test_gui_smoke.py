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
    import arxps_tab as AT
    import depth_tab as DT
except ImportError:
    TT = AT = DT = None

import synthetic as S  # noqa: E402
import vamas as VM  # noqa: E402


def _make_app(testcase):
    app = A.App()
    app.withdraw()
    testcase.addCleanup(app.destroy)
    return app


@unittest.skipUnless(HAVE_DISPLAY and TT is not None, "needs a display and numpy + matplotlib")
class ThickogramTabSmoke(unittest.TestCase):
    def test_window_has_the_five_tabs_and_the_suite_title(self):
        app = _make_app(self)
        self.assertEqual(app.title(), "XPS Thickness Suite")
        nb = app.nametowidget(app.calc_tab.winfo_parent())
        self.assertEqual([nb.tab(i, "text") for i in range(nb.index("end"))],
                         ["Calculator & simulation", "Thickogram (different core levels)", "ARXPS (tilt series)",
                          "Depth profile (ARXPS)", "Method / References"])

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


@unittest.skipUnless(HAVE_DISPLAY and AT is not None, "needs a display and numpy + matplotlib")
class ARXPSTabSmoke(unittest.TestCase):
    """The ARXPS tab driven with synthetic VAMAS tilt series (see synthetic.py)."""

    def _write(self, blocks):
        import tempfile
        d = tempfile.TemporaryDirectory()
        self.addCleanup(d.cleanup)
        path = os.path.join(d.name, "series.vms")
        VM.write(path, blocks)
        return path

    def _states_tab(self, blocks):
        app = _make_app(self)
        tab = app.arxps
        self.assertTrue(tab.load_file(self._write(blocks)))
        tab.sub_preset.set("Si (metal/elemental)")
        tab._apply_preset("sub")
        tab.ov_preset.set("SiO2")
        tab._apply_preset("ov")
        return app, tab

    def test_empty_tab_shows_a_prompt_and_no_result(self):
        app = _make_app(self)
        self.assertEqual(app.arxps.readout_var.get(), "-.-- nm")
        self.assertEqual(app.arxps.results, {})

    def test_presets_fill_lambda_and_rinf_like_the_calculator(self):
        app, tab = self._states_tab(S.two_state_series())
        self.assertAlmostEqual(float(tab.v["lam"].get()), 3.89)
        self.assertAlmostEqual(float(tab.v["rinf"].get()), 0.5547, places=3)

    def test_two_states_calculation_matches_the_library_and_fills_the_table(self):
        import arxps as AR
        app, tab = self._states_tab(S.two_state_series(d_nm=1.5))
        tab.auto_windows()
        self.assertTrue(tab.v["env_lo"].get() and tab.v["split"].get())
        tab.v["max_angle"].set("80")
        tab.calculate()
        self.assertEqual(len(tab.results), 1)
        res = tab.results[None]
        self.assertAlmostEqual(res.d_nm, 1.5, delta=0.03)
        self.assertEqual(tab.readout_var.get(), f"{res.d_nm:.2f} nm")
        self.assertEqual(len(tab.table.get_children()), 6)
        self.assertEqual(tab.badge.cget("text"), "The angles agree: x falls on a line through the origin")
        self.assertEqual(len(tab.figs["fit"][0].axes), 1)
        self.assertGreater(len(tab.figs["spectra"][0].axes[0].lines), 6)

    def test_high_angles_are_excluded_and_say_so(self):
        app, tab = self._states_tab(S.two_state_series())
        tab.calculate()                               # windows are found automatically when blank
        self.assertIn("75 deg", tab.badge.cget("text"))
        self.assertEqual(sum(r.used for r in tab.results[None].rows), 5)

    def test_two_peaks_mode(self):
        app = _make_app(self)
        tab = app.arxps
        tab.load_file(self._write(S.two_peak_series(d_nm=1.5)))
        tab.mode.set(AT.MODE_PEAKS)
        tab._mode_changed()
        tab.v["region_sub"].set("Si2p Scan")
        tab.v["region_ov"].set("O1s Scan")
        tab._region_changed()
        for key, value in {"lam": "3.0", "so": "2.93", "ss": "0.82", "eo": "532", "es": "99.3"}.items():
            tab.v[key].set(value)
        tab.v["max_angle"].set("80")
        tab.calculate()
        self.assertAlmostEqual(tab.results[None].d_nm, 1.5, delta=0.02)
        self.assertAlmostEqual(tab.results[None].b, (1486.68 - 532) / (1486.68 - 99.3), places=9)

    def test_several_azimuths_give_a_mean_a_selector_and_a_polar_plot(self):
        app, tab = self._states_tab(S.two_state_series(azimuths=(0.0, 90.0), d_by_azimuth={0.0: 1.2, 90.0: 1.8}))
        self.assertEqual(tab.v["azim_src"].get(), "Experimental variable 2")
        tab.v["max_angle"].set("80")
        tab.calculate()
        self.assertEqual(sorted(tab.results), [0.0, 90.0])
        self.assertAlmostEqual(tab.summary[0], 1.5, delta=0.05)
        self.assertEqual(tab.readout_var.get(), f"{tab.summary[0]:.2f} nm")
        self.assertEqual(list(tab.view_cb.cget("values")), ["0", "90"])
        self.assertTrue(tab.view_cb.winfo_manager())
        self.assertEqual(tab.figs["azimuth"][0].axes[0].name, "polar")
        tab.v["view_az"].set("90")
        tab._view_changed()
        self.assertEqual(tab.table.item(tab.table.get_children()[0], "values")[0], "90")

    def test_a_single_azimuth_says_there_is_nothing_to_plot(self):
        app, tab = self._states_tab(S.two_state_series())
        tab.calculate()
        texts = [t.get_text() for t in tab.figs["azimuth"][0].axes[0].texts]
        self.assertTrue(any("single azimuth" in t for t in texts))

    def test_errors_become_messages(self):
        app = _make_app(self)
        tab = app.arxps
        with mock.patch.object(AT.messagebox, "showerror") as shown:
            tab.calculate()                                   # nothing loaded
            self.assertIn("Open a VAMAS file", shown.call_args[0][1])
        with mock.patch.object(AT.messagebox, "showerror") as shown:
            path = os.path.join(os.path.dirname(__file__), "test_arxps.py")       # not a VAMAS file
            self.assertFalse(tab.load_file(path))
            self.assertTrue(shown.called)
        self.assertIsNone(tab.file)

    def test_a_missing_overlayer_number_is_reported(self):
        app, tab = self._states_tab(S.two_state_series())
        tab.v["rinf"].set("")
        with mock.patch.object(AT.messagebox, "showerror") as shown:
            tab.calculate()
        self.assertIn("R∞", shown.call_args[0][1])

    def test_a_file_with_a_casaxps_fit_defaults_to_the_fitted_components(self):
        app, tab = self._states_tab(S.two_state_series(fitted=True))
        self.assertEqual(tab.area_source.get(), AT.SRC_FITS)
        self.assertEqual(tab.setup_nb.tab(tab.wins_page, "text"), "3  Components")
        self.assertIn("disabled", tab.per_second_cb.state())
        self.assertIn("Si El", tab.fit_info.cget("text"))
        self.assertEqual({k: v.get() for (sec, k), v in tab.roles.items()}, {"i0": "Substrate", "i1": "Overlayer"})
        tab.calculate()                                   # no windows are needed
        self.assertAlmostEqual(tab.results[None].d_nm, 1.5, places=4)
        self.assertEqual(tab.readout_var.get(), "1.50 nm")
        self.assertIn("stored CasaXPS fits", tab.details.get("1.0", "end"))
        self.assertEqual(len(tab.table.get_children()), 6)
        spectra = tab.figs["spectra"][0].axes[0]
        self.assertGreaterEqual(len(spectra.lines), 6 + 2)               # six spectra and the two state markers

    def test_the_assignment_can_be_edited(self):
        app, tab = self._states_tab(S.two_state_series(fitted=True))
        tab.roles[("main", "i1")].set("Ignore")
        with mock.patch.object(AT.messagebox, "showerror") as shown:
            tab.calculate()
        self.assertIn("at least one overlayer", shown.call_args[0][1])
        tab._refresh_components()                         # Auto-assign puts the guess back
        self.assertEqual(tab.roles[("main", "i1")].get(), "Overlayer")

    def test_switching_back_to_windows_works_on_the_same_file(self):
        app, tab = self._states_tab(S.two_state_series(fitted=True))
        tab.area_source.set(AT.SRC_WINDOWS)
        tab._source_changed()
        self.assertEqual(tab.setup_nb.tab(tab.wins_page, "text"), "3  Windows")
        self.assertNotIn("disabled", tab.per_second_cb.state())
        tab.v["max_angle"].set("80")
        tab.calculate()
        self.assertAlmostEqual(tab.results[None].d_nm, 1.5, delta=0.03)
        self.assertNotIn("CasaXPS", tab.details.get("1.0", "end"))

    def test_a_file_without_a_fit_stays_on_windows_and_the_fitted_choice_explains(self):
        app, tab = self._states_tab(S.two_state_series())
        self.assertEqual(tab.area_source.get(), AT.SRC_WINDOWS)
        tab.area_source.set(AT.SRC_FITS)
        tab._source_changed()
        self.assertIn("no CasaXPS fit", tab.comp_info.cget("text"))
        with mock.patch.object(AT.messagebox, "showerror") as shown:
            tab.calculate()
        self.assertTrue(shown.called)

    def test_two_peaks_from_fits_prefills_the_rsfs_and_binding_energies(self):
        app = _make_app(self)
        tab = app.arxps
        tab.load_file(self._write(S.two_peak_series(d_nm=1.5, fitted=True)))
        tab.mode.set(AT.MODE_PEAKS)
        tab._mode_changed()
        tab.v["region_sub"].set("Si2p Scan")
        tab.v["region_ov"].set("O1s Scan")
        tab._region_changed()
        self.assertEqual(tab.area_source.get(), AT.SRC_FITS)
        self.assertEqual((tab.v["so"].get(), tab.v["ss"].get()), ("2.93", "0.82"))
        self.assertEqual((tab.v["eo"].get(), tab.v["es"].get()), ("532.00", "99.30"))
        tab.v["lam"].set("3.0")
        tab.v["max_angle"].set("80")
        tab.calculate()
        self.assertAlmostEqual(tab.results[None].d_nm, 1.5, places=4)
        self.assertEqual(len(tab.figs["spectra"][0].axes), 2)

    def test_changing_the_region_clears_the_windows(self):
        app, tab = self._states_tab(S.two_state_series())
        tab.auto_windows()
        tab._region_changed()
        self.assertEqual(tab.v["env_lo"].get(), "")


@unittest.skipUnless(HAVE_DISPLAY and DT is not None, "needs a display and numpy + matplotlib")
class DepthTabSmoke(unittest.TestCase):
    """The depth-profile tab (Surrey Trend Plot) driven with synthetic VAMAS files."""

    def _write(self, blocks):
        import tempfile
        d = tempfile.TemporaryDirectory()
        self.addCleanup(d.cleanup)
        path = os.path.join(d.name, "series.vms")
        VM.write(path, blocks)
        return path

    def _tab(self, blocks):
        app = _make_app(self)
        tab = app.depth
        self.assertTrue(tab.load_file(self._write(blocks)))
        return app, tab

    def _table(self, tab):
        return {tab.trend_table.item(i, "values")[0]: tab.trend_table.item(i, "values")
                for i in tab.trend_table.get_children()}

    def test_empty_tab_shows_a_prompt(self):
        app = _make_app(self)
        self.assertEqual(app.depth.species, [])
        self.assertEqual(len(app.depth.trend_table.get_children()), 0)

    def test_a_fitted_file_gives_one_species_per_state_and_the_exact_trend(self):
        app, tab = self._tab(S.two_state_series(d_nm=1.5, fitted=True))
        self.assertEqual(tab.v["source"].get(), DT.SRC_AUTO)
        self.assertEqual([sp.name for sp in tab.species], ["Si El", "Si Ox"])
        self.assertAlmostEqual(tab.species[0].lam_nm, 3.89, delta=0.05)          # TPP-2M default for an oxide-like matrix
        rows = self._table(tab)
        self.assertAlmostEqual(float(rows["Si El"][2]), -1.5 / tab.species[0].lam_nm, delta=0.003)
        self.assertAlmostEqual(float(rows["Si El"][5]), 1.5, delta=0.01)         # apparent depth = -lambda x slope
        self.assertEqual(len(tab.trend_fig.axes), 1)
        self.assertEqual(len(tab.trend_fig.axes[0].get_legend().get_texts()), 2)

    def test_a_file_without_fits_uses_whole_regions(self):
        app, tab = self._tab(S.two_state_series())
        self.assertEqual(tab.v["source"].get(), DT.SRC_REGIONS)
        self.assertEqual([sp.name for sp in tab.species], ["Si 2p"])

    def test_the_fitted_choice_without_fits_explains_itself(self):
        app, tab = self._tab(S.two_state_series())
        tab.v["source"].set(DT.SRC_FITS)
        tab.reload()
        self.assertEqual(tab.species, [])
        self.assertIn("CasaXPS fit", tab.status.get())

    def test_unticking_a_species_and_editing_lambda_update_the_table(self):
        app, tab = self._tab(S.two_state_series(d_nm=1.5, fitted=True))
        tab.rows[1]["use"].set(False)
        tab.refresh()
        self.assertEqual(list(self._table(tab)), ["Si El"])
        tab.rows[0]["lam"].set("2.0")
        tab.refresh()
        self.assertAlmostEqual(float(self._table(tab)["Si El"][5]), 2.0 * 1.5 / 3.89, delta=0.01)

    def test_the_axis_and_normalise_options_redraw(self):
        app, tab = self._tab(S.two_state_series(fitted=True))
        tab.x_mode.set("theta")
        tab.normalise.set(False)
        tab._draw()
        ax = tab.trend_fig.axes[0]
        self.assertIn("Emission angle", ax.get_xlabel())
        self.assertIn("ln( I )", ax.get_ylabel())

    def test_stratification_places_the_substrate_state_below_the_overlayer_state(self):
        app, tab = self._tab(S.two_state_series(d_nm=1.5, fitted=True))
        self.assertEqual((tab.strat_a.get(), tab.strat_b.get()), ("0", "60"))        # 75 is above the 70 degree limit
        rows = {tab.strat_table.item(i, "values")[0]: tab.strat_table.item(i, "values")
                for i in tab.strat_table.get_children()}
        self.assertEqual(set(rows), {"Si El", "Si Ox"})
        self.assertGreater(float(rows["Si El"][1]), float(rows["Si Ox"][1]))         # the substrate is deeper
        self.assertAlmostEqual(sum(float(r[3]) for r in rows.values()), 100.0, delta=0.2)
        self.assertEqual(len(tab.strat_fig.axes), 1)

    def test_the_stratification_angles_and_modes_can_be_changed(self):
        app, tab = self._tab(S.two_state_series(fitted=True))
        tab.strat_b.set("45")
        tab.strat_all.set(True)
        tab._draw_strat()
        vals = [tab.strat_table.item(i, "values") for i in tab.strat_table.get_children()]
        self.assertEqual(len(vals), 2)
        self.assertEqual(vals[0][4], "-")                                             # no single ratio in the fit mode

    def test_the_geometry_option_divides_the_intensities_by_cos_theta(self):
        app, tab = self._tab(S.two_state_series(fitted=True))
        before = tab.species[0].intensity.copy()
        tab.v["geometry"].set(DT.GEO_COS)
        tab.reload()
        import numpy as np
        want = before / np.cos(np.radians(tab.species[0].angles))
        self.assertTrue(np.allclose(tab.species[0].intensity, want))

    def _layers_tab(self):
        blocks = S.layered_series([0.8, 1.5], [dict(name="C 1s", be=285.0, layer=0), dict(name="O 1s", be=532.0, layer=1),
                                               dict(name="Al 2p", be=73.0, layer=-1)])
        app, tab = self._tab(blocks)
        tab.v["max_angle"].set("80")
        for sp in tab.species:
            sp.lam_nm = 3.0
        for row in tab.rows:
            row["lam"].set("3.0")
        tab.transmission.set(False)
        tab.reload()
        for row in tab.rows:
            row["lam"].set("3.0")
        return app, tab

    def test_layers_default_assignment_follows_depth_and_the_fit_recovers_the_stack(self):
        app, tab = self._layers_tab()
        self.assertEqual([v.get() for v in tab.layer_vars], ["1", "2", DT.SUB_LABEL])
        tab.fit_layers_now()
        self.assertAlmostEqual(float(tab.layer_d[0].get()), 0.8, delta=0.1)
        self.assertAlmostEqual(float(tab.layer_d[1].get()), 1.5, delta=0.15)
        self.assertEqual(len(tab.layer_fig.axes), 2)                      # curves and the layer stack
        rows = [tab.layer_table.item(i, "values") for i in tab.layer_table.get_children()]
        self.assertEqual([r[0] for r in rows], ["1", "2", "Sub"])
        self.assertEqual(rows[0][1], "C 1s")

    def test_murdoch_panel_fits_and_the_eal_button_uses_the_kinetic_energy(self):
        app, tab = self._layers_tab()
        self.assertEqual([g["name"] for g in tab.mur_grid], ["C 1s", "O 1s", "Al 2p"])
        self.assertEqual([g["order"].get() for g in tab.mur_grid], ["1", "2", "3"])
        tab.mur_log.set("-4")
        tab._mur_fit_now()
        self.assertIsNotNone(tab.mur_fit)
        self.assertEqual(len(tab.mur_table.get_children()), 3)
        tab._mur_scurve_now()
        self.assertIsNotNone(tab.mur_scurve)
        tab._mur_ensemble_now()
        self.assertIsNotNone(tab.mur_ens)
        tab._mur_eal()
        self.assertIn("kinetic energy", tab.mur_note.get())
        tab.mur_grid[1]["order"].set("1")                                  # two peaks in the same place: refused, no crash
        tab._mur_fit_now()
        self.assertIn("different position", tab.mur_note.get())

    def test_murdoch_use_boxes_switch_peaks_on_and_off(self):
        app, tab = self._layers_tab()
        tab.mur_log.set("-4")
        self.assertTrue(all(g["use"].get() for g in tab.mur_grid))
        tab.mur_grid[0]["use"].set(False)                                  # untick C 1s
        tab._mur_use_changed()
        self.assertEqual(sorted(g["order"].get() for g in tab.mur_grid if g["use"].get()), ["1", "2"])
        tab._mur_fit_now()
        self.assertEqual(tab.mur_fit.names, ["O 1s", "Al 2p"])
        tab.mur_grid[1]["use"].set(False)                                  # only one left
        tab._mur_use_changed()
        tab._mur_fit_now()
        self.assertIn("at least two", tab.mur_note.get())
        tab.mur_grid[0]["use"].set(True)
        tab.mur_grid[1]["use"].set(True)
        tab._mur_use_changed()
        tab._mur_fit_now()
        self.assertEqual(len(tab.mur_fit.names), 3)

    def test_scanning_the_ranges_fills_the_table_and_editing_a_thickness_clears_them(self):
        app, tab = self._layers_tab()
        tab.scan_layers_now()
        rows = [tab.layer_table.item(i, "values") for i in tab.layer_table.get_children()]
        self.assertIn(" to ", rows[0][3])
        self.assertIn(" to ", rows[1][3])
        tab.layer_d[0].set("2.0")                                          # the user moves a thickness: the ranges are stale
        self.assertIsNone(tab.layer_fit)

    def test_changing_the_number_of_layers_rebuilds_the_assignment(self):
        app, tab = self._layers_tab()
        tab.n_layers.set("3")
        tab._default_layers()
        self.assertEqual([v.get() for v in tab.layer_vars], ["1", "2", DT.SUB_LABEL])
        self.assertIn("3", tab._layer_labels())
        tab.n_layers.set("1")
        tab._default_layers()
        self.assertEqual([v.get() for v in tab.layer_vars], ["1", "1", DT.SUB_LABEL])

    def test_fit_without_a_species_in_any_layer_says_so(self):
        app, tab = self._layers_tab()
        for v in tab.layer_vars:
            v.set(DT.SUB_LABEL)
        tab.fit_layers_now()
        self.assertIn("at least one species to a layer", tab.layer_note.get())

    def test_fits_where_available_mixes_fitted_states_and_whole_regions(self):
        blocks = S.two_state_series(fitted=True) + S.layered_series([0.8], [dict(name="C 1s", be=285.0, layer=0)])
        app, tab = self._tab(blocks)
        self.assertEqual(tab.v["source"].get(), DT.SRC_AUTO)
        self.assertEqual(sorted(sp.name for sp in tab.species), ["C 1s", "Si El", "Si Ox"])
        self.assertEqual(sorted(sp.source for sp in tab.species), ["fit", "fit", "region"])

    def _stack(self, noise=False):
        blocks = S.layered_series([0.8, 1.5], [dict(name="C 1s", be=285.0, layer=0), dict(name="O 1s", be=532.0, layer=1),
                                               dict(name="Al 2p", be=73.0, layer=-1)])
        app, tab = self._tab(blocks)
        tab.v["max_angle"].set("80")
        tab.transmission.set(False)
        tab.reload()
        for row in tab.rows:
            row["lam"].set("3.0")
        tab.active()
        return app, tab

    def test_the_profile_panel_gives_one_profile_per_species_and_a_table(self):
        app, tab = self._stack()
        tab._calc_profile()
        self.assertEqual([p[0].name for p in tab.profiles], ["C 1s", "O 1s", "Al 2p"])
        rows = [tab.prof_table.item(i, "values") for i in tab.prof_table.get_children()]
        self.assertEqual([r[0] for r in rows], ["C 1s", "O 1s", "Al 2p"])
        self.assertEqual(len(tab.prof_fig.axes), 2)                        # profiles, and measured against predicted
        self.assertLess(float(rows[0][2]), 5.0)

    def test_non_negative_bulk_constant_and_fraction_options(self):
        app, tab = self._stack()
        tab.prof_nonneg.set(True)
        tab.prof_bulk[2][0].set(tab.BULK_LABELS["constant"])
        tab.prof_bulk[0][0].set(tab.BULK_LABELS["zero"])
        tab.prof_bulk[0][1].set("4")
        tab.prof_frac.set(True)
        tab._calc_profile()
        by_name = {p[0].name: p for p in tab.profiles}
        import numpy as np
        for sp, f, bulk, zc in tab.profiles:
            self.assertGreaterEqual(float(f.c.min()), 0.0)
        c_zero = by_name["C 1s"][1]
        self.assertTrue(np.all(c_zero.c[c_zero.z_nm > 4.0] == 0.0))
        al = by_name["Al 2p"][1]
        self.assertTrue(np.allclose(al.c[al.z_nm > 8.0], al.c[al.z_nm > 8.0][0]))
        self.assertIn("Atomic fraction", tab.prof_fig.axes[0].get_ylabel())

    def test_a_fixed_alpha_changes_the_smoothness(self):
        app, tab = self._stack()
        tab._calc_profile()
        auto = tab.profiles[0][1].alpha
        tab.prof_auto.set(False)
        tab.prof_logalpha.set("-6")
        tab._calc_profile()
        self.assertAlmostEqual(tab.profiles[0][1].alpha, 1e-6, places=12)
        self.assertNotAlmostEqual(auto, 1e-6)

    def test_bad_profile_settings_are_reported_not_raised(self):
        app, tab = self._stack()
        tab.prof_n.set("3")
        tab._calc_profile()
        self.assertEqual(tab.profiles, [])
        self.assertIn("between 5 and 60", tab.prof_note.get())

    def test_the_attenuation_length_factor_scales_the_starting_lambdas(self):
        app, tab = self._tab(S.two_state_series(fitted=True))
        before = [sp.lam_nm for sp in tab.species]
        tab.v["al_factor"].set("0.8")
        tab.reload()
        for b, sp in zip(before, tab.species):
            self.assertAlmostEqual(sp.lam_nm, 0.8 * b, delta=0.01)      # the table shows two decimals
        tab.v["al_factor"].set("-1")
        with mock.patch.object(DT.messagebox, "showerror"):
            tab.reload()
        self.assertEqual(tab.species, [])

    def test_with_fits_the_substrate_state_is_the_default_substrate_and_a_mixed_region_starts_unticked(self):
        import numpy as np
        blocks = S.two_state_series(fitted=True) + S.layered_series([0.8, 9.0], [dict(name="Mix 2s", be=150.0, layer=1)])
        app, tab = self._tab(blocks)
        names = [sp.name for sp in tab.species]
        sub = [n for n, v in zip(names, tab.layer_vars) if v.get() == DT.SUB_LABEL]
        self.assertEqual(sub, ["Si El"])
        self.assertTrue(tab.rows[names.index("Si El")]["use"].get())

    def test_it_can_take_the_file_from_the_arxps_tab(self):
        app = _make_app(self)
        app.arxps.load_file(self._write(S.two_state_series(fitted=True)))
        app.depth.use_arxps_file()
        self.assertIsNotNone(app.depth.file)
        self.assertEqual([sp.name for sp in app.depth.species], ["Si El", "Si Ox"])

    def test_a_bad_file_is_refused_with_a_message(self):
        app = _make_app(self)
        with mock.patch.object(DT.messagebox, "showerror") as shown:
            self.assertFalse(app.depth.load_file(os.path.join(os.path.dirname(__file__), "test_arxps.py")))
        self.assertTrue(shown.called)
        self.assertIsNone(app.depth.file)


@unittest.skipUnless(HAVE_DISPLAY, "needs a display")
class WithoutPlottingPackages(unittest.TestCase):
    def test_calculator_still_works_and_the_tab_shows_an_install_hint(self):
        with mock.patch.object(A, "HAVE_THICKOGRAM", False), mock.patch.object(A, "HAVE_ARXPS", False), \
                mock.patch.object(A, "HAVE_DEPTH", False), mock.patch.object(A, "HAVE_PLOT", False):
            app = _make_app(self)
            self.assertFalse(hasattr(app, "thickogram"))
            self.assertFalse(hasattr(app, "arxps"))
            self.assertFalse(hasattr(app, "depth"))
            texts = [w.cget("text") for w in app.thick_tab.winfo_children() if w.winfo_class() == "TLabel"]
            self.assertTrue(any("pip install numpy matplotlib" in t for t in texts))
            texts = [w.cget("text") for w in app.arxps_tab.winfo_children() if w.winfo_class() == "TLabel"]
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
