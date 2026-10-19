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
        self.assertEqual(tab.v["source"].get(), DT.SRC_FITS)
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
