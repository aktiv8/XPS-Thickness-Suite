"""
The "Thickogram" tab of the XPS Thickness Suite: Cumpson's method for an overlayer and a substrate measured from
DIFFERENT core levels, with the classic nomogram drawn (embedded, with the corrected red line).

Maths and drawing live in thickogram.py; this module is only the Tk/matplotlib glue. It is a ttk.Frame that app.py
puts into a notebook tab. Needs numpy + matplotlib (app.py shows an install hint instead if they are missing).

Angle convention: THIS tab uses theta from the surface NORMAL (0 = normal emission), as in Cumpson's paper and the
original script. The Calculator tab uses theta from the SURFACE; "Copy from Calculator" converts (90 - theta).
"""
import tkinter as tk
from tkinter import ttk, messagebox

from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

import thickogram as T
import topofactor as TFAC
from gui_common import ACCENT, BADGE, CARD, MUTED, Tooltip, _to_float

SHORT_SHAPE = {"cylinder": "cylinder", "sphere": "sphere", "sphere_shadowed": "shadowed sphere"}
SHORT_MODE = {"eq6": "Eq. (6)", "constant": "Constant average"}       # the compact layout has less room

# key, label, hint (the hints are the original script's)
FIELDS = [
    ("Eo", "Overlayer binding energy Eo (eV)", "e.g. 707 for Fe 2p"),
    ("Es", "Substrate binding energy Es (eV)", "e.g. 100 for Si 2p"),
    ("Io", "Overlayer peak area Io", "Peak area as measured in your software for the overlayer"),
    ("Is", "Substrate peak area Is", "Peak area as measured in your software for the substrate"),
    ("So", "Overlayer RSF So", "Sensitivity factor of the overlayer peak (optionally multiplied by the atomic fraction "
                               "of the element, as in Eq. 10 of Shard et al. 2009)"),
    ("Ss", "Substrate RSF Ss", "Sensitivity factor of the substrate peak (optionally multiplied by the atomic "
                               "fraction of the element, as in Eq. 10 of Shard et al. 2009)"),
    ("Lambda", "Lambda L (nm)", "Attenuation length of the overlayer's photoelectrons in the overlayer"),
    ("Theta", "Theta from the surface NORMAL (deg)", "0 = normal emission (Cumpson's convention)"),
    ("Esource", "X-ray source energy (eV)", "e.g. 1486.6 (Al Ka) or 1253.6 (Mg Ka)"),
]

# Demo numbers only (not a real measurement): Fe 2p overlayer on Si 2p, Al Ka
EXAMPLE = {"Eo": "707", "Es": "100", "Io": "20000", "Is": "6000", "So": "10.8", "Ss": "0.82", "Lambda": "1.6",
           "Theta": "0", "Esource": "1486.6"}


class ThickogramTab(ttk.Frame):
    def __init__(self, parent, app, substrates, overlayers, compact=False):
        super().__init__(parent)
        self.app, self.substrates, self.overlayers, self.compact = app, substrates, overlayers, compact
        self.vars = {key: tk.StringVar() for key, _, _ in FIELDS}
        self.vars["Esource"].set("1486.6")
        self.vars["Theta"].set("0")
        self.result = None          # planar Thickogram result
        self.corrected = None       # topofactor.CorrectedResult, or None for a flat sample
        self.entries = {}
        self.shape_var = tk.StringVar(value=TFAC.SHAPES["flat"].label)
        self.mode_var = tk.StringVar(value=self._mode_label("eq6"))
        self.shape_hint = tk.StringVar()

        self.columnconfigure(0, weight=0)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)
        left = ttk.Frame(self)
        left.grid(row=0, column=0, sticky="ns", padx=(0, 8), pady=(8, 4))
        left.rowconfigure(1, weight=1)
        self._build_inputs(left).grid(row=0, column=0, sticky="ew")
        self._build_results(left).grid(row=1, column=0, sticky="nsew", pady=(8, 0))
        self._build_plot().grid(row=0, column=1, sticky="nsew", pady=(8, 4))
        self._draw()

    # --- inputs -----------------------------------------------------------------------------
    def _build_inputs(self, parent):
        card = ttk.LabelFrame(parent, text="  Inputs  ", style="Card.TLabelframe", padding=(10, 4, 10, 8))
        row = 0

        presets = ttk.Frame(card, style="Card.TFrame")
        presets.grid(row=row, column=0, columnspan=3, sticky="ew", pady=(0, 6))
        self.ov_preset = tk.StringVar()
        self.sub_preset = tk.StringVar()
        for col, (text, var, choices, kind) in enumerate((
                ("Overlayer preset:", self.ov_preset, [k for k in self.overlayers if k != "None"], "ov"),
                ("Substrate preset:", self.sub_preset, list(self.substrates), "sub"))):
            if self.compact:     # both presets on one row to save height on laptops
                ttk.Label(presets, text=text.replace(" preset", ""), style="Card.TLabel").grid(
                    row=0, column=2 * col, sticky="w", pady=1)
                combo = ttk.Combobox(presets, textvariable=var, values=choices, width=14, state="readonly")
                combo.grid(row=0, column=2 * col + 1, padx=(4, 10), pady=1, sticky="w")
            else:
                ttk.Label(presets, text=text, style="Card.TLabel").grid(row=col, column=0, sticky="w", pady=1)
                combo = ttk.Combobox(presets, textvariable=var, values=choices, width=24, state="readonly")
                combo.grid(row=col, column=1, padx=6, pady=1, sticky="w")
            combo.bind("<<ComboboxSelected>>", lambda e, k=kind: self._apply_preset(k))
            Tooltip(combo, "Optional: fills the binding energy (and, for the overlayer, lambda) from the "
                           "Calculator's reference table. RSFs and peak areas are always typed in.")
        row += 1

        for key, label, hint in FIELDS:
            ttk.Label(card, text=label, style="Card.TLabel").grid(row=row, column=0, sticky="w", pady=1)
            entry = ttk.Entry(card, textvariable=self.vars[key], width=11)
            entry.grid(row=row, column=1, padx=6, pady=1)
            entry.bind("<Return>", lambda e: self.calculate())
            Tooltip(entry, hint)
            self.entries[key] = entry
            extra = ttk.Frame(card, style="Card.TFrame")
            extra.grid(row=row, column=2, sticky="w")
            if key == "Esource":
                for text, ev in (("Al Kα", "1486.6"), ("Mg Kα", "1253.6")):
                    ttk.Button(extra, text=text, width=6, command=lambda v=ev: self.vars["Esource"].set(v)).pack(
                        side="left", padx=(0, 3))
            elif key == "Theta":
                self.theta_note = ttk.Label(extra, text="", style="CardMuted.TLabel")
                self.theta_note.pack(side="left")
                self.vars["Theta"].trace_add("write", lambda *a: self._update_theta_note())
                self._update_theta_note()
            elif not self.compact:
                ttk.Label(extra, text=hint, style="CardMuted.TLabel", wraplength=250, justify="left").pack(side="left")
            row += 1

        # --- sample morphology (Topofactor, Shard, Wang & Spencer 2009) ---------------------------------------
        if not self.compact:
            ttk.Separator(card, orient="horizontal").grid(row=row, column=0, columnspan=3, sticky="ew", pady=(6, 4))
            row += 1
        morph = ttk.Frame(card, style="Card.TFrame")
        morph.grid(row=row, column=0, columnspan=3, sticky="ew")
        morph.columnconfigure(1, weight=1)
        ttk.Label(morph, text="Sample shape:", style="Card.TLabel").grid(row=0, column=0, sticky="w", pady=1)
        shape_combo = ttk.Combobox(morph, textvariable=self.shape_var, state="readonly", width=1,
                                   values=[s.label for s in TFAC.SHAPES.values()])
        shape_combo.grid(row=0, column=1, sticky="ew", padx=(6, 0), pady=1)
        shape_combo.bind("<<ComboboxSelected>>", lambda e: self._morphology_changed())
        Tooltip(shape_combo, "Correct the planar (flat-sample) thickness for a fibre or particle using a Topofactor "
                             "(Shard, Wang & Spencer, Surf. Interface Anal. 2009). Analyse with the sample's mean "
                             "normal facing the analyser (theta = 0).")
        ttk.Label(morph, text="Topofactor:", style="Card.TLabel").grid(row=1, column=0, sticky="w", pady=1)
        mode_combo = ttk.Combobox(morph, textvariable=self.mode_var, state="readonly", width=1,
                                  values=[self._mode_label(m) for m in TFAC.MODES])
        mode_combo.grid(row=1, column=1, sticky="ew", padx=(6, 0), pady=1)
        mode_combo.bind("<<ComboboxSelected>>", lambda e: self._morphology_changed())
        Tooltip(mode_combo, "Eq. (6) of the paper depends on the thickness and on B (accurate to ~1% against its "
                            "calculations); the constant average values (0.79 cylinder, 0.67 sphere, 0.71 shadowed "
                            "sphere) are within ~10% for typical samples.")
        row += 1
        if not self.compact:
            ttk.Label(card, textvariable=self.shape_hint, style="CardMuted.TLabel", wraplength=560,
                      justify="left").grid(row=row, column=0, columnspan=3, sticky="w", pady=(2, 0))
            row += 1
        self._update_shape_hint()

        buttons = ttk.Frame(card, style="Card.TFrame")
        buttons.grid(row=row, column=0, columnspan=3, sticky="ew", pady=(8, 0))
        example_btn = ttk.Button(buttons, text="Load example", command=self.load_example)
        example_btn.pack(side="left")
        copy_btn = ttk.Button(buttons, text="Copy from Calculator", command=self.copy_from_calculator)
        copy_btn.pack(side="left", padx=6)
        clear_btn = ttk.Button(buttons, text="Clear", command=self.clear)
        clear_btn.pack(side="left")
        Tooltip(example_btn, "Fill in demo numbers (Fe 2p overlayer on Si 2p) - not a real measurement")
        Tooltip(copy_btn, "Take the take-off angle (converted to the surface-normal convention), the binding "
                          "energies and lambda from the Calculator tab. Peak areas and RSFs are not copied.")
        return card

    # --- morphology helpers ------------------------------------------------------------------------
    def _mode_label(self, key):
        return SHORT_MODE[key] if self.compact else TFAC.MODE_LABELS[key]

    def _shape_key(self):
        for key, shape in TFAC.SHAPES.items():
            if shape.label == self.shape_var.get():
                return key
        return "flat"

    def _mode_key(self):
        for key in TFAC.MODES:
            if self._mode_label(key) == self.mode_var.get():
                return key
        return "eq6"

    def _update_shape_hint(self):
        key = self._shape_key()
        shape = TFAC.SHAPES[key]
        if key == "flat":
            self.shape_hint.set("Choose a fibre or particle to correct for its curvature (Topofactor). "
                                + TFAC.NANOSCALE_NOTE)
        else:
            self.shape_hint.set(f"{shape.hint} Analyse with the sample's mean normal facing the analyser (theta = 0). "
                                + TFAC.NANOSCALE_NOTE)

    def _morphology_changed(self):
        self._update_shape_hint()
        if self.result is not None:
            self._present()

    def _update_theta_note(self):
        th = _to_float(self.vars["Theta"].get())
        self.theta_note.configure(text=(f"= {T.surface_from_normal(th):g}° from the surface" if th is not None
                                        else ""))

    def _apply_preset(self, kind):
        comps, var = (self.overlayers, self.ov_preset) if kind == "ov" else (self.substrates, self.sub_preset)
        comp = comps.get(var.get())
        if comp is None:
            return
        be_key = "Eo" if kind == "ov" else "Es"
        if comp.be_ev:
            self.vars[be_key].set(f"{comp.be_ev:g}")
            note = f"{var.get()}: binding energy {comp.be_ev:g} eV ({comp.core_level})"
        else:
            self.vars[be_key].set("")
            note = f"'{var.get()}' is an Auger preset and has no binding energy - type the BE of your peak"
        if kind == "ov":
            self.vars["Lambda"].set(f"{comp.imfp_nm:g}")
            note += f"; lambda = {comp.imfp_nm:g} nm (an IMFP: attenuation lengths are ~10% shorter)"
        self._set_badge("idle", note)

    def load_example(self):
        for key, value in EXAMPLE.items():
            self.vars[key].set(value)
        self.ov_preset.set("")
        self.sub_preset.set("")
        self._set_badge("idle", "Demo numbers loaded (Fe 2p overlayer on Si 2p, Al Kα). Press Calculate & plot.")

    def clear(self):
        for key in self.vars:
            self.vars[key].set("")
        self.vars["Esource"].set("1486.6")
        self.vars["Theta"].set("0")
        self.ov_preset.set("")
        self.sub_preset.set("")
        self.shape_var.set(TFAC.SHAPES["flat"].label)
        self.mode_var.set(self._mode_label("eq6"))
        self._update_shape_hint()

    def copy_from_calculator(self):
        app, copied = self.app, []
        theta = _to_float(app.theta_var.get())
        if theta is not None:
            self.vars["Theta"].set(f"{T.normal_from_surface(theta):g}")
            copied.append(f"angle ({theta:g}° from the surface = {T.normal_from_surface(theta):g}° from the normal)")
        sub_be = _to_float(app.substrate_row.be_var.get())
        if sub_be:
            self.vars["Es"].set(f"{sub_be:g}")
            copied.append("substrate BE")
        if app.substrate_row.choice_var.get() in self.substrates:
            self.sub_preset.set(app.substrate_row.choice_var.get())
        for row in app.overlayer_rows:
            if row.active_var.get() and row.choice_var.get() != "None":
                ov_be, lam = _to_float(row.be_var.get()), _to_float(row.imfp_var.get())
                if ov_be:
                    self.vars["Eo"].set(f"{ov_be:g}")
                    copied.append("overlayer BE")
                if lam:
                    self.vars["Lambda"].set(f"{lam:g}")
                    copied.append("lambda")
                if row.choice_var.get() in self.overlayers:
                    self.ov_preset.set(row.choice_var.get())
                break
        self._set_badge("idle", "Copied from the Calculator: " + (", ".join(copied) if copied else "nothing to copy")
                        + ". Enter the peak areas and RSFs.")

    # --- results -------------------------------------------------------------------------------
    def _build_results(self, parent):
        card = ttk.LabelFrame(parent, text="  Calculation  ", style="Card.TLabelframe",
                              padding=(10, 2, 10, 4) if self.compact else 10)
        card.columnconfigure(0, weight=1)
        card.rowconfigure(3, weight=1)
        self.calc_btn = ttk.Button(card, text="Calculate & plot", style="Primary.TButton", command=self.calculate)
        self.calc_btn.grid(row=0, column=0, sticky="ew")
        Tooltip(self.calc_btn, "Solve for the thickness and mark it on the nomogram (Enter in any box does the same)")

        readout = ttk.Frame(card, style="Card.TFrame")
        readout.grid(row=1, column=0, sticky="ew", pady=(8, 0))
        self.readout_var = tk.StringVar(value="-.-- nm")
        self.readout_sub_var = tk.StringVar(value="Thickness will appear here")
        self.readout_lbl = tk.Label(readout, textvariable=self.readout_var, bg=CARD, fg="#b8bec7", anchor="w",
                                    font=("Segoe UI", 20 if self.compact else 30, "bold"))
        self.readout_lbl.pack(anchor="w")
        self.readout_sub = tk.Label(readout, textvariable=self.readout_sub_var, bg=CARD, fg=MUTED, anchor="w",
                                    font=("Segoe UI", 9), wraplength=340, justify="left")
        self.readout_sub.pack(anchor="w")
        self.readout_tip = Tooltip(self.readout_sub, "")          # the full details, shown on hover (see _present)

        self.badge = tk.Label(card, text="", font=("Segoe UI", 9, "bold"), padx=8, pady=4, anchor="w", justify="left",
                              wraplength=330)
        self.badge.grid(row=2, column=0, sticky="ew", pady=(8, 8))
        self._set_badge("idle", "No calculation yet")

        frame = ttk.Frame(card, style="Card.TFrame")
        frame.grid(row=3, column=0, sticky="nsew")
        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(0, weight=1)
        self.details = tk.Text(frame, height=6, width=46, wrap="word", state="disabled", relief="flat", bg="#f8fafc",
                               font=("Consolas", 8 if self.compact else 9), highlightthickness=1,
                               highlightbackground="#d1d5db", padx=6, pady=4)
        scroll = ttk.Scrollbar(frame, orient="vertical", command=self.details.yview)
        self.details.configure(yscrollcommand=scroll.set)
        self.details.grid(row=0, column=0, sticky="nsew")
        scroll.grid(row=0, column=1, sticky="ns")
        if self.compact:     # short screens: the readout + badge carry the result, the details are on hover
            frame.grid_remove()
        return card

    def _set_badge(self, kind, text):
        bg, fg = BADGE[kind]
        self.badge.configure(text=text, bg=bg, fg=fg)

    def calculate(self):
        try:
            values = {}
            for key, label, _ in FIELDS:
                try:
                    values[key] = float(self.vars[key].get())
                except ValueError:
                    raise ValueError(f"Please enter a number for '{label}'.")
            res = T.thickness(io=values["Io"], is_=values["Is"], so=values["So"], ss=values["Ss"],
                              be_o=values["Eo"], be_s=values["Es"], source_ev=values["Esource"],
                              lam_nm=values["Lambda"], theta_normal_deg=values["Theta"])
        except ValueError as e:
            messagebox.showerror("Input error", str(e))
            return
        self.result = res
        self._present()

    def _present(self):
        """Show the planar result, or the Topofactor-corrected one when a curved shape is selected."""
        res = self.result
        key, mode = self._shape_key(), self._mode_key()
        corr = TFAC.apply_topofactor(res, key, mode) if key != "flat" else None
        self.corrected = corr
        warnings = list(res.warnings) + (corr.warnings if corr else [])

        if corr is None:
            value_nm = res.t_nm
            self.readout_sub_var.set(f"{res.t_nm * 10:.1f} Å     (x = t/(L cosθ) = {res.x:.3f})")
            ok_text = "Within the printed nomogram scales"
            head = [f"==> THICKNESS t = {res.t_nm:.3f} nm  ({res.t_nm * 10:.2f} Angstrom)"]
        else:
            value_nm = corr.t_nm
            if self.compact:
                self.readout_sub_var.set(f"{corr.t_nm * 10:.1f} Å = {corr.tf:.3f} × planar {res.t_nm:.2f} nm")
            else:
                self.readout_sub_var.set(f"{corr.t_nm * 10:.1f} Å  corrected for a {SHORT_SHAPE[key]}: Topofactor "
                                         f"{corr.tf:.3f} × planar t = {res.t_nm:.2f} nm")
            ok_text = f"Topofactor applied ({SHORT_SHAPE[key]}, {corr.regime})"
            head = [f"==> CORRECTED THICKNESS T = {corr.t_nm:.3f} nm  ({corr.t_nm * 10:.2f} Angstrom)",
                    f"    = Topofactor {corr.tf:.4f} x planar (Thickogram) t = {res.t_nm:.3f} nm",
                    f"    shape: {TFAC.SHAPES[key].label}; {TFAC.MODE_LABELS[mode]}; regime: {corr.regime}",
                    f"    t/L = {corr.t_over_l:.3f}, B (paper) = (KEo/KEs)^0.75 = {corr.b_paper:.3f}"]
        self.readout_var.set(f"{value_nm:.2f} nm")
        self.readout_lbl.configure(fg=ACCENT)
        if warnings:
            self._set_badge("warn", " ".join(warnings))
        else:
            self._set_badge("ok", ok_text)

        lines = head + [
            "",
            f"Planar Thickogram: x = {res.x:.5f}   (t = x * L * cos(theta) = {res.x:.5f} * {res.lam_nm:g} nm * "
            f"cos({res.theta_normal_deg:g} deg))",
            f"A = Io/(Is*Ro) = {res.a:.5g}   with Ro = So/Ss = {res.ro:.5g}",
            f"B = KEo/KEs    = {res.b:.5g}   (KEo = {res.ke_o:.1f} eV, KEs = {res.ke_s:.1f} eV)",
            "",
            "On the nomogram the red line runs from A on the left scale to B; where it crosses the green curve "
            "is the answer (X = 2x).",
        ]
        if corr is not None:
            lines += [""] + ["Note: " + n for n in corr.notes]
        text = "\n".join(lines)
        self.details.configure(state="normal")
        self.details.delete("1.0", "end")
        self.details.insert("1.0", text)
        self.details.configure(state="disabled")
        self.readout_tip.text = text
        self._draw()

    # --- plot -----------------------------------------------------------------------------------
    def _build_plot(self):
        card = ttk.LabelFrame(self, text="  Thickogram (nomogram)  ", style="Card.TLabelframe", padding=(6, 4, 6, 4))
        self.fig = Figure(figsize=(8.5, 6.5), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=card)
        toolbar = NavigationToolbar2Tk(self.canvas, card, pack_toolbar=False)
        toolbar.set_message = lambda s: None          # a coordinate read-out means nothing on a nomogram
        for child in [toolbar] + list(toolbar.winfo_children()):
            try:
                child.configure(bg=CARD)
            except tk.TclError:
                pass
        toolbar.pack(side="bottom", fill="x")
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        return card

    def _draw(self):
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        annotation = None
        if self.corrected is not None:
            c = self.corrected
            annotation = (f"planar t = {c.planar_t_nm:.2f} nm\nTopofactor {c.tf:.3f} "
                          f"({SHORT_SHAPE[c.shape]}, {c.regime})\nT = {c.t_nm:.2f} nm")
        T.draw_nomogram(ax, self.result, annotation=annotation)
        self.fig.subplots_adjust(left=0.03, right=0.98, top=0.94, bottom=0.07)
        self.canvas.draw_idle()
