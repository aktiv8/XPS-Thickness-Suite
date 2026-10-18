"""
The "ARXPS" tab of the XPS Thickness Suite: overlayer thickness from a tilt series (VAMAS), from the intensity ratio
against 1/cos(theta), with the azimuth read from the file when it is there.

Maths in arxps.py, file reading in vamas.py; this module is only the Tk/matplotlib glue. It is a ttk.Frame that app.py
puts into a notebook tab, and needs numpy + matplotlib (app.py shows an install hint instead if they are missing).

Angle convention: theta is from the surface NORMAL (0 = normal emission) and equals the stage tilt plus an offset.
"""
import math
import os
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

import arxps as AR
import vamas as VM
from gui_common import ACCENT, BADGE, CARD, MUTED, Tooltip, _to_float

SUB_COLOUR, OV_COLOUR = "#2563eb", "#d97706"
NO_AZIMUTH = "None (single azimuth)"
MODE_STATES, MODE_PEAKS = 1, 2


def _az_label(az):
    return "-" if az is None else f"{az:g}"


def _az_sort(keys):
    return sorted(keys, key=lambda a: (a is not None, a if a is not None else 0.0))


class ARXPSTab(ttk.Frame):
    def __init__(self, parent, app, substrates, overlayers, volume_density, compact=False):
        super().__init__(parent)
        self.app, self.substrates, self.overlayers = app, substrates, overlayers
        self.volume_density, self.compact = volume_density, compact
        self.file = None
        self.results = {}                   # azimuth -> ARXPSResult
        self.summary = None                 # (mean thickness, std between azimuths or None)
        self.mode = tk.IntVar(value=MODE_STATES)
        self.per_second = tk.BooleanVar(value=True)
        self.transmission = tk.BooleanVar(value=True)
        keys = ("region", "region_sub", "region_ov", "angle_src", "azim_src", "view_az", "offset", "max_angle", "lam",
                "rinf", "so", "ss", "eo", "es", "hv", "env_lo", "env_hi", "split", "s_lo", "s_hi", "o_lo", "o_hi")
        self.v = {k: tk.StringVar() for k in keys}
        self.v["offset"].set("0")
        self.v["max_angle"].set(f"{AR.DEFAULT_MAX_ANGLE:g}")
        self.v["so"].set("1")
        self.v["ss"].set("1")
        self.sub_preset, self.ov_preset = tk.StringVar(), tk.StringVar()
        self.file_label = tk.StringVar(value="No file loaded")

        self.columnconfigure(0, weight=0)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)
        left = ttk.Frame(self)
        left.grid(row=0, column=0, sticky="ns", padx=(0, 8), pady=(8, 4))
        left.rowconfigure(1, weight=1)
        self._build_setup(left).grid(row=0, column=0, sticky="ew")
        self._build_results(left).grid(row=1, column=0, sticky="nsew", pady=(8, 0))
        self._build_right().grid(row=0, column=1, sticky="nsew", pady=(8, 4))
        self._mode_changed()
        self._draw()

    # --- setup pages ------------------------------------------------------------------------------
    def _entry(self, parent, key, row, label, hint, col=0, width=9, span=1):
        ttk.Label(parent, text=label, style="Card.TLabel").grid(row=row, column=col, sticky="w", pady=2)
        e = ttk.Entry(parent, textvariable=self.v[key], width=width)
        e.grid(row=row, column=col + 1, columnspan=span, sticky="w", padx=6, pady=2)
        Tooltip(e, hint)
        return e

    def _combo(self, parent, var, row, label, hint, values=(), col=0, width=26):
        lbl = ttk.Label(parent, text=label, style="Card.TLabel")
        lbl.grid(row=row, column=col, sticky="w", pady=2)
        c = ttk.Combobox(parent, textvariable=var, values=list(values), state="readonly", width=width)
        c.grid(row=row, column=col + 1, columnspan=3, sticky="w", padx=6, pady=2)
        Tooltip(c, hint)
        return lbl, c

    def _build_setup(self, parent):
        card = ttk.LabelFrame(parent, text="  Setup  ", style="Card.TLabelframe", padding=(6, 2, 6, 6))
        nb = ttk.Notebook(card, style="Sub.TNotebook")
        nb.pack(fill="both", expand=True)
        data = ttk.Frame(nb, style="Card.TFrame", padding=(8, 6))
        mats = ttk.Frame(nb, style="Card.TFrame", padding=(8, 6))
        wins = ttk.Frame(nb, style="Card.TFrame", padding=(8, 6))
        nb.add(data, text="1  Data")
        nb.add(mats, text="2  Materials")
        nb.add(wins, text="3  Windows")
        self.setup_nb = nb
        self._build_data_page(data)
        self._build_materials_page(mats)
        self._build_windows_page(wins)
        return card

    def _build_data_page(self, f):
        top = ttk.Frame(f, style="Card.TFrame")
        top.grid(row=0, column=0, columnspan=4, sticky="ew")
        btn = ttk.Button(top, text="Open VAMAS file…", command=self.open_file)
        btn.pack(side="left")
        Tooltip(btn, "Load a VAMAS (.vms) file with a tilt series: one block per region and tilt. Thermo Avantage "
                     "(K-Alpha) exports are supported; the stage tilt is read from experimental variable 1.")
        ttk.Label(top, textvariable=self.file_label, style="CardMuted.TLabel", wraplength=230,
                  justify="left").pack(side="left", padx=8)

        ttk.Radiobutton(f, text="Two states of one peak (Si 2p: Si + SiO₂)", value=MODE_STATES, variable=self.mode,
                        style="Card.TRadiobutton", command=self._mode_changed).grid(
            row=1, column=0, columnspan=4, sticky="w", pady=(6, 0))
        ttk.Radiobutton(f, text="Two different peaks (e.g. O 1s over Si 2p)", value=MODE_PEAKS, variable=self.mode,
                        style="Card.TRadiobutton", command=self._mode_changed).grid(
            row=2, column=0, columnspan=4, sticky="w")

        self.region_lbl, self.region_cb = self._combo(f, self.v["region"], 3, "Region:",
                                                      "The block (region) holding both chemical states.")
        self.region_cb.bind("<<ComboboxSelected>>", lambda e: self._region_changed())
        self.rsub_lbl, self.rsub_cb = self._combo(f, self.v["region_sub"], 4, "Substrate peak:",
                                                  "The region holding the substrate peak.")
        self.rov_lbl, self.rov_cb = self._combo(f, self.v["region_ov"], 5, "Overlayer peak:",
                                                "The region holding the overlayer peak.")
        for cb in (self.rsub_cb, self.rov_cb):
            cb.bind("<<ComboboxSelected>>", lambda e: self._region_changed())

        ttk.Separator(f, orient="horizontal").grid(row=6, column=0, columnspan=4, sticky="ew", pady=5)
        _, self.angle_cb = self._combo(f, self.v["angle_src"], 7, "Angle from:",
                                       "Where the tilt comes from. The K-Alpha writes the stage tilt into experimental "
                                       "variable 1, which the file only calls 'Exp Variable 1'.", width=24)
        self._entry(f, "offset", 8, "Offset (°):", "Added to the tilt to give the emission angle from the surface "
                    "normal. 0 for a stage tilt measured from normal emission (K-Alpha).")
        _, self.azim_cb = self._combo(f, self.v["azim_src"], 9, "Azimuth from:",
                                      "Optional: the field holding the sample azimuth. If the file has several azimuths "
                                      "each is fitted separately and the thickness is shown against azimuth.", width=24)
        self.view_lbl, self.view_cb = self._combo(f, self.v["view_az"], 10, "Show azimuth:",
                                                  "Which azimuth the Spectra plot and the table show.", width=10)
        self.view_cb.bind("<<ComboboxSelected>>", lambda e: self._view_changed())
        self._entry(f, "max_angle", 11, "Max angle (°):", "Angles above this are drawn but left out of the fit: the "
                    "simple exponential model fails at grazing emission (elastic scattering, roughness).")
        opts = ttk.Frame(f, style="Card.TFrame")
        opts.grid(row=12, column=0, columnspan=4, sticky="w", pady=(4, 0))
        c1 = ttk.Checkbutton(opts, text="Counts per second", variable=self.per_second, style="Card.TCheckbutton")
        c1.pack(side="left")
        c2 = ttk.Checkbutton(opts, text="Transmission-correct", variable=self.transmission, style="Card.TCheckbutton")
        c2.pack(side="left", padx=(10, 0))
        Tooltip(c1, "Divide counts by scans × dwell. The K-Alpha stores raw counts summed over the scans, so this is "
                    "needed when the two peaks have different scan counts. It cancels for two states of one peak.")
        Tooltip(c2, "Divide by the Transmission column of the file. Needed when the two peaks are at different "
                    "kinetic energies; it cancels (to ~1%) for two states of one peak.")
        self.angle_cb.bind("<<ComboboxSelected>>", lambda e: self._clear_windows())

    def _build_materials_page(self, f):
        self.sub_lbl, self.sub_cb = self._combo(f, self.sub_preset, 0, "Substrate:", "Optional preset (from the "
                                                "Calculator's table): fills the binding energy, and Rinf with the "
                                                "overlayer.", list(self.substrates), width=24)
        self.ov_lbl, self.ov_cb = self._combo(f, self.ov_preset, 1, "Overlayer:", "Optional preset: fills the binding "
                                              "energy and lambda, and Rinf with the substrate.",
                                              [k for k in self.overlayers if k != "None"], width=24)
        self.sub_cb.bind("<<ComboboxSelected>>", lambda e: self._apply_preset("sub"))
        self.ov_cb.bind("<<ComboboxSelected>>", lambda e: self._apply_preset("ov"))
        self._entry(f, "lam", 2, "λ overlayer (nm):", "Attenuation length of the overlayer photoelectrons in the "
                    "overlayer. The thickness is slope × lambda. The presets hold IMFPs (attenuation lengths are ~10% "
                    "shorter).", width=9, span=2)

        self.states_frame = ttk.Frame(f, style="Card.TFrame")
        self.states_frame.grid(row=3, column=0, columnspan=4, sticky="w")
        self._entry(self.states_frame, "rinf", 0, "R∞:", "R∞ = (RSF·N·λ)overlayer / (RSF·N·λ)substrate: the "
                    "overlayer/substrate area ratio of infinitely thick layers. Filled from the two presets, or type "
                    "it (e.g. from a reference sample).", width=9)
        ttk.Label(self.states_frame, text="from the presets, or type it", style="CardMuted.TLabel").grid(
            row=0, column=2, sticky="w")

        self.peaks_frame = ttk.Frame(f, style="Card.TFrame")
        self.peaks_frame.grid(row=4, column=0, columnspan=4, sticky="w")
        pf = self.peaks_frame
        self._entry(pf, "so", 0, "RSF overlayer:", "Relative sensitivity factor of the overlayer peak (Ro = So/Ss).")
        self._entry(pf, "ss", 1, "RSF substrate:", "Relative sensitivity factor of the substrate peak.")
        self._entry(pf, "eo", 2, "BE overlayer (eV):", "Binding energy of the overlayer peak (sets B = KEo/KEs).")
        self._entry(pf, "es", 3, "BE substrate (eV):", "Binding energy of the substrate peak.")
        self._entry(pf, "hv", 4, "Source (eV):", "X-ray photon energy; read from the file.")
        ttk.Label(f, text="Presets only fill in numbers: change anything.", style="CardMuted.TLabel").grid(
            row=5, column=0, columnspan=4, sticky="w", pady=(6, 0))

    def _build_windows_page(self, f):
        self.win_states = ttk.Frame(f, style="Card.TFrame")
        self.win_states.grid(row=0, column=0, sticky="w")
        w = self.win_states
        self._entry(w, "env_lo", 0, "Envelope from (eV):", "Low-BE end of the region integrated for both states "
                    "(one Shirley background across the whole envelope).")
        self._entry(w, "env_hi", 1, "Envelope to (eV):", "High-BE end of the envelope.")
        self._entry(w, "split", 2, "Split at (eV):", "BE between the two states: the area below it is the substrate "
                    "state, above it the overlayer state. Auto-detect puts it at the valley.")

        self.win_peaks = ttk.Frame(f, style="Card.TFrame")
        self.win_peaks.grid(row=0, column=0, sticky="w")
        p = self.win_peaks
        for r, (tag, name) in enumerate((("s", "Substrate"), ("o", "Overlayer"))):
            ttk.Label(p, text=f"{name} peak", style="RowHead.TLabel").grid(row=3 * r, column=0, columnspan=2, sticky="w",
                                                                           pady=(4 if r else 0, 0))
            self._entry(p, f"{tag}_lo", 3 * r + 1, "  window from (eV):", "Low-BE end of the window (Shirley background "
                        "between the window ends).")
            self._entry(p, f"{tag}_hi", 3 * r + 2, "  window to (eV):", "High-BE end of the window.")

        btns = ttk.Frame(f, style="Card.TFrame")
        btns.grid(row=1, column=0, sticky="w", pady=(8, 0))
        auto = ttk.Button(btns, text="Auto-detect windows", command=self.auto_windows)
        auto.pack(side="left")
        Tooltip(auto, "Find the peak(s) near the preset (or typed) binding energies in the loaded data and set the "
                      "windows. Always check them against the Spectra plot.")
        self.win_note = tk.Label(f, text="", bg=CARD, fg=MUTED, font=("Segoe UI", 8), justify="left", wraplength=330,
                                 anchor="w")
        self.win_note.grid(row=2, column=0, sticky="w", pady=(6, 0))

    # --- mode / selections --------------------------------------------------------------------------
    def _mode_changed(self):
        states = self.mode.get() == MODE_STATES
        for w in (self.region_lbl, self.region_cb):
            (w.grid if states else w.grid_remove)()
        for w in (self.rsub_lbl, self.rsub_cb, self.rov_lbl, self.rov_cb):
            (w.grid_remove if states else w.grid)()
        (self.states_frame.grid if states else self.states_frame.grid_remove)()
        (self.peaks_frame.grid_remove if states else self.peaks_frame.grid)()
        (self.win_states.grid if states else self.win_states.grid_remove)()
        (self.win_peaks.grid_remove if states else self.win_peaks.grid)()
        self._clear_windows()
        self._fill_regions()

    def _fill_regions(self):
        names = self.file.regions() if self.file else []
        for cb in (self.region_cb, self.rsub_cb, self.rov_cb):
            cb.configure(values=names)
        if not names:
            return
        spectra = [n for n in names if "survey" not in n.lower()] or names
        if self.mode.get() == MODE_STATES:
            if self.v["region"].get() not in names:
                self.v["region"].set(spectra[0])
        else:
            if self.v["region_sub"].get() not in names:
                self.v["region_sub"].set(spectra[0])
            if self.v["region_ov"].get() not in names:
                self.v["region_ov"].set(spectra[min(1, len(spectra) - 1)])
        self._region_changed()

    def _regions_in_use(self):
        if self.mode.get() == MODE_STATES:
            return [self.v["region"].get()]
        return [self.v["region_sub"].get(), self.v["region_ov"].get()]

    def _region_changed(self):
        if not self.file:
            return
        regs = [r for r in self._regions_in_use() if r]
        if not regs:
            return
        present = [s for s in VM.ANGLE_SOURCES if all(self.file.sources_present(r).count(s) for r in regs)]
        labels = [VM.ANGLE_SOURCES[s] for s in present]
        self.angle_cb.configure(values=labels or list(VM.ANGLE_SOURCES.values()))
        self.azim_cb.configure(values=[NO_AZIMUTH] + labels)
        angle, azim = self.file.default_sources(regs[0])
        self.v["angle_src"].set(VM.ANGLE_SOURCES[angle] if angle else "")
        self.v["azim_src"].set(VM.ANGLE_SOURCES[azim] if azim and azim in present else NO_AZIMUTH)
        blocks = self.file.blocks_of(regs[0])
        if blocks and blocks[0].hv_ev and not self.v["hv"].get():
            self.v["hv"].set(f"{blocks[0].hv_ev:g}")
        self._clear_windows()
        self._set_badge("idle", f"{self.file_label.get()}: choose the materials, then Calculate & plot.")

    def _source_key(self, label):
        for key, text in VM.ANGLE_SOURCES.items():
            if text == label:
                return key
        return None

    def _clear_windows(self):
        for k in ("env_lo", "env_hi", "split", "s_lo", "s_hi", "o_lo", "o_hi"):
            self.v[k].set("")
        if hasattr(self, "win_note"):
            self.win_note.configure(text="")

    def _apply_preset(self, kind):
        comps, var = (self.overlayers, self.ov_preset) if kind == "ov" else (self.substrates, self.sub_preset)
        comp = comps.get(var.get())
        if comp is None:
            return
        if comp.be_ev:
            self.v["eo" if kind == "ov" else "es"].set(f"{comp.be_ev:g}")
        if kind == "ov":
            self.v["lam"].set(f"{comp.imfp_nm:g}")
            if comp.rsf:
                self.v["so"].set(f"{comp.rsf:g}")
        elif comp.rsf:
            self.v["ss"].set(f"{comp.rsf:g}")
        sub, ov = self.substrates.get(self.sub_preset.get()), self.overlayers.get(self.ov_preset.get())
        if sub is not None and ov is not None:
            try:
                n_s = self.volume_density(sub.density_g_cm3, sub.molar_mass_g_mol, sub.atoms_per_formula)
                n_o = self.volume_density(ov.density_g_cm3, ov.molar_mass_g_mol, ov.atoms_per_formula)
                rinf = (ov.rsf * n_o * ov.imfp_nm) / (sub.rsf * n_s * sub.imfp_nm)
                self.v["rinf"].set(f"{rinf:.4g}")
            except (ValueError, ZeroDivisionError):
                pass
        self._clear_windows()
        self._set_badge("idle", f"{var.get()}: preset values filled in (edit anything). Check the windows, then calculate.")

    # --- file ---------------------------------------------------------------------------------------
    def open_file(self):
        path = filedialog.askopenfilename(title="Open a VAMAS tilt series",
                                          filetypes=[("VAMAS files", "*.vms *.vamas"), ("All files", "*.*")])
        if path:
            self.load_file(path)

    def load_file(self, path):
        try:
            vf = VM.read(path)
        except VM.VamasError as exc:
            messagebox.showerror("Cannot read the file", str(exc))
            return False
        if not vf.blocks:
            messagebox.showerror("Cannot read the file", "The file holds no spectra.")
            return False
        self.file = vf
        self.results, self.summary = {}, None
        self.v["hv"].set("")
        self.file_label.set(f"{os.path.basename(path)}: {len(vf.blocks)} blocks, regions {', '.join(vf.regions())}")
        self._fill_regions()
        self._draw()
        return True

    # --- windows ------------------------------------------------------------------------------------
    def _reference_series(self, region):
        """The first azimuth's blocks of a region as [(angle, block), ...], lowest angle first (None if no angle)."""
        angle_src = self._source_key(self.v["angle_src"].get())
        if angle_src is None:
            raise ValueError("Choose where the angle comes from (Data page) first.")
        series = self.file.series(region, angle_src, self._source_key(self.v["azim_src"].get()))
        if not series:
            raise ValueError(f"No block of '{region}' has an angle in the chosen field.")
        return series[_az_sort(series)[0]]

    def _spectra_of(self, region):
        ps, tr = self.per_second.get(), self.transmission.get()
        return [(b.be_ev, AR.intensity(b, ps, tr)) for _, b in self._reference_series(region)]

    def auto_windows(self):
        if not self.file:
            messagebox.showinfo("No data", "Open a VAMAS file first.")
            return
        try:
            if self.mode.get() == MODE_STATES:
                sub = self.substrates.get(self.sub_preset.get())
                ov = self.overlayers.get(self.ov_preset.get())
                if sub is None or ov is None or not sub.be_ev or not ov.be_ev:
                    messagebox.showinfo("Need the materials", "Choose the substrate and overlayer presets on the Materials "
                                        "page (they give the expected binding energies) so the two states can be found.")
                    return
                # The preset BEs are only a starting point: charging moves the measured peaks by a few tenths of an eV or more.
                w = AR.auto_two_state_windows(self._spectra_of(self.v["region"].get()), sub.be_ev, ov.be_ev)
                for key, val in (("env_lo", w.env_lo), ("env_hi", w.env_hi), ("split", w.split)):
                    self.v[key].set(f"{val:.2f}")
                self.win_note.configure(text=f"Substrate state found at {w.sub_pos:.2f} eV and overlayer state at "
                                             f"{w.ov_pos:.2f} eV (as they change with angle); split at {w.split:.2f} eV. "
                                             "Check the shaded areas on the Spectra plot.")
            else:
                for tag, region, be_key in (("s", self.v["region_sub"].get(), "es"), ("o", self.v["region_ov"].get(), "eo")):
                    nominal = _to_float(self.v[be_key].get())
                    if not nominal:
                        messagebox.showinfo("Need the binding energies", "Enter the binding energies of both peaks on "
                                            "the Materials page so the peaks can be found.")
                        return
                    w = AR.auto_peak_window(self._spectra_of(region), nominal)
                    self.v[f"{tag}_lo"].set(f"{w.lo:.2f}")
                    self.v[f"{tag}_hi"].set(f"{w.hi:.2f}")
                self.win_note.configure(text="Peaks found near the binding energies you entered. Check the shaded areas "
                                             "on the Spectra plot after calculating.")
        except ValueError as exc:
            messagebox.showerror("Auto-detect failed", str(exc))

    # --- results card -------------------------------------------------------------------------------
    def _build_results(self, parent):
        card = ttk.LabelFrame(parent, text="  Calculation  ", style="Card.TLabelframe",
                              padding=(10, 2, 10, 4) if self.compact else 10)
        card.columnconfigure(0, weight=1)
        card.rowconfigure(3, weight=1)
        self.calc_btn = ttk.Button(card, text="Calculate & plot", style="Primary.TButton", command=self.calculate)
        self.calc_btn.grid(row=0, column=0, sticky="ew")
        Tooltip(self.calc_btn, "Integrate the peaks at every angle, solve for the thickness and plot x against 1/cos θ")
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
        self.readout_tip = Tooltip(self.readout_sub, "")
        self.badge = tk.Label(card, text="", font=("Segoe UI", 9, "bold"), padx=8, pady=4, anchor="w", justify="left",
                              wraplength=330)
        self.badge.grid(row=2, column=0, sticky="ew", pady=(8, 8))
        self._set_badge("idle", "Open a VAMAS tilt series to start")
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
        if self.compact:
            frame.grid_remove()
        return card

    def _set_badge(self, kind, text):
        bg, fg = BADGE[kind]
        self.badge.configure(text=text, bg=bg, fg=fg)

    # --- calculation --------------------------------------------------------------------------------
    def _num(self, key, label, positive=False):
        x = _to_float(self.v[key].get())
        if x is None or not math.isfinite(x):
            raise ValueError(f"Please enter a number for '{label}'.")
        if positive and x <= 0:
            raise ValueError(f"'{label}' must be greater than zero.")
        return x

    def calculate(self):
        try:
            self._calculate()
        except (ValueError, VM.VamasError) as exc:
            messagebox.showerror("Input error", str(exc))

    def _calculate(self):
        if not self.file:
            raise ValueError("Open a VAMAS file first.")
        angle_src = self._source_key(self.v["angle_src"].get())
        if angle_src is None:
            raise ValueError("Choose where the angle comes from (Data page): this file has no usable angle field.")
        az_src = self._source_key(self.v["azim_src"].get())
        offset = self._num("offset", "Offset")
        max_angle = self._num("max_angle", "Max angle", positive=True)
        lam = self._num("lam", "λ overlayer", positive=True)
        ps, tr = self.per_second.get(), self.transmission.get()
        results = {}
        if self.mode.get() == MODE_STATES:
            region = self.v["region"].get()
            rinf = self._num("rinf", "R∞", positive=True)
            if not self.v["env_lo"].get():
                self.auto_windows()
            w = AR.TwoStateWindows(self._num("env_lo", "Envelope from"), self._num("env_hi", "Envelope to"),
                                   self._num("split", "Split at"))
            series = self.file.series(region, angle_src, az_src)
            if not series:
                raise ValueError("No block of that region has an angle in the chosen field.")
            first = self._reference_series(region)[0][1]          # every azimuth is aligned to this spectrum
            reference = (first.be_ev, AR.intensity(first, ps, tr))
            for az, rows in series.items():
                results[az] = AR.analyse_two_states(rows, w, lam, rinf, offset, max_angle, ps, tr, reference=reference)
        else:
            reg_s, reg_o = self.v["region_sub"].get(), self.v["region_ov"].get()
            if reg_s == reg_o:
                raise ValueError("Choose two different regions for the substrate and overlayer peaks.")
            if not self.v["s_lo"].get():
                self.auto_windows()
            sw = AR.Window(self._num("s_lo", "Substrate window from"), self._num("s_hi", "Substrate window to"))
            ow = AR.Window(self._num("o_lo", "Overlayer window from"), self._num("o_hi", "Overlayer window to"))
            ser_s, ser_o = self.file.series(reg_s, angle_src, az_src), self.file.series(reg_o, angle_src, az_src)
            if not ser_s or not ser_o:
                raise ValueError("No block of those regions has an angle in the chosen field.")
            refs = []
            for region in (reg_s, reg_o):
                first = self._reference_series(region)[0][1]
                refs.append((first.be_ev, AR.intensity(first, ps, tr)))
            for az in ser_s:
                if az in ser_o:
                    results[az] = AR.analyse_two_peaks(
                        ser_s[az], ser_o[az], sw, ow, self._num("so", "RSF overlayer", True),
                        self._num("ss", "RSF substrate", True), self._num("eo", "BE overlayer"),
                        self._num("es", "BE substrate"), self._num("hv", "Source", True), lam, offset, max_angle, ps, tr,
                        sub_reference=refs[0], ov_reference=refs[1])
            if not results:
                raise ValueError("The two regions have no azimuth in common.")
        for az, res in results.items():
            res.azimuth = az
        self.results = results
        self.summary = AR.summarise_azimuths(results)
        keys = _az_sort(results)
        self.v["view_az"].set(_az_label(keys[0]))
        self.view_cb.configure(values=[_az_label(k) for k in keys])
        for w in (self.view_lbl, self.view_cb):
            (w.grid if len(keys) > 1 else w.grid_remove)()
        self._present()

    def _view_result(self):
        for az, res in self.results.items():
            if _az_label(az) == self.v["view_az"].get():
                return res
        return next(iter(self.results.values())) if self.results else None

    def _view_changed(self):
        self._draw()
        self._fill_table()

    def _present(self):
        res = self._view_result()
        mean, spread = self.summary
        multi = len(self.results) > 1
        warnings = []
        for az, r in self.results.items():
            warnings += [(f"[azimuth {az:g}°] " if multi and az is not None else "") + w for w in r.warnings]
        self.readout_var.set(f"{mean:.2f} nm")
        self.readout_lbl.configure(fg=ACCENT)
        if multi:
            sub = f"{mean * 10:.1f} Å: mean of {len(self.results)} azimuths" + (f", spread {spread:.2f} nm" if spread else "")
        else:
            se = f" ± {res.d_se_nm:.2f} nm" if res.d_se_nm is not None else ""
            sub = f"{mean * 10:.1f} Å{se} (slope fit through the origin, {len(res.used)} angles)"
        self.readout_sub_var.set(sub)
        if warnings:
            self._set_badge("warn", " ".join(warnings))
        else:
            self._set_badge("ok", "The angles agree: x falls on a line through the origin")

        mode_txt = ("Two states of one peak: x = ln(1 + R/R∞)" if res.mode == 1
                    else "Two different peaks: A = (1 - e^-x) e^(Bx), solved for x")
        lines = [f"==> THICKNESS d = {mean:.3f} nm  ({mean * 10:.2f} Angstrom)"
                 + (f"   (mean over azimuths, spread {spread:.3f} nm)" if multi and spread else ""), "", mode_txt]
        for az in _az_sort(self.results):
            r = self.results[az]
            head = f"azimuth {az:g} deg: " if az is not None else ""
            lines.append(f"{head}slope k = d/lambda = {r.k:.4f}" + (f" ± {r.k_se:.4f}" if r.k_se is not None else "")
                         + f"  -> d = {r.d_nm:.3f} nm (lambda = {r.lam_nm:g} nm)")
            if r.d_mean_nm is not None:
                lines.append(f"  per-angle mean {r.d_mean_nm:.3f} nm"
                             + (f", spread {r.d_spread_nm:.3f} nm" if r.d_spread_nm is not None else ""))
            if r.line is not None:
                lines.append(f"  free intercept fit: x = {r.line[0]:+.3f} + {r.line[1]:.4f}/cos θ  (d = {r.d_free_nm:.3f} nm)")
            if r.rinf_free is not None:
                if r.rinf_free["edge"]:
                    lines.append("  R∞-free fit: not determined by these angles (the ratio is too close to linear in 1/cos θ)")
                else:
                    lines.append(f"  R∞-free fit: R∞ = {r.rinf_free['rinf']:.3g} (preset {r.rinf:.3g}), "
                                 f"d = {r.d_rinf_free_nm:.3f} nm")
            if r.rinf is not None:
                lines.append(f"  R∞ used = {r.rinf:.4g}")
            if r.b is not None:
                lines.append(f"  B = KEo/KEs = {r.b:.4g}")
        lines += ["", "Angles from the surface normal; windows follow the peak at each angle (Shirley background)."]
        text = "\n".join(lines)
        self.details.configure(state="normal")
        self.details.delete("1.0", "end")
        self.details.insert("1.0", text)
        self.details.configure(state="disabled")
        self.readout_tip.text = text
        self._draw()
        self._fill_table()

    # --- right-hand side ----------------------------------------------------------------------------
    def _build_right(self):
        frame = ttk.Frame(self)
        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(0, weight=1)
        nb = ttk.Notebook(frame, style="Sub.TNotebook")
        nb.grid(row=0, column=0, sticky="nsew")
        self.plot_nb = nb
        self.figs = {}
        for key, title in (("spectra", "Spectra"), ("fit", "x vs 1/cos θ"), ("azimuth", "Azimuth")):
            page = ttk.Frame(nb, style="Card.TFrame", padding=4)
            nb.add(page, text=title)
            fig = Figure(figsize=(8, 5), dpi=100)
            canvas = FigureCanvasTkAgg(fig, master=page)
            toolbar = NavigationToolbar2Tk(canvas, page, pack_toolbar=False)
            toolbar.configure(width=250)
            for child in [toolbar] + list(toolbar.winfo_children()):
                try:
                    child.configure(bg=CARD)
                except tk.TclError:
                    pass
            toolbar.pack(side="bottom", fill="x")
            canvas.get_tk_widget().pack(fill="both", expand=True)
            self.figs[key] = (fig, canvas)
        cols = ("az", "tilt", "theta", "isub", "iov", "ratio", "x", "d", "shift", "used")
        heads = ("Azimuth", "Tilt", "θ", "I substrate", "I overlayer", "R or A", "x", "d (nm)", "Shift (eV)", "Used")
        widths = (60, 50, 50, 90, 90, 70, 60, 70, 70, 150)
        self.table = ttk.Treeview(frame, columns=cols, show="headings", height=5 if self.compact else 7)
        for c, h, w in zip(cols, heads, widths):
            self.table.heading(c, text=h)
            self.table.column(c, width=w, anchor="e" if c != "used" else "w", stretch=(c == "used"))
        self.table.grid(row=1, column=0, sticky="ew", pady=(6, 0))
        return frame

    def _fill_table(self):
        self.table.delete(*self.table.get_children())
        res = self._view_result()
        if res is None:
            return
        for r in res.rows:
            self.table.insert("", "end", values=(
                _az_label(res.azimuth), f"{r.tilt:g}", f"{r.angle_deg:g}", f"{r.i_sub:.4g}", f"{r.i_ov:.4g}",
                "-" if r.ratio is None else f"{r.ratio:.4f}", "-" if r.x is None else f"{r.x:.4f}",
                "-" if r.d_nm is None else f"{r.d_nm:.3f}", f"{r.shift_ev:+.2f}", "yes" if r.used else r.note))

    # --- plots --------------------------------------------------------------------------------------
    def _placeholder(self, fig, text):
        fig.clear()
        ax = fig.add_subplot(111)
        ax.axis("off")
        ax.text(0.5, 0.5, text, ha="center", va="center", color=MUTED, fontsize=11, wrap=True)

    def _draw(self):
        res = self._view_result() if self.results else None
        for key, (fig, canvas) in self.figs.items():
            if res is None:
                self._placeholder(fig, "Open a VAMAS tilt series, set the materials and windows,\nthen press "
                                       "Calculate & plot.")
            elif key == "spectra":
                self._draw_spectra(fig, res)
            elif key == "fit":
                self._draw_fit(fig)
            else:
                self._draw_azimuth(fig)
            canvas.draw_idle()

    def _angle_colours(self, res):
        from matplotlib import cm
        angles = [r.angle_deg for r in res.rows]
        lo, hi = min(angles), max(angles)
        return {id(r): cm.viridis(0.0 if hi == lo else 0.9 * (r.angle_deg - lo) / (hi - lo)) for r in res.rows}

    def _draw_spectra(self, fig, res):
        fig.clear()
        colours = self._angle_colours(res)
        panels = ([("sub", "Both states, normalised to the substrate-state area")] if res.mode == 1
                  else [("sub", "Substrate peak"), ("ov", "Overlayer peak")])
        axes = fig.subplots(1, len(panels), squeeze=False)[0]
        for ax, (which, title) in zip(axes, panels):
            first = None
            for r in res.rows:
                ex = r.sub if which == "sub" else r.ov
                if ex is None or not ex.area_sub > 0:
                    continue
                y = ex.y / ex.area_sub
                m = (ex.grid >= ex.lo - 1.5) & (ex.grid <= ex.hi + 1.5)
                style = "-" if r.used else ":"
                ax.plot(ex.grid[m], y[m], style, color=colours[id(r)], lw=1.2, label=f"{r.tilt:g}°")
                if first is None and r.used:
                    first = (ex, y)
            if first is not None:
                ex, y = first
                inside = (ex.grid >= ex.lo) & (ex.grid <= ex.hi)
                bg = ex.bg / ex.area_sub
                ax.plot(ex.grid[inside], bg[inside], "--", color="#6b7280", lw=1, label="Shirley")
                if ex.split is not None:
                    left, right = inside & (ex.grid <= ex.split), inside & (ex.grid >= ex.split)
                    ax.fill_between(ex.grid[left], bg[left], y[left], color=SUB_COLOUR, alpha=0.25, lw=0)
                    ax.fill_between(ex.grid[right], bg[right], y[right], color=OV_COLOUR, alpha=0.25, lw=0)
                    ax.axvline(ex.split, color="#6b7280", lw=0.8, ls=":")
                else:
                    ax.fill_between(ex.grid[inside], bg[inside], y[inside], color=SUB_COLOUR if which == "sub" else OV_COLOUR,
                                    alpha=0.25, lw=0)
            ax.invert_xaxis()
            ax.set_xlabel("Binding energy (eV), aligned to the reference spectrum")
            ax.set_ylabel("Intensity (normalised)")
            ax.set_title(title, fontsize=10)
            ax.legend(title="tilt", fontsize=8, title_fontsize=8, loc="upper left")
        if res.mode == 1:
            axes[0].text(0.99, 0.97, "blue: substrate state   orange: overlayer state", transform=axes[0].transAxes,
                         ha="right", va="top", fontsize=8, color=MUTED)
        fig.suptitle("Check that the shaded areas cover each state at every angle (dotted = excluded angles)",
                     fontsize=9, color=MUTED)
        fig.subplots_adjust(left=0.08, right=0.98, top=0.88, bottom=0.14, wspace=0.25)

    def _draw_fit(self, fig):
        fig.clear()
        ax = fig.add_subplot(111)
        keys = _az_sort(self.results)
        markers = "osD^v<>"
        umax = 1.0
        for i, az in enumerate(keys):
            res = self.results[az]
            colour = "#1f6feb" if len(keys) == 1 else None
            used, excl = res.used, [r for r in res.rows if not r.used and r.x is not None]
            lab = "" if az is None else f"azimuth {az:g}°"
            pts = ax.plot([r.u for r in used], [r.x for r in used], markers[i % len(markers)], ms=7, color=colour,
                          label=(lab + " " if lab else "") + "used")[0]
            if excl:
                ax.plot([r.u for r in excl], [r.x for r in excl], markers[i % len(markers)], ms=7, mfc="none",
                        color="#9ca3af", label="excluded")
            umax = max([umax] + [r.u for r in res.rows if r.x is not None])
            us = np.linspace(0, umax * 1.05, 50)
            ax.plot(us, res.k * us, "-", color=pts.get_color(), lw=1.5,
                    label=f"through origin: d = {res.d_nm:.2f} nm")
            if res.line is not None and i == 0 and len(keys) == 1:
                ax.plot(us, res.line[0] + res.line[1] * us, "--", color="#6b7280", lw=1,
                        label=f"free intercept ({res.line[0]:+.2f}): d = {res.d_free_nm:.2f} nm")
        ax.set_xlim(0, umax * 1.08)
        ax.set_ylim(bottom=0)
        ax.set_xlabel("1 / cos θ")
        res0 = self.results[keys[0]]
        ax.set_ylabel("x = ln(1 + R/R∞)" if res0.mode == 1 else "x from the Thickogram equation")
        ax.set_title("x = d / (λ cos θ): a straight line through the origin means a uniform overlayer", fontsize=10)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, loc="upper left")
        fig.subplots_adjust(left=0.09, right=0.98, top=0.92, bottom=0.12)

    def _draw_azimuth(self, fig):
        fig.clear()
        keys = [k for k in _az_sort(self.results) if k is not None]
        if len(keys) < 2:
            fig.add_subplot(111).axis("off")
            fig.axes[0].text(0.5, 0.5, "This data set has a single azimuth, so there is nothing to plot here.\n"
                                       "If your file stores the azimuth, choose its field under 'Azimuth from' (Data page).",
                             ha="center", va="center", color=MUTED, fontsize=11)
            return
        ax = fig.add_subplot(111, projection="polar")
        theta = np.radians(keys + [keys[0]])
        d = [self.results[k].d_nm for k in keys]
        ax.plot(theta, d + [d[0]], "o-", color=ACCENT)
        mean, spread = self.summary
        ax.plot(np.linspace(0, 2 * np.pi, 200), np.full(200, mean), "--", color="#6b7280", lw=1,
                label=f"mean {mean:.2f} nm" + (f", spread {spread:.2f}" if spread else ""))
        ax.set_title("Thickness (nm) against sample azimuth", fontsize=10, pad=14)
        ax.set_rlim(bottom=0)
        ax.legend(loc="lower right", bbox_to_anchor=(1.15, -0.05), fontsize=8)
        fig.subplots_adjust(left=0.08, right=0.9, top=0.88, bottom=0.1)
