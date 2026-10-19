"""
The "Depth profile (ARXPS)" tab of the XPS Thickness Suite: what a tilt series says about the DEPTH distribution of
several species at once. Panels (a notebook on the right): the Surrey Trend Plot, and (added in later steps)
Stratification, Layers and Profile.

Maths in depthprofile.py, VAMAS reading in vamas.py, CasaXPS fits in casa.py; this module is the Tk/matplotlib glue
and needs numpy + matplotlib (app.py shows an install hint instead if they are missing). Inspired by the NPL ARCtick
package (P. J. Cumpson, 1995) but an independent implementation.

Angle convention: theta is from the surface NORMAL (0 = normal emission) and equals the stage tilt plus an offset.
"""
import math
import os
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import numpy as np
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

import arxps as AR
import depthlayers as DL
import depthprofile as DP
import depthreg as DR
import vamas as VM
from gui_common import ACCENT, BADGE, CARD, MUTED, Tooltip, _to_float

NO_AZIMUTH = "None (single azimuth)"
SRC_FITS = "Fitted states (CasaXPS)"
SRC_REGIONS = "Whole regions (Shirley)"
SRC_AUTO = "Fits where available, else regions"
MAX_LAYERS = 4
SUB_LABEL = "Sub"
GEO_NONE = "none (X-ray spot sets the area)"
GEO_COS = "divide by cos θ (analyser sets it)"
# a colour-blind-safe cycle (Okabe and Ito) with a few extras
PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000",
           "#7F3C8D", "#11A579", "#3969AC", "#E73F74"]


def _az_label(az):
    return "-" if az is None else f"{az:g}"


def _fmt(v, spec="{:.3f}", none="-"):
    return none if v is None else spec.format(v)


class DepthTab(ttk.Frame):
    def __init__(self, parent, app, lam_default=None, compact=False):
        super().__init__(parent)
        self.app, self.compact = app, compact
        self.lam_default = lam_default or (lambda ke: 2.0)
        self.file = None
        self.species = []                    # depthprofile.Species, as loaded
        self.rows = []                       # per species: dict(use, lam, rsf) tk variables
        self.v = {k: tk.StringVar() for k in ("angle_src", "azim_src", "azimuth", "offset", "max_angle", "source", "geometry", "al_factor")}
        self.v["offset"].set("0")
        self.v["max_angle"].set(f"{DP.DEFAULT_MAX_ANGLE:g}")
        self.v["source"].set(SRC_AUTO)
        self.v["geometry"].set(GEO_NONE)
        self.v["al_factor"].set("1.0")
        self.transmission = tk.BooleanVar(value=True)
        self.x_mode = tk.StringVar(value="sec")
        self.normalise = tk.BooleanVar(value=True)
        self.show_fits = tk.BooleanVar(value=True)
        self.strat_a, self.strat_b = tk.StringVar(), tk.StringVar()
        self.strat_all = tk.BooleanVar(value=False)
        self.n_layers = tk.StringVar(value="2")
        self.layer_d = [tk.StringVar(value="1.0") for _ in range(MAX_LAYERS)]
        self.layer_vars = []                 # per species: StringVar holding "1".."K" or SUB_LABEL
        self.layer_fit = None                # depthlayers.LayerFit of the last fit / scan
        self.layer_note = tk.StringVar(value="")
        self._setting = False                # True while the program (not the user) changes the thickness boxes
        self.prof_depth, self.prof_n = tk.StringVar(value="12"), tk.StringVar(value="20")
        self.prof_auto, self.prof_logalpha = tk.BooleanVar(value=True), tk.StringVar(value="3.0")
        self.prof_nonneg, self.prof_frac = tk.BooleanVar(value=False), tk.BooleanVar(value=False)
        self.prof_bands = tk.BooleanVar(value=True)
        self.prof_note = tk.StringVar(value="")
        self.prof_bulk = []                  # per species: (mode StringVar, z_c StringVar)
        self.profiles = []                   # depthreg.ProfileFit of the last calculation
        self.file_label = tk.StringVar(value="No file loaded")
        self.status = tk.StringVar(value="Open a VAMAS tilt series (or use the one in the ARXPS tab).")
        self._after = None

        self.columnconfigure(0, weight=0)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)
        left = ttk.Frame(self)
        left.grid(row=0, column=0, sticky="ns", padx=(0, 8), pady=(8, 4))
        left.rowconfigure(1, weight=1)
        self._build_data(left).grid(row=0, column=0, sticky="ew")
        self._build_species(left).grid(row=1, column=0, sticky="nsew", pady=(8, 0))
        self._build_right().grid(row=0, column=1, sticky="nsew", pady=(8, 4))
        self.bind("<Destroy>", lambda e: self._cancel_after() if e.widget is self else None)
        self._draw()

    def _cancel_after(self):
        if self._after is not None:
            try:
                self.after_cancel(self._after)
            except tk.TclError:
                pass
            self._after = None

    # --- data card ------------------------------------------------------------------------------------
    def _combo(self, parent, var, row, label, hint, values=(), width=24):
        ttk.Label(parent, text=label, style="Card.TLabel").grid(row=row, column=0, sticky="w", pady=2)
        c = ttk.Combobox(parent, textvariable=var, values=list(values), state="readonly", width=width)
        c.grid(row=row, column=1, columnspan=3, sticky="w", padx=6, pady=2)
        Tooltip(c, hint)
        return c

    def _entry(self, parent, key, row, label, hint, col=0, width=6):
        ttk.Label(parent, text=label, style="Card.TLabel").grid(row=row, column=col, sticky="w", pady=2)
        e = ttk.Entry(parent, textvariable=self.v[key], width=width)
        e.grid(row=row, column=col + 1, sticky="w", padx=6, pady=2)
        Tooltip(e, hint)
        return e

    def _build_data(self, parent):
        card = ttk.LabelFrame(parent, text="  Data  ", style="Card.TLabelframe", padding=(8, 4, 8, 6))
        top = ttk.Frame(card, style="Card.TFrame")
        top.grid(row=0, column=0, columnspan=4, sticky="ew")
        b1 = ttk.Button(top, text="Open VAMAS file…", command=self.open_file)
        b1.pack(side="left")
        b2 = ttk.Button(top, text="Use ARXPS tab's file", command=self.use_arxps_file)
        b2.pack(side="left", padx=6)
        Tooltip(b1, "Load a VAMAS (.vms) file with a tilt series: one block per region and tilt.")
        Tooltip(b2, "Take the file that is already open on the ARXPS tab, with its angle settings.")
        ttk.Label(card, textvariable=self.file_label, style="CardMuted.TLabel", wraplength=330, justify="left").grid(
            row=1, column=0, columnspan=4, sticky="w", pady=(2, 4))
        self.angle_cb = self._combo(card, self.v["angle_src"], 2, "Angle from:", "Where the tilt comes from (the "
                                    "K-Alpha stores it in experimental variable 1).")
        self._entry(card, "offset", 3, "Offset (°):", "Added to the tilt to give the emission angle from the surface "
                    "normal.")
        self._entry(card, "max_angle", 3, "Max angle (°):", "Angles above this are drawn but left out of fits: the "
                    "simple exponential model fails at grazing emission.", col=2)
        self.azim_cb = self._combo(card, self.v["azim_src"], 4, "Azimuth from:", "Optional: the field holding the "
                                   "azimuth. One azimuth is analysed at a time.")
        self.az_cb = self._combo(card, self.v["azimuth"], 5, "Azimuth:", "Which azimuth to analyse.", width=10)
        self.az_cb.bind("<<ComboboxSelected>>", lambda e: self.reload())
        self.src_cb = self._combo(card, self.v["source"], 6, "Species from:", "Fitted states: every chemical state of "
                                  "the CasaXPS fits stored in the file. Whole regions: each region integrated with a "
                                  "Shirley background (any file). Fits where available: the states of regions that "
                                  "have a fit, whole regions for the rest.", [SRC_AUTO, SRC_FITS, SRC_REGIONS])
        for cb in (self.angle_cb, self.azim_cb, self.src_cb):
            cb.bind("<<ComboboxSelected>>", lambda e: self.reload())
        self.geo_cb = self._combo(card, self.v["geometry"], 7, "Geometry:", "Any change of the analysed area with tilt "
                                  "scales every intensity by the same factor. 'none' assumes the X-ray spot defines "
                                  "the area. If all species fall with angle much faster than their depths can explain "
                                  "(a species placed deeper than the whole sample is thick), try dividing by cos θ.",
                                  [GEO_NONE, GEO_COS])
        self.geo_cb.bind("<<ComboboxSelected>>", lambda e: self.reload())
        ck = ttk.Checkbutton(card, text="Transmission-correct", variable=self.transmission, style="Card.TCheckbutton",
                             command=self.reload)
        ck.grid(row=8, column=0, columnspan=2, sticky="w", pady=(2, 0))
        Tooltip(ck, "Divide by the transmission function at each peak's kinetic energy so peaks at different energies "
                    "can be compared.")
        ttk.Button(card, text="Reload species", command=self.reload).grid(row=8, column=2, columnspan=2, sticky="e")
        self._entry(card, "al_factor", 9, "AL / IMFP:", "Starting attenuation lengths are TPP-2M inelastic mean free paths. "
                    "Attenuation lengths, which include some elastic scattering, are typically 10 to 30 percent shorter "
                    "(Cumpson and Seah 1997): set 0.8 for instance, then press Reload species. You can still edit each "
                    "lambda in the table.")
        return card

    # --- species card ---------------------------------------------------------------------------------
    def _build_species(self, parent):
        card = ttk.LabelFrame(parent, text="  Species  ", style="Card.TLabelframe", padding=(8, 4, 8, 6))
        card.rowconfigure(1, weight=1)
        card.columnconfigure(0, weight=1)
        ttk.Label(card, textvariable=self.status, style="CardMuted.TLabel", wraplength=330, justify="left").grid(
            row=0, column=0, columnspan=2, sticky="w")
        holder = ttk.Frame(card, style="Card.TFrame")
        holder.grid(row=1, column=0, sticky="nsew", pady=(4, 0))
        holder.rowconfigure(0, weight=1)
        holder.columnconfigure(0, weight=1)
        self.sp_canvas = tk.Canvas(holder, bg=CARD, highlightthickness=0, height=150 if self.compact else 220, width=340)
        scroll = ttk.Scrollbar(holder, orient="vertical", command=self.sp_canvas.yview)
        self.sp_canvas.configure(yscrollcommand=scroll.set)
        self.sp_canvas.grid(row=0, column=0, sticky="nsew")
        scroll.grid(row=0, column=1, sticky="ns")
        self.sp_frame = ttk.Frame(self.sp_canvas, style="Card.TFrame")
        self.sp_canvas.create_window((0, 0), window=self.sp_frame, anchor="nw")
        self.sp_frame.bind("<Configure>", lambda e: self.sp_canvas.configure(scrollregion=self.sp_canvas.bbox("all")))
        return card

    def _fill_species(self):
        for w in self.sp_frame.winfo_children():
            w.destroy()
        self.rows = []
        for c, text in enumerate(("Use", "Species", "λ (nm)", "RSF", "KE (eV)")):
            ttk.Label(self.sp_frame, text=text, style="ColHead.TLabel").grid(row=0, column=c, padx=3, sticky="w")
        for i, sp in enumerate(self.species, start=1):
            use = tk.BooleanVar(value=True)
            lam = tk.StringVar(value=f"{sp.lam_nm:.2f}")
            rsf = tk.StringVar(value=f"{sp.rsf:g}")
            ttk.Checkbutton(self.sp_frame, variable=use, style="Card.TCheckbutton", command=self.refresh).grid(
                row=i, column=0, padx=3)
            tk.Label(self.sp_frame, text="●", fg=PALETTE[(i - 1) % len(PALETTE)], bg=CARD).grid(row=i, column=1, sticky="w")
            ttk.Label(self.sp_frame, text=sp.name, style="Card.TLabel").grid(row=i, column=1, sticky="w", padx=(16, 0))
            for col, var, hint in ((2, lam, "Attenuation length of this species' photoelectrons (nm). The default is a "
                                            "TPP-2M estimate for an oxide-like matrix: replace it with a value for your "
                                            "material."),
                                   (3, rsf, "Relative sensitivity factor (the fit's own RSF when it has one).")):
                e = ttk.Entry(self.sp_frame, textvariable=var, width=6)
                e.grid(row=i, column=col, padx=3, pady=1)
                Tooltip(e, hint)
                var.trace_add("write", lambda *a: self._schedule())
            ttk.Label(self.sp_frame, text=_fmt(sp.ke_ev, "{:.0f}"), style="CardMuted.TLabel").grid(row=i, column=4, padx=3)
            self.rows.append(dict(use=use, lam=lam, rsf=rsf))

    # --- loading ------------------------------------------------------------------------------------------
    def open_file(self):
        path = filedialog.askopenfilename(title="Open a VAMAS tilt series",
                                          filetypes=[("VAMAS files", "*.vms *.vamas"), ("All files", "*.*")])
        if path:
            self.load_file(path)

    def use_arxps_file(self):
        arx = getattr(self.app, "arxps", None)
        if arx is None or arx.file is None:
            messagebox.showinfo("No file", "Open a VAMAS file on the ARXPS tab first.")
            return
        self._set_file(arx.file, os.path.basename(arx.file.path or "file"))

    def load_file(self, path):
        try:
            vf = VM.read(path)
        except VM.VamasError as exc:
            messagebox.showerror("Cannot read the file", str(exc))
            return False
        if not vf.blocks:
            messagebox.showerror("Cannot read the file", "The file holds no spectra.")
            return False
        self._set_file(vf, os.path.basename(path))
        return True

    def _set_file(self, vf, name):
        self.file = vf
        self.file_label.set(f"{name}: {len(vf.blocks)} blocks, regions {', '.join(vf.regions())}")
        angle, azim = vf.default_sources()
        present = vf.sources_present()
        labels = [VM.ANGLE_SOURCES[s] for s in present]
        self.angle_cb.configure(values=labels or list(VM.ANGLE_SOURCES.values()))
        self.azim_cb.configure(values=[NO_AZIMUTH] + labels)
        self.v["angle_src"].set(VM.ANGLE_SOURCES[angle] if angle else "")
        self.v["azim_src"].set(VM.ANGLE_SOURCES[azim] if azim and azim in present else NO_AZIMUTH)
        has_fits = any(AR.fit_regions(self._series(r)) for r in self._spectral_regions() if self._series(r))
        self.v["source"].set(SRC_AUTO if has_fits else SRC_REGIONS)
        self.reload()

    def _key(self, label):
        for key, text in VM.ANGLE_SOURCES.items():
            if text == label:
                return key
        return None

    def _spectral_regions(self):
        names = self.file.regions() if self.file else []
        return [n for n in names if "survey" not in n.lower()] or names

    def _series(self, region):
        """The blocks of a region at the chosen azimuth as [(tilt, block), ...] (empty if there is no angle)."""
        angle_src = self._key(self.v["angle_src"].get())
        if angle_src is None or not self.file:
            return []
        series = self.file.series(region, angle_src, self._key(self.v["azim_src"].get()))
        if not series:
            return []
        want = self.v["azimuth"].get()
        keys = sorted(series, key=lambda a: (a is not None, a if a is not None else 0.0))
        self.az_cb.configure(values=[_az_label(k) for k in keys])
        for k in keys:
            if _az_label(k) == want:
                return series[k]
        self.v["azimuth"].set(_az_label(keys[0]))
        return series[keys[0]]

    def reload(self):
        """(Re)build the species from the file and the current settings."""
        if not self.file:
            return
        try:
            offset = self._number("offset", "Offset")
            factor = self._number("al_factor", "AL / IMFP")
            if factor <= 0:
                raise ValueError("AL / IMFP must be positive.")
            tr = self.transmission.get()
            species = []
            for region in self._spectral_regions():
                series = self._series(region)
                if len(series) < 2:
                    continue
                source = self.v["source"].get()
                fit_names = AR.fit_regions(series) if source in (SRC_FITS, SRC_AUTO) else []
                if fit_names:
                    for name in fit_names:
                        species += DP.species_from_fits(series, offset, tr, name)
                elif source != SRC_FITS:
                    species.append(DP.species_from_region(series, offset, True, tr))
            if not species:
                raise ValueError("No species found: " + ("none of the regions holds a CasaXPS fit (try Whole regions)"
                                                         if self.v["source"].get() == SRC_FITS else
                                                         "check the angle source (Data card)") + ".")
        except ValueError as exc:
            self.species, self.rows = [], []
            self._fill_species()
            self.status.set(str(exc))
            self._draw()
            return
        for sp in species:
            sp.lam_nm = (float(self.lam_default(sp.ke_ev)) if sp.ke_ev else 2.0) * factor
        species = DP.apply_geometry(species, 1.0 if self.v["geometry"].get() == GEO_COS else 0.0)
        self.species = species
        self._fill_angle_choices()
        self._fill_species()
        self._default_layers()
        self._default_bulk()
        self.status.set(f"{len(species)} species at {len(species[0].angles)} angles; set each λ for your material.")
        self.refresh()

    def _number(self, key, label):
        x = _to_float(self.v[key].get())
        if x is None or not math.isfinite(x):
            raise ValueError(f"Please enter a number for '{label}'.")
        return x

    def _schedule(self):
        if self._after is not None:
            self.after_cancel(self._after)
        self._after = self.after(250, self.refresh)

    def active(self):
        """The species currently ticked, with the typed λ and RSF applied."""
        out = []
        for sp, row in zip(self.species, self.rows):
            if not row["use"].get():
                continue
            lam, rsf = _to_float(row["lam"].get()), _to_float(row["rsf"].get())
            if lam and lam > 0:
                sp.lam_nm = lam
            if rsf and rsf > 0:
                sp.rsf = rsf
            out.append(sp)
        return out

    def refresh(self):
        self._after = None
        self._draw()
        self._draw_strat()
        self._draw_layers()

    # --- right-hand side ------------------------------------------------------------------------------
    def _build_right(self):
        frame = ttk.Frame(self)
        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(0, weight=1)
        nb = ttk.Notebook(frame, style="Sub.TNotebook")
        nb.grid(row=0, column=0, sticky="nsew")
        self.plot_nb = nb
        self._build_trend_page(nb)
        self._build_strat_page(nb)
        self._build_layers_page(nb)
        self._build_profile_page(nb)
        return frame

    def _canvas(self, parent):
        fig = Figure(figsize=(8, 5), dpi=100)
        canvas = FigureCanvasTkAgg(fig, master=parent)
        toolbar = NavigationToolbar2Tk(canvas, parent, pack_toolbar=False)
        toolbar.configure(width=250)
        for child in [toolbar] + list(toolbar.winfo_children()):
            try:
                child.configure(bg=CARD)
            except tk.TclError:
                pass
        toolbar.pack(side="bottom", fill="x")
        canvas.get_tk_widget().pack(fill="both", expand=True)
        return fig, canvas

    def _build_trend_page(self, nb):
        page = ttk.Frame(nb, style="Card.TFrame", padding=4)
        nb.add(page, text="Surrey Trend Plot")
        bar = ttk.Frame(page, style="Card.TFrame")
        bar.pack(side="top", fill="x")
        for text, val in (("ln I against sec θ", "sec"), ("ln I against θ", "theta")):
            ttk.Radiobutton(bar, text=text, value=val, variable=self.x_mode, style="Card.TRadiobutton",
                            command=self._draw).pack(side="left", padx=(0, 10))
        for text, var, tip in (("Normalise to the first angle", self.normalise,
                                "Divide each curve by its value at the first angle, so parallel lines are easy to see."),
                               ("Fit lines", self.show_fits, "Straight-line fit of ln I against sec θ for each species.")):
            c = ttk.Checkbutton(bar, text=text, variable=var, style="Card.TCheckbutton", command=self._draw)
            c.pack(side="left", padx=(0, 10))
            Tooltip(c, tip)
        lower = ttk.Frame(page, style="Card.TFrame")
        lower.pack(side="bottom", fill="x", pady=(6, 0))
        cols = ("sp", "ke", "slope", "se", "rms", "depth", "n", "similar", "note")
        heads = ("Species", "KE (eV)", "Slope", "± se", "Scatter", "Depth (nm)", "Points", "Parallel to", "Note")
        widths = (90, 60, 60, 60, 60, 75, 50, 130, 260)
        self.trend_table = ttk.Treeview(lower, columns=cols, show="headings", height=4 if self.compact else 6)
        for c, h, w in zip(cols, heads, widths):
            self.trend_table.heading(c, text=h)
            self.trend_table.column(c, width=w, anchor="w" if c in ("sp", "similar", "note") else "e",
                                    stretch=(c == "note"))
        self.trend_table.pack(fill="x")
        self.trend_fig, self.trend_canvas = self._canvas(page)

    # --- stratification -----------------------------------------------------------------------------------
    def _fill_angle_choices(self):
        angles = sorted({float(a) for sp in self.species for a in sp.angles})
        labels = [f"{a:g}" for a in angles]
        for cb in (self.strat_a_cb, self.strat_b_cb):
            cb.configure(values=labels)
        if labels:
            self.strat_a.set(labels[0])
            try:
                limit = self._number("max_angle", "Max angle")
            except ValueError:
                limit = DP.DEFAULT_MAX_ANGLE
            ok = [a for a in angles if a <= limit and a != angles[0]] or angles[1:] or angles
            self.strat_b.set(f"{min(ok, key=lambda a: abs(a - 67.5)):g}")

    def _build_strat_page(self, nb):
        page = ttk.Frame(nb, style="Card.TFrame", padding=4)
        nb.add(page, text="Stratification")
        bar = ttk.Frame(page, style="Card.TFrame")
        bar.pack(side="top", fill="x")
        ttk.Label(bar, text="Angles:", style="Card.TLabel").pack(side="left")
        self.strat_a_cb = ttk.Combobox(bar, textvariable=self.strat_a, state="readonly", width=6)
        self.strat_b_cb = ttk.Combobox(bar, textvariable=self.strat_b, state="readonly", width=6)
        self.strat_a_cb.pack(side="left", padx=4)
        ttk.Label(bar, text="and", style="Card.TLabel").pack(side="left")
        self.strat_b_cb.pack(side="left", padx=4)
        for cb in (self.strat_a_cb, self.strat_b_cb):
            cb.bind("<<ComboboxSelected>>", lambda e: self._draw_strat())
        Tooltip(self.strat_a_cb, "The near-normal angle. Two measurements are enough: normal emission and one oblique "
                                 "angle (about 67.5 degrees is typical).")
        c = ttk.Checkbutton(bar, text="Use every angle (fit instead)", variable=self.strat_all,
                            style="Card.TCheckbutton", command=self._draw_strat)
        c.pack(side="left", padx=14)
        Tooltip(c, "Fit ln I against sec θ over all angles up to the maximum instead of using just two.")
        lower = ttk.Frame(page, style="Card.TFrame")
        lower.pack(side="bottom", fill="x", pady=(6, 0))
        cols = ("sp", "z", "se", "amount", "ratio", "curv", "note")
        heads = ("Species", "Depth (nm)", "± se", "Amount (%)", "I(b)/I(a)", "Curvature t", "Note")
        widths = (90, 80, 65, 80, 75, 85, 340)
        self.strat_table = ttk.Treeview(lower, columns=cols, show="headings", height=4 if self.compact else 6)
        for col, h, w in zip(cols, heads, widths):
            self.strat_table.heading(col, text=h)
            self.strat_table.column(col, width=w, anchor="w" if col in ("sp", "note") else "e", stretch=(col == "note"))
        self.strat_table.pack(fill="x")
        self.strat_fig, self.strat_canvas = self._canvas(page)

    def _draw_strat(self):
        fig = self.strat_fig
        fig.clear()
        self.strat_table.delete(*self.strat_table.get_children())
        species = self.active() if self.species else []
        if not species:
            ax = fig.add_subplot(111)
            ax.axis("off")
            ax.text(0.5, 0.5, "Each species is placed at its average depth from just two angles (normal and "
                              "oblique).\nBar position = depth, bar width = relative amount. Misleading if a species "
                              "sits in two layers.", ha="center", va="center", color=MUTED, fontsize=11)
            self.strat_canvas.draw_idle()
            return
        try:
            max_angle = self._number("max_angle", "Max angle")
            a, b = float(self.strat_a.get()), float(self.strat_b.get())
        except ValueError:
            max_angle, a, b = DP.DEFAULT_MAX_ANGLE, None, None
        rows = DP.stratification(species, a, b, self.strat_all.get(), max_angle)
        ax = fig.add_subplot(111)
        valid = [r for r in rows if r.z_nm is not None and r.amount_pct is not None]
        span = max([abs(r.z_nm) for r in valid] + [1.0])
        for sp, r in zip(species, rows):
            if r.z_nm is None or r.amount_pct is None:
                continue
            colour = PALETTE[self.species.index(sp) % len(PALETTE)]
            width = max(0.5 * span * r.amount_pct / 100.0, 0.02 * span)
            ax.add_patch(Rectangle((r.z_nm - width / 2, 0), width, 1.0, color=colour, alpha=0.55, lw=0))
            ax.text(r.z_nm, 1.03, sp.name, rotation=90, ha="center", va="bottom", fontsize=8, color=colour)
            if r.z_se_nm:
                ax.errorbar([r.z_nm], [0.5], xerr=[r.z_se_nm], color="black", capsize=3, lw=1)
        ax.set_xlim(min(-0.05 * span, min([r.z_nm for r in valid] + [0]) - 0.1 * span), 1.15 * span)
        ax.set_ylim(0, 1.45)
        ax.set_yticks([])
        ax.set_xlabel("Average depth below the surface (nm)")
        ax.set_title("Stratification: bar position = average depth, bar width = relative amount", fontsize=10)
        ax.grid(alpha=0.3, axis="x")
        fig.subplots_adjust(left=0.05, right=0.98, top=0.9, bottom=0.12)
        self.strat_canvas.draw_idle()
        for r in rows:
            self.strat_table.insert("", "end", values=(
                r.name, _fmt(r.z_nm, "{:.2f}"), _fmt(r.z_se_nm, "{:.2f}"), _fmt(r.amount_pct, "{:.1f}"),
                _fmt(r.ratio, "{:.3f}"), _fmt(r.curved_t, "{:.1f}"), r.note))

    # --- layers (Gries) -----------------------------------------------------------------------------------
    def _k(self):
        try:
            return max(1, min(MAX_LAYERS, int(self.n_layers.get())))
        except ValueError:
            return 2

    def _layer_labels(self):
        return [str(i + 1) for i in range(self._k())] + [SUB_LABEL]

    def _default_layers(self):
        """First guess: layer of each species by average depth, deepest in the substrate."""
        self.layer_fit = None
        self.layer_vars = []
        if not self.species:
            return
        try:
            guess = DL.default_assignment(self.species, self._k(), self._max_angle())
        except ValueError:
            guess = [0] * len(self.species)
        labels = self._layer_labels()
        self.layer_vars = [tk.StringVar(value=SUB_LABEL if g == DL.SUBSTRATE_LAYER else labels[min(g, self._k() - 1)])
                           for g in guess]
        for v in self.layer_vars:
            v.trace_add("write", lambda *a: self._layers_changed())
        self._rebuild_layer_grid()

    def _max_angle(self):
        try:
            return self._number("max_angle", "Max angle")
        except ValueError:
            return DP.DEFAULT_MAX_ANGLE

    def _build_layers_page(self, nb):
        page = ttk.Frame(nb, style="Card.TFrame", padding=4)
        nb.add(page, text="Layers (Gries)")
        bar = ttk.Frame(page, style="Card.TFrame")
        bar.pack(side="top", fill="x")
        ttk.Label(bar, text="Layers above the substrate:", style="Card.TLabel").pack(side="left")
        cb = ttk.Combobox(bar, textvariable=self.n_layers, values=[str(i) for i in range(1, MAX_LAYERS + 1)],
                          state="readonly", width=3)
        cb.pack(side="left", padx=4)
        cb.bind("<<ComboboxSelected>>", lambda e: self._default_layers() or self._draw_layers())
        for text, cmd, tip in (("Auto-assign", self._auto_assign, "Put each species in a layer by its average depth: the "
                                "shallowest in layer 1, the deepest in the substrate. Change any row afterwards."),
                               ("Fit thicknesses", self.fit_layers_now, "Least-squares thicknesses from the shape of each "
                                "species' intensity against angle (no sensitivity factors needed)."),
                               ("Scan ranges", self.scan_layers_now, "For each layer, the range of thickness that fits "
                                "within the noise (the others re-fitted). A wide range means the data cannot say.")):
            b = ttk.Button(bar, text=text, command=cmd)
            b.pack(side="left", padx=4)
            Tooltip(b, tip)
        ttk.Label(bar, textvariable=self.layer_note, style="CardMuted.TLabel", wraplength=420, justify="left").pack(
            side="left", padx=8)
        mid = ttk.Frame(page, style="Card.TFrame")
        mid.pack(side="top", fill="x", pady=(4, 0))
        self.layer_grid = ttk.Frame(mid, style="Card.TFrame")
        self.layer_grid.pack(side="left", anchor="n")
        lower = ttk.Frame(page, style="Card.TFrame")
        lower.pack(side="bottom", fill="x", pady=(6, 0))
        cols = ("layer", "species", "d", "range", "note")
        heads = ("Layer", "Species", "Thickness (nm)", "Range that fits (nm)", "Note")
        widths = (60, 220, 100, 140, 330)
        self.layer_table = ttk.Treeview(lower, columns=cols, show="headings", height=4)
        for col, h, w in zip(cols, heads, widths):
            self.layer_table.heading(col, text=h)
            self.layer_table.column(col, width=w, anchor="w" if col in ("species", "note") else "e",
                                    stretch=(col == "note"))
        self.layer_table.pack(fill="x")
        self.layer_fig, self.layer_canvas = self._canvas(page)

    def _rebuild_layer_grid(self):
        self._setting = True
        try:
            self._build_layer_grid()
        finally:
            self._setting = False

    def _build_layer_grid(self):
        for w in self.layer_grid.winfo_children():
            w.destroy()
        k = self._k()
        ttk.Label(self.layer_grid, text="Thickness (nm)", style="ColHead.TLabel").grid(row=0, column=0, columnspan=3,
                                                                                       sticky="w")
        for i in range(k):
            ttk.Label(self.layer_grid, text=f"Layer {i + 1}", style="Card.TLabel").grid(row=i + 1, column=0, sticky="w")
            sc = ttk.Scale(self.layer_grid, from_=0.0, to=10.0, orient="horizontal", length=130,
                           command=lambda v, j=i: self._scale_moved(j, v))
            sc.set(_to_float(self.layer_d[i].get()) or 1.0)
            sc.grid(row=i + 1, column=1, padx=4)
            e = ttk.Entry(self.layer_grid, textvariable=self.layer_d[i], width=6)
            e.grid(row=i + 1, column=2)
            self.layer_d[i].trace_add("write", lambda *a: self._layers_changed())
            Tooltip(e, "Layer thickness in nm: type it or drag the slider and watch the curves.")
        ttk.Label(self.layer_grid, text="  ", style="Card.TLabel").grid(row=0, column=3)
        ttk.Label(self.layer_grid, text="Species in layer", style="ColHead.TLabel").grid(row=0, column=4, columnspan=2,
                                                                                         sticky="w")
        for i, (sp, var) in enumerate(zip(self.species, self.layer_vars)):
            ttk.Label(self.layer_grid, text=sp.name, style="Card.TLabel").grid(row=i + 1, column=4, sticky="w", padx=(0, 6))
            ttk.Combobox(self.layer_grid, textvariable=var, values=self._layer_labels(), state="readonly",
                         width=5).grid(row=i + 1, column=5, pady=1)

    def _scale_moved(self, k, value):
        if self._setting:
            return
        self.layer_d[k].set(f"{float(value):.2f}")

    def _auto_assign(self):
        self._default_layers()
        self._draw_layers()

    def _layers_changed(self):
        """The user changed a thickness or an assignment: any stored fit and its ranges no longer apply."""
        if not self._setting:
            self.layer_fit = None
        self._schedule()

    def _set_thicknesses(self, values):
        self._setting = True
        try:
            for i, v in enumerate(values):
                self.layer_d[i].set(f"{v:.2f}")
        finally:
            self._setting = False

    def _assignment(self, species):
        """Layer index (0-based, or SUBSTRATE_LAYER) of each of the given species."""
        index = {id(sp): v for sp, v in zip(self.species, self.layer_vars)}
        out = []
        for sp in species:
            label = index[id(sp)].get() if id(sp) in index else SUB_LABEL
            out.append(DL.SUBSTRATE_LAYER if label == SUB_LABEL else min(int(label) - 1, self._k() - 1))
        return out

    def _thicknesses(self):
        d = []
        for i in range(self._k()):
            v = _to_float(self.layer_d[i].get())
            d.append(max(v, 0.0) if v is not None else 0.0)
        return d

    def fit_layers_now(self):
        species = self.active()
        if not species:
            return
        try:
            fit = DL.fit_layers(species, self._assignment(species), self._k(), d0=self._thicknesses(),
                                max_angle=self._max_angle())
        except ValueError as exc:
            self.layer_note.set(str(exc))
            return
        self.layer_fit = fit
        self._set_thicknesses(fit.d_nm)
        self.layer_note.set(f"Fitted: chi-squared per point {_fmt(fit.chi2_red, '{:.2f}')}. " + " ".join(fit.warnings))
        self._rebuild_layer_grid()
        self._draw_layers()

    def scan_layers_now(self):
        species = self.active()
        if not species:
            return
        self.fit_layers_now()
        if self.layer_fit is None:
            return
        self.layer_note.set("Scanning the thickness ranges…")
        self.update_idletasks()
        DL.scan_ranges(self.layer_fit, species, self._k(), max_angle=self._max_angle())
        self.layer_note.set(f"Fitted: chi-squared per point {_fmt(self.layer_fit.chi2_red, '{:.2f}')}; the ranges are "
                            "where chi-squared stays within the noise.")
        self._draw_layers()

    def _draw_layers(self):
        fig = self.layer_fig
        fig.clear()
        self.layer_table.delete(*self.layer_table.get_children())
        species = self.active() if self.species else []
        if not species or not self.layer_vars:
            ax = fig.add_subplot(111)
            ax.axis("off")
            ax.text(0.5, 0.5, "Put each species in a layer, then drag the thickness sliders until the predicted curves "
                              "follow the points,\nor press Fit thicknesses. Only the shape of each peak against angle "
                              "is used: no sensitivity factors.", ha="center", va="center", color=MUTED, fontsize=11)
            self.layer_canvas.draw_idle()
            return
        assignment = self._assignment(species)
        d = self._thicknesses()
        max_angle = self._max_angle()
        ax1, ax2 = fig.subplots(1, 2, gridspec_kw={"width_ratios": [3, 1.4]})
        th = np.linspace(0.0, 80.0, 81)
        for sp, layer in zip(species, assignment):
            colour = PALETTE[self.species.index(sp) % len(PALETTE)]
            meas = DL.measured_normalised(sp, max_angle)
            if meas is None:
                continue
            angles, y, sig, _ = meas
            ax1.errorbar(angles, y, yerr=sig, fmt="o", color=colour, ms=5, capsize=2, label=sp.name)
            lam = [sp.lam_nm] * len(d)
            curve = DL.layer_intensity(th, d, layer, lam)
            ref = DL.layer_intensity([angles[0]], d, layer, lam)[0]
            if ref > 1e-12:
                ax1.plot(th, curve / ref, "-", color=colour, lw=1.2)
        ax1.set_xlabel("Emission angle θ (°)")
        ax1.set_ylabel("Intensity / intensity at the first angle")
        ax1.set_title("Measured (points) and predicted (lines)", fontsize=10)
        ax1.grid(alpha=0.3)
        ax1.legend(fontsize=8, loc="best")
        # the layer stack
        names = {k: [] for k in list(range(self._k())) + [DL.SUBSTRATE_LAYER]}
        for sp, layer in zip(species, assignment):
            names[layer].append(sp.name)
        top, sub_h = 0.0, max(0.3 * max(sum(d), 1.0), 0.5)
        cmap = ["#cfe3f6", "#fde2c4", "#d4efdf", "#f1d6e6"]
        for k, thick in enumerate(d):
            ax2.add_patch(Rectangle((0, top), 1, thick, facecolor=cmap[k % 4], edgecolor="black", lw=0.8))
            label = f"{k + 1}: {thick:.2f} nm" + ("\n" + ", ".join(names[k]) if names[k] else "")
            ax2.text(0.5, top + thick / 2, label, ha="center", va="center", fontsize=8)
            top += thick
        ax2.add_patch(Rectangle((0, top), 1, sub_h, facecolor="#e5e7eb", edgecolor="black", lw=0.8, hatch="//"))
        ax2.text(0.5, top + sub_h / 2, "Substrate\n" + ", ".join(names[DL.SUBSTRATE_LAYER]), ha="center", va="center",
                 fontsize=8, bbox=dict(facecolor="white", alpha=0.8, lw=0))
        ax2.set_xlim(0, 1)
        ax2.set_ylim(top + sub_h, 0)
        ax2.set_xticks([])
        ax2.set_ylabel("Depth (nm)")
        ax2.set_title("Layer stack", fontsize=10)
        fig.subplots_adjust(left=0.07, right=0.98, top=0.9, bottom=0.12, wspace=0.25)
        self.layer_canvas.draw_idle()
        fit = self.layer_fit
        for k, thick in enumerate(d):
            rng = fit.ranges[k] if fit is not None and k < len(fit.ranges) and fit.ranges[k] else None
            note = ""
            if rng is not None and rng[1] - rng[0] > 2.0 * max(thick, 0.2):
                note = "poorly determined: the data allow a wide range"
            self.layer_table.insert("", "end", values=(
                k + 1, ", ".join(names[k]) or "-", f"{thick:.2f}",
                f"{rng[0]:.2f} to {rng[1]:.2f}" if rng else "-", note))
        chi, npts = DL.gries_chi2(d, species, assignment, max_angle)
        self.layer_table.insert("", "end", values=("Sub", ", ".join(names[DL.SUBSTRATE_LAYER]) or "-", "-", "-",
                                                    f"chi-squared per point {chi / npts:.2f}" if npts else ""))

    # --- profile (regularisation) -------------------------------------------------------------------------
    BULK_LABELS = {"free": "free", "constant": "constant below z", "zero": "zero below z"}

    def _default_bulk(self):
        self.profiles = []
        self.prof_bulk = [(tk.StringVar(value=self.BULK_LABELS["free"]), tk.StringVar(value="8")) for _ in self.species]
        self._rebuild_bulk_grid()

    def _build_profile_page(self, nb):
        page = ttk.Frame(nb, style="Card.TFrame", padding=4)
        nb.add(page, text="Profile (regularised)")
        bar = ttk.Frame(page, style="Card.TFrame")
        bar.pack(side="top", fill="x")
        for text, var, width, tip in (("Max depth (nm):", self.prof_depth, 5, "How deep the profile is allowed to go."),
                                      ("Depth points:", self.prof_n, 4, "Number of depth slabs (5 to 60).")):
            ttk.Label(bar, text=text, style="Card.TLabel").pack(side="left", padx=(0, 2))
            e = ttk.Entry(bar, textvariable=var, width=width)
            e.pack(side="left", padx=(0, 8))
            Tooltip(e, tip)
        ck = ttk.Checkbutton(bar, text="Auto α", variable=self.prof_auto, style="Card.TCheckbutton",
                             command=self._calc_profile)
        ck.pack(side="left")
        Tooltip(ck, "Choose the smoothing strength from the data: the smoothest profile that still fits within the "
                    "counting noise. Untick to set log10 of α by hand.")
        e = ttk.Entry(bar, textvariable=self.prof_logalpha, width=5)
        e.pack(side="left", padx=(2, 8))
        Tooltip(e, "log10 of the smoothing strength (relative units) used when Auto α is off: larger is smoother.")
        for text, var, tip in (("Non-negative", self.prof_nonneg, "Force every concentration to be zero or positive."),
                               ("Atomic fraction", self.prof_frac, "Show the species as shares that add up to 1 at each "
                                                                  "depth."),
                               ("Noise bands", self.prof_bands, "Draw ±1 standard deviation from the counting noise "
                                                                "(it does not include the smoothing bias).")):
            c = ttk.Checkbutton(bar, text=text, variable=var, style="Card.TCheckbutton", command=self._calc_profile)
            c.pack(side="left", padx=(0, 6))
            Tooltip(c, tip)
        b = ttk.Button(bar, text="Calculate profile", command=self._calc_profile)
        b.pack(side="left", padx=6)
        ttk.Label(page, textvariable=self.prof_note, style="CardMuted.TLabel", wraplength=900, justify="left").pack(
            side="top", fill="x", pady=(2, 0))
        mid = ttk.Frame(page, style="Card.TFrame")
        mid.pack(side="top", fill="x", pady=(2, 0))
        self.bulk_grid = ttk.Frame(mid, style="Card.TFrame")
        self.bulk_grid.pack(side="left", anchor="n")
        lower = ttk.Frame(page, style="Card.TFrame")
        lower.pack(side="bottom", fill="x", pady=(6, 0))
        cols = ("sp", "alpha", "chi", "bulk", "peak", "note")
        heads = ("Species", "log10 α", "χ² per point", "Bulk", "Peak conc.", "Note")
        widths = (90, 70, 90, 140, 90, 420)
        self.prof_table = ttk.Treeview(lower, columns=cols, show="headings", height=4)
        for col, h, w in zip(cols, heads, widths):
            self.prof_table.heading(col, text=h)
            self.prof_table.column(col, width=w, anchor="w" if col in ("sp", "bulk", "note") else "e",
                                   stretch=(col == "note"))
        self.prof_table.pack(fill="x")
        self.prof_fig, self.prof_canvas = self._canvas(page)

    def _rebuild_bulk_grid(self):
        for w in self.bulk_grid.winfo_children():
            w.destroy()
        ttk.Label(self.bulk_grid, text="Prior knowledge: where a species is constant or absent (nm)",
                  style="ColHead.TLabel").grid(row=0, column=0, columnspan=7, sticky="w")
        for i, (sp, (mode, zc)) in enumerate(zip(self.species, self.prof_bulk)):
            r, c = 1 + i // 2, 3 * (i % 2)
            ttk.Label(self.bulk_grid, text=sp.name, style="Card.TLabel").grid(row=r, column=c, sticky="w", padx=(0, 4))
            cb = ttk.Combobox(self.bulk_grid, textvariable=mode, values=list(self.BULK_LABELS.values()), state="readonly",
                              width=15)
            cb.grid(row=r, column=c + 1, padx=2, pady=1)
            Tooltip(cb, "free: no assumption. constant below z: the species has one concentration below depth z (a bulk "
                        "species). zero below z: the species is not present below depth z (a surface species).")
            ttk.Entry(self.bulk_grid, textvariable=zc, width=5).grid(row=r, column=c + 2, padx=(2, 14))

    def _bulk_of(self, sp):
        i = self.species.index(sp)
        mode, zc = self.prof_bulk[i] if i < len(self.prof_bulk) else (None, None)
        label = mode.get() if mode else self.BULK_LABELS["free"]
        key = next((k for k, v in self.BULK_LABELS.items() if v == label), "free")
        return key, (_to_float(zc.get()) if zc else None)

    def _calc_profile(self):
        species = self.active() if self.species else []
        self.profiles = []
        if not species:
            self._draw_profile()
            return
        try:
            depth = float(self.prof_depth.get())
            n = int(float(self.prof_n.get()))
            if not (5 <= n <= 60) or depth <= 0:
                raise ValueError("Max depth must be positive and the number of depth points between 5 and 60.")
        except ValueError as exc:
            self.prof_note.set(str(exc) if "Max depth" in str(exc) else "Enter numbers for the depth and the points.")
            self._draw_profile()
            return
        alpha = None
        if not self.prof_auto.get():
            la = _to_float(self.prof_logalpha.get())
            alpha = 10.0 ** la if la is not None else None
        notes = []
        for sp in species:
            bulk, zc = self._bulk_of(sp)
            try:
                self.profiles.append((sp, DR.fit_profile(sp, depth, n, alpha, self.prof_nonneg.get(), bulk, zc,
                                                         self._max_angle()), bulk, zc))
            except ValueError as exc:
                notes.append(str(exc))
        self.prof_note.set(" ".join(notes))
        self._draw_profile()

    def _draw_profile(self):
        fig = self.prof_fig
        fig.clear()
        self.prof_table.delete(*self.prof_table.get_children())
        if not self.profiles:
            ax = fig.add_subplot(111)
            ax.axis("off")
            ax.text(0.5, 0.5, "Press Calculate profile: each species gets a smooth concentration-depth profile that "
                              "reproduces its\nintensities at the measured angles. Prior knowledge (constant or absent "
                              "below a depth, non-negative) sharpens it.\nOnly a few angles are measured, so the depth "
                              "resolution is limited: read the noise bands and the notes.", ha="center", va="center",
                    color=MUTED, fontsize=11)
            self.prof_canvas.draw_idle()
            return
        ax1, ax2 = fig.subplots(1, 2, gridspec_kw={"width_ratios": [3, 2]})
        fits = [p[1] for p in self.profiles]
        frac = DR.atomic_fraction(fits) if self.prof_frac.get() and len(fits) > 1 else None
        for i, (sp, f, bulk, zc) in enumerate(self.profiles):
            colour = PALETTE[self.species.index(sp) % len(PALETTE)]
            c = frac[i] if frac is not None else f.c / max(float(np.max(np.abs([p[1].c for p in self.profiles]))), 1e-30)
            ax1.plot(f.z_nm, c, "-", color=colour, lw=1.6, label=sp.name)
            if self.prof_bands.get() and f.c_se is not None and frac is None:
                scale = max(float(np.max(np.abs([p[1].c for p in self.profiles]))), 1e-30)
                ax1.fill_between(f.z_nm, c - f.c_se / scale, c + f.c_se / scale, color=colour, alpha=0.18, lw=0)
            ax2.errorbar(f.angles, f.measured, yerr=f.sigma, fmt="o", color=colour, ms=4, capsize=2)
            ax2.plot(f.angles, f.predicted, "-", color=colour, lw=1.2)
            self.prof_table.insert("", "end", values=(
                sp.name, f"{math.log10(f.alpha):.1f}", _fmt(f.chi2_red, "{:.2f}"), self.BULK_LABELS[bulk] + (
                    f" {zc:g}" if bulk != "free" and zc is not None else ""),
                f"{float(np.max(f.c)):.3g}", "; ".join(n for n in f.notes if "resolution" not in n)))
        ax1.set_xlabel("Depth below the surface (nm)")
        ax1.set_ylabel("Atomic fraction" if frac is not None else "Concentration (relative to the largest)")
        ax1.set_title("Regularised depth profiles", fontsize=10)
        ax1.grid(alpha=0.3)
        ax1.axhline(0, color="black", lw=0.6)
        ax1.legend(fontsize=8, loc="best")
        ax2.set_xlabel("Emission angle θ (°)")
        ax2.set_ylabel("Intensity / sensitivity factor")
        ax2.set_title("Measured (points) and from the profile (lines)", fontsize=10)
        ax2.grid(alpha=0.3)
        fig.subplots_adjust(left=0.07, right=0.98, top=0.9, bottom=0.12, wspace=0.28)
        self.prof_canvas.draw_idle()
        res = [p[1].notes[-1] for p in self.profiles if p[1].notes]
        if res and not self.prof_note.get():
            self.prof_note.set("Note: " + res[0] + ".")

    # --- trend plot ---------------------------------------------------------------------------------------
    def _draw(self):
        fig = self.trend_fig
        fig.clear()
        species = self.active() if self.species else []
        self.trend_table.delete(*self.trend_table.get_children())
        if not species:
            ax = fig.add_subplot(111)
            ax.axis("off")
            ax.text(0.5, 0.5, "Open a VAMAS tilt series: every species is drawn as ln(intensity) against angle.\n"
                              "Species whose lines run parallel (and have similar attenuation lengths) may be "
                              "distributed together.", ha="center", va="center", color=MUTED, fontsize=11)
            self.trend_canvas.draw_idle()
            return
        try:
            max_angle = self._number("max_angle", "Max angle")
        except ValueError:
            max_angle = DP.DEFAULT_MAX_ANGLE
        rows = DP.trend_table(species, max_angle)
        ax = fig.add_subplot(111)
        use_sec = self.x_mode.get() == "sec"
        for i, (sp, tr) in enumerate(zip(species, rows)):
            colour = PALETTE[self.species.index(sp) % len(PALETTE)]
            y = np.asarray(sp.intensity, dtype=float)
            ok = y > 0
            if self.normalise.get():
                y = DP.normalised(sp)
            x = sp.sec if use_sec else sp.angles
            ly = np.log(np.where(y > 0, y, np.nan))
            err = None
            if sp.sigma is not None:
                rel = np.where(sp.intensity > 0, sp.sigma / np.where(sp.intensity > 0, sp.intensity, 1.0), np.nan)
                err = rel
            inside = sp.angles <= max_angle
            ax.errorbar(x[ok & inside], ly[ok & inside], yerr=None if err is None else err[ok & inside], fmt="o",
                        color=colour, ms=5, capsize=2, label=sp.name)
            if (~inside & ok).any():
                ax.plot(x[ok & ~inside], ly[ok & ~inside], "o", mfc="none", color=colour, ms=5)
            if self.show_fits.get() and tr.slope is not None:
                xs = np.linspace(float(np.min(x[ok & inside])), float(np.max(x[ok & inside])), 20)
                if use_sec:
                    line = tr.intercept + tr.slope * xs
                else:
                    line = tr.intercept + tr.slope / np.cos(np.radians(xs))
                if self.normalise.get():
                    first = np.flatnonzero(sp.intensity > 0)[0]
                    line = line - math.log(sp.intensity[first])
                ax.plot(xs, line, "-", color=colour, lw=1)
        ax.set_xlabel("sec θ  (1 / cos θ)" if use_sec else "Emission angle θ (°)")
        ax.set_ylabel("ln( I / I at first angle )" if self.normalise.get() else "ln( I )  (I in counts/s·eV)")
        ax.set_title("Surrey Trend Plot: parallel lines may mean species that are distributed together", fontsize=10)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, loc="best")
        fig.subplots_adjust(left=0.09, right=0.98, top=0.92, bottom=0.12)
        self.trend_canvas.draw_idle()
        for r in rows:
            self.trend_table.insert("", "end", values=(
                r.name, _fmt(r.ke_ev, "{:.0f}"), _fmt(r.slope), _fmt(r.slope_se), _fmt(r.rms), _fmt(r.depth_nm, "{:.2f}"),
                r.n, ", ".join(r.similar) or "-", r.note))
