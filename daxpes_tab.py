"""
The "DAXPES" tab of the XPS Thickness Suite: thickness and coverage of a thin film from one core level measured with a
soft and a hard X-ray source (dual-anode XPS, B. J. Murdoch, Surf. Interface Anal. 58, 125 (2026), peak-area method).

Maths in daxpes.py, peak areas from arxps.py, file reading in vamas.py; this module is only the Tk/matplotlib glue. It is
a ttk.Frame that app.py puts into a notebook tab and needs numpy + matplotlib. Angle convention: theta is from the
surface NORMAL (0 = normal emission). Every area can also be typed, which overrides the file.
"""
import math
import os
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

import arxps as AR
import daxpes as DX
import depthprofile as DP
import tougaard as TGD
import vamas as VM
from gui_common import BADGE, CARD, MUTED, Tooltip, _to_float

SOFT_COLOUR, HARD_COLOUR = "#2563eb", "#d97706"
SLOTS = (("samp_soft", "Film, soft source"), ("samp_hard", "Film, hard source"),
         ("ref_soft", "Thick reference, soft"), ("ref_hard", "Thick reference, hard"))
NO_REGION = "(none)"
CUSTOM = "custom (type the values)"
METHOD_AREA = "Peak areas (Eq. 6)"
METHOD_BG = "Inelastic background B1 (Eqs. 7-8)"
AL_SEAH = "Seah 2012 universal curve (S3)"
AL_TPP = "TPP-2M IMFP x AL/IMFP"


class _Slot:
    """One input spectrum: a VAMAS file, the region used, and the area read from it."""

    def __init__(self, key, title):
        self.key, self.title = key, title
        self.file = None
        self.path = tk.StringVar(value="no file")
        self.region = tk.StringVar(value=NO_REGION)
        self.combo = None
        self.area = self.sigma = None
        self.extraction = None
        self.block = None


class DaxpesTab(ttk.Frame):
    def __init__(self, parent, app, tpp2m, compact=False):
        super().__init__(parent)
        self.app, self.tpp2m, self.compact = app, tpp2m, compact
        self.slots = {key: _Slot(key, title) for key, title in SLOTS}
        self.result = None
        self.align = tk.BooleanVar(value=True)
        self.material = tk.StringVar(value="TiO2")
        self.al_method = tk.StringVar(value=AL_SEAH)
        self.method = tk.StringVar(value=METHOD_AREA)
        self.cross = tk.StringVar(value=next(iter(TGD.CROSS_SECTIONS)))
        self.b1_fits = {}                      # "soft"/"hard" -> tougaard.B1Fit of the last background calculation
        keys = ("be", "w_lo", "w_hi", "theta", "lam_soft", "lam_hard", "hv_soft", "hv_hard", "factor", "i_soft", "i_hard",
                "inf_soft", "inf_hard", "rho", "c", "d", "delta_e", "b1_soft", "b1_hard")
        self.v = {k: tk.StringVar() for k in keys}
        self.v["theta"].set("0")
        self.v["factor"].set("0.9")
        self.v["delta_e"].set(f"{TGD.DEFAULT_DELTA_E:g}")
        self._cross_changed()
        self.v["hv_soft"].set(f"{DX.SOURCES['Al Kα']:g}")
        self.v["hv_hard"].set(f"{DX.SOURCES['Ag Lα']:g}")

        self.columnconfigure(0, weight=0)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)
        left = ttk.Frame(self)
        left.grid(row=0, column=0, sticky="ns", padx=(0, 8), pady=(8, 4))
        left.rowconfigure(1, weight=1)
        self._build_setup(left).grid(row=0, column=0, sticky="ew")
        self._build_results(left).grid(row=1, column=0, sticky="nsew", pady=(8, 0))
        self._build_right().grid(row=0, column=1, sticky="nsew", pady=(8, 4))
        self._draw()

    # --- setup ------------------------------------------------------------------------------------------
    def _entry(self, parent, key, row, label, hint, col=0, width=9):
        ttk.Label(parent, text=label, style="Card.TLabel").grid(row=row, column=col, sticky="w", pady=2)
        e = ttk.Entry(parent, textvariable=self.v[key], width=width)
        e.grid(row=row, column=col + 1, sticky="w", padx=6, pady=2)
        Tooltip(e, hint)
        return e

    def _build_setup(self, parent):
        card = ttk.LabelFrame(parent, text="  Data  ", style="Card.TLabelframe", padding=(10, 2, 10, 4) if self.compact else 10)
        top = ttk.Frame(card, style="Card.TFrame")
        top.grid(row=0, column=0, sticky="ew")
        ttk.Label(top, text="Method:", style="Card.TLabel").pack(side="left")
        mcb = ttk.Combobox(top, textvariable=self.method, values=[METHOD_AREA, METHOD_BG], state="readonly", width=32)
        mcb.pack(side="left", padx=6)
        mcb.bind("<<ComboboxSelected>>", lambda e: self._method_changed())
        Tooltip(mcb, "Peak areas: the film and thick-reference peak areas of each source (Eq. 6), which needs a reference "
                     "for each source. Inelastic background: the Tougaard background parameter B1 of the film spectrum of "
                     "each source (Eqs. 7-8), which needs no reference sample.")
        self.setup_nb = nb = ttk.Notebook(card, style="Sub.TNotebook")
        nb.grid(row=1, column=0, sticky="ew", pady=(6, 0))
        pages = []
        for title in ("Spectra", "Lengths", "Method"):
            page = ttk.Frame(nb, style="Card.TFrame", padding=6)
            nb.add(page, text=title)
            pages.append(page)
        self._build_spectra_page(pages[0])
        self._build_lengths_page(pages[1])
        self._build_method_page(pages[2])
        self._method_changed()
        return card

    def _build_spectra_page(self, card):
        r = 0
        ttk.Label(card, text="One core level, measured with a soft and a hard X-ray source. Spectra come from VAMAS "
                             "files (or type the areas on the Method page).", style="CardMuted.TLabel", wraplength=360,
                  justify="left").grid(row=r, column=0, columnspan=4, sticky="w")
        for key, title in SLOTS:
            r += 1
            slot = self.slots[key]
            b = ttk.Button(card, text=title + "…", command=lambda s=slot: self.open_file(s), width=22)
            b.grid(row=r, column=0, sticky="w", pady=1)
            Tooltip(b, "Open a VAMAS file holding this spectrum. A file may hold both the film and the reference regions "
                       "(open it twice).")
            slot.combo = ttk.Combobox(card, textvariable=slot.region, values=[NO_REGION], state="readonly", width=16)
            slot.combo.grid(row=r, column=1, columnspan=2, sticky="w", padx=4)
            slot.combo.bind("<<ComboboxSelected>>", lambda e: self._clear_result())
            ttk.Label(card, textvariable=slot.path, style="CardMuted.TLabel", wraplength=110).grid(row=r, column=3, sticky="w")
        r += 1
        self._entry(card, "be", r, "Peak BE (eV):", "Binding energy of the peak, to centre the window. Blank: the tallest "
                    "peak of the film spectrum.")
        self._entry(card, "w_lo", r, "Window from:", "Lower binding energy of the peak window (the Shirley background of "
                    "the area method, and the range for the centroid and area of the background method).", col=2, width=7)
        r += 1
        self._entry(card, "w_hi", r, "to:", "Upper binding energy of the window.", col=2, width=7)
        b = ttk.Button(card, text="Auto window", command=self._auto_clicked)
        b.grid(row=r, column=0, sticky="w", pady=2)
        Tooltip(b, "Set the window from the film spectrum of the soft source.")
        r += 1
        c = ttk.Checkbutton(card, text="Follow small energy shifts between spectra (charging)", variable=self.align,
                            style="Card.TCheckbutton")
        c.grid(row=r, column=0, columnspan=4, sticky="w")
        Tooltip(c, "Area method: shift each spectrum onto the film spectrum of the soft source before integrating, so "
                   "one window serves all four.")

    def _build_lengths_page(self, card):
        r = 0
        self._entry(card, "theta", r, "θ from normal (°):", "Take-off angle from the surface NORMAL, the same for both "
                    "sources.")
        r += 1
        ttk.Label(card, text="Attenuation lengths (nm)", style="ColHead.TLabel").grid(row=r, column=0, columnspan=4, sticky="w")
        r += 1
        self._entry(card, "lam_soft", r, "λ soft:", "Attenuation length of the photoelectrons excited by the SOFT source in "
                    "the film.")
        self._entry(card, "lam_hard", r, "λ hard:", "Attenuation length for the HARD source (must be the longer one).",
                    col=2, width=7)
        r += 1
        self._entry(card, "hv_soft", r, "hν soft (eV):", "Soft photon energy; read from the file when it is opened.")
        self._entry(card, "hv_hard", r, "hν hard:", "Hard photon energy (Ag Lα 2984.3, Ga Kα 9252).", col=2, width=7)
        r += 1
        cb = ttk.Combobox(card, textvariable=self.material, values=list(DX.FILM_MATERIALS) + [CUSTOM], state="readonly",
                          width=12)
        cb.grid(row=r, column=0, sticky="w", pady=2)
        Tooltip(cb, "Film material for the starting attenuation lengths (handbook densities, compositions and heats of "
                    "formation; the paper's own values are not used).")
        mcb = ttk.Combobox(card, textvariable=self.al_method, values=[AL_SEAH, AL_TPP], state="readonly", width=28)
        mcb.grid(row=r, column=1, columnspan=3, sticky="w", padx=4)
        Tooltip(mcb, "Seah (2012, Eq. 11): an attenuation length from the kinetic energy, the mean atomic number, the "
                     "atomic size and the heat of formation, RMS 8 % against computed values. TPP-2M: the inelastic mean "
                     "free path times the AL/IMFP factor.")
        r += 1
        self._entry(card, "factor", r, "AL/IMFP:", "Only for the TPP-2M method: the attenuation length as a fraction of "
                    "the IMFP (1 = the IMFP, which is too long; 0.85-0.95 matches Murdoch's Table 1).", col=0, width=6)
        b = ttk.Button(card, text="Fill λ", command=self.fill_lambdas)
        b.grid(row=r, column=3, sticky="w")
        Tooltip(b, "Estimate both attenuation lengths at the kinetic energy each source gives the peak. Replace them "
                   "with values you trust; the thickness scales with them.")

    def _build_method_page(self, card):
        r = 0
        ttk.Label(card, text="Inelastic background (B1) method", style="ColHead.TLabel").grid(row=r, column=0, columnspan=4,
                                                                                            sticky="w")
        r += 1
        cb = ttk.Combobox(card, textvariable=self.cross, values=list(TGD.CROSS_SECTIONS) + [CUSTOM], state="readonly",
                          width=34)
        cb.grid(row=r, column=0, columnspan=4, sticky="w", pady=2)
        cb.bind("<<ComboboxSelected>>", lambda e: self._cross_changed())
        Tooltip(cb, "The parameters C and D of Tougaard's universal inelastic cross-section K(T) = T/((C-T²)²+D T²) "
                    "(Eq. 2): Murdoch used C = -1643, D = 0 for metal oxides and C = 325, D = 96 for silicon.")
        r += 1
        self._entry(card, "c", r, "C (eV²):", "Cross-section parameter C of Eq. 2 (eV²).", width=8)
        self._entry(card, "d", r, "D (eV²):", "Cross-section parameter D of Eq. 2 (eV²).", col=2, width=7)
        r += 1
        self._entry(card, "delta_e", r, "ΔE (eV):", "B1 is the value that makes the primary spectrum zero this far above "
                    "the peak in binding energy (15-30 eV; 30 in the paper). The spectrum must extend at least that far, "
                    "and should start on a flat zero level below the peak.", width=8)
        r += 1
        self._entry(card, "b1_soft", r, "B1 soft (eV²):", "Type B1 instead of fitting it from the spectrum (for example "
                    "from other software).", width=8)
        self._entry(card, "b1_hard", r, "B1 hard:", "B1 of the hard-source spectrum.", col=2, width=7)
        r += 1
        ttk.Separator(card).grid(row=r, column=0, columnspan=4, sticky="ew", pady=4)
        r += 1
        ttk.Label(card, text="Peak areas typed here override the files (blank = from the spectra); a typed ratio "
                             "overrides the reference spectra", style="ColHead.TLabel", wraplength=360,
                  justify="left").grid(row=r, column=0, columnspan=4, sticky="w")
        r += 1
        self._entry(card, "i_soft", r, "I film soft:", "Peak area of the film with the soft source.")
        self._entry(card, "i_hard", r, "I film hard:", "Peak area of the film with the hard source.", col=2, width=7)
        r += 1
        self._entry(card, "inf_soft", r, "I∞ soft:", "Peak area of the thick reference with the soft source.")
        self._entry(card, "inf_hard", r, "I∞ hard:", "Peak area of the thick reference with the hard source.", col=2, width=7)
        r += 1
        self._entry(card, "rho", r, "or I∞H/I∞X:", "Ratio of the thick-reference areas, hard over soft, if you have no "
                    "reference spectra (area method). Then only the thickness comes out, not the coverage.")

    def _cross_changed(self):
        cd = TGD.CROSS_SECTIONS.get(self.cross.get())
        if cd is not None:
            self.v["c"].set(f"{cd[0]:g}")
            self.v["d"].set(f"{cd[1]:g}")

    def _method_changed(self):
        self._clear_result()

    def _build_results(self, parent):
        card = ttk.LabelFrame(parent, text="  Calculation  ", style="Card.TLabelframe",
                              padding=(10, 2, 10, 4) if self.compact else 10)
        card.columnconfigure(0, weight=1)
        card.rowconfigure(3, weight=1)
        self.calc_btn = ttk.Button(card, text="Calculate & plot", style="Primary.TButton", command=self.calculate)
        self.calc_btn.grid(row=0, column=0, sticky="ew")
        Tooltip(self.calc_btn, "Integrate the peaks, solve Murdoch's Eq. 6 for the thickness and find the coverage")
        readout = ttk.Frame(card, style="Card.TFrame")
        readout.grid(row=1, column=0, sticky="ew", pady=(8, 0))
        self.readout_var = tk.StringVar(value="-.-- nm")
        self.readout_sub_var = tk.StringVar(value="DAXPES thickness will appear here")
        tk.Label(readout, textvariable=self.readout_var, bg=CARD, fg="#b8bec7", anchor="w",
                 font=("Segoe UI", 20 if self.compact else 30, "bold")).pack(anchor="w")
        tk.Label(readout, textvariable=self.readout_sub_var, bg=CARD, fg=MUTED, anchor="w", font=("Segoe UI", 9),
                 wraplength=340, justify="left").pack(anchor="w")
        self.badge = tk.Label(card, text="", font=("Segoe UI", 9, "bold"), padx=8, pady=4, anchor="w", justify="left",
                              wraplength=330)
        self.badge.grid(row=2, column=0, sticky="ew", pady=(8, 8))
        self._set_badge("idle", "Open the spectra (or type the areas), then calculate")
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

    def _build_right(self):
        frame = ttk.Frame(self)
        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(0, weight=1)
        nb = ttk.Notebook(frame, style="Sub.TNotebook")
        nb.grid(row=0, column=0, sticky="nsew")
        self.figs = {}
        for key, title in (("result", "Thickness and coverage"), ("spectra", "Spectra")):
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
        return frame

    def _set_badge(self, kind, text):
        bg, fg = BADGE[kind]
        self.badge.configure(text=text, bg=bg, fg=fg)

    # --- files and windows --------------------------------------------------------------------------------
    def open_file(self, slot):
        path = filedialog.askopenfilename(title=f"Open the spectrum: {slot.title}",
                                          filetypes=[("VAMAS files", "*.vms *.vamas"), ("All files", "*.*")])
        if path:
            self.load_file(slot.key, path)

    def load_file(self, key, path, region=None):
        slot = self.slots[key]
        try:
            vf = VM.read(path)
        except VM.VamasError as exc:
            messagebox.showerror("Cannot read the file", str(exc))
            return False
        if not vf.blocks:
            messagebox.showerror("Cannot read the file", "The file holds no spectra.")
            return False
        slot.file = vf
        slot.path.set(os.path.basename(path))
        regions = vf.regions()
        slot.combo.configure(values=regions)
        if region in regions:
            slot.region.set(region)
        else:                                       # film files: first region; reference slots: a region named "ref"
            pick = next((r for r in regions if ("ref" in r.lower()) == key.startswith("ref")), regions[0])
            slot.region.set(pick)
        hv = vf.blocks_of(slot.region.get())[0].hv_ev
        if hv and key.endswith("soft"):
            self.v["hv_soft"].set(f"{hv:g}")
        elif hv and key.endswith("hard"):
            self.v["hv_hard"].set(f"{hv:g}")
        self._clear_result()
        return True

    def _block(self, slot):
        if slot.file is None or slot.region.get() == NO_REGION:
            return None
        blocks = slot.file.blocks_of(slot.region.get())
        return blocks[0] if blocks else None

    def _spectrum(self, blk):
        return np.asarray(blk.be_ev, dtype=float), AR.intensity(blk, True, True)

    def auto_window(self):
        blk = self._block(self.slots["samp_soft"]) or self._block(self.slots["samp_hard"])
        if blk is None:
            raise ValueError("Open a film spectrum first.")
        be, y = self._spectrum(blk)
        nominal = _to_float(self.v["be"].get())
        if nominal is None:
            nominal = float(be[int(np.argmax(AR.smooth(y)))])
            self.v["be"].set(f"{nominal:.2f}")
        win = AR.auto_peak_window([(be, y)], nominal, below_ev=4.0, above_ev=4.0)
        self.v["w_lo"].set(f"{win.lo:.2f}")
        self.v["w_hi"].set(f"{win.hi:.2f}")

    def _auto_clicked(self):
        try:
            self.auto_window()
        except ValueError as exc:
            messagebox.showerror("Input error", str(exc))

    def fill_lambdas(self):
        try:
            material = self.material.get()
            if material == CUSTOM:
                raise ValueError("Pick a film material for the TPP-2M estimate, or type the attenuation lengths.")
            be = self._num("be", "Peak BE", positive=True)
            factor = self._num("factor", "AL/IMFP", positive=True) if self.al_method.get() == AL_TPP else 1.0
            for key, hv_key in (("lam_soft", "hv_soft"), ("lam_hard", "hv_hard")):
                ke = self._num(hv_key, "Photon energy", positive=True) - be
                if ke <= 50:
                    raise ValueError("The photon energy must exceed the binding energy by more than 50 eV.")
                lam = (DX.seah_start_nm(material, ke) if self.al_method.get() == AL_SEAH
                       else DX.tpp2m_start_nm(self.tpp2m, material, ke, factor))
                self.v[key].set(f"{lam:.3g}")
        except ValueError as exc:
            messagebox.showerror("Input error", str(exc))

    def _clear_result(self):
        self.result = None
        if hasattr(self, "badge"):
            self._set_badge("idle", "Inputs changed: press Calculate & plot")

    # --- calculation --------------------------------------------------------------------------------------
    def _num(self, key, label, positive=False):
        x = _to_float(self.v[key].get())
        if x is None or not math.isfinite(x):
            raise ValueError(f"Please enter a number for '{label}'.")
        if positive and x <= 0:
            raise ValueError(f"'{label}' must be greater than zero.")
        return x

    def _optional(self, key, label):
        text = self.v[key].get().strip()
        return None if not text else self._num(key, label, positive=True)

    def calculate(self):
        try:
            self._calculate()
        except (ValueError, VM.VamasError) as exc:
            messagebox.showerror("Input error", str(exc))

    def _extract(self, slot, win, reference):
        """Area and counting sigma of one spectrum (None, None when the slot is empty)."""
        blk = self._block(slot)
        slot.block, slot.extraction = blk, None
        if blk is None:
            return None, None
        be, y = self._spectrum(blk)
        ex = AR.single_peak_area(be, y, win, reference=reference, align=self.align.get() and reference is not None)
        slot.extraction = ex
        tf = AR._transmission_factor(blk, blk.hv_ev - 0.5 * (win.lo + win.hi))
        sig = DP.area_sigma(blk, blk.hv_ev - win.hi, blk.hv_ev - win.lo, tf)
        return ex.area_sub, sig

    def _calculate(self):
        theta = self._num("theta", "θ")
        if not 0.0 <= theta < 90.0:
            raise ValueError("θ must be from 0 up to (not including) 90 degrees.")
        lam_s, lam_h = self._num("lam_soft", "λ soft", True), self._num("lam_hard", "λ hard", True)
        if self.method.get() == METHOD_BG:
            return self._calculate_background(theta, lam_s, lam_h)
        if not self.v["w_lo"].get() or not self.v["w_hi"].get():
            if any(self._block(s) is not None for s in self.slots.values()):
                self.auto_window()
        have_file = any(self._block(s) is not None for s in self.slots.values())
        if have_file:
            win = AR.Window(self._num("w_lo", "Window from"), self._num("w_hi", "Window to"))
            if win.hi <= win.lo:
                raise ValueError("The window must run from a lower to a higher binding energy.")
            ref_blk = self._block(self.slots["samp_soft"])
            reference = self._spectrum(ref_blk) if ref_blk is not None else None
        for key, _ in SLOTS:
            slot = self.slots[key]
            slot.area = slot.sigma = None
            if have_file:
                slot.area, slot.sigma = self._extract(slot, win, None if key == "samp_soft" else reference)
        typed = {"samp_soft": self._optional("i_soft", "I film soft"), "samp_hard": self._optional("i_hard", "I film hard"),
                 "ref_soft": self._optional("inf_soft", "I∞ soft"), "ref_hard": self._optional("inf_hard", "I∞ hard")}
        vals, sigs = {}, {}
        for key, _ in SLOTS:
            slot = self.slots[key]
            if typed[key] is not None:
                vals[key], sigs[key] = typed[key], None
            else:
                vals[key], sigs[key] = slot.area, slot.sigma
        if vals["samp_soft"] is None or vals["samp_hard"] is None:
            raise ValueError("Give the film peak area for both sources: open the spectra or type 'I film soft' and "
                             "'I film hard'.")
        rho = self._optional("rho", "I∞H/I∞X")
        typed_inf = typed["ref_soft"] is not None and typed["ref_hard"] is not None
        if rho is not None and not typed_inf:             # a typed ratio overrides reference spectra
            vals["ref_soft"] = vals["ref_hard"] = None
        both_inf = vals["ref_soft"] is not None and vals["ref_hard"] is not None
        if not both_inf and rho is None:
            raise ValueError("Give the thick-reference area of both sources, or their ratio I∞H/I∞X: without it the "
                             "two sources cannot be compared.")
        res = DX.analyse(vals["samp_soft"], vals["samp_hard"], lam_s, lam_h, theta,
                         iinf_soft=vals["ref_soft"] if both_inf else None,
                         iinf_hard=vals["ref_hard"] if both_inf else None, rho=None if both_inf else rho,
                         sig_soft=sigs["samp_soft"], sig_hard=sigs["samp_hard"],
                         sig_inf_soft=sigs["ref_soft"] if both_inf else None,
                         sig_inf_hard=sigs["ref_hard"] if both_inf else None)
        self.result = (res, lam_s, lam_h, theta, vals if both_inf else None)
        self._present(res)
        self._draw()

    def _calculate_background(self, theta, lam_s, lam_h):
        """Murdoch's background method: B1 of the film spectrum of each source, then Eqs. 7 and 8."""
        c, d = self._num("c", "C"), self._num("d", "D")
        if d < 0:
            raise ValueError("D must not be negative.")
        delta_e = self._num("delta_e", "ΔE", positive=True)
        typed = {"soft": self._optional("b1_soft", "B1 soft"), "hard": self._optional("b1_hard", "B1 hard")}
        for slot in self.slots.values():
            slot.extraction = None
        self.b1_fits = {}
        any_file = any(self._block(sl) is not None for sl in self.slots.values())
        if any_file and (not self.v["w_lo"].get() or not self.v["w_hi"].get()):
            self.auto_window()
        win = None
        if self.v["w_lo"].get() and self.v["w_hi"].get():
            win = (self._num("w_lo", "Window from"), self._num("w_hi", "Window to"))
            if win[1] <= win[0]:
                raise ValueError("The window must run from a lower to a higher binding energy.")
        b1, sig = {}, {}
        for src in ("soft", "hard"):
            if typed[src] is not None:
                b1[src], sig[src] = typed[src], None
                continue
            blk = self._block(self.slots["samp_" + src])
            if blk is None:
                raise ValueError(f"Open the film spectrum measured with the {src} source (or type 'B1 {src}').")
            if win is None:
                raise ValueError("Set the peak window (Auto window) first.")
            be, y = self._spectrum(blk)
            fit = TGD.fit_b1(be, y, win, c, d, delta_e=delta_e)
            tf = AR._transmission_factor(blk, blk.hv_ev - 0.5 * (win[0] + win[1]))
            scale = 1.0 / ((blk.scans or 1) * (blk.dwell_s or 1.0)) / max(tf, 1e-12)
            b1[src] = fit.b1
            sig[src] = TGD.b1_sigma(be, blk.counts, scale, win, c, d, delta_e=delta_e, draws=40)
            self.b1_fits[src] = fit
        # the coverage needs the peak areas of the film and a thick reference (files or typed); it is optional here
        vals = {}
        typed_a = {"samp_soft": self._optional("i_soft", "I film soft"), "samp_hard": self._optional("i_hard", "I film hard"),
                   "ref_soft": self._optional("inf_soft", "I∞ soft"), "ref_hard": self._optional("inf_hard", "I∞ hard")}
        ref_blk = self._block(self.slots["samp_soft"])
        reference = self._spectrum(ref_blk) if ref_blk is not None else None
        for key, _ in SLOTS:
            if typed_a[key] is not None:
                vals[key] = typed_a[key]
            elif win is not None and self._block(self.slots[key]) is not None:
                try:
                    vals[key] = self._extract(self.slots[key], AR.Window(*win), None if key == "samp_soft" else reference)[0]
                except ValueError:
                    vals[key] = None
        res = DX.analyse_background(b1["soft"], b1["hard"], lam_s, lam_h, theta, sig.get("soft"), sig.get("hard"),
                                    i_soft=vals.get("samp_soft"), iinf_soft=vals.get("ref_soft"),
                                    i_hard=vals.get("samp_hard"), iinf_hard=vals.get("ref_hard"))
        self.result = (res, lam_s, lam_h, theta, None)
        self._present(res)
        self._draw()

    def _present_background(self, res):
        err = f" ± {res.d_se_nm:.2f}" if res.d_se_nm else ""
        self.readout_var.set(f"{res.d_nm:.2f} nm")
        self.readout_sub_var.set(
            f"background method, Eq. 8{err}; B1 {res.b1_soft:.0f} / {res.b1_hard:.0f} eV²" +
            (f"; coverage {100 * res.coverage:.0f} %" if res.coverage is not None else "; no coverage (no reference areas)"))
        self._set_badge("warn" if res.warnings else "ok", res.warnings[0] if res.warnings else
                        "The B1 ratio lies inside the range the two sources can resolve")
        lines = [f"==> THICKNESS d = {res.d_nm:.3f} nm  ({res.d_nm * 10:.2f} Angstrom)" + (f" ± {res.d_se_nm:.3f}" if res.d_se_nm else ""),
                 "", "Murdoch (2026) background method: B1 of each source's Tougaard background,",
                 f"  B1 soft = {res.b1_soft:.1f} eV², B1 hard = {res.b1_hard:.1f} eV²,  B1_H/B1_X = {res.q:.4f}",
                 f"  Eq. 7:  d* = {res.d_star_nm:.3f} nm  (d* = d / (2 cos θ) when elastic scattering is ignored)",
                 f"  Eq. 8:  d = 3 λ_hard cos θ (1 - exp(-d*/(λ_hard cos θ))) = {res.d_nm:.3f} nm",
                 f"  simple thin-film reading 2 d* cos θ = {res.d_simple_nm:.3f} nm (no 3λ cap)"]
        if res.coverage is not None:
            lines.append(f"coverage = (I/I∞) / (1 - exp(-d/λcosθ)) = {100 * res.coverage:.1f} %")
        if res.t_soft_nm is not None:
            lines.append(f"single-source peak-area thickness, soft: {res.t_soft_nm:.3f} nm" +
                         (f", hard: {res.t_hard_nm:.3f} nm" if res.t_hard_nm is not None else ""))
        lines.append(f"elastic limit 3 λ_hard cos θ = {res.limit_nm:.2f} nm")
        for key, text in (("ratio", "λ_hard ±10 % (the λ ratio)"), ("common", "both λ ±10 %")):
            band = res.lam_band.get(key)
            if band is None:
                lines.append(f"if {text}: no solution for one of the two")
            else:
                cov = (f", coverage {100 * band[2]:.0f} to {100 * band[3]:.0f} %" if band[2] is not None else "")
                lines.append(f"if {text}: d = {band[0]:.2f} to {band[1]:.2f} nm" + cov)
        lines += [""] + [f"! {w}" for w in res.warnings]
        lines += ["", "B1 does not depend on the coverage, the flux or the cross-section of either source; the uncertainty "
                      "is the counting-statistics part only (Monte Carlo)."]
        self.details.configure(state="normal")
        self.details.delete("1.0", "end")
        self.details.insert("1.0", "\n".join(lines))
        self.details.configure(state="disabled")

    def _present(self, res):
        if res.method == "background":
            return self._present_background(res)
        err = f" ± {res.d_se_nm:.2f}" if res.d_se_nm else ""
        self.readout_var.set(f"{res.d_nm:.2f} nm")
        parts = [f"DAXPES thickness{err}", "independent of the coverage"]
        if res.coverage is not None:
            parts = [f"coverage {100 * res.coverage:.0f} %" + (f" ± {100 * res.coverage_se:.0f}" if res.coverage_se else "")]
        self.readout_sub_var.set("; ".join(parts) if res.coverage is not None
                                 else "DAXPES thickness, independent of the coverage (no reference spectra: no coverage)")
        self._set_badge("warn" if res.warnings else "ok", res.warnings[0] if res.warnings else
                        "Both sources agree on a film inside the resolvable range")
        lines = [f"==> THICKNESS d = {res.d_nm:.3f} nm  ({res.d_nm * 10:.2f} Angstrom)" + (f" ± {res.d_se_nm:.3f}" if res.d_se_nm else ""),
                 "", "Murdoch (2026) Eq. 6: the coverage cancels in the ratio of the two sources,",
                 f"  (I_H/I∞H)/(I_X/I∞X) = {res.q:.4f}"]
        if res.t_soft_nm is not None:
            lines.append(f"single-source thickness, soft: {res.t_soft_nm:.3f} nm" +
                         (f", hard: {res.t_hard_nm:.3f} nm" if res.t_hard_nm is not None else ""))
        if res.coverage is not None:
            lines.append(f"coverage = (I/I∞) / (1 - exp(-d/λcosθ)) = {100 * res.coverage:.1f} %" +
                         (f" ± {100 * res.coverage_se:.1f}" if res.coverage_se else ""))
        if res.coverage_paper is not None:
            lines.append(f"coverage as the paper's thickness ratio (hard / DAXPES): {100 * res.coverage_paper:.1f} %")
        lines.append(f"elastic limit 3 λ_hard cos θ = {res.limit_nm:.2f} nm")
        for key, text in (("ratio", "λ_hard ±10 % (the λ ratio)"), ("common", "both λ ±10 %")):
            band = res.lam_band.get(key)
            if band is None:
                lines.append(f"if {text}: no solution for one of the two")
            else:
                cov = (f", coverage {100 * band[2]:.0f} to {100 * band[3]:.0f} %" if band[2] is not None else "")
                lines.append(f"if {text}: d = {band[0]:.2f} to {band[1]:.2f} nm" + cov)
        lines += [""] + [f"! {w}" for w in res.warnings]
        lines += ["", "Both sources must see the SAME analysed area; the areas of each source are normalised only to "
                      "the thick reference of the same source."]
        self.details.configure(state="normal")
        self.details.delete("1.0", "end")
        self.details.insert("1.0", "\n".join(lines))
        self.details.configure(state="disabled")

    # --- plots --------------------------------------------------------------------------------------------
    def _placeholder(self, fig, text):
        fig.clear()
        ax = fig.add_subplot(111)
        ax.axis("off")
        ax.text(0.5, 0.5, text, ha="center", va="center", color=MUTED, fontsize=11)

    def _draw(self):
        res = self.result
        for key, (fig, canvas) in self.figs.items():
            if res is None:
                self._placeholder(fig, "Open the soft- and hard-source spectra of the film and of a thick reference\n"
                                       "(or type the areas), set the attenuation lengths, then press Calculate & plot.")
            elif key == "result":
                (self._draw_result_bg if res[0].method == "background" else self._draw_result)(fig, res)
            else:
                (self._draw_spectra_bg if res[0].method == "background" else self._draw_spectra)(fig)
            canvas.draw_idle()

    def _draw_result(self, fig, packed):
        res, lam_s, lam_h, theta, vals = packed
        fig.clear()
        ax1, ax2 = fig.subplots(1, 2)
        dmax = max(2.5 * res.d_nm, 3.5 * lam_s * math.cos(math.radians(theta)), 1.0)
        d = np.linspace(1e-3 * dmax, dmax, 400)
        fs = np.array([DX.film_factor(x, lam_s, theta) for x in d])
        fh = np.array([DX.film_factor(x, lam_h, theta) for x in d])
        ax1.plot(d, fs, color=SOFT_COLOUR, label=f"soft (λ {lam_s:.3g} nm), complete film")
        ax1.plot(d, fh, color=HARD_COLOUR, label=f"hard (λ {lam_h:.3g} nm), complete film")
        if vals is not None:
            rs, rh = vals["samp_soft"] / vals["ref_soft"], vals["samp_hard"] / vals["ref_hard"]
            ax1.axhline(rs, color=SOFT_COLOUR, ls=":", lw=1)
            ax1.axhline(rh, color=HARD_COLOUR, ls=":", lw=1)
            ax1.plot([res.d_nm], [rs], "o", color=SOFT_COLOUR)
            ax1.plot([res.d_nm], [rh], "o", color=HARD_COLOUR)
            ax1.plot(d, res.coverage * fs, color=SOFT_COLOUR, lw=0.8, ls="--")
            ax1.plot(d, res.coverage * fh, color=HARD_COLOUR, lw=0.8, ls="--")
        ax1.axvline(res.d_nm, color="black", lw=1)
        ax1.axvline(res.limit_nm, color="#6b7280", lw=0.8, ls="-.")
        ax1.set_xlabel("Film thickness (nm)")
        ax1.set_ylabel("Peak area / thick-reference area")
        ax1.set_title("Each source (solid: complete film; dashed: at the found coverage)", fontsize=9)
        ax1.set_ylim(0, 1.05)
        ax1.grid(alpha=0.3)
        ax1.legend(fontsize=7, loc="lower right")
        q = np.divide(fh, fs, out=np.full_like(fh, np.nan), where=fs > 0)
        ax2.plot(d, q, color="#374151")
        ax2.axhline(res.q, color="#dc2626", lw=1)
        ax2.axvline(res.d_nm, color="black", lw=1)
        ax2.axhline(lam_s / lam_h, color="#6b7280", lw=0.8, ls=":")
        ax2.axvline(res.limit_nm, color="#6b7280", lw=0.8, ls="-.")
        ax2.set_ylim(0, 1.05)
        ax2.set_xlabel("Film thickness (nm)")
        ax2.set_ylabel("(I_H/I∞H) / (I_X/I∞X)")
        ax2.set_title("Eq. 6: the ratio fixes d whatever the coverage", fontsize=9)
        ax2.grid(alpha=0.3)
        ax2.text(0.98, 0.04, "dotted: thin-film limit λ_soft/λ_hard\ndash-dot: elastic limit 3λ_hard cosθ", transform=ax2.transAxes,
                 ha="right", va="bottom", fontsize=7, color=MUTED)
        fig.subplots_adjust(left=0.08, right=0.98, top=0.9, bottom=0.12, wspace=0.28)

    def _draw_result_bg(self, fig, packed):
        res, lam_s, lam_h, theta, _ = packed
        c = math.cos(math.radians(theta))
        fig.clear()
        ax1, ax2 = fig.subplots(1, 2)
        dmax = max(2.0 * res.d_star_nm, 3.5 * lam_h, 1.0)
        ds = np.linspace(1e-3 * dmax, dmax, 400)
        q = (1 - np.exp(-ds / lam_h)) / (1 - np.exp(-ds / lam_s))
        ax1.plot(ds, q, color="#374151")
        ax1.axhline(res.q, color="#dc2626", lw=1)
        ax1.axvline(res.d_star_nm, color="black", lw=1)
        ax1.axhline(lam_s / lam_h, color="#6b7280", lw=0.8, ls=":")
        ax1.set_ylim(0, 1.05)
        ax1.set_xlabel("d* (nm)")
        ax1.set_ylabel("B1(hard) / B1(soft)")
        ax1.set_title("Eq. 7: the B1 ratio fixes d*", fontsize=9)
        ax1.grid(alpha=0.3)
        eq8 = 3.0 * lam_h * c * (1 - np.exp(-ds / (lam_h * c)))
        ax2.plot(ds, eq8, color=HARD_COLOUR, label="Eq. 8: 3 λ_hard cosθ (1 - exp(-d*/(λ_hard cosθ)))")
        ax2.plot(ds, 2 * ds * c, color="#6b7280", ls="--", lw=1, label="2 d* cosθ (thin film, no cap)")
        ax2.axhline(res.limit_nm, color="#6b7280", lw=0.8, ls="-.")
        ax2.plot([res.d_star_nm], [res.d_nm], "o", color="black")
        ax2.set_ylim(0, max(1.2 * res.limit_nm, 1.2 * res.d_nm))
        ax2.set_xlabel("d* (nm)")
        ax2.set_ylabel("Film thickness d (nm)")
        ax2.set_title("Eq. 8: d from d* (capped at 3λ)", fontsize=9)
        ax2.grid(alpha=0.3)
        ax2.legend(fontsize=7, loc="lower right")
        fig.subplots_adjust(left=0.08, right=0.98, top=0.9, bottom=0.12, wspace=0.28)

    def _draw_spectra_bg(self, fig):
        fig.clear()
        if not self.b1_fits:
            self._placeholder(fig, "B1 values were typed in: there are no spectra to show.")
            return
        axes = fig.subplots(2, len(self.b1_fits), squeeze=False)
        for col, (src, fit) in enumerate(self.b1_fits.items()):
            colour = SOFT_COLOUR if src == "soft" else HARD_COLOUR
            for row in range(2):
                ax = axes[row][col]
                ax.plot(fit.be, fit.j, color=colour, lw=0.7, label="spectrum")
                ax.plot(fit.be, fit.background, color="black", lw=1.0, label=f"B1 = {fit.b1:.0f} eV²")
                ax.plot(fit.be, fit.f, color="#16a34a", lw=0.7, label="primary F")
                ax.axvline(fit.centroid_be, color="#6b7280", lw=0.7, ls=":")
                ax.axvline(fit.eval_be, color="#dc2626", lw=0.8, ls="--")
                ax.set_xlim(fit.be[-1], fit.be[0])
                ax.grid(alpha=0.3)
                ax.tick_params(labelsize=7)
                if row == 1:
                    top = max(3.0 * float(np.interp(fit.eval_be, fit.be, fit.background)), 1e-9)
                    ax.set_ylim(-0.5 * top, 2 * top)
                    ax.set_xlabel("Binding energy (eV)", fontsize=8)
                else:
                    ax.set_title(f"{src} source: F = 0 at {fit.eval_be:.1f} eV", fontsize=8)
                    ax.legend(fontsize=6, loc="upper left")
        fig.subplots_adjust(left=0.09, right=0.99, top=0.92, bottom=0.1, wspace=0.25, hspace=0.18)

    def _draw_spectra(self, fig):
        fig.clear()
        have = [(key, title) for key, title in SLOTS if self.slots[key].extraction is not None]
        if not have:
            self._placeholder(fig, "Areas were typed in: there are no spectra to show.")
            return
        axes = fig.subplots(1, len(have), squeeze=False)[0]
        for ax, (key, title) in zip(axes, have):
            ex = self.slots[key].extraction
            colour = SOFT_COLOUR if key.endswith("soft") else HARD_COLOUR
            ax.plot(ex.grid, ex.y, color=colour, lw=0.9)
            ax.plot(ex.grid, ex.bg, color="black", lw=0.8)
            ax.axvspan(ex.lo, ex.hi, color=colour, alpha=0.08)
            ax.set_xlim(ex.grid[-1], ex.grid[0])
            ax.set_title(f"{title}\narea {ex.area_sub:.4g}", fontsize=8)
            ax.set_xlabel("Binding energy (eV)", fontsize=8)
            ax.tick_params(labelsize=7)
            ax.grid(alpha=0.3)
        fig.subplots_adjust(left=0.07, right=0.99, top=0.84, bottom=0.14, wspace=0.3)
