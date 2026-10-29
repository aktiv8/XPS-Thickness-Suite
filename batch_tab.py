"""
The "Batch (VAMAS)" tab of the XPS Thickness Suite: a Strohmeier thickness for every sample of a CasaXPS-fitted VAMAS
file (e.g. Al 2p fitted as metal + oxide, one block per sample).

Maths in batch.py, file reading in vamas.py / casa.py; this module is only the Tk/matplotlib glue. It is a ttk.Frame
that app.py puts into a notebook tab, and needs numpy + matplotlib (app.py shows an install hint if they are missing).
"""
import csv
import os
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

import batch as B
import vamas as VM
from gui_common import ACCENT, BADGE, CARD, MUTED, Tooltip, _to_float

ROLE_TEXT = {B.SUBSTRATE: "Substrate", B.OVERLAYER: "Overlayer", B.IGNORE: "Ignore"}
TEXT_ROLE = {v: k for k, v in ROLE_TEXT.items()}
FIXED_ANGLE = "One fixed angle"
OK_COLOUR, WARN_COLOUR = "#2563eb", "#d97706"


def _norm(text):
    return "".join(str(text).lower().split())


class BatchTab(ttk.Frame):
    def __init__(self, parent, app, substrates, overlayers, volume_density, compact=False):
        super().__init__(parent)
        self.app, self.substrates, self.overlayers = app, substrates, overlayers
        self.volume_density, self.compact = volume_density, compact
        self.file = None
        self.samples, self.states, self.rows = [], [], []
        self.role_vars, self.mat_vars = {}, {}           # state key -> StringVar
        self.region = tk.StringVar()
        self.theta = tk.StringVar(value="90")
        self.angle_src = tk.StringVar(value=FIXED_ANGLE)
        self.offset = tk.StringVar(value="0")
        self._src_by_label = {}
        self.sub_preset = tk.StringVar()
        self.file_label = tk.StringVar(value="No file loaded")
        self.info = tk.StringVar(value="")
        self.summary = tk.StringVar(value="")

        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)
        left = ttk.LabelFrame(self, text="  Setup  ", style="Card.TLabelframe", padding=(8, 4, 8, 8))
        left.grid(row=0, column=0, sticky="ns", padx=(0, 8), pady=(8, 4))
        self._build_setup(left)
        self._build_right().grid(row=0, column=1, sticky="nsew", pady=(8, 4))
        self._draw()

    # --- setup ----------------------------------------------------------------------------------------
    def _build_setup(self, f):
        f.columnconfigure(0, weight=1)
        top = ttk.Frame(f, style="Card.TFrame")
        top.grid(row=0, column=0, sticky="ew")
        btn = ttk.Button(top, text="Open VAMAS file…", command=self.open_file)
        btn.pack(side="left")
        Tooltip(btn, "A VAMAS (.vms) file saved from CasaXPS with the peak fits in it: one block per sample, the same "
                     "core level in each, fitted with the same components (e.g. Al 2p Metal + Al 2p Oxide).")
        ttk.Label(top, textvariable=self.file_label, style="CardMuted.TLabel", wraplength=240,
                  justify="left").pack(side="left", padx=8)

        grid = ttk.Frame(f, style="Card.TFrame")
        grid.grid(row=1, column=0, sticky="ew", pady=(8, 0))
        ttk.Label(grid, text="Region:", style="Card.TLabel").grid(row=0, column=0, sticky="w", pady=2)
        self.region_cb = ttk.Combobox(grid, textvariable=self.region, state="readonly", width=24)
        self.region_cb.grid(row=0, column=1, sticky="w", padx=6, pady=2)
        self.region_cb.bind("<<ComboboxSelected>>", lambda e: self._region_changed())
        Tooltip(self.region_cb, "The core level to analyse. Only regions with a stored CasaXPS fit are listed.")
        ttk.Label(grid, text="Take-off angle (°):", style="Card.TLabel").grid(row=1, column=0, sticky="w", pady=2)
        e = ttk.Entry(grid, textvariable=self.theta, width=7)
        e.grid(row=1, column=1, sticky="w", padx=6, pady=2)
        self.theta_entry = e
        Tooltip(e, "Measured from the SURFACE, 90 = normal emission (the K-Alpha at zero tilt).")
        ttk.Label(grid, text="Angle from:", style="Card.TLabel").grid(row=2, column=0, sticky="w", pady=2)
        self.angle_cb = ttk.Combobox(grid, textvariable=self.angle_src, values=[FIXED_ANGLE], state="readonly",
                                     width=24)
        self.angle_cb.grid(row=2, column=1, sticky="w", padx=6, pady=2)
        self.angle_cb.bind("<<ComboboxSelected>>", lambda e: self._angle_changed())
        Tooltip(self.angle_cb, "One fixed take-off angle for every sample, or the stage tilt stored in each block (a "
                               "tilt series: one Strohmeier thickness per tilt). The K-Alpha keeps the tilt in "
                               "experimental variable 1.")
        ttk.Label(grid, text="Tilt offset (°):", style="Card.TLabel").grid(row=3, column=0, sticky="w", pady=2)
        self.offset_entry = ttk.Entry(grid, textvariable=self.offset, width=7)
        self.offset_entry.grid(row=3, column=1, sticky="w", padx=6, pady=2)
        Tooltip(self.offset_entry, "Added to the tilt to give the emission angle from the surface normal; the take-off "
                                   "angle is 90 minus that. 0 for a tilt measured from normal emission.")
        ttk.Label(grid, text="Substrate:", style="Card.TLabel").grid(row=4, column=0, sticky="w", pady=2)
        cb = ttk.Combobox(grid, textvariable=self.sub_preset, values=list(self.substrates), state="readonly", width=24)
        cb.grid(row=4, column=1, sticky="w", padx=6, pady=2)
        Tooltip(cb, "Density, molar mass and IMFP of the substrate (the Calculator presets).")

        ttk.Label(f, text="Fitted states", style="Card.TLabel").grid(row=2, column=0, sticky="w", pady=(10, 0))
        self.states_frame = ttk.Frame(f, style="Card.TFrame")
        self.states_frame.grid(row=3, column=0, sticky="ew")
        ttk.Label(f, textvariable=self.info, style="CardMuted.TLabel", wraplength=330, justify="left").grid(
            row=4, column=0, sticky="w", pady=(6, 0))

        self.calc_btn = ttk.Button(f, text="Calculate all", style="Primary.TButton", command=self.calculate)
        self.calc_btn.grid(row=5, column=0, sticky="w", pady=(10, 0))
        self.badge = tk.Label(f, text="", bg=CARD, fg=MUTED, anchor="w", justify="left", wraplength=330,
                              padx=6, pady=3)
        self.badge.grid(row=6, column=0, sticky="ew", pady=(8, 0))
        ttk.Label(f, textvariable=self.summary, style="CardMuted.TLabel", wraplength=330, justify="left").grid(
            row=7, column=0, sticky="w", pady=(4, 0))
        ttk.Label(f, style="CardMuted.TLabel", wraplength=330, justify="left",
                  text="Areas are used as stored in the file. They belong to one core level, so the transmission "
                       "function cancels in the ratio.").grid(row=8, column=0, sticky="w", pady=(8, 0))

    def _build_right(self):
        frame = ttk.Frame(self)
        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(0, weight=1)
        page = ttk.Frame(frame, style="Card.TFrame", padding=4)
        page.grid(row=0, column=0, sticky="nsew")
        self.fig = Figure(figsize=(8, 4), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=page)
        toolbar = NavigationToolbar2Tk(self.canvas, page, pack_toolbar=False)
        toolbar.configure(width=250)
        for child in [toolbar] + list(toolbar.winfo_children()):
            try:
                child.configure(bg=CARD)
            except tk.TclError:
                pass
        toolbar.pack(side="bottom", fill="x")
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        cols = ("sample", "theta", "isub", "iov", "ratio", "d", "da", "ok")
        heads = ("Sample", "θ (°)", "I substrate", "I overlayer", "I ov / I sub", "d (nm)", "d (Å)", "Validity / note")
        widths = (220, 55, 90, 90, 80, 70, 70, 170)
        self.table = ttk.Treeview(frame, columns=cols, show="headings", height=6 if self.compact else 9)
        for c, h, w in zip(cols, heads, widths):
            self.table.heading(c, text=h)
            self.table.column(c, width=w, anchor="w" if c in ("sample", "ok") else "e", stretch=(c == "sample"))
        self.table.grid(row=1, column=0, sticky="ew", pady=(6, 0))
        sb = ttk.Scrollbar(frame, orient="vertical", command=self.table.yview)
        sb.grid(row=1, column=1, sticky="ns", pady=(6, 0))
        self.table.configure(yscrollcommand=sb.set)
        bar = ttk.Frame(frame)
        bar.grid(row=2, column=0, sticky="w", pady=(4, 0))
        ttk.Button(bar, text="Export CSV…", command=self.export_csv).pack(side="left")
        ttk.Button(bar, text="Copy table", command=self.copy_table).pack(side="left", padx=6)
        return frame

    # --- file -----------------------------------------------------------------------------------------
    def open_file(self):
        path = filedialog.askopenfilename(title="Open a fitted VAMAS file",
                                          filetypes=[("VAMAS files", "*.vms *.vamas"), ("All files", "*.*")])
        if path:
            self.load_file(path)

    def load_file(self, path):
        try:
            vf = VM.read(path)
        except VM.VamasError as exc:
            messagebox.showerror("Cannot read the file", str(exc))
            return False
        regions = B.fitted_regions(vf)
        if not regions:
            messagebox.showerror("No fits found", "No block of this file carries a CasaXPS peak fit. Save the file "
                                 "from CasaXPS with the fits in place.")
            return False
        self.file, self.rows = vf, []
        self.file_label.set(f"{os.path.basename(path)}: {len(vf.blocks)} blocks")
        self.region_cb.configure(values=regions)
        self.region.set(regions[0])
        self._region_changed()
        return True

    def _region_changed(self):
        self.rows = []
        self.samples, skipped = B.collect(self.file, self.region.get())
        self.states = B.fit_states(self.samples)
        note = f"{len(self.samples)} samples"
        if skipped:
            note += f"; {skipped} block(s) without a fit left out"
        self.info.set(note)
        for w in self.states_frame.winfo_children():
            w.destroy()
        self.role_vars, self.mat_vars = {}, {}
        roles = B.guess_roles(self.states)
        sub_default = ""
        for r, st in enumerate(self.states):
            ttk.Label(self.states_frame, text=f"{st.name} ({st.be:.1f} eV)", style="Card.TLabel").grid(
                row=r, column=0, sticky="w", pady=2)
            rv = tk.StringVar(value=ROLE_TEXT[roles[st.key]])
            cb = ttk.Combobox(self.states_frame, textvariable=rv, values=list(ROLE_TEXT.values()), state="readonly",
                              width=10)
            cb.grid(row=r, column=1, padx=4, pady=2)
            self.role_vars[st.key] = rv
            mv = tk.StringVar(value=self._guess_material(st.name, self.overlayers)
                              or self._guess_material(self.region.get(), self.overlayers))
            mcb = ttk.Combobox(self.states_frame, textvariable=mv, state="readonly", width=16,
                               values=[k for k, v in self.overlayers.items() if v is not None])
            mcb.grid(row=r, column=2, padx=4, pady=2)
            Tooltip(mcb, "Material of this state when it is an overlayer (density, molar mass, IMFP).")
            self.mat_vars[st.key] = mv
            if roles[st.key] == B.SUBSTRATE:
                sub_default = (self._guess_material(st.name, self.substrates)
                               or self._guess_material(self.region.get(), self.substrates))
        self.sub_preset.set(sub_default)
        self._fill_angle_sources()
        self._fill_table()
        self._draw()
        self._set_badge("ok" if self.samples else "warn", "Check the roles and materials, then Calculate all."
                        if self.samples else "No fitted samples in this region.")

    def _fill_angle_sources(self):
        """Offer the tilt sources that the region's blocks actually hold; a varying one is chosen by default."""
        region = self.region.get()
        present = self.file.sources_present(region) if self.file else []
        self._src_by_label = {VM.ANGLE_SOURCES[s]: s for s in present}
        self.angle_cb.configure(values=[FIXED_ANGLE] + list(self._src_by_label))
        best = self.file.default_sources(region)[0] if self.file else None
        varies = best is not None and len({b.value(best) for b in self.file.blocks_of(region)}) > 1
        self.angle_src.set(VM.ANGLE_SOURCES[best] if varies else FIXED_ANGLE)
        self._angle_changed()

    def _angle_changed(self):
        fixed = self.angle_src.get() == FIXED_ANGLE
        self.theta_entry.configure(state="normal" if fixed else "disabled")
        self.offset_entry.configure(state="disabled" if fixed else "normal")

    @staticmethod
    def _guess_material(state_name, presets):
        """First preset whose core level ('Al 2p', 'Si 2p3/2' -> 'Si 2p') appears in the state name ('Al 2p Metal')."""
        n = _norm(state_name)
        for key, comp in presets.items():
            level = _norm(comp.core_level).split("/")[0].rstrip("0123456789") if comp is not None else ""
            level = level or ("" if comp is None else _norm(comp.core_level))
            if level and level in n:
                return key
        return ""

    # --- calculation ----------------------------------------------------------------------------------
    def _set_badge(self, kind, text):
        bg, fg = BADGE[kind]
        self.badge.configure(text=text, bg=bg, fg=fg)

    def calculate(self):
        if not self.samples:
            self._set_badge("warn", "Open a fitted VAMAS file first.")
            return
        if self.angle_src.get() == FIXED_ANGLE:
            theta = _to_float(self.theta.get())
            if theta is None or not (0 < theta <= 90):
                messagebox.showerror("Input error", "The take-off angle must be between 0 and 90 degrees.")
                return
        else:
            offset = _to_float(self.offset.get())
            if offset is None:
                messagebox.showerror("Input error", "The tilt offset must be a number.")
                return
            theta = B.takeoff_angles(self.file, self.samples, self._src_by_label[self.angle_src.get()], offset)
        sub = self.substrates.get(self.sub_preset.get())
        if sub is None:
            messagebox.showerror("Input error", "Choose a substrate material.")
            return
        roles = {k: TEXT_ROLE[v.get()] for k, v in self.role_vars.items()}
        mats = {k: self.overlayers.get(v.get()) for k, v in self.mat_vars.items()}
        try:
            self.rows = B.compute(self.samples, roles, sub, mats, theta, self.volume_density)
        except ValueError as exc:
            messagebox.showerror("Input error", str(exc))
            return
        done = [r for r in self.rows if r.d_nm is not None]
        bad = len(self.rows) - len(done)
        over = sum(1 for r in done if not r.within_validity)
        if done:
            ds = [r.d_nm for r in done]
            self.summary.set(f"d from {min(ds):.2f} to {max(ds):.2f} nm; mean {sum(ds) / len(ds):.2f} nm "
                             f"({len(done)} samples)")
        else:
            self.summary.set("")
        if bad:
            self._set_badge("bad", f"{bad} sample(s) could not be calculated (see the table).")
        elif over:
            self._set_badge("warn", f"{over} sample(s) exceed the ~3λ validity limit: treat with caution.")
        else:
            self._set_badge("ok", "All samples within the validity limit (~3λ).")
        self._fill_table()
        self._draw()

    def _fill_table(self):
        self.table.delete(*self.table.get_children())
        if not self.rows:
            for s in self.samples:
                self.table.insert("", "end", values=(s.sample_id, "", "", "", "", "", "", ""))
            return
        for r in self.rows:
            if r.d_nm is None:
                self.table.insert("", "end", values=(r.sample.sample_id, self._theta_text(r),
                                                     f"{r.i_sub:.6g}" if r.i_sub else "-",
                                                     f"{r.i_ov:.6g}" if r.i_ov else "-", "-", "-", "-", r.error))
                continue
            note = "within ~3λ" if r.within_validity else f"EXCEEDS ~3λ ({r.validity_limit_nm:.1f} nm)"
            self.table.insert("", "end", values=(r.sample.sample_id, self._theta_text(r), f"{r.i_sub:.6g}",
                                                 f"{r.i_ov:.6g}",
                                                 f"{r.i_ov / r.i_sub:.4f}", f"{r.d_nm:.3f}", f"{r.d_nm * 10:.2f}",
                                                 note))

    @staticmethod
    def _theta_text(r):
        return "-" if r.theta_deg is None else f"{r.theta_deg:g}"

    def _table_lines(self):
        head = ["Sample", "Take-off angle (deg)", "I substrate", "I overlayer", "I ov / I sub", "d (nm)", "d (Angstrom)", "Validity / note"]
        out = [head]
        for r in self.rows:
            if r.d_nm is None:
                out.append([r.sample.sample_id, self._theta_text(r), r.i_sub, r.i_ov, "", "", "", r.error])
            else:
                out.append([r.sample.sample_id, self._theta_text(r), r.i_sub, r.i_ov, round(r.i_ov / r.i_sub, 5),
                            round(r.d_nm, 4),
                            round(r.d_nm * 10, 3), "within 3 lambda" if r.within_validity else "exceeds 3 lambda"])
        return out

    def export_csv(self):
        if not self.rows:
            messagebox.showinfo("No data", "Run Calculate all first.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv")],
                                            initialfile="batch_thickness.csv")
        if path:
            with open(path, "w", newline="", encoding="utf-8") as fh:
                csv.writer(fh).writerows(self._table_lines())
            messagebox.showinfo("Exported", f"Saved {len(self.rows)} row(s) to:\n{path}")

    def copy_table(self):
        if self.rows:
            self.clipboard_clear()
            self.clipboard_append("\n".join("\t".join(str(c) for c in line) for line in self._table_lines()))

    # --- plot -----------------------------------------------------------------------------------------
    def _draw(self):
        fig = self.fig
        fig.clear()
        ax = fig.add_subplot(111)
        done = [r for r in self.rows if r.d_nm is not None]
        if not done:
            ax.axis("off")
            ax.text(0.5, 0.5, "Open a CasaXPS-fitted VAMAS file, check the roles and press Calculate all.",
                    ha="center", va="center", color=MUTED, transform=ax.transAxes, wrap=True)
        else:
            xs = list(range(len(done)))
            ax.bar(xs, [r.d_nm for r in done],
                   color=[OK_COLOUR if r.within_validity else WARN_COLOUR for r in done])
            ax.set_xticks(xs)
            tilted = self.angle_src.get() != FIXED_ANGLE
            one_sample = len({r.sample.sample_id for r in done}) == 1
            ax.set_xticklabels([(f"θ = {r.theta_deg:g}°" if one_sample else f"{r.sample.sample_id}  θ={r.theta_deg:g}°")
                                if tilted else r.sample.sample_id for r in done],
                               rotation=60 if not (tilted and one_sample) else 0, ha="right" if not (tilted and one_sample)
                               else "center", fontsize=7 if not (tilted and one_sample) else 9)
            ax.set_ylabel("Overlayer thickness d (nm)")
            ax.grid(axis="y", alpha=0.3)
            fig.subplots_adjust(bottom=0.38, left=0.09, right=0.98, top=0.95)
        self.canvas.draw_idle()
