# XPS Thickness Suite

Tkinter tools for estimating the thickness of thin overlayers (oxides, hydroxides,
nitrides, passivation layers) from XPS peak areas, in one window with one tab per method:

| Tab | Method | Use it when |
|---|---|---|
| **Calculator & simulation** | Strohmeier equation (B. R. Strohmeier, *Surf. Interface Anal.* 1990, **15**, 51; [HarwellXPS Tech Note #1](https://subsite.harwellxps.uk/wp-content/uploads/2018/02/HarwellXPS_TechNote_01.pdf)) | The overlayer and the substrate are two chemical states of the **same core level** (metal vs. oxide components of Al 2p, Cu LMM, Si 2p ...). Uses densities and IMFPs; includes a simulated spectrum. |
| **Thickogram (different core levels)** | Cumpson's thickogram (P. J. Cumpson, *Surf. Interface Anal.* 2000, **29**, 403) | The overlayer and the substrate are **different peaks** (e.g. a Fe 2p overlayer on Si 2p). Uses RSFs and kinetic-energy scaling; draws the classic nomogram. |
| **Method / References** | | Equations, assumptions and caveats for both. |

For equal kinetic energies the thickogram reduces exactly to the Strohmeier equation, and the tests check that.

## Running it

The thickness calculator needs no third-party packages (only the Python
standard library / Tkinter, which ships with standard Windows Python installs).
The **Simulated spectrum** pane and the **Thickogram** tab need numpy and matplotlib
(without them those parts show an install hint and the Calculator still works):

```bash
pip install numpy matplotlib   # optional
python app.py
```

## Calculator tab: what it does

- Pick a **substrate** (Cu, Al, Si, or custom) and up to 4 **overlayer**
  components (CuO, Cu2O, Cu(OH)2, Al2O3, SiO2, Si3N4, SiON, or custom) —
  tick the ones present and enter their fitted peak areas.
- Density, molar mass, atoms/formula unit, typical binding energy and IMFP
  are pre-filled from a small reference table and are all editable.
- Enter the photoelectron take-off angle and the substrate ("metal") peak
  area, then click **Calculate thickness**.
- Handles both single-component (metal vs. one oxide) and mixed-layer
  (e.g. oxide + hydroxide present together) cases automatically, using
  intensity-weighted average density/IMFP for the mixed case, matching the
  two forms of the equation given in the tech note.
- An optional RSF (relative sensitivity factor) field on each row
  generalises the equation for cases mixing different core levels/elements;
  leave RSF = 1.0 for the standard same-element case, where it cancels.
- Flags whether the result is within the method's validity limit
  (~3x the overlayer IMFP).
- Results can be exported to CSV via **Export history to CSV...**.
- **One window:** inputs on the left, the calculation (big blue
  **Calculate thickness** button, a large thickness readout and a green/amber/red
  validity badge) on the right, and the simulated spectrum underneath. Drag the
  divider between them to trade height. On short screens (laptops) it switches
  to a compact layout; pressing Enter in any input box also calculates.
- **Simulated spectrum (estimate)** pane: after a calculation it draws what that
  thickness would look like for the substrate and oxide core-level peaks
  (spin-orbit doublets, Shirley-type background) next to a to-scale layer
  diagram, clearly labelled as an estimate. Use **Load my spectrum...**
  (2-column text/CSV of binding energy, intensity) to overlay your own data,
  the thickness box/slider for what-if values, the BE shift for charge
  correction and the FWHM boxes to match your resolution. The figure can be
  saved as PNG from its toolbar. Peak splitting/FWHM are typical defaults.
- **What the simulation leaves out is stated on the figure:** for simplicity it
  does not draw shake-up satellites, plasmon/loss features, Auger lines or
  asymmetric metal tails. For Cu 2p3/2 with CuO or Cu(OH)₂ the Cu(II) shake-up
  satellite region is shaded and labelled "not simulated" (the satellites carry
  roughly 35–39 % of the Cu(II) area, so the drawn Cu(II) main peak is
  overstated). Auger presets (e.g. Cu LMM) are not simulated.

See the **Method / References** tab inside the app for the full equations,
derivation notes, and important caveats about the reference density/IMFP
defaults (in particular: IMFP is source/instrument-dependent and film
density can differ from bulk values — verify both for rigorous work).

Several presets come straight from the original papers, and the tests pin them
to the published numbers:

| System | Source | Inputs |
|---|---|---|
| Al / Al₂O₃ (Al Kα and Mg Kα) | Strohmeier, *Surf. Interface Anal.* 1990, **15**, 51 | IMFP 2.8 / 2.6 nm (Al Kα), 2.4 / 2.2 nm (Mg Kα), Nm/Nox = 1.5 |
| Cu / Cu₂O / CuO / Cu(OH)₂ | Biesinger, *Surf. Interface Anal.* 2017, **49**, 1325 | 8.96 / 6.00 / 6.31 / 3.37 g/cm³; IMFP 1.54 / 1.77 / 1.76 / 2.19 nm |
| W / WO₃ | Carlson & McGuire, *J. Electron Spectrosc.* 1972/73, **1**, 161 | 19.3 / 6.8 g/cm³; IMFP 1.28 / 2.63 nm; 45° take-off |

The Al preset also reproduces M. C. Biesinger's Al oxide/hydroxide spreadsheet to
0.1 Å. The Cu IMFPs are for Cu LMM Auger electrons (KE ≈ 918 eV), as in the
paper, so use the `(LMM)` presets with Cu LMM areas. Separate `(2p3/2)` presets
(IMFP 1.01 / 1.23 / 1.22 / 1.51 nm for Cu / Cu₂O / CuO / Cu(OH)₂ at KE ≈ 554 eV) are
**computed** with TPP-2M (`tpp2m_imfp_nm()` in `app.py`), not published values; the same
code reproduces Biesinger's oxide IMFPs at 918 eV to 0.2 %. Don't mix the two sets —
thickness scales almost linearly with the IMFPs.

Si, SiO₂, Si₃N₄, Ti and TiO₂ take density and band gap from the Thermo Avantage
`Densities.txt` (2.33/2.196/3.44/4.54/4.26 g/cm³), with IMFPs **computed** by TPP-2M
(3.09/3.89/3.00/2.19/2.17 nm at Al Kα). Note the 3.44 g/cm³ Si₃N₄ density is high for
thin amorphous films. SiON is still a rough placeholder.

More pairs: Mg/MgO, Ge/GeO₂, Zn/ZnO (Zn LMM line), Zr/ZrO₂, Hf/HfO₂, Sn/SnO₂, SiC (→ SiO₂), GaAs and GaP
(Ga 3d → Ga₂O₃), CdTe (Te 3d5/2 → TeO₂), InP (In 3d5/2 → In₂O₃; small shift) and PbS (S 2p → PbSO₄). Densities/band gaps come from the Avantage file where it has
them, otherwise standard handbook values (not from the user's files); IMFPs are computed with TPP-2M.
Their Nm/Nox ratios agree with the Pilling–Bedworth ratios (Mg 0.79, Zr 1.55, Hf 1.62, Sn 1.34, Zn 1.58).
Not included: ZnS, diamond and graphite (see CLAUDE.md for why).

`BAND_GAPS_EV` in `app.py` is a reference table of band gaps (used only for TPP-2M IMFP
estimates): the Avantage file's values, Biesinger's Cu oxide gaps, and four user-supplied
additions (InP 1.27, GaP 2.25, CdTe 1.44, PbS 0.37 eV). Existing entries are never overwritten.

## Thickogram tab

Cumpson's method for an overlayer and a substrate measured from **different core levels**.
Enter the two binding energies, peak areas and RSFs, the overlayer's attenuation length λ,
the take-off angle and the X-ray source energy; the tab solves for the thickness and marks it
on the thickogram (red line from the intensity-ratio scale through the energy-ratio scale; its
crossing with the green curve is the answer). The nomogram is the original one, embedded
(zoom, pan, save as PNG from its toolbar).

- With x = t/(λ cos θ): `A = Io/(Is·Ro) = 2·exp((B^0.75 − 0.5)x)·sinh(x/2)`, `Ro = So/Ss`,
  `B = KEo/KEs`, `t = x·λ·cos θ`.
- **Angle convention:** this tab follows Cumpson (θ from the surface *normal*, 0 = normal
  emission); the Calculator tab uses θ from the *surface* (90 = normal). The tab shows the
  conversion live and **Copy from Calculator** applies it (and copies the binding energies and λ;
  peak areas and RSFs are not copied).
- **Load example** fills demo numbers; the optional preset dropdowns fill binding energy and λ
  from the Calculator's reference table (an IMFP is ~10 % longer than an attenuation length).
- The printed scales cover B = 0.4–3.0 and x = 0.1–6; outside that the answer is still
  calculated, with a warning that the red line leaves the printed scales.
- Fixed relative to the original script: the red construction line now ends on the B scale
  (the original ended it at the thickness position, its own "redline plots correctly" to-do),
  inputs are validated, the solver cannot hang, and the plot is embedded instead of opening a
  separate pyplot window.

## Tests

```bash
python -m unittest discover -s tests -v
```

`tests/test_regression.py` checks published numbers from the papers above and Biesinger's Al
spreadsheet (65/35 → 35.8 Å, 23/77 → 9.7 Å, 85/15 → 61.2 Å); `test_thickogram.py` checks the
thickogram solver against the original script's results, the B = 1 equivalence with the
Strohmeier equation and the red-line geometry; `test_spectrum_sim.py` covers the spectrum
simulation. The numpy tests skip themselves if numpy is missing, and `test_gui_smoke.py` (drives
the real window) skips itself when no display is available.

## Files

- `app.py` — the application window, the Calculator tab, presets and the Strohmeier calculation.
- `gui_common.py` — shared GUI palette, ttk theme/styles and tooltips.
- `spectrum_sim.py` — simulated-spectrum maths (numpy only, no GUI).
- `thickogram.py` — thickogram maths and nomogram drawing (numpy only, no Tk).
- `thickogram_tab.py` — the Thickogram tab.
- `legacy/thickogram_v1_3.py` — the original standalone Thickogram script, unmodified, for provenance.
- `tests/` — regression, spectrum, thickogram and GUI smoke tests.
- `NOTICE` — attribution and licensing notes.

## Credits

Thickogram: Peter J. Cumpson, *Surf. Interface Anal.* **29**, 403–406 (2000); the program is a
Python/Tk conversion of the GlowScript program by Prof. Jong Wan Lee (Hallym University), via
the standalone `aktiv8/Thickogram` repository. See `NOTICE`.
