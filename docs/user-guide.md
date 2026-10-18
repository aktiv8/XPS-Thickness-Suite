# User guide

The XPS Thickness Suite estimates the thickness of thin overlayers (oxides, hydroxides, nitrides, passivation layers,
contamination films) from XPS peak areas. It has one window with a tab for each method.

Contents: [Starting the app](#starting-the-app), [The Calculator tab](#the-calculator-tab),
[The Thickogram tab](#the-thickogram-tab), [Which tab should I use?](#which-tab-should-i-use),
[Troubleshooting](#troubleshooting).

## Starting the app

You need Python 3 with Tkinter, which comes with the standard Windows installer. The calculator itself needs nothing
else. The simulated spectrum, the Thickogram tab and the typeset equations in this documentation also need numpy and
matplotlib:

```
pip install numpy matplotlib
python app.py
```

The window sizes itself to your screen. On a short screen, such as a 1366 by 768 laptop, it switches to a compact
layout with smaller fonts and tighter rows. Hover over most buttons and boxes for a short tip.

## The Calculator tab

![The Calculator tab: inputs on the left, the calculation on the right, the simulated spectrum below](images/calculator.png)

Use this tab when the overlayer and the substrate are two chemical states of the **same core level**, for example
metallic Al and aluminium oxide in the Al 2p spectrum. It applies the Strohmeier equation described in
[Methods](methods.md#strohmeier-equation-calculator-tab).

### Inputs

1. Choose the **substrate** from the first drop-down, for example Al (metal). The row fills in with reference values.
2. For each **overlayer component** present, pick it from the drop-down and tick its box. Tick several for a mixed
   layer such as oxide plus hydroxide. Untick or choose None to leave a component out.
3. Type the fitted **peak areas** (the last column). Enter the substrate area on the substrate row and the area of each
   overlayer component on its own row. Any consistent unit works, because only ratios matter.
4. Set the **take-off angle** $\theta$ measured **from the surface**: 90 means normal emission.

The other columns are pre-filled but editable:

| Column | Meaning |
|---|---|
| Density, Molar mass | Used to work out the volume density of the metal atoms in each phase |
| Atoms/formula | Metal atoms per formula unit, for example 2 for Al2O3 |
| Core level, BE | The peak the row refers to, and its binding energy (also used to position the simulated spectrum) |
| IMFP | Inelastic mean free path of the photoelectrons in that phase, in nm |
| RSF | Relative sensitivity factor; leave at 1 when both areas come from the same core level |
| Peak area | The fitted area of the peak |

Every preset value is a default to check, not an authority: the density of a thin film can differ from the bulk, and
the IMFP depends on your X-ray source. Details of where each value comes from are in
[Presets and data](presets.md#presets-and-data). In particular, **use the Cu (LMM) presets with areas from the Cu LMM
Auger fit, and the Cu (2p3/2) presets with areas from the Cu 2p3/2 fit**; mixing the sets changes the answer by about
30 percent.

### Calculation

Press the blue **Calculate thickness** button (or press Enter in any input box). The right-hand card shows

- the **thickness** in nm, with the same value in angstrom;
- a **badge**: green when the thickness is within the method's validity limit of about three IMFPs, amber when it is
  close, red when it exceeds the limit and the result may be unreliable;
- the details: the weighted density and IMFP of the overlayer, and the density ratio.

**Clear peak areas** empties the area boxes and keeps the presets. **Export history to CSV** saves every
calculation from the session.

### Simulated spectrum

Below the inputs, the app draws what that thickness looks like in the core-level spectrum, next to a to-scale picture
of the layer. This is an **estimate** for visual comparison, not a fit, and the figure says so.

- **Thickness** box and slider: try other thicknesses (the what-if values); **Reset** returns to the calculated one.
- **BE shift**: move the simulated peaks to match your charge correction.
- **FWHM**: the width of the metal and of the oxide peaks; set them to match your instrument.
- **Load my spectrum** overlays a two-column text or CSV file (binding energy, intensity). It is scaled to the
  height of the simulation, so it compares shapes and positions, not absolute counts.
- The toolbar zooms, pans and saves the figure as an image. Drag the divider above the plot to give it more room.

The note under the plot states what is **not** drawn: satellites, plasmon and loss features, Auger lines and
asymmetric metal tails. For Cu 2p3/2 with CuO or Cu(OH)2 the shake-up satellite region is shaded and marked "not
simulated". Auger presets, such as Cu LMM, are not simulated at all.

![A Cu 2p3/2 simulation with the satellite region shaded and the note under the plot](images/calculator-cu.png)

### Worked example: Al2O3 on Al

Choose Al (metal) as the substrate and Al2O3 as the overlayer. With Al 2p areas of 35 for the metal and 65 for the
oxide at a take-off angle of 90 degrees, the calculator gives a thickness of **3.58 nm (35.8 angstrom)**. This is
the same number as Biesinger's Al oxide spreadsheet, which uses Strohmeier's values.

### Worked example: oxide and hydroxide on Cu

Choose Cu (metal, LMM) and tick Cu2O (LMM) and Cu(OH)2 (LMM). A Cu LMM fit of 25 percent metal, 63 percent Cu2O and
13 percent Cu(OH)2 at 90 degrees gives **3.22 nm**, matching the roughly 3.2 nm that Biesinger reports for this
spectrum [5].

## The Thickogram tab

![The Thickogram tab with a cylinder Topofactor applied](images/thickogram.png)

Use this tab when the overlayer and the substrate give **different peaks**, for example a Fe 2p overlayer on a Si 2p
substrate, or carbon on silica. It applies Cumpson's Thickogram and, if the sample is a fibre or a particle, a
Topofactor. The method is in [Methods](methods.md#thickogram-thickogram-tab).

### Inputs

| Field | What to enter |
|---|---|
| Overlayer / substrate binding energy | Binding energy of each peak in eV, for example 707 for Fe 2p and 100 for Si 2p |
| Overlayer / substrate peak area | The fitted peak areas as measured in your software |
| Overlayer / substrate RSF | The sensitivity factor of each peak (see below) |
| Lambda L | Attenuation length of the **overlayer's** photoelectrons in the overlayer, in nm |
| Theta | Take-off angle from the surface **normal**: 0 means normal emission |
| X-ray source energy | 1486.6 eV for Al K-alpha, 1253.6 eV for Mg K-alpha; use the buttons |

Things to watch:

- **The angle convention differs from the Calculator tab.** Here 0 degrees is normal emission and the angle is
  measured from the surface normal; on the Calculator tab 90 degrees is normal emission. The tab shows the
  equivalent surface angle next to the box.
- The **RSFs** are typed in. If you want to follow Cumpson and Shard exactly, enter each sensitivity factor already
  multiplied by the atomic fraction of the element in the pure material.
- The optional **preset** drop-downs fill the binding energy (and lambda for the overlayer) from the Calculator's
  reference table. Presets give an IMFP, which is about 10 percent longer than an attenuation length.
- **Load example** fills demonstration numbers. **Copy from Calculator** takes the angle (converted), binding energies
  and lambda from the Calculator tab, but not peak areas or RSFs.

Press **Calculate and plot**. The big number is the thickness; the nomogram marks the solution with a red line from the
intensity ratio on the left scale through the energy ratio on the right-hand scales, and the thickness is where it
crosses the green curve.

### Sample shape and Topofactors

By default the sample is treated as flat. For a **fibre or a particle** choose its shape under *Sample shape*; the tab
then applies a Topofactor [4], the factor $\mathrm{TF} = T/t$ that converts the apparent planar thickness $t$ into the true
conformal thickness $T$. The big readout shows the corrected value, with the planar value and the factor beside it,
and the nomogram is annotated with all three.

| Shape | Use it for | Typical factor |
|---|---|---|
| Cylinder (fibre) | A fibre lying with its axis perpendicular to the plane of source, sample and analyser | about 0.79 |
| Sphere, X-rays reach all of it | A particle with the X-ray source collinear with the analyser | about 0.67 |
| Sphere, X-ray shadowed | A particle with the source at 45 to 65 degrees to the analyser | about 0.71 |

*Topofactor* chooses between **Eq. (6)**, which depends on the thickness and on the energy ratio, and the **constant
average** values in the table, which are good to about 10 percent.

Before relying on a Topofactor, check that:

- the sample's mean normal faces the analyser, with the take-off angle set to 0 (the tab warns otherwise);
- the fibres or particles are **much larger than the attenuation length**. Topofactors do **not** apply to nanoparticles
  or nanofibres;
- the overlayer is uniform and conformal, and the two peaks are cleanly separated.

### Worked example: a coating on a fibre

Take carbon (C 1s at 285 eV) as the overlayer on a silica fibre (Si 2p at 99.5 eV), with Al K-alpha X-rays, peak areas of
800 and 3000, sensitivity factors of 1.0 and 3.0, an attenuation length of 3.0 nm and normal emission. The planar
Thickogram gives **1.85 nm**. Choosing Cylinder (fibre) with Eq. (6) applies a Topofactor of **0.736** and gives a
corrected thickness of **1.36 nm**.

## Which tab should I use?

| Your situation | Use |
|---|---|
| Oxide on its own metal, one core level, flat sample | Calculator tab |
| Overlayer and substrate are different elements or different peaks, flat sample | Thickogram tab, shape Flat |
| Coating on a micro-scale fibre or particle | Thickogram tab with a Topofactor |
| Nanoparticles or nanofibres | Neither: the Topofactor model is not valid at this scale |

For equal kinetic energies the Thickogram reduces exactly to the Strohmeier equation, so the two tabs agree where
their assumptions overlap.

## Troubleshooting

- **The Thickogram tab or the spectrum pane says to install numpy and matplotlib.** Run `pip install numpy matplotlib`
  and restart the app. The Calculator's thickness calculation still works without them.
- **The equations in this documentation appear as plain text.** The same packages are needed for typeset equations.
- **A badge turns amber or red.** Read it: it names the problem, for example a thickness beyond the validity limit, a
  binding-energy ratio outside the nomogram scales, or a Topofactor applied off normal emission.
- **A preset says it has no binding energy.** Auger presets, such as Cu LMM, carry no binding energy; type the
  binding energy of your own peak.
- **The numbers look about 30 percent off for copper.** Check that you did not mix the LMM and 2p3/2 presets.
