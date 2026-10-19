# User guide

The XPS Thickness Suite estimates the thickness of thin overlayers (oxides, hydroxides, nitrides, passivation layers,
contamination films) from XPS peak areas. It has one window with a tab for each method.

Contents: [Starting the app](#starting-the-app), [The Calculator tab](#the-calculator-tab),
[The Thickogram tab](#the-thickogram-tab), [The ARXPS tab](#the-arxps-tab),
[The Depth profile tab](#the-depth-profile-tab),
[Which tab should I use?](#which-tab-should-i-use), [Troubleshooting](#troubleshooting).

## Starting the app

You need Python 3 with Tkinter, which comes with the standard Windows installer. The calculator itself needs nothing
else. The simulated spectrum, the Thickogram and ARXPS tabs and the typeset equations in this documentation also need
numpy and matplotlib:

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

## The ARXPS tab

![The ARXPS tab with a Si 2p tilt series: spectra at six angles, the integration regions and the table of results](images/arxps.png)

Use this tab when you have measured the same spot at **several emission angles** (a tilt series) and want the
overlayer thickness from the slope of a line instead of from a single ratio. It reads VAMAS files, finds and aligns the
peaks, subtracts a Shirley background, integrates the areas and fits the thickness. The method is in
[Methods](methods.md#angle-resolved-xps-arxps-tab).

### Workflow

The setup has three pages. Work through them from left to right and press **Calculate and plot**.

**1 Data.** Press *Open VAMAS file* and choose a `.vms` file. Then choose:

| Field | What to enter |
|---|---|
| Mode | *Two states of one peak* when the overlayer and substrate are two chemical states in one region, for example Si and SiO2 in Si 2p. *Two different peaks* when they are separate regions, for example O 1s over Si 2p |
| Region | The block (region) to analyse; in the two-peak mode one region for each peak |
| Angle from | Where the tilt comes from. The K-Alpha stores the stage tilt in experimental variable 1, the default |
| Offset | Added to the tilt to give the emission angle from the surface normal (0 for the K-Alpha) |
| Azimuth from | Optional: the field that holds the sample azimuth. If the file has several azimuths each is fitted separately |
| Max angle | Angles above this are drawn but left out of the fit (default 70 degrees) |
| Counts per second, Transmission-correct | Normalise the intensities; they cancel for two states of one peak but matter for two different peaks |
| Areas from | *Windows + Shirley* integrates the spectra. *Fitted components* uses a peak fit that CasaXPS stored in the file; it is chosen automatically when the region has one |

**2 Materials.** Optionally choose the substrate and overlayer presets, the same ones as on the Calculator tab: they fill
the binding energies, lambda for the overlayer and, for two states of one peak, the ratio $R_\infty$ of the two layers
when infinitely thick. Everything can be edited. For two different peaks enter the sensitivity factors, the binding
energies and lambda yourself; the photon energy is read from the file.

**3 Windows.** Press *Auto-detect windows*. For two states of one peak this sets the envelope that holds both
peaks and the energy at which it is split; for two peaks it sets a window round each. The detection uses how the spectrum
changes with angle, so it also works when the oxide peak is the taller one. **Always check the result on the Spectra
plot**: the blue and orange areas should cover each state at every angle.

**3 Components.** When the areas come from a CasaXPS fit this page replaces the windows. It lists the fitted states
(components that share an INDEX, such as the two halves of a doublet, count as one) with each state's share of the fitted
area at the lowest and highest angle, and a selector: *Substrate*, *Overlayer* or *Ignore*. The tab fills it in by
the trend, a falling share meaning substrate and a rising share overlayer; change any row, or press *Auto-assign* to
start again.

### Reading the results

- The big number is the thickness from the slope of the line through the origin, with the standard error of the fit.
- **Spectra** shows every angle, aligned and normalised, with the Shirley background and the two integrated areas. Dotted
  lines are angles left out of the fit.
- **x vs 1/cos theta** shows the per-angle reduced thickness $x$ against $1/\cos\theta$. For a uniform overlayer the points
  lie on a straight line through the origin; the dashed line is a free-intercept fit and its offset is the check.
- **Azimuth** appears as a polar plot when the file has more than one azimuth. The mean thickness is shown, with the
  spread between azimuths.
- The table lists, for each angle, the two areas, the ratio, $x$, the thickness from that angle alone and the energy
  shift that was removed.
- An amber badge names any problem: an intercept away from the origin, angles left out, fewer than three angles, a
  peak split that cannot separate overlapping states, or a thickness beyond three attenuation lengths.

![The same data as x against 1/cos theta with the fitted line](images/arxps-fit.png)

### Worked example: a synthetic Si 2p tilt series

The file `docs/sample/si-native-oxide-arxps.vms` is a **synthetic** series (not a measurement) for a 1.5 nm layer of SiO2 on
silicon, with the peaks drifting by 0.3 eV with tilt, Poisson noise and tilts of 0 to 75 degrees. Open it, choose the
region Si2p Scan, the presets *Si (metal/elemental)* and *SiO2* (lambda 3.89 nm, $R_\infty$ = 0.555) and press
*Auto-detect windows*: the envelope is 96.5 to 106.1 eV, split at 100.8 eV. The result is **1.48 nm**, within 0.02 nm of
the value that was put in. The 75 degree point is left out by the 70 degree limit, and the free-intercept fit passes
through the origin to within its error.

### Using a peak fit from CasaXPS

![The ARXPS tab on a file with a stored CasaXPS fit: the Components page and the measured spectra with the fitted positions](images/arxps-fitted.png)

If you have fitted the region in CasaXPS and saved the file as VAMAS, open that file. The tab finds the fit, sets
*Areas from* to *Fitted components* and fills the list on the Components page. The areas are taken exactly as CasaXPS
stored them, so nothing is re-integrated and no window or background has to be chosen. This is the way to analyse
states that overlap too much for a window split, such as the oxide and substrate components of Ga 3d on GaAs.

- For **two states of one peak** pick the region; the substrate and overlayer states are sums of the states you assign.
- For **two different peaks** choose the region holding the substrate fit and the region holding the overlayer fit
  (each must be fitted); the sensitivity factors and binding energies are filled in from the fit and can be edited.
  Keep *Transmission-correct* ticked, because the stored areas are raw counts.
- If a block holds several fitted regions, such as Ga 3d and As 3d in one scan, a selector on the Components page
  chooses which to use.
- On the Spectra page the measured spectra are drawn with a dashed line at each fitted component position; the peak
  shapes are not redrawn.

### Worked example: a synthetic file with a stored fit

The file `docs/sample/si-native-oxide-arxps-fitted.vms` holds the same synthetic series as above with a CasaXPS-style fit
stored in each block (two doublets, Si El and Si Ox). Open it, choose the region Si2p Scan and the presets *Si (metal/elemental)*
and *SiO2*: the tab selects the fitted components, makes Si El the substrate and Si Ox the overlayer, and gives
**1.50 nm** after *Calculate and plot*, with no windows to set.

### Tips

- Use a series with at least three angles spanning a good range of $1/\cos\theta$. Angles between about 0 and 60 degrees
  are the most reliable.
- Two states less than about 2 eV apart (for example the oxide and the substrate component of Ga 3d in GaAs) overlap
  too much for a window split. The tab warns, and the thickness is then only a rough estimate; fit the peaks in
  CasaXPS and use the fitted components instead.
- The thickness uses the attenuation length you give, so its accuracy is limited by that number and by $R_\infty$. The
  *R-infinity-free* line in the details shows what the data alone prefer.
- Real data often show an intercept above zero: contamination, roughness and elastic scattering all bend the line. The
  badge says when it is more than the scatter allows.

## The Depth profile tab

![The Depth profile tab: the Surrey Trend Plot of a fitted synthetic Si 2p series](images/depth-trend.png)

Use this tab when you want to know where several species sit in depth, not just the thickness of one overlayer. It reads the
same VAMAS tilt series as the ARXPS tab. The method is in [Methods](methods.md#depth-profiles-depth-profile-tab).

1. Press *Open VAMAS file* or *Use ARXPS tab's file* (which also takes the file you already loaded there).
2. Check *Angle from*, *Offset* and *Max angle* as on the ARXPS tab.
3. Choose *Species from*: *Fitted states* uses every chemical state of the CasaXPS fits in the file (chosen automatically if
   there are any); *Whole regions* integrates each region with a Shirley background and works with any file.
4. In the Species list untick anything you do not want and set each attenuation length *lambda* for your material. The
   starting values are only a guide.

### Surrey Trend Plot

Each species is drawn as ln(intensity) against sec theta (or against theta), with error bars from counting statistics when
they are known and a straight-line fit. Species whose lines run parallel, and whose kinetic energies are similar, may be
distributed together. The table gives each line's slope, its error, the scatter of the points, an apparent depth
(minus lambda times the slope), and which other species run parallel. Use it as a first look before choosing a model: it
tells you how many distinct depth distributions the data hold and how noisy they are.

## Which tab should I use?

| Your situation | Use |
|---|---|
| Oxide on its own metal, one core level, flat sample | Calculator tab |
| Overlayer and substrate are different elements or different peaks, flat sample | Thickogram tab, shape Flat |
| Coating on a micro-scale fibre or particle | Thickogram tab with a Topofactor |
| A tilt series (several emission angles) in a VAMAS file | ARXPS tab |
| Where several species sit in depth, from a tilt series | Depth profile tab |
| Nanoparticles or nanofibres | Neither: the Topofactor model is not valid at this scale |

For equal kinetic energies the Thickogram reduces exactly to the Strohmeier equation, so the two tabs agree where
their assumptions overlap.

## Troubleshooting

- **The Thickogram tab or the spectrum pane says to install numpy and matplotlib.** Run `pip install numpy matplotlib`
  and restart the app. The Calculator's thickness calculation still works without them.
- **The equations in this documentation appear as plain text.** The same packages are needed for typeset equations.
- **A badge turns amber or red.** Read it: it names the problem, for example a thickness beyond the validity limit, a
  binding-energy ratio outside the nomogram scales, or a Topofactor applied off normal emission.
- **Auto-detect windows fails or the shaded areas look wrong on the ARXPS tab.** Check that the presets (or binding
  energies) are close to your peaks and that the lowest angle is first; set the envelope and split by hand on the
  Windows page if the two states are not resolved.
- **The ARXPS thickness changes a lot when I move the split.** The two states overlap; see the tips under the ARXPS tab.
- **The Components list is empty or says there is no fit.** The file has to be saved from CasaXPS with the fit in
  place; a file that CasaXPS only opened, or the original instrument file, carries no fit. Switch *Areas from* to
  *Windows + Shirley* to integrate the spectra instead.
- **A state is assigned to the wrong side.** The guess uses how the share of each state changes with angle, which
  fails for a state that hardly changes; set it by hand in the list.
- **A preset says it has no binding energy.** Auger presets, such as Cu LMM, carry no binding energy; type the
  binding energy of your own peak.
- **The numbers look about 30 percent off for copper.** Check that you did not mix the LMM and 2p3/2 presets.
