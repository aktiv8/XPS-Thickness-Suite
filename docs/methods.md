# Methods

This page gives the equations behind each tab, the assumptions they rest on, and where they come from. Numbers in
square brackets, such as [3], refer to the [reference list](references.md#references).

Contents: [Overview](#overview), [Strohmeier equation](#strohmeier-equation-calculator-tab),
[Thickogram](#thickogram-thickogram-tab), [Topofactors](#topofactors-fibres-and-particles),
[Angle-resolved XPS](#angle-resolved-xps-arxps-tab), [Depth profiles](#depth-profiles-depth-profile-tab), [Dual-anode XPS](#dual-anode-xps-daxpes-tab),
[Spectrum simulation](#spectrum-simulation),
[IMFPs from TPP-2M](#imfps-from-tpp-2m).

## Overview

All the methods assume a **uniform overlayer of constant thickness on a semi-infinite substrate**, with the
overlayer and substrate signals separable by peak fitting. They differ in which signals are compared.

| Method | Signals compared | Needs | Valid for |
|---|---|---|---|
| Strohmeier (Calculator tab) | Two chemical states of the same core level, for example metallic Al and Al oxide in Al 2p | Densities, IMFPs | Flat samples, thickness up to about 3 IMFPs |
| Thickogram (Thickogram tab) | Two different peaks, for example Fe 2p overlayer and Si 2p substrate | RSFs, kinetic energies, attenuation length | Flat samples |
| Topofactor (Thickogram tab) | Correction applied to the Thickogram result | Sample shape | Fibres and particles much larger than the attenuation length |
| Depth profile (Depth profile tab) | Several species followed against emission angle | Attenuation lengths, a tilt series | Flat samples, angles up to about 70 degrees |
| ARXPS (ARXPS tab) | The same signals as above, measured at several emission angles from a VAMAS tilt series | Presets or RSFs, lambda, a tilt series | Flat samples, emission angles up to about 70 degrees |
| DAXPES (DAXPES tab) | One core level measured with a soft and a hard X-ray source, film and thick reference | Two attenuation lengths, a reference for each source | Films that may cover only part of the area; the same analysed area for both sources |

## Strohmeier equation (Calculator tab)

For an oxide (or hydroxide, or nitride) on its own metal, the metal and oxide components of one core level are
measured in a single spectrum. With $I_m$ and $I_{ox}$ their fitted peak areas, the thickness $d$ follows from [1]:

```math
d = \lambda_{ox}\,\sin\theta\;\ln\left[\frac{N_m\,\lambda_m\,I_{ox}}{N_{ox}\,\lambda_{ox}\,I_m}+1\right]
```

Here $\lambda_m$ and $\lambda_{ox}$ are the inelastic mean free paths (IMFPs) of the photoelectrons in the metal and
in the oxide, $N_m$ and $N_{ox}$ are the volume densities of the *metal atoms* in each phase, and $\theta$ is the
take-off angle measured **from the sample surface** (90 degrees is normal emission). The equation is the solved
form of the intensity ratio for a uniform overlayer:

```math
\frac{I_m}{I_{ox}} = \frac{N_m\,\lambda_m}{N_{ox}\,\lambda_{ox}}\;\frac{\exp\left[-d/(\lambda_{ox}\sin\theta)\right]}{1-\exp\left[-d/(\lambda_{ox}\sin\theta)\right]}
```

### Volume density

The volume density of metal atoms in a compound is

```math
N = \frac{n\,\rho}{M}
```

with $\rho$ the density, $M$ the molar mass and $n$ the number of metal atoms per formula unit (2 for Cu2O and
Al2O3, 3 for Si3N4). Only the ratio $N_m/N_{ox}$ enters the equation, so the units cancel.

### Mixed overlayers

When several overlayer components are present, for example oxide and hydroxide, their areas are summed and the
overlayer density and IMFP are replaced by intensity-weighted averages [2]:

```math
N_{ave}=\frac{\sum_i N_i I_i}{\sum_i I_i},\qquad
\lambda_{ave}=\frac{\sum_i \lambda_i I_i}{\sum_i I_i}
```

```math
d = \lambda_{ave}\,\sin\theta\;\ln\left[\frac{N_m\,\lambda_m\,\sum_i I_i}{N_{ave}\,\lambda_{ave}\,I_m}+1\right]
```

### Sensitivity factors

The Calculator also accepts a relative sensitivity factor (RSF) for the substrate and each overlayer row. For the
usual same-core-level case the RSFs are identical and cancel, so they default to 1. If you combine intensities
from different core levels, the general form is

```math
d = \lambda_{ov}\,\sin\theta\;\ln\left[\frac{R_{sub}\,N_{sub}\,\lambda_{sub}\,I_{ov}}{R_{ov}\,N_{ov}\,\lambda_{ov}\,I_{sub}}+1\right]
```

For different peaks the [Thickogram](#thickogram-thickogram-tab) is the better-founded method.

### Validity

- The method assumes a uniform, continuous overlayer of constant composition.
- It is reliable up to roughly three IMFPs of the overlayer, $d \lesssim 3\lambda$, which is also about the
  sampling depth of XPS. Beyond that the metal signal vanishes into the noise. The app flags results above this
  limit.
- Elastic scattering is neglected, and the IMFP is used in place of the (slightly shorter) attenuation length.
- Strohmeier's own estimate of the accuracy is 10 to 20 percent in absolute terms, with a precision of about
  2 angstrom for repeated measurements [1].

## Thickogram (Thickogram tab)

Cumpson's Thickogram [3] gives the thickness of an overlayer from one peak of the overlayer and one peak of the
substrate, which need not be the same core level or even the same element. It also assumes a flat sample.

Let $x = t/(L\,\cos\theta)$ be the thickness in units of the overlayer attenuation length $L$ along the emission
direction; here $\theta$ is measured **from the surface normal**, so $\theta = 0$ is normal emission. The ratio of
the normalised signals is

```math
\frac{I_o/I_o^{\infty}}{I_s/I_s^{\infty}} = \frac{1-e^{-x}}{e^{-B x}}
```

where $I_o^{\infty}$ and $I_s^{\infty}$ are the intensities from infinitely thick samples of pure overlayer and pure
substrate. Substrate electrons are attenuated in the overlayer with a length $L/B$, and Cumpson takes the energy
dependence of the attenuation length as

```math
B = \left(\frac{E_o}{E_s}\right)^{0.75}
```

where $E_o$ and $E_s$ are the kinetic energies of the overlayer and substrate photoelectrons, each equal to the
X-ray energy minus the binding energy. Defining the measured ratio

```math
A = \frac{I_o}{I_s\,R_o},\qquad R_o = \frac{I_o^{\infty}}{I_s^{\infty}} = \frac{X_o S_o}{X_s S_s}
```

where $S$ is the sensitivity factor of the peak and $X$ the atomic fraction of the element in the pure material
(Eq. 10 of [4]), the equation becomes

```math
A = (1-e^{-x})\,e^{B x} = 2\,\exp\left[\left(B-\frac{1}{2}\right)x\right]\sinh\frac{x}{2}
```

and the thickness is

```math
t = x\,L\,\cos\theta
```

The equation has exactly one root $x$ for any positive $A$ and $B$, found numerically. The thickogram nomogram on
the tab is the graphical solution: draw a line from $A$ on the left scale to $B$ on the right-hand scales, and
read the answer where it crosses the green curve.

### Practical notes

- The tab asks for $S_o$ and $S_s$ as typed values, so you can enter the sensitivity factors already multiplied by
  the atomic fractions $X_o$ and $X_s$ to follow Eq. 10 of [4] exactly.
- For equal kinetic energies $B = 1$ and the equation reduces to $x = \ln(1+A)$. This is exactly the Strohmeier
  equation with $N_m = N_{ox}$ and $\lambda_m = \lambda_{ox}$, and the test suite checks it.
- The nomogram scales cover $0.4 \leq E_o/E_s \leq 3$ and $0.1 \leq x \leq 6$. Outside that range the thickness is
  still calculated but the red line leaves the printed scales.
- **Angle conventions differ between the tabs.** The Thickogram uses $\theta$ from the surface normal and the
  Calculator uses $\theta$ from the surface, so that $\theta_{surface} = 90^{\circ} - \theta_{normal}$. The
  Thickogram tab shows the conversion and its Copy from Calculator button applies it.
- The red construction line in the original Thickogram script ended at the thickness position instead of on the
  $B$ scale. The suite corrects this, so the line now meets the green curve exactly at the answer.

![The Thickogram tab with a cylinder Topofactor applied](images/thickogram.png)

## Topofactors: fibres and particles

The Thickogram treats the sample as flat. For a fibre or a particle analysed with its mean surface normal pointing
at the analyser, the planar equation overestimates the overlayer thickness, because much of the surface is steeply
inclined and its overlayer is viewed along a longer path. Shard, Wang and Spencer [4] introduced the **Topofactor**,
the ratio of the true conformal thickness $T$ to the apparent planar thickness $t$ given by the Thickogram:

```math
\mathrm{TF} = \frac{T}{t}
```

so the corrected thickness is $T = \mathrm{TF}\cdot t$. Analysis is done at $\theta = 0$, the planar thickness is
calculated from the Thickogram, and the Topofactor is applied to the result.

### Where the Topofactor comes from

For a surface of known shape, the local thickness along the emission direction is $T/\cos\alpha$, where $\alpha$ is
the local take-off angle, and the signals are summed over the analysed area. The overlayer to substrate ratio is
then (Eq. 2 of [4])

```math
\frac{I_o/I_o^{\infty}}{I_s/I_s^{\infty}} = \frac{\int\left[1-\exp\left(-T/\cos\alpha\right)\right]dA}{\int\exp\left(-B\,T/\cos\alpha\right)dA}
```

Equating this to the planar Thickogram ratio gives $t$ for a given $T$, and hence $\mathrm{TF} = T/t$. For an
infinitesimally thin overlayer this reduces to $\mathrm{TF} = 1/\langle\sec\alpha\rangle$, which equals the ratio of the
projected to the true analysed area:

```math
\mathrm{TF}_{t\to 0} = \frac{1}{\langle\sec\alpha\rangle} = \frac{A_{proj}}{A_{true}}
```

This gives $2/\pi \approx 0.64$ for a cylinder (fibre axis perpendicular to the analyser), $0.50$ for a sphere whose
whole surface is illuminated, and $0.57 \pm 0.01$ for a sphere with X-ray shadowing at a source-analyser angle of
45 to 65 degrees. For thicker overlayers the numerical Topofactor rises with thickness.

### Equation 6 of the paper

For $0.1 < t < 3$, with $t$ the planar thickness in units of the attenuation length $L$, the paper's numerical results
are reproduced to about 1 percent by

```math
\mathrm{TF} \approx c_0 + \left(c_1 + c_2\ln B\right)t^{k} + c_3\,e^{-c_4 B}
```

with $B = (E_o/E_s)^{0.75}$ as above, and the shape parameters of Table 1 of [4]:

| Shape | c0 | c1 | c2 | c3 | c4 | k | Thin limit | Constant |
|---|---|---|---|---|---|---|---|---|
| Cylinder (fibre) | 0.46 | 0.25 | 0.09 | 0.16 | 1.1 | 0.22 | 0.64 | 0.79 |
| Sphere, X-rays reach all of it | 0.40 | 0.20 | 0.09 | 0.14 | 1.5 | 0.33 | 0.50 | 0.67 |
| Sphere, X-ray shadowed | 0.46 | 0.19 | 0.08 | 0.13 | 1.4 | 0.30 | 0.57 | 0.71 |

The last column is the single average value the paper recommends when a quick answer within about 10 percent will
do (0.79 for cylinders, 0.67 for spheres, calculated at $t = 1.5$ and $B = 1$). The tab offers both: **Eq. (6)**, which
follows the thickness and $B$, and **Constant average**.

![Topofactors against planar thickness, from Eq. 6](images/topofactor-curves.png)

### How the tab applies it

1. The tab solves the planar Thickogram for $t$ in nanometres, then forms $t/L$.
2. It converts the energy ratio to the paper's $B = (E_o/E_s)^{0.75}$. The details box lists the plain energy ratio
   as B and the exponentiated value as B (paper).
3. Below $t/L = 0.1$ it uses the thin-film limit, between 0.1 and 3 it uses Eq. 6, and above 3 it holds the factor at
   its value for 3 and warns. The factor therefore changes slightly at $t/L = 0.1$, where the limit and Eq. 6 meet.
4. The corrected thickness $T = \mathrm{TF}\cdot t$ becomes the headline result; the planar value is kept beside it.

### Conditions and limits

- **Mean surface normal towards the analyser** ($\theta = 0$). The paper recommends this geometry and advises against
  the so-called magic angle. The tab warns if $\theta$ is not 0.
- **Orientation.** Fibres should lie with their axis perpendicular to the plane containing the source, sample and
  analyser, to minimise X-ray shadowing. The shadowed-sphere row applies when half the sphere is in shadow.
- **The overlayer is uniform, conformal and defect free**, and the overlayer and substrate each have a unique peak.
- **Size.** Topofactors are for objects much larger than the electron attenuation length, from micrometre-scale
  particles and fibres upwards. They are **not valid for nanoparticles and nanofibres**, whose radii are comparable
  to $L$; for those the paper gives a separate model (its Eqs. 7 to 9) that the suite does not implement.
- **Accuracy.** Eq. 6 is best when $B$ is close to 1; the paper notes that $B$ usually lies between 0.7 and 1.4, and
  the tab warns outside that range. Experiments on polystyrene-coated silica capillaries gave a mean Topofactor of
  about 0.79, within 10 percent of the cylinder value.
- **Inconsistency in the paper.** Table 1 is captioned as parameters for Eqn (7) but belongs to Eqn (6). The sphere
  value 0.67 is quoted at $t = 1.5$ in the text and at $T = 1.5$ in Table 2; Eq. 6 gives 0.66 at $t = 1.5$.

## Angle-resolved XPS (ARXPS tab)

Measuring the same spot at several emission angles $\theta$ (from the surface normal) changes the path length of the
photoelectrons through the overlayer to $d/\cos\theta$. The overlayer to substrate ratio therefore grows with angle in
a way that depends only on the thickness, so a tilt series gives the thickness from the **slope** of a line rather than
from a single ratio [12]. Departures from a straight line show up as a warning instead of a silent error.

### From ratios to a line

At each angle the tab turns the measured ratio into the reduced thickness $x_i = d/(\lambda\cos\theta_i)$, using the
same equations as the other tabs. For two states of one core level (the Strohmeier case, Calculator tab) with
$R_i = I_{ov}/I_{sub}$ at angle $i$:

```math
x_i = \ln\left(1+\frac{R_i}{R_\infty}\right),\qquad R_\infty = \frac{S_{ov}\,N_{ov}\,\lambda_{ov}}{S_{sub}\,N_{sub}\,\lambda_{sub}}
```

where $R_\infty$ is the ratio for infinitely thick layers, built from the presets exactly as on the Calculator tab. For two
different peaks (the Thickogram case) $x_i$ is the root of $A_i = (1-e^{-x_i})\,e^{B x_i}$ with $A_i$ and $B$ as in the
[Thickogram](#thickogram-thickogram-tab) section; the first form is the special case $B = 1$. In both cases

```math
x_i = k\,u_i,\qquad u_i = \frac{1}{\cos\theta_i},\qquad k = \frac{d}{\lambda}
```

so a plot of $x_i$ against $1/\cos\theta$ is a straight line through the origin whose slope is $d/\lambda$, with
$\lambda$ the attenuation length of the overlayer electrons in the overlayer. The tab fits the line through the origin
by least squares and reports the thickness and its standard error:

```math
k = \frac{\sum_i u_i x_i}{\sum_i u_i^{2}},\qquad d = k\,\lambda,\qquad \sigma_k = \sqrt{\frac{\sum_i \left(x_i - k\,u_i\right)^{2}}{(n-1)\sum_i u_i^{2}}}
```

The error covers the scatter between angles only, not the uncertainty in $\lambda$ or $R_\infty$, which are usually larger.
Two cross-checks are shown with it:

- A **free-intercept fit** $x = c + k\,u$. A line that does not pass through the origin means the windows, $R_\infty$ or
  the angle offset are wrong, or the layer is not uniform (for example contamination on top, or roughness). The tab
  warns when the intercept is clearly non-zero.
- For two states of one peak, an **$R_\infty$-free fit** of $R_i = R_\infty\left(e^{k u_i}-1\right)$ for both $k$ and
  $R_\infty$. It shows how much the answer depends on the densities, IMFPs and sensitivity factors in the presets.
  When $R$ is almost proportional to $1/\cos\theta$ the data cannot separate the two parameters, and the tab says so.

### Partial coverage (two states of one peak)

The equations above assume a complete, uniform layer. If the overlayer covers only a fraction $f$ of the surface and the rest
is bare substrate, the ratio for equal attenuation lengths becomes [23] (after Flinn and McIntyre and McKeown and co-workers)

```math
R = R_\infty\,\frac{f\left(1-e^{-x}\right)}{(1-f)+f\,e^{-x}}\quad\Longleftrightarrow\quad
e^{-x} = \frac{R_\infty f - R(1-f)}{f\,(R + R_\infty)},\qquad x = \frac{d}{\lambda\cos\theta}
```

which reduces to $\ln(1+R/R_\infty)$ at $f=1$. Because $x$ is still proportional to $1/\cos\theta$, the through-origin fit and
its error are unchanged; only the conversion from $R$ to $x$ differs. However thick the patches, $R$ cannot exceed
$R_\infty f/(1-f)$, and a larger ratio is refused for that coverage. Treating a patchy layer as complete gives a thickness
that is too small (for $f=0.7$ the example in the tests returns 1.04 nm for a true 2.0 nm) and a line that does not pass through the
origin. $f$ is usually not known; McKeown and co-workers varied it until the thickness no longer changed with angle, and the
*Coverage f* box does the same when it is set to `fit`: it chooses the $f$ (0.05 to 1) that minimises the relative scatter of
$x_i/u_i$ and reports that scatter next to the one for a complete layer. With few angles, or data that a complete layer
already fits, $f$ is poorly determined and the tab warns. The $R_\infty$-free fit assumes a complete layer, so it is left out when
$f \ne 1$, and the model applies only to two states of one peak (equal attenuation length), not to the two-peaks mode.

### Peak areas

By default the areas come from the spectra themselves, with no fitting of line shapes. (The alternative, using a
peak fit you made in CasaXPS, is described [below](#peak-areas-from-fitted-components).)

1. **Alignment.** Charging and the flood gun move the peaks by a few tenths of an electronvolt between angles. Each
   spectrum is shifted to match the previous (already aligned) one by minimising the squared difference between the
   two line shapes inside the envelope, each scaled to unit height, then refined below the step size. Matching shapes
   rather than peak maxima works when the taller peak is the overlayer, and for shoulders. The shifts are listed in
   the table.
2. **Background.** A Shirley background [15] is computed over the envelope (two states) or each window (two peaks), with
   the ends set by the mean of the first and last three points:

```math
B(E) = I_{low} + \left(I_{high}-I_{low}\right)\frac{\int_{E_{low}}^{E}\left[I(E^{\prime})-B(E^{\prime})\right]dE^{\prime}}{\int_{E_{low}}^{E_{high}}\left[I(E^{\prime})-B(E^{\prime})\right]dE^{\prime}}
```

3. **Areas.** The background-subtracted signal is integrated. For two states of one peak a single background is used
   across both, and the signal is split at one energy between them (the valley between the peaks).
4. **Automatic windows.** The substrate state is found as the place where the normalised spectrum falls from the lowest
   to the highest angle, and the overlayer state where it rises, near the binding energies of the presets. This tells
   them apart even when the oxide peak is the taller one.

Counts are divided by scans times dwell time (the K-Alpha stores raw counts summed over the scans) and by the
transmission function of the file. Both cancel for two states of one peak but matter when two different peaks are
compared.

### Peak areas from fitted components

If you have already fitted the spectra in CasaXPS [17] and saved the file as VAMAS, the tab can use that fit instead of
integrating windows. CasaXPS stores the fit in each block's comment lines: for every component its name, line shape,
area, width, position and an INDEX that links components belonging together, such as the two spin-orbit partners of a
doublet. The tab reads these lines with a parser ported from the author's eXPoSe SpectraDeck project [18].

- **Areas are used as stored.** They are counts per second times eV as CasaXPS computed them. Nothing is re-integrated,
  no peak shape is redrawn and no background is recomputed, so the result is exactly what your fit says.
- **States.** Components that share an INDEX form one chemical state and their areas add up; a component with no
  shared INDEX is a state on its own. Sensitivity factors of single components are never used (a linked doublet partner
  often stores zero).
- **Which state is the substrate.** The tab compares each state's share of the total fitted area at the lowest and
  the highest angle:

```math
\Delta_s = \frac{A_s(\theta_{max})}{\sum_j A_j(\theta_{max})} - \frac{A_s(\theta_{min})}{\sum_j A_j(\theta_{min})}
```

  A state whose share falls by more than 0.02 is taken as substrate, one whose share rises by more than 0.02 as
  overlayer, and one that hardly moves is left unassigned. This works whichever state is the larger one. Every
  assignment can be changed in the Components list; several states can belong to the overlayer, for example the
  sub-oxides of silicon.
- **Ratios.** With $A_{sub}$ and $A_{ov}$ the sums of the assigned areas at one angle, $R_i = A_{ov}/A_{sub}$ goes into
  the equations above exactly as a ratio of integrated windows does, so the line fit, the intercept check and the
  azimuth handling are the same.
- **Two different peaks.** The substrate states are taken from one fitted region and the overlayer states from another,
  matched by tilt. The stored areas are raw counts, so with *Transmission-correct* each component area is divided by the
  transmission function of its block at the kinetic energy of that component; the sensitivity factors are filled in
  from the fit and can be edited. For two states of one peak the correction cancels to within about 1 percent.
- **Charge correction.** A `Calib` line in the comments, which CasaXPS writes when you shift the energy scale, changes
  the displayed binding energies but not the areas, so it does not change the thickness.

Because the fit already separates overlapping components, this route does not suffer from the window-split limit
below. The thickness is only as good as the fit, so check it in CasaXPS: a poorly constrained fit gives poorly
constrained areas.

### Limits

- **Overlapping states.** A window split cannot separate peaks that are less than about 2 eV apart or that merge,
  as the oxide and substrate components of Ga 3d on an oxidised GaAs wafer do. The tab warns. The tab does not fit peaks
  itself, but it can use a fit made in CasaXPS (see above), which does not have this limit.
- **Grazing angles.** Above about 70 degrees elastic scattering, roughness and surface contamination make the simple
  exponential unreliable, and the effective attenuation length differs from the IMFP [13]. Those angles are drawn
  but excluded from the fit (the limit is editable).
- **Non-uniform layers.** Islands, roughness and an over-layer of adventitious carbon all bend the line. The ratio of
  two states of one peak is largely unaffected by a carbon layer on top, because both are attenuated alike, but the
  thickness and the intercept check then describe the oxide only.
- **Angle convention.** The emission angle is the stage tilt plus an offset (default 0) and is measured from the
  surface normal, as on the Thickogram tab. The K-Alpha lens axis is along the normal at zero tilt; for an instrument
  where it is not, set the offset. A tilt series is also not strictly depth profiling: it needs the sample to be
  flat and to stay at the analysis position as it tilts. ISO 14701 [14] describes the method for silicon oxide on silicon.

### Azimuth

If the file records the sample azimuth (rotation about the surface normal) as well as the tilt, the tab reads both,
groups the blocks by azimuth, fits each azimuth separately and reports the mean thickness and the spread between
azimuths. A plot of thickness against azimuth is drawn in polar form. The thickness of a uniform layer does not depend
on azimuth, so a systematic variation points to something else: a crystalline substrate, where photoelectron diffraction
and forward focusing change the intensity ratios with azimuth and polar angle, a structured surface, or a sample
that moves off the analysis spot as it rotates. The suite reports the variation and does not correct for it.

### VAMAS files

The tab reads VAMAS files (ISO 14976, the format described in [16]) with experiment mode NORM and a regular energy
axis, which is what the Thermo K-Alpha writes. Unspecified values are stored as 1e37 and are treated as missing. The
file does not say what its experimental variables mean, so the tab lets you choose where the tilt and the azimuth come
from; it defaults to the first experimental variable for the tilt, which is where the K-Alpha stores it. Binding
energy is taken as the photon energy minus the kinetic energy. Files saved from CasaXPS keep the fit in the block
comments, which the tab reads as described under peak areas from fitted components; positions in those lines are
kinetic energies.

## Depth profiles (Depth profile tab)

The ARXPS tab answers one question, the thickness of one overlayer. The Depth profile tab asks where **several species**
sit in depth. It follows the ideas of the NPL ARCtick package [12], implemented independently from the published physics.

### Species

A *species* is anything whose intensity can be followed against emission angle. Two sources are offered:

- **Fitted states.** Every chemical state of every CasaXPS fit stored in the VAMAS file (see
  [peak areas from fitted components](#peak-areas-from-fitted-components)). The stored areas are used as they are, and
  each component area is divided by the block's transmission function at its kinetic energy when *Transmission-correct*
  is ticked.
- **Whole regions.** Each region integrated with a Shirley background over the part of the spectrum that holds the
  peaks (the stretch where the smoothed spectrum stays above a tenth of its range, plus 1.5 eV each side), each spectrum
  aligned to the previous one. A Shirley background over a wide region with a long loss tail can fail to converge, so the
  Shirley routine now solves for the background directly and falls back to a straight line when no sensible solution exists.

Every species has an attenuation length $\lambda$ for its photoelectrons. The starting value is a TPP-2M estimate for an
oxide-like matrix, which is only a guide: replace it for your material. Counting statistics give each area an uncertainty,
$\sigma = \Delta E\,\sqrt{\sum_i N_i}\,/\,(n_{scans}\,t_{dwell})$ over the integrated points.

### Surrey Trend Plot

For a thin layer of a species at depth $z$ the intensity falls with angle as

```math
I(\theta) = I(0)\,\exp\left(-\frac{z}{\lambda\cos\theta}\right),\qquad \ln I = \ln I(0) - \frac{z}{\lambda}\,\sec\theta
```

so a plot of $\ln I$ against $\sec\theta$ is a straight line whose slope is $-z/\lambda$. The Surrey Trend Plot, suggested by
J. Castle and J. Watts of the University of Surrey and included in ARCtick, plots exactly this for all species at once. It is
exploratory, not a depth profile:

- **Parallel lines** mean species that fall off together. If their attenuation lengths are similar they may be
  distributed together in depth. The table lists, for each species, the others whose slope differs by less than twice the
  combined standard error and whose kinetic energy is within a factor of 0.7 to 1.4.
- **The scatter** of the points about each line (the rms residual of $\ln I$) shows the precision of the measurements.
- **The slope** gives an apparent depth $z_{app} = -\lambda\times\mathrm{slope}$. It is exact for a thin layer; for a species
  spread over a thickness it is a weighted mean depth, and a line that curves says the species is not in a single thin layer
  (an overlayer, whose intensity is $1 - e^{-x}$, curves upwards).
- The comparison is less useful when kinetic energies, and hence attenuation lengths, differ a lot. The table says so.

Curves can be normalised to their first angle so that parallel lines lie on top of each other. Angles above the maximum
angle are drawn as open symbols and left out of the fits.

### Geometry

An absolute intensity also carries any change of the analysed area with tilt, a factor common to every species. Ratios of
species at one angle (as on the ARXPS tab) do not need it, and parallel lines on the Trend Plot stay parallel, but depths and
amounts do. Cumpson [12] distinguishes two instrument types. If the **X-ray beam** defines the area, it does not change with
tilt and nothing needs correcting (the default). If the **analyser** defines it (its field of view), the area on the sample is
$A_0\sec\theta$, so the raw intensity rises as $\sec\theta$; the *Geometry* option then multiplies every intensity by
$\cos\theta$. This is the convention of the ARCtick spreadsheet too (its Stratification, AR/SR/XPS and Regularisation
sheets all differ between the two settings by $\sec\theta$). If every species, the substrate included, falls with angle much
faster than a thin layer can explain (a species placed deeper than the whole sample is thick), the analysed area probably
changes with tilt and this setting is worth checking. On a K-Alpha (beam-defined) tilt series of a silicon wafer the oxygen
1s signal is placed 2.6 nm deep without a correction, outside a 1.7 nm oxide, and 0.6 nm deep when every intensity is divided by
$\cos\theta$. That third option (*divide by cos θ (empirical: tilted beam)*) is the opposite of the analyser correction and is
not from ARCtick or the paper: it is kept because it removed an impossible result on that instrument, perhaps because the beam's
flux per unit area falls on tilt, and it should be judged against your own instrument.

### Stratification

A thin layer of areal amount $n$ at depth $z$ gives $I(\theta) = S\,n\,\exp(-z\sec\theta/\lambda)$, with $S$ the sensitivity
factor. Two angles $a$ and $b$ therefore give the depth and the amount:

```math
R = \frac{I(\theta_b)}{I(\theta_a)},\qquad z = -\frac{\lambda\,\ln R}{\sec\theta_b-\sec\theta_a},\qquad n = \frac{I(\theta_a)\,\exp\left(z\sec\theta_a/\lambda\right)}{S}
```

The amount is the intensity corrected for the attenuation by the material above, divided by the sensitivity factor; the tab
reports each species' share of the total. Measurements at only normal emission and one oblique angle (about 67.5 degrees is
usual) are enough. The uncertainty of the depth comes from the counting statistics of the two areas:

```math
\sigma_z = \frac{\lambda}{\sec\theta_b-\sec\theta_a}\sqrt{\left(\frac{\sigma_a}{I_a}\right)^2+\left(\frac{\sigma_b}{I_b}\right)^2}
```

With *Use every angle* the line $\ln I = \ln(S\,n) - (z/\lambda)\sec\theta$ is fitted instead, which gives the same depth and
amount for a thin layer and an error from the scatter. The result is drawn as bars: position is the average depth and width
the relative amount.

Limits: a species that is spread over a thickness, or sits in two layers, is reported at one average depth. The tab tests the
curvature of $\ln I$ against $\sec\theta$ (the t statistic of the quadratic term) and warns when it is large. A species whose
intensity rises with angle lies at the surface and gets a negative depth with a note. A substrate that extends to infinite
depth is not a thin layer: its signal is reported at its information depth and the curvature warning usually fires.

### Relative depth plot

The relative depth plot (RDP) is the simplest angle-resolved analysis: it says **which species lie nearer the surface**, and
nothing else [22, 23]. It has no depth scale and no thickness. For each species it compares the signal at a high (near-grazing)
and a low (near-normal) angle; the nearer the surface a species sits, the more it gains with angle, so the larger its value.
Plotted with the surface at the top, the species fall in the order of the layers (a self-assembled monolayer on silver, for
example, puts Ag 3d lowest, then S 2p, then the carbon and oxygen of the chain [22]). The tab offers two forms:

```math
\mathrm{RDP}_{\mathrm{pub}} = \frac{\ln p(\theta_{high})}{\ln p(\theta_{low})},\qquad
\mathrm{RDP}_{\mathrm{ratio}} = \ln\frac{I(\theta_{high})}{I(\theta_{low})} = -\frac{(\sec\theta_{high}-\sec\theta_{low})\,z}{\lambda}
```

The first is the definition in the papers, a ratio of logarithms; a logarithm needs a dimensionless argument, so $p$ is the
species' atomic percentage (its $I/S$ as a share of all the species shown) at each angle, which also removes any factor common
to every species, such as a change of the analysed area with tilt. It is undefined where a share is about 1 % or less, because the
logarithm is then zero or negative, and the tab says so. The second is unit-free, is the Stratification depth in units of $\lambda$
with the sign reversed, and so always orders the species as Stratification does; a factor common to every species only shifts
all the values together. Uncertainties come from the counting statistics when they are known (for the first form through the
shares, so the other species enter).

Limits, from the papers: it works only if the sample really is layered; a species present in two layers is averaged, which the
curvature warning of Stratification also reports; and species with widely separated kinetic energies have different attenuation
lengths, which the plain index ignores, so the tab flags large kinetic-energy differences and Stratification (which scales by
$\lambda$) should be used to compare them. Brundle and co-workers [22] note Cumpson's point that it is neither a depth
profile nor a measure of layer thickness. The angles are those chosen on the Stratification page.

### Layers (Gries signal ratios)

For a stack of layers $1\ldots K$ with thicknesses $d_k$ on a substrate, the intensity of a species that sits in layer $j$
(unit concentration, attenuation length $\lambda_{j}$ in each layer) is attenuated by every layer above it:

```math
I_j(\theta) = \exp\left(-\sum_{m<j}\frac{d_m}{\lambda_m\cos\theta}\right)\left[1-\exp\left(-\frac{d_j}{\lambda_j\cos\theta}\right)\right]
```

and a substrate species gives $\exp\left(-\sum_{m=1}^{K} d_m/(\lambda_m\cos\theta)\right)$. W. Gries [19] noted that if every
measured curve and every prediction is divided by its own value at normal emission, **no sensitivity factors are needed**: only
the shape of each peak's intensity against angle is compared, and the thicknesses are adjusted until the predicted curves follow
the points. The tab draws the predicted curves over the measured ones with a slider for each thickness, so you can do that by
eye, and can also fit them: the thicknesses minimise

```math
\chi^2 = \sum_i\sum_{\theta}\left(\frac{y_i(\theta)-p_i(\theta)}{\sigma_i(\theta)}\right)^2
```

with $y_i$ the measured and $p_i$ the predicted normalised curve, using a bounded simplex search started from several
points. The uncertainties come from the counting statistics (3 percent of the value when they are unknown).

**Range of validity.** For each layer the thickness is held at each value on a grid while the others are re-fitted, and the range
where $\chi^2$ stays within one unit of its minimum is reported (scaled up when the reduced $\chi^2$ is above 1, so
unrealistically small errors do not give unrealistically narrow ranges). A wide range means the data cannot say: information fades
with depth, and the deeper layers are always less well determined than the shallow ones.

Limits and cautions:

- Each species must belong to **one** layer (or the substrate). A region that mixes several states, such as a whole Si 2p
  spectrum holding silicon and its oxide, does not belong in the stack: use the fitted states.
- The model assumes uniform, flat layers and a known geometry. A tilt-dependent change of the analysed area (see Geometry
  above) acts on every curve; the cosine option is a guess, and a poor fit (a reduced $\chi^2$ far above 1) is the signal that
  something is missing from the model. On a K-Alpha tilt series of a silicon wafer with a carbon layer, an oxide and silicon the
  stack did not fit (reduced $\chi^2$ in the thousands, even with the geometry exponent free), so for a single
  oxide thickness the ratio of two states on the ARXPS tab, which cancels all common factors, is the more robust route.
- The tab uses one attenuation length per species in every layer. In reality it depends on the material it passes
  through, so using each layer's own value would change the result.

### Regularised depth profile

When the species are not in discrete layers, the concentration can be a smooth function of depth. For one species with
sensitivity factor $S$ and attenuation length $\lambda$, divide the depth into slabs $n$ with concentrations $c_n$; the
intensity at emission angle $\theta_m$ is

```math
\frac{I(\theta_m)}{S} = \sum_n c_n\,A_{mn},\qquad A_{mn} = \lambda\cos\theta_m\left[e^{-z_n^{-}/(\lambda\cos\theta_m)}-e^{-z_n^{+}/(\lambda\cos\theta_m)}\right]
```

where the slab runs from depth $z_n^{-}$ to $z_n^{+}$. A handful of angles cannot determine many slabs, so the profile is found
by Tikhonov regularisation: the weighted misfit plus a penalty on the roughness of the profile (its second difference) is
minimised,

```math
\min_{c}\ \sum_m\left(\frac{\sum_n A_{mn}c_n - b_m}{\sigma_m}\right)^{2} + \alpha\sum_n\left(c_{n-1}-2c_n+c_{n+1}\right)^{2}
```

with $b_m = I(\theta_m)/S$ and $\sigma_m$ from the counting statistics. A very large $\alpha$ gives a straight line and a very
small one reproduces the noise as wild oscillations. By default $\alpha$ is chosen from the data as the **smoothest profile
whose misfit is no larger than the noise allows** (the discrepancy principle); generalised cross-validation is available in
the code but under-smooths badly with so few angles, which is why it is not the default. The result carries a depth
resolution of roughly $\lambda\cos\theta_{max}$, and the tab reports it.

Prior knowledge sharpens the profile, and it is the most valuable input because the angle-resolved data hold little
information:

- **Non-negative.** Concentrations are forced to be zero or positive (a Lawson-Hanson active-set solution). A profile that
  needs negative concentrations to fit is flagged: either an assumption (attenuation length, geometry, which species) is
  wrong, or the measurement has lost the information. Looking at the unconstrained profile first is a useful check.
- **Constant, or zero, below a depth.** A bulk species such as the metal under an oxide has one concentration below some
  depth, and a surface species has none: the slabs deeper than that depth are merged into one parameter, or removed.

The noise bands drawn on the plot are the covariance of the linear estimator,

```math
\mathrm{cov}(c) = G\,G^{T},\qquad G = \left(A_w^{T}A_w + \alpha L^{T}L\right)^{-1}A_w^{T},\qquad A_w = A/\sigma
```

with $L$ the second-difference matrix. They show how far counting noise alone moves the profile; they do **not** include the
bias that smoothing always adds, so a narrow band does not mean an accurate profile. The species are solved one at a time;
*Atomic fraction* rescales them to add to 1 at each depth for display only.

### Three-point profile (the ARCtick method)

The *Method* box on the same tab also offers the method of the ARCtick spreadsheet [12], which treats the problem differently:
instead of smoothing a fine grid, it uses a model with only as many parameters as the data can support. Cumpson shows that the
depth resolution is limited by the signal-to-noise ratio, not by the number of angles, and that five angles of 5% precision hold
about three degrees of freedom. The highest generalised frequency the data can resolve, for $n$ angles of relative uncertainty
$\sigma/I$, is

```math
\omega_{max} = \frac{1}{\pi}\cosh^{-1}\!\left[\pi (n-1)\left(\frac{I}{\sigma}\right)^{2}\right]
```

(2.93 for $n=5$ and 5%, the value in the spreadsheet's K27 cell), which fixes three knots spaced geometrically,

```math
z_1 = \lambda\cos\theta_{max},\qquad z_2 = z_1 e^{\pi/\omega_{max}},\qquad z_3 = z_2 e^{\pi/\omega_{max}}
```

The profile (Paynter's model) is $c_1$ down to $z_1$, straight lines to $c_2$ at $z_2$ and $c_3$ at $z_3$, then constant at $c_3$.
The intensity is linear in the three concentrations, so they come from weighted linear least squares with no smoothing
parameter; the covariance $(A_w^T A_w)^{-1}$ gives the noise bands. If a species is *zero below z*, the concentration is zero
below the third knot instead of continuing at $c_3$. *Auto alpha*, *Non-negative* and the number of depth points do not apply.
Where the two methods differ is a smooth profile (Tikhonov) against a three-number description (three-point); both agree where
the data are informative, and any feature narrower than the resolution is not in the data under either.

*Not identical to the spreadsheet.* The knots, the model and the least-squares solution follow Cumpson's paper and the
spreadsheet's formulas, but the spreadsheet's integrals are in a VBA module that could not be read, so its exact treatment of
the bulk tail and of its depth-scale slider is inferred from the paper, and the Suite has no equivalent of the slider.

### Profile (Murdoch)

Murdoch and McCulloch [20] fit the **apparent atomic concentrations** rather than the intensities. For peak $j$ at emission
angle $\theta_k$ the measured value is $X_{jk}$, the share of peak $j$ in the sum of all the peak areas divided by their
sensitivity factors. Any factor common to every peak (flux, analysed area, the tilt geometry) cancels in that share, which
is the method's main attraction. With the sample divided into slabs $i$ at depths $z_i$ and concentrations $c_{ij}$ in at%,

```math
I_{jk} = \sum_i c_{ij}\,e^{-z_i/(\lambda_j\cos\theta_k)},\qquad I^{\infty}_{jk} = \sum_i e^{-z_i/(\lambda_j\cos\theta_k)},\qquad X^{calc}_{jk} = \frac{I_{jk}/I^{\infty}_{jk}}{\sum_l I_{lk}/I^{\infty}_{lk}}
```

where $\lambda_j$ is the attenuation length of peak $j$ and $I^{\infty}$ the signal of a uniform sample. The profile minimises

```math
Q = \chi^{2} + \alpha R,\qquad \chi^{2} = \sum_{j,k}\left(\frac{X^{calc}_{jk}-X_{jk}}{\sigma}\right)^{2},\qquad R = \sum_{i,j}\frac{(c_{i+1,j}-c_{ij})^{2}}{\Delta}
```

with the concentrations in at%, the slab thickness $\Delta$ in angstrom and $\sigma = 1$ for fractions, the convention in which the
paper's $\alpha$ of about $10^{-6}$ to $10^{-4}$ makes sense. The penalty is the squared concentration **gradient**, which unlike
a penalty on the squared concentration does not push the substrate towards zero. Here the concentrations are written as a
softmax over the peaks in each slab, so they add to exactly 100 at% and stay positive, and Q is minimised with its analytic
gradient (an L-BFGS iteration from the initial guess).

The result depends on the **initial guess**, a layer stack with guessed thicknesses, which is how prior knowledge enters.
*Initial-guess ensemble* repeats the fit from guesses whose thicknesses vary by up to 50 percent and plots the mean and
spread. The **S-curve** is $\chi^{2}$ against $\log_{10}\alpha$; the suggested $\alpha$ is the largest one whose $\chi^{2}$ is
within 25 percent of the minimum (a rule of thumb that automates reading the curve, not part of the paper). The *Murdoch EAL*
button sets every attenuation length from the kinetic energy alone, $143/E^{2} + 0.054\sqrt{E}$ nm, the semi-quantitative
approximation in the authors' spreadsheet (after Seah and Dench [21]); use material specific lengths where you have them.

Differences from the authors' Excel implementation: their spreadsheet optimises free amplitudes per slab and rescales to 100
at% afterwards. That is not well posed, because the apparent concentrations do not depend on a common scale, so the
penalty (which falls with the scale) drives all amplitudes to zero when they are free to move; the closed form above does not
have that degeneracy. Fed the spreadsheet's example (WS2 on Al2O3, five angles) with its attenuation lengths, the forward model
reproduces its stored solution's chi-squared (0.0026 against 0.0023), and the closed fit recovers the same stack: carbon on
top, a sulfur and tungsten layer about 1 nm thick, then oxygen and aluminium. Five angles cannot resolve a profile
below about $\lambda\cos\theta_{max}$, so features narrower than that are the regulariser and the guess talking.

### ARXPS MEMS (Rossi)

The maximum-entropy method (MEM) of Scorciapino, Navarra, Elsener and Rossi [27] finds the depth profile that satisfies the
angle-resolved data but contains the least structure: the entropy is the regulariser. It is built on the same model as the Profile
(Murdoch) page, apparent concentrations with a closure over the species, but replaces the roughness penalty by an entropy and lets the
attenuation lengths follow the composition. The sample is cut into layers of thickness $t$ (1 Å in the paper) at depths $z_i=i\,t$
with the atomic fraction $n_{j,i}$ of species $j$ in layer $i$ and $\sum_j n_{j,i}=1$. The intensity of species $j$ at emission angle
$\theta$ is $I_j=k_j\sum_i n_{j,i}\,e^{-z_i/(\lambda_{j,i}\cos\theta)}$, and that of the pure species $I_j^{\infty}=k_j\sum_i
e^{-z_i/(\lambda_{j,i}\cos\theta)}$, so the apparent concentration, which is what the sensitivity-corrected intensities give, is

```math
X_j(\theta)=\frac{I_j/I_j^{\infty}}{\sum_l I_l/I_l^{\infty}}
```

independent of the $k_j$ and of any factor common to all peaks. With $m_{j,i}$ the initial estimate of the profile the entropy,
the misfit and the probability function that is maximised are

```math
S=\sum_{j,i}\left[n_{j,i}-m_{j,i}-n_{j,i}\ln\frac{n_{j,i}}{m_{j,i}}\right],\qquad
C=\chi^2=\sum_{j,k}\frac{\left(X^{calc}_{j,k}-X^{obs}_{j,k}\right)^2}{\sigma_{j,k}^2},\qquad Q=\alpha S-\frac{C}{2}
```

$S$ is zero at $n=m$ and negative elsewhere, so a large $\alpha$ keeps the profile at the start model (and does not fit the data), while
a small $\alpha$ fits the data, noise included. The paper lets $\alpha$ fall from 100 towards 1 over successive calculation cycles
(0.1 only for simple structures), each started from the previous solution. Here $Q$ is maximised by minimising $C/2-\alpha S$ over
softmax variables, so that every layer sums to 1 and no fraction is negative, with an analytic gradient and the L-BFGS minimiser of the
Murdoch page; the cycles stop at the first (largest) $\alpha$ whose $\chi^2$ per point is within the noise if *Stop at χ²/N ≤ 1* is ticked.
$\sigma$ is a relative error of every apparent concentration (6 % is about the standard deviation of the paper's ±10 % uniform error).

**Attenuation lengths that follow the composition.** The IMFP of the electrons of species $j$ in layer $i$ is the thickness-weighted mean
of their IMFPs in all the materials they cross (the paper's Eqs. 4-7), so $\lambda_{j,i}=\frac{1}{i+1}\sum_{l=0}^{i}\sum_m n_{m,l}\,\Lambda_{j,m}$,
with $\Lambda_{j,m}$ the IMFP of species $j$'s electrons in the material of species $m$. It depends on the profile, so it is recomputed
from each cycle's solution until no length changes by more than 1 %. In a synthetic test this matters: with lengths that follow the
composition, a fit that keeps them fixed misplaces a 1.6 nm overlayer by 9 % and its composition by 19 percentage points, while the
iteration recovers both to within 3 % and 0.4 percentage points. $\Lambda_{j,m}$ comes from the **G-1 equation** of Gries [28], which
the paper uses through NIST SRD 71 [29]:

```math
\lambda_{G1}\,[\mathrm{nm}]=k_1\,\frac{V_a}{Z^*}\,\frac{E}{\lg E-k_2},\quad V_a=\frac{\sum_i q_iM_i}{\rho\sum_i q_i},\quad
Z^*=\frac{\sum_i q_i\sqrt{Z_i}}{\sum_i q_i}
```

for a material of formula $\mathrm{A}_{q_1}\mathrm{B}_{q_2}\dots$, density $\rho$ (g cm$^{-3}$) and an electron of kinetic energy $E$ (eV); $(k_1,k_2)$ is
(0.0020, 1.30) for the 3d elements Ti-Cu, (0.0019, 1.35) for 4d Zr-Ag, (0.0019, 1.45) for 5d Hf-Au, (0.0014, 1.10) for the other
elements, (0.0018, 1.00) for organic and (0.0019, 1.30) for inorganic compounds. The tab asks for a formula and a density per
species (presets for common materials, or type them) and reproduces the G-1 values of Gries' Table 3 for water. G-1 was fitted to
200-2000 eV (the tab says when an energy lies outside), and NIST cautions that it can deviate by up to about 50 % for some materials
and that TPP-2M is generally preferred; the alternative on the page uses each species' own $\lambda$ from the species table scaled
to the electron energy as $E^{0.75}$.

**The starting profile is the prior.** The paper stresses that the reconstruction is only as good as the starting profile (the
version of its algorithm with free layers failed from a flat start on structures with seven or eight components). Its protocol
finds the start model with a simulator and refines it with a four-parameter pseudo-Gaussian version of the algorithm before
the free-layer version used here; this page implements the free-layer version only, with the start model entered as a table of
layers (thickness and composition in %, the last row being the bulk, by default one species per layer in the order of the
Stratification depths) and used as the entropy prior $m$. Layers below the *Free to* depth are held at the bulk composition,
as in the paper's constrained bulk. On the paper's synthetic structures (7 angles from 0 to 60 degrees, which the text does not
give; the IMFPs, layers and compositions are those of its Table 1) with start models whose layer thicknesses are 10 % off, the
overlayer thickness is recovered to within 12 % and the compositions of the upper layers to within 16 percentage points,
with or without ±10 % random error, as the paper reports (thickness ±10-20 %, composition ±10-30 %); a 0.6 nm pure layer between
two thicker ones is at the limit of what the fit resolves. From a start 30 % off the interfaces follow the start model rather than
the data, which the start-model ensemble (mean and spread over 8 starts) makes visible. The page also draws the paper's relative
depth plot, the ratio of the model's apparent concentration at the largest angle to that at the smallest, and reports each species'
curve area (nm of that species), peak depth and the interfaces where the dominant species changes.

### Attenuation lengths and the range of validity

Every method here depends on the attenuation lengths. The starting values are TPP-2M inelastic mean free paths for an
oxide-like matrix, which is only a guide. Attenuation lengths, which include some of the effect of elastic scattering, are
typically 10 to 30 percent shorter than inelastic mean free paths [13]; the *AL / IMFP* box scales the starting values (0.8,
for example) and each length can be edited. Empirical estimates are good to about 6 to 20 percent, so the thickness and
depth scales inherit that uncertainty.

The simple exponential model is reliable up to emission angles of about 60 to 70 degrees; beyond that elastic scattering,
roughness and surface contamination matter (the *Max angle* setting leaves those angles out of fits), and a crystalline
substrate adds photoelectron diffraction that changes the intensities with polar and azimuthal angle. The analysis assumes
flat, uniform samples.

## Dual-anode XPS (DAXPES tab)

A single X-ray source cannot tell a thin film that covers the whole analysed area from a thicker one that covers only part of
it: both lower the film signal. Murdoch [24] separates the two by measuring the **same core level with a soft and a hard
source** (for example Al Kα, 1486.7 eV, and Ag Lα, 2984.3 eV, or Ga Kα, 9.25 keV). The hard-source electrons have a longer
attenuation length, so the film is seen through a different fraction of its thickness. For a film of thickness $d$ covering a
fraction $c$ of the analysed area, with $I_\infty$ the area of a thick reference of the same material measured with the same
source, each source gives

```math
\frac{I}{I_\infty} = c\left(1-e^{-d/(\lambda\cos\theta)}\right)\qquad\Longleftrightarrow\qquad
t = -\lambda\cos\theta\,\ln\left(1-\frac{I}{I_\infty}\right)\quad(c=1)
```

where $\lambda$ is the attenuation length of the electrons of that source, $\theta$ is the take-off angle from the surface normal
and $t$ is the thickness a single source would report if the film were complete. The coverage $c$ is a common factor, so it
cancels in the ratio of the two sources (X = soft, H = hard), which fixes the thickness on its own [24]:

```math
\frac{I_H/I_{\infty,H}}{I_X/I_{\infty,X}} = \frac{1-e^{-d/(\lambda_H\cos\theta)}}{1-e^{-d/(\lambda_X\cos\theta)}},
\qquad c = \frac{I_X/I_{\infty,X}}{1-e^{-d/(\lambda_X\cos\theta)}}
```

The ratio on the left rises monotonically from $\lambda_X/\lambda_H$ for a vanishingly thin film to 1 for a film thicker than the
soft source can see, so the method is sensitive only in between, and the tab solves the equation for $d$ by bisection and
refuses a ratio outside that range, saying why. The coverage is what the single-source areas leave over; the paper also gives
it as the ratio of the hard-source single-source thickness to the DAXPES thickness (the tab shows both, and they agree for thin
films). Only the ratio $I_{\infty,H}/I_{\infty,X}$ is needed for the thickness: it can be measured on a thick reference of the
same material with both sources or typed in, but without a reference for each source there is no coverage.

Assumptions and limits, from the paper. Both sources must see **the same analysed area**, which is easy on an instrument that
defines the area by an aperture or by imaging, and not otherwise (an analysis area set by the X-ray spot changes with the
source, and then the single-source areas are not comparable). The thickness cannot exceed about $3\lambda_H\cos\theta$, the
elastic-scattering limit of the hard-source electrons, and above about $3\lambda_X\cos\theta$ the soft source is saturated and
the ratio is close to 1, so the answer rests on a small difference; very thin films are difficult for the opposite reason, because
the sources differ little. The attenuation lengths matter, and the tab offers two starting points (use values you trust, because the thickness scales
with them). The default is Seah's universal curve [26], which needs only the kinetic energy $E$ of the electrons (hard-source
electrons are faster, so $\lambda_H > \lambda_X$), the mean atomic number $Z$ of the film, the mean atomic size $a$ and a small
term $W$ for strongly bonded compounds:

```math
\lambda = \frac{5.8 + 0.0041\,Z^{1.7} + 0.088\,E^{0.93}}{Z^{0.38}\,(1-W)}\,a^{1.82}\ \mathrm{nm},\qquad
a^{3} = \frac{10^{21}M}{\rho N_A (g+h)}\ \mathrm{nm}^3
```

with $M$, $\rho$ and $g+h$ the molar mass, density and atoms per formula unit, $W = 0.06H$ with $H$ the heat of formation in eV per
atom (0 for an element) and $E$ in eV; it agrees with attenuation lengths computed from optical data to an RMS of 8 % between
100 eV and 30 keV for emission up to 65 degrees [26]. The other is the TPP-2M IMFP times an AL/IMFP factor (0.9 by default). The
densities, compositions and heats of formation of the few built-in films (TiO$_2$, Ga$_2$O$_3$, SnO$_2$, SiO$_2$, Si) are handbook
values, not from the papers. Because the answer rests on the ratio of the two sources, an error in the *ratio* $\lambda_H/\lambda_X$ moves the thickness
and the coverage far more than an error in both lengths together (which scales $d$ by the same factor and leaves the coverage
unchanged), so the details box lists the thickness and coverage for $\lambda_H$ alone and for both lengths changed by $\pm 10$ %;
for the synthetic example a 10 % change of $\lambda_H$ alone moves 3.0 nm to 2.1-3.7 nm and the coverage from 70 % to 55 %,
while changing both by 10 % gives 2.7-3.3 nm and the same 60 %. Against the lengths in Murdoch's Table 1 the Seah curve is 18 % below to 6 % above, and
TPP-2M times 0.9 within 15 %, which is the scatter to expect from either; for another material, type the lengths in. The quoted error is the counting-statistics part only: on the synthetic sample the choice of window (half-width 3 to 7 eV) and the Shirley end points move the area-method thickness by about 6 % (2.8 to 3.2 nm for a true 3.0 nm), more than the 0.05 nm statistical error, so repeat the calculation with a few windows before quoting a precision. The auto window is the envelope of the signal near the peak energy, padded by 3 eV and never narrower than 6 eV either side of it (a doublet such as Ti 2p has its second component 5-6 eV away). The peak areas come from the same windows and
Shirley background as the ARXPS tab, with one window in binding energy for all four spectra; they may also be typed in.

### Background method

The paper's second method needs **no reference sample**: instead of peak areas it uses the Tougaard inelastic-background
parameter $B_1$ of the film spectrum of each source [24] (the idea of the background-based thickness of Murdoch and co-workers
[25]). The peak and its tail are described by Tougaard's equation, with $J$ the measured spectrum, $F$ the primary spectrum, $E$
the kinetic energy and $T$ the energy lost,

```math
F(E) = J(E) - B_1\int_E^{E_{max}} J(E')\,K(E'-E)\,dE',\qquad K(T) = \frac{T}{(C-T^2)^2 + D\,T^2}
```

with $C=-1643$ eV$^2$ and $D=0$ for metal oxides and $C=325$ eV$^2$, $D=96$ eV$^2$ for silicon. $B_1$ (in eV$^2$) is the one number
that makes $F$ zero at a point $\Delta E$ (30 eV in the paper; 15-30 is allowed) below the centroid of the peak in kinetic
energy, that is above it in binding energy. It is a measure of the fraction of the photoelectrons that have lost about one
quantum of energy, so it grows with the film thickness and is smaller for the faster electrons of the hard source. In the
model of the paper $B_1/B_0 = 1-e^{-d/(2\lambda\cos\theta)}$, with $B_0$ the value for a thick sample, and $B_0$ cancels in the ratio
of the two sources, which gives $d^*$ (equal to $d/(2\cos\theta)$ when elastic scattering is ignored) and then the thickness:

```math
\frac{B_1^{H}}{B_1^{X}} = \frac{1-e^{-d^*/\lambda_H}}{1-e^{-d^*/\lambda_X}},\qquad
d = 3\lambda_H\cos\theta\left(1-e^{-d^*/(\lambda_H\cos\theta)}\right)
```

The first equation has the form of Eq. 6 with $\theta=0$, so the tab solves it with the same bisection and refuses a ratio outside
$(\lambda_X/\lambda_H, 1)$. The second converts $d^*$ to a thickness and caps it at the elastic limit $3\lambda_H\cos\theta$. It is
the equation of the paper, and it reproduces the paper's Table 2 (sputtered Si on polystyrene, $\theta=0$, B$_1$ values from its
Fig. 9: 2.9, 13.8 and 15.6 nm against the published 2.7, 13.8 and 15.7 nm). Its small-$d^*$ behaviour is worth knowing: it gives about
$3d^*$ where the simple thin-film relation from the $B_1$ model is $2d^*\cos\theta$, so the tab also prints $2d^*\cos\theta$ (the
two readings cross near $d^*\approx 0.86\,\lambda_H$), and the paper reports that the method underestimates ultra-thin films (below about
3 nm) in its simulations. How $\cos\theta$ enters Eqs. 7 and 8 off normal emission could not be verified from the paper (its
Fig. 4 caption prints $d^*=2d\cos\theta$ against the text's $d/2\cos\theta$), so the tab uses the equations as printed and the
check against Table 2 holds only at normal emission.

The tab fits $B_1$ from the spectrum itself: the zero level is the mean of the first points at the low-binding-energy end
(which must be flat and free of peaks), the spectrum is smoothed over five points for the fit, the centroid is found
self-consistently from $F$, and $B_1$ is the ratio of the (smoothed, zero-level-corrected) signal at $E_c+\Delta E$ to the
integral of Eq. 1 there. The spectrum must therefore run from a flat level before the peak to at least $\Delta E$ above the
centroid, and the peak window sets the range used for the centroid and the area. The uncertainty is the counting-statistics part
(a Monte Carlo over Poisson draws of every count); it does not include the choice of $C$ and $D$, of $\Delta E$, or of the
background range, which move $B_1$ much more. Because $B_1$ is relative to the peak area, it does not depend on the coverage, the
flux or the cross-section of either source, but the *coverage* does need a peak area and a thick-reference area (the same
window and Shirley area as the area method, files or typed), from $c=(I/I_\infty)/(1-e^{-d/(\lambda\cos\theta)})$.
Only elemental information comes out of the background: chemical states are not separated.

### Both methods together: the attenuation-length ratio

Both DAXPES methods rest on the **ratio** $\rho=\lambda_H/\lambda_X$ of the two attenuation lengths far more than on either length
alone, and this is where the universal curves disagree most. For Ti 2p in TiO$_2$ Seah's curve [26] gives 1.83 and 3.99 nm
($\rho=2.18$) and Murdoch's Table 1 [24] 1.87 and 3.77 nm ($\rho=2.02$): the soft lengths agree to 2 % and the hard ones differ by
6 %. A true 3.0 nm film (60 % coverage) then reads 3.5 nm (+17 %) by peak areas and 4.0 nm (+35 %) by the background (the model
reading $2d^*\cos\theta$), because $B_1$ behaves like a film of half the thickness, nearer the thin-film limit $\lambda_X/\lambda_H$
where the ratio is most sensitive. The background method removes the need for references and for equal flux, not the dependence
on $\rho$.

What the two methods can do together is **measure** $\rho$. With $u = d/(\lambda_X\cos\theta)$ the paper's models give

```math
q_A = \frac{1-e^{-u/\rho}}{1-e^{-u}},\qquad q_B = \frac{1-e^{-u/(2\rho)}}{1-e^{-u/2}}
```

for the peak-area ratio and the $B_1$ ratio: two equations for the two unknowns $u$ and $\rho$, well conditioned (condition number
about 10), so a film measured both ways fixes its thickness and its own length ratio with no assumed $\rho$; only $\lambda_X$ is
needed to turn $u$ into nanometres. In the example above the exact ratios return $\rho=2.016$ and $d=3.0\times1.83/1.87=2.93$ nm
(-2 %) from Seah's soft length, against +17 % from fixing the ratio. With realistic noise the scatter is the price: 1 % on both
ratios gives $\rho$ to about 5 % and $d$ to 10 %, 4 % on the $B_1$ ratio gives $\rho$ to 10-19 % and $d$ to 30 %. The *Both methods*
choice therefore fits $(u,\rho)$ to the two ratios with a Gaussian prior on $\rho$ centred on the lengths you entered (width
set in the Method page, 10 % by default), reports the thickness and the length ratio the data prefer next to the one you entered, and
says when the data pull away from the prior. The quoted errors are the posterior errors (data and prior) from the Jacobian, counting
statistics only. It rests on the model of Eq. 5 (elastic scattering and the factor 2 are approximations) and on equal analysed areas,
so read the recovered $\rho$ as a diagnostic.

The tab also compares the two thicknesses like for like, the peak-area $d$ against $2d^*\cos\theta$ (Eq. 8's capped value is a
different, empirical mapping), and flags a disagreement beyond two combined errors or 15 % with the likely causes: the length
ratio, unequal analysed areas, or the $C$, $D$ and $\Delta E$ of the $B_1$ fit. The last is checked in a table: on the synthetic
sample the single-point $B_1$ ratio at $\Delta E$ = 15, 20 and 30 eV is 0.656, 0.725 and 0.664 (thickness 6.3, 7.8 and 6.5 nm by Eq. 8
for a true 6 nm film) because the hard-source tail is noisy at one point, and the silicon parameters instead of the oxide ones move it
to 7.6 nm. The optional **band fit** (least squares over 15 to 45 eV above the peak, which is the same number for a spectrum that obeys Eq. 1)
is about three times less noisy (hard-source $B_1$ scatter 9 against 27 eV$^2$), gives a ratio of 0.682 (true 0.687) and 6.9 nm, and is a
variant of the paper's definition, so it is off by default.

## Spectrum simulation

The Simulated spectrum pane draws what a given thickness would look like. It is an **estimate, not a fit**.

The thickness is turned into relative peak areas by solving the Strohmeier equation for the overlayer to substrate
ratio $q$, where $R$ are the RSFs:

```math
q = \frac{I_{ov}}{I_{sub}} = \frac{R_{ov}N_{ov}\lambda_{ov}}{R_{sub}N_{sub}\lambda_{sub}}\left[\exp\left(\frac{d}{\lambda_{ov}\sin\theta}\right)-1\right]
```

The substrate then carries a fraction $1/(1+q)$ of the total area and the overlayer $q/(1+q)$, shared between
overlayer components in proportion to the areas you entered. Simulating at the calculated thickness therefore
reproduces your entered area ratio exactly.

Each component is a pseudo-Voigt line of the given width, with a Gaussian fraction $(1-\eta)$ and a Lorentzian
fraction $\eta$ of equal full width at half maximum, and $\eta = 0.3$:

```math
(1-\eta)\,G(E) + \eta\,L(E)
```

Core levels with spin-orbit splitting are drawn as a doublet with the standard area ratio (2:1 for p, 3:2 for d, 4:3
for f). The background is a self-consistent Shirley step, rising across the peak in proportion to the peak area
accumulated at lower binding energy:

```math
B(E) = B_0 + \Delta B\,\frac{\int_{E_{min}}^{E} P(E')\,dE'}{\int_{E_{min}}^{E_{max}} P(E')\,dE'}
```

with $\Delta B$ equal to 15 percent of the peak height. Splitting, branching ratio and width are typical literature
values, not measured ones, and the widths are editable.

**What the simulation leaves out.** For simplicity it does not draw shake-up or shake-off satellites, plasmon or loss
features, Auger lines, or the asymmetric tails of metal peaks. The most important case is Cu 2p3/2: the Cu(II)
satellites of CuO and Cu(OH)2 carry roughly 35 to 39 percent of the Cu(II) area (main to satellite ratios of 1.89
and 1.57, from [5]), so the drawn Cu(II) main peak is overstated. The figure says so, and shades the satellite
region. If you use Cu 2p3/2 areas for a thickness, include the satellite in the Cu(II) area. Auger presets, such as
the Cu LMM line, are not simulated.

## IMFPs from TPP-2M

Where the suite *computes* an IMFP instead of quoting a published value, it uses the TPP-2M predictive equation of
Tanuma, Powell and Penn [6]. For a kinetic energy $E$ in eV it gives the IMFP in angstrom:

```math
\lambda = \frac{E}{E_p^{2}\left[\beta\ln(\gamma E) - C/E + D/E^{2}\right]}
```

with the parameters

```math
\beta = -0.10 + \frac{0.944}{\sqrt{E_p^{2}+E_g^{2}}} + 0.069\,\rho^{0.1}
```

```math
\gamma = 0.191\,\rho^{-0.5},\qquad C = 1.97-0.91\,U,\qquad D = 53.4-20.8\,U
```

```math
E_p = 28.8\sqrt{\frac{N_v\rho}{M}},\qquad U = \frac{N_v\rho}{M} = \frac{E_p^{2}}{829.4}
```

where $\rho$ is the density in g/cm3, $M$ the molar mass, $N_v$ the number of valence electrons per formula unit and
$E_g$ the band gap in eV. The suite divides by 10 to give nanometres. The same code reproduces the oxide IMFPs
that Biesinger quotes for Cu LMM electrons [5] to within 0.2 percent. Which presets use computed values, and the
inputs behind them, are listed in [Presets and data](presets.md#presets-and-data).
