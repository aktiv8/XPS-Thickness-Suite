"""
Angle-resolved XPS: peak areas from a tilt series and the overlayer thickness from intensity against 1/cos(theta).
Pure maths (numpy only, no Tk), so it is easy to test.

Angles are measured from the surface NORMAL (0 = normal emission), the Thickogram convention. For a stage-tilt series
the emission angle is the stage tilt plus an offset (default 0).

The model is the one the other tabs use. With x = d / (lambda cos(theta)) and u = 1 / cos(theta):

    two states of one core level:  R / Rinf = e^x - 1,          x = ln(1 + R / Rinf)       (Strohmeier)
    two different peaks:           A = (1 - e^-x) e^(B x),       x = thickogram.solve_x(A, B)

(the first is the B = 1 case of the second). In both cases x_i = k u_i with k = d / lambda, so the per-angle x_i fall
on a straight line through the origin when plotted against 1/cos(theta). The slope gives the thickness.

Peak areas come from one of two places:

* the spectra themselves: every spectrum is first shifted to match the reference (the lowest angle) by comparing line
  shapes, because charging and the flood gun move the peaks by a few tenths of an eV between angles; then a Shirley
  background is subtracted and the signal is integrated over fixed windows;
* the peak fits stored in a VAMAS file saved by CasaXPS (see casa.py): the substrate and overlayer states are sums of the
  stored component areas, which are used as they are. The states are told apart by how their share of the total area
  changes with angle (falling = substrate, rising = overlayer) and the assignment can be edited.
"""
import math
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

import casa
import thickogram as TG

_trapezoid = getattr(np, "trapezoid", None) or np.trapz

SHIRLEY_END_POINTS = 3          # points averaged at each end of a window to fix the background end values
SEARCH_PAD_EV = 1.0             # how far a spectrum may have moved from the reference (charging, flood gun)
MIN_STATE_SEPARATION_EV = 2.0   # closer than this, two states overlap too much for a window split
TRANSMISSION_SCALE = 1.0e4      # corrected intensities are counts/s / (T / this), so they stay near count rates
DEFAULT_MAX_ANGLE = 70.0        # beyond this the simple exponential model is unreliable (elastic scattering, roughness)


# ---------------------------------------------------------------------------------------------------------------
# Intensities, background, peak finding
# ---------------------------------------------------------------------------------------------------------------

def intensity(block, per_second: bool = True, transmission: bool = True) -> np.ndarray:
    """
    Block intensity on its own (ascending BE) axis. The K-Alpha stores raw counts summed over the scans (the noise is
    sqrt(counts)), so `per_second` divides by scans x dwell, which is needed to compare regions acquired with different
    scan counts. `transmission` divides by the Transmission column, which is needed to compare peaks at different kinetic
    energies; for two states of one core level both are constant and cancel in the ratio.
    """
    y = np.asarray(block.counts, dtype=float).copy()
    if per_second and block.scans and block.dwell_s:
        y /= block.scans * block.dwell_s
    if transmission and block.transmission is not None:
        t = np.asarray(block.transmission, dtype=float)
        if np.all(t > 0):
            y /= t / TRANSMISSION_SCALE
    return y


def smooth(y: np.ndarray, n: int = 5) -> np.ndarray:
    """Moving average with edge padding (for locating peaks and valleys only, never for areas)."""
    if n <= 1 or len(y) < n:
        return np.asarray(y, dtype=float)
    pad = n // 2
    p = np.concatenate((np.full(pad, y[0]), y, np.full(pad, y[-1])))
    return np.convolve(p, np.ones(n) / n, mode="valid")[:len(y)]


def shirley(y: np.ndarray, dx: float, n_end: int = SHIRLEY_END_POINTS, max_iter: int = 200,
            tol: float = 1e-8) -> np.ndarray:
    """
    Iterative Shirley (Proctor-Sherwood) background on an ASCENDING binding-energy grid: the background rises from the
    mean of the first n_end points to the mean of the last n_end, in proportion to the signal area at lower BE.
    """
    y = np.asarray(y, dtype=float)
    low, high = float(np.mean(y[:n_end])), float(np.mean(y[-n_end:]))
    bg = np.full_like(y, low)
    for _ in range(max_iter):
        sig = y - bg
        cum = np.concatenate(([0.0], np.cumsum(0.5 * (sig[1:] + sig[:-1]) * dx)))
        total = cum[-1]
        new = low + (high - low) * (cum / total if abs(total) > 0 else 0.0)
        done = np.max(np.abs(new - bg)) <= tol * max(1.0, abs(high - low))
        bg = new
        if done:
            break
    return bg


def shift_onto(be: np.ndarray, y: np.ndarray, shift: float, grid: np.ndarray) -> np.ndarray:
    """The spectrum moved by -shift eV (a feature found at E + shift lands at E) and sampled on `grid`."""
    return np.interp(grid + shift, be, y)


def xcorr_shift(ref_be: np.ndarray, ref_y: np.ndarray, be: np.ndarray, y: np.ndarray, lo: float, hi: float,
                pad: float = SEARCH_PAD_EV, step: float = 0.01) -> float:
    """
    How far (eV) the spectrum (be, y) has moved from the reference, found by matching the line shapes inside [lo, hi]
    (each scaled to 0-1): the shift s that minimises the squared difference between y(E + s) and the reference. Works
    whichever peak is tallest and for shoulders, which a peak-maximum search does not. Searched within +/- pad eV and
    refined to better than `step` by a parabola.
    """
    g = ref_be[(ref_be >= lo) & (ref_be <= hi)]
    if len(g) < 7:
        raise ValueError(f"The region {lo:g}-{hi:g} eV holds too few points to align the spectra.")
    r = np.interp(g, ref_be, ref_y)
    r = (r - r.min()) / max(float(np.ptp(r)), 1e-12)
    shifts = np.arange(-pad, pad + step / 2, step)
    cost = np.empty(len(shifts))
    for i, s in enumerate(shifts):
        v = np.interp(g + s, be, y)
        v = (v - v.min()) / max(float(np.ptp(v)), 1e-12)
        cost[i] = float(np.sum((v - r) ** 2))
    j = int(np.argmin(cost))
    if 0 < j < len(shifts) - 1:
        c0, c1, c2 = cost[j - 1], cost[j], cost[j + 1]
        denom = c0 - 2 * c1 + c2
        if denom > 0:
            return float(shifts[j] + 0.5 * step * (c0 - c2) / denom)
    return float(shifts[j])


# ---------------------------------------------------------------------------------------------------------------
# Windows and areas
# ---------------------------------------------------------------------------------------------------------------

@dataclass
class Window:
    """One peak: the BE window [lo, hi] (in the reference spectrum's scale) that gets a Shirley background."""
    lo: float
    hi: float


@dataclass
class TwoStateWindows:
    """Substrate and overlayer states of ONE region: shared envelope and Shirley background, split at `split`."""
    env_lo: float
    env_hi: float
    split: float
    sub_pos: Optional[float] = None       # where the substrate / overlayer state was found (information only)
    ov_pos: Optional[float] = None


@dataclass
class Extraction:
    """What was integrated in one spectrum, kept for plotting."""
    grid: np.ndarray
    y: np.ndarray               # aligned onto the reference scale
    bg: np.ndarray              # Shirley background (nan outside the envelope)
    lo: float
    hi: float
    split: Optional[float]
    shift_ev: float
    area_sub: float
    area_ov: Optional[float]
    peak_sub: Optional[float] = None      # BE of the signal maximum in each part (tells resolved states from merged ones)
    peak_ov: Optional[float] = None


def _prepare(be, y, reference, lo, hi, grid, align):
    shift = 0.0
    if align and reference is not None:
        shift = xcorr_shift(reference[0], reference[1], be, y, lo, hi)
    g = np.asarray(grid if grid is not None else be, dtype=float)
    return shift, g, shift_onto(be, y, shift, g)


def single_peak_area(be: np.ndarray, y: np.ndarray, win: Window, reference=None, grid: Optional[np.ndarray] = None,
                     align: bool = True) -> Extraction:
    """
    Shirley-corrected area of one peak. With a `reference` (be, y) the spectrum is first shifted to match it, so the
    window follows the peak from angle to angle.
    """
    shift, g, ya = _prepare(be, y, reference, win.lo, win.hi, grid, align)
    m = (g >= win.lo) & (g <= win.hi)
    if m.sum() < 2 * SHIRLEY_END_POINTS + 1:
        raise ValueError(f"The window {win.lo:g}-{win.hi:g} eV holds too few points.")
    x, yy = g[m], ya[m]
    bg = shirley(yy, float(np.mean(np.diff(x))))
    full = np.full_like(g, np.nan)
    full[m] = bg
    return Extraction(g, ya, full, win.lo, win.hi, None, shift, float(_trapezoid(yy - bg, x)), None)


def two_state_areas(be: np.ndarray, y: np.ndarray, w: TwoStateWindows, reference=None,
                    grid: Optional[np.ndarray] = None, align: bool = True) -> Extraction:
    """One Shirley background across the whole envelope, then the signal is split at `split` into the two states."""
    shift, g, ya = _prepare(be, y, reference, w.env_lo, w.env_hi, grid, align)
    m = (g >= w.env_lo) & (g <= w.env_hi)
    if m.sum() < 2 * SHIRLEY_END_POINTS + 1:
        raise ValueError(f"The envelope {w.env_lo:g}-{w.env_hi:g} eV holds too few points.")
    if not (w.env_lo < w.split < w.env_hi):
        raise ValueError("The split energy must lie inside the envelope.")
    x, yy = g[m], ya[m]
    bg = shirley(yy, float(np.mean(np.diff(x))))
    sig = yy - bg
    ys = float(np.interp(w.split, x, sig))                 # put the split exactly on the grid
    left = np.concatenate((x[x < w.split], [w.split]))
    right = np.concatenate(([w.split], x[x > w.split]))
    sl, sr = np.concatenate((sig[x < w.split], [ys])), np.concatenate(([ys], sig[x > w.split]))
    full = np.full_like(g, np.nan)
    full[m] = bg
    return Extraction(g, ya, full, w.env_lo, w.env_hi, w.split, shift, float(_trapezoid(sl, left)),
                      float(_trapezoid(sr, right)), float(left[int(np.argmax(smooth(sl)))]),
                      float(right[int(np.argmax(smooth(sr)))]))


def _aligned_unit_area(spectra, lo, hi):
    """
    All spectra shifted onto the first one's scale, baseline-subtracted and scaled to unit area inside [lo, hi]. Each
    spectrum is matched to the previous (already aligned) one, so the composition change between neighbours is small.
    """
    be0, y0 = spectra[0]
    m = (be0 >= lo) & (be0 <= hi)
    prev = (be0, y0)
    out = []
    for be, y in spectra:
        ya = shift_onto(be, y, xcorr_shift(prev[0], prev[1], be, y, lo, hi), be0)
        prev = (be0, ya)
        v = ya - np.min(ya[m])
        out.append(v / max(float(_trapezoid(v[m], be0[m])), 1e-12))
    return be0, out


def auto_two_state_windows(spectra: Sequence[Tuple[np.ndarray, np.ndarray]], sub_nominal: float, ov_nominal: float,
                           below_ev: float = 2.5, above_ev: float = 3.0, half_range: float = 2.5) -> TwoStateWindows:
    """
    Windows for a two-state region from the data. `spectra` must run from the lowest angle to the highest. The spectra
    are aligned, and the substrate and overlayer states are told apart by how they change with angle: the substrate
    state is where the normalised spectrum FALLS from the lowest to the highest angle, the overlayer state where it
    RISES. This works when the overlayer peak is the taller one (a thick oxide), unlike 'take the tallest peak'.
    Each is searched within +/- half_range eV of its nominal BE (charging moves the peaks by a few tenths of an eV or
    more). The envelope runs below_ev below the substrate state to above_ev above the overlayer state (clipped to the
    data). The split is the minimum of the summed spectra between the states, or half-way if there is no valley.
    """
    if len(spectra) < 2:
        raise ValueError("Auto-detection needs at least two angles.")
    be0 = spectra[0][0]
    lo = max(float(be0[0]), min(sub_nominal, ov_nominal) - half_range - below_ev)
    hi = min(float(be0[-1]), max(sub_nominal, ov_nominal) + half_range + above_ev)
    _, norm = _aligned_unit_area(list(spectra), lo, hi)
    diff = smooth(norm[-1] - norm[0], 7)
    ms = (be0 >= sub_nominal - half_range) & (be0 <= sub_nominal + half_range)
    mo = (be0 >= ov_nominal - half_range) & (be0 <= ov_nominal + half_range)
    if ms.sum() < 3 or mo.sum() < 3:
        raise ValueError("The data do not cover the expected binding energies of the two states.")
    sub_pos = float(be0[ms][int(np.argmin(diff[ms]))])
    ov_pos = float(be0[mo][int(np.argmax(diff[mo]))])
    if ov_pos <= sub_pos:
        raise ValueError("Could not tell the two states apart (the overlayer state did not grow with angle). "
                         "Set the windows by hand.")
    total = np.sum(norm, axis=0)
    gap = ov_pos - sub_pos
    m = (be0 >= sub_pos + 0.3 * gap) & (be0 <= ov_pos - 0.25 * gap)
    sm = smooth(total, 5)[m]
    j = int(np.argmin(sm)) if len(sm) else 0
    split = float(be0[m][j]) if len(sm) > 2 and 0 < j < len(sm) - 1 else 0.5 * (sub_pos + ov_pos)
    return TwoStateWindows(max(float(be0[0]), sub_pos - below_ev), min(float(be0[-1]), ov_pos + above_ev), split,
                           sub_pos, ov_pos)


def auto_peak_window(spectra: Sequence[Tuple[np.ndarray, np.ndarray]], nominal: float, below_ev: float = 4.0,
                     above_ev: float = 4.0, half_range: float = 2.5) -> Window:
    """A window around the tallest peak, in the tallest spectrum, found near `nominal` BE."""
    best, pos = -math.inf, None
    for be, y in spectra:
        m = (be >= nominal - half_range) & (be <= nominal + half_range)
        if m.sum() < 3:
            continue
        sm = smooth(y)[m]
        if float(sm.max()) > best:
            best, pos = float(sm.max()), float(be[m][int(np.argmax(sm))])
    if pos is None:
        raise ValueError("The data do not cover that binding energy.")
    be0 = spectra[0][0]
    return Window(max(float(be0[0]), pos - below_ev), min(float(be0[-1]), pos + above_ev))


# ---------------------------------------------------------------------------------------------------------------
# Fits
# ---------------------------------------------------------------------------------------------------------------

def fit_through_origin(u: Sequence[float], x: Sequence[float]) -> Tuple[float, Optional[float]]:
    """Least-squares slope of x = k u and its standard error (None with a single point)."""
    u, x = np.asarray(u, dtype=float), np.asarray(x, dtype=float)
    suu = float(np.sum(u * u))
    k = float(np.sum(u * x)) / suu
    if len(u) < 2:
        return k, None
    s = math.sqrt(float(np.sum((x - k * u) ** 2)) / (len(u) - 1))
    return k, s / math.sqrt(suu)


def fit_line(u: Sequence[float], x: Sequence[float]) -> Optional[Tuple[float, float, float, float]]:
    """Ordinary least squares x = c + k u: (c, k, se_c, se_k). None with fewer than 3 points or no spread in u."""
    u, x = np.asarray(u, dtype=float), np.asarray(x, dtype=float)
    n = len(u)
    sxx = float(np.sum((u - u.mean()) ** 2))
    if n < 3 or sxx <= 0:
        return None
    k = float(np.sum((u - u.mean()) * (x - x.mean()))) / sxx
    c = float(x.mean() - k * u.mean())
    s2 = float(np.sum((x - c - k * u) ** 2)) / (n - 2)
    return c, k, math.sqrt(s2 * (1.0 / n + u.mean() ** 2 / sxx)), math.sqrt(s2 / sxx)


def fit_rinf_free(u: Sequence[float], r: Sequence[float]) -> Optional[dict]:
    """
    Fit R_i = Rinf (e^(k u_i) - 1) for BOTH k and Rinf (relative errors, so every angle counts alike). For a trial k the
    best Rinf is linear; k is found by a log-spaced scan refined by golden-section search. Returns None unless there are
    at least 3 positive ratios; `edge` is True when the data cannot pin the pair down (best k on the scan limit, or nearly linear in u).
    """
    u, r = np.asarray(u, dtype=float), np.asarray(r, dtype=float)
    if len(u) < 3 or np.any(r <= 0):
        return None

    def solve(k):
        g = np.expm1(k * u) / r
        rinf = float(np.sum(g)) / float(np.sum(g * g))
        return float(np.sum((1.0 - rinf * g) ** 2)), rinf

    ks = np.logspace(-3, 1.3, 500)
    costs = [solve(k)[0] for k in ks]
    j = int(np.argmin(costs))
    lo, hi = ks[max(j - 1, 0)], ks[min(j + 1, len(ks) - 1)]
    phi = (math.sqrt(5) - 1) / 2
    a, b = lo, hi
    for _ in range(80):
        c1, c2 = b - phi * (b - a), a + phi * (b - a)
        if solve(c1)[0] < solve(c2)[0]:
            b = c2
        else:
            a = c1
    k = 0.5 * (a + b)
    cost, rinf = solve(k)
    # Not determined when the best k is on the scan limit or the exponential is nearly linear over the angles used
    # (k u_max below 0.5): R_inf and k are then interchangeable and any pair will do.
    weak = k * float(np.max(u)) < 0.5
    return {"k": k, "rinf": rinf, "cost": cost, "edge": j in (0, len(ks) - 1) or weak}


# ---------------------------------------------------------------------------------------------------------------
# Series analysis
# ---------------------------------------------------------------------------------------------------------------

@dataclass
class FitView:
    """What the Spectra page needs to show one angle of a fitted region (peak shapes are not redrawn)."""
    block: object
    shift_ev: float                   # drift of the substrate state relative to the first angle
    norm: float                       # raw stored area used to normalise the drawn spectrum
    marks: List[Tuple[str, float, str]] = field(default_factory=list)      # (state name, BE, role)


@dataclass
class AngleRow:
    angle_deg: float                  # emission angle from the surface normal
    tilt: float                       # as stored in the file
    i_sub: float
    i_ov: float
    ratio: Optional[float]            # R = I_ov / I_sub (mode 1) or A (mode 2)
    x: Optional[float] = None         # d / (lambda cos theta) from this angle alone
    d_nm: Optional[float] = None
    used: bool = False
    note: str = ""
    shift_ev: float = 0.0
    sub: Optional[Extraction] = None
    ov: Optional[Extraction] = None
    fit_sub: Optional[FitView] = None          # set when the areas came from stored CasaXPS fits
    fit_ov: Optional[FitView] = None

    @property
    def u(self) -> float:
        return 1.0 / math.cos(math.radians(self.angle_deg))


@dataclass
class ARXPSResult:
    mode: int
    lam_nm: float
    rows: List[AngleRow]
    rinf: Optional[float] = None
    b: Optional[float] = None
    k: Optional[float] = None                 # d / lambda, slope through the origin
    k_se: Optional[float] = None
    d_nm: Optional[float] = None
    d_se_nm: Optional[float] = None
    line: Optional[Tuple[float, float, float, float]] = None      # free-intercept fit (c, k, se_c, se_k)
    rinf_free: Optional[dict] = None
    azimuth: Optional[float] = None
    warnings: List[str] = field(default_factory=list)

    @property
    def used(self) -> List[AngleRow]:
        return [r for r in self.rows if r.used]

    @property
    def d_mean_nm(self) -> Optional[float]:
        d = [r.d_nm for r in self.used if r.d_nm is not None]
        return float(np.mean(d)) if d else None

    @property
    def d_spread_nm(self) -> Optional[float]:
        d = [r.d_nm for r in self.used if r.d_nm is not None]
        return float(np.std(d, ddof=1)) if len(d) > 1 else None

    @property
    def d_free_nm(self) -> Optional[float]:
        return self.line[1] * self.lam_nm if self.line else None

    @property
    def d_rinf_free_nm(self) -> Optional[float]:
        return self.rinf_free["k"] * self.lam_nm if self.rinf_free else None


def fit_series(rows: List[AngleRow], lam_nm: float, solve_x: Callable[[float], float], max_angle: float,
               mode: int, rinf: Optional[float] = None, b: Optional[float] = None) -> ARXPSResult:
    """
    From rows with areas and ratios to the thickness: per-angle x, the through-origin fit of x against 1/cos(theta),
    the free-intercept check and (mode 1) the Rinf-free fit. `solve_x` turns one ratio into x or raises ValueError.
    """
    if not math.isfinite(lam_nm) or lam_nm <= 0:
        raise ValueError("Lambda (attenuation length) must be a positive number of nm.")
    res = ARXPSResult(mode=mode, lam_nm=lam_nm, rows=rows, rinf=rinf, b=b)
    for r in rows:
        if not (0.0 <= r.angle_deg < 90.0):
            r.note = "angle outside 0-90 degrees"
            continue
        if r.angle_deg > max_angle:
            r.note = f"above the {max_angle:g} deg limit"
            if r.ratio is not None and r.ratio > 0:
                try:                                           # still shown, never fitted
                    r.x = solve_x(r.ratio)
                    r.d_nm = r.x * lam_nm / r.u
                except ValueError:
                    pass
            continue
        if r.ratio is None or not (r.ratio > 0) or not math.isfinite(r.ratio):
            r.note = "area is zero or negative"
            continue
        try:
            r.x = solve_x(r.ratio)
        except ValueError as exc:
            r.note = str(exc)
            continue
        r.d_nm = r.x * lam_nm / r.u
        r.used = True

    used = res.used
    if not used:
        raise ValueError("No usable angle: check the windows (an area is zero or negative) and the maximum angle.")
    u = [r.u for r in used]
    xs = [r.x for r in used]
    res.k, res.k_se = fit_through_origin(u, xs)
    res.d_nm = res.k * lam_nm
    res.d_se_nm = res.k_se * lam_nm if res.k_se is not None else None
    res.line = fit_line(u, xs)
    if mode == 1:
        res.rinf_free = fit_rinf_free(u, [r.ratio for r in used])

    skipped = [r for r in rows if not r.used and r.note and "limit" not in r.note]
    for r in skipped:
        res.warnings.append(f"{r.tilt:g} deg skipped: {r.note}.")
    over = [r for r in rows if "limit" in r.note]
    if over:
        res.warnings.append("Excluded from the fit (above the maximum angle): "
                            + ", ".join(f"{r.tilt:g} deg" for r in over) + ".")
    if len(used) < 3:
        res.warnings.append(f"Only {len(used)} usable angle(s): the line, its error and the intercept check need 3 or more.")
    elif max(u) - min(u) < 0.3:
        res.warnings.append("The angles span a very small range of 1/cos(theta); the fit is poorly determined.")
    if res.line is not None:
        c, _, se_c, _ = res.line
        if abs(c) > max(3.0 * se_c, 0.05):
            res.warnings.append(f"The free-intercept fit does not pass through the origin (intercept {c:+.3f} in x): "
                                "the windows, Rinf or the angle offset may be wrong, or the layer is not uniform.")
    if res.d_nm > 3.0 * lam_nm:
        res.warnings.append(f"The thickness exceeds 3 lambda ({3.0 * lam_nm:.2f} nm): the substrate signal is small and the "
                            "result unreliable.")
    shifts = [abs(r.shift_ev) for r in rows if r.shift_ev]
    if shifts and max(shifts) > 0.5:
        res.warnings.append(f"Peaks moved by up to {max(shifts):.2f} eV between angles (charging); windows followed them.")
    return res


def _emission(tilt: float, offset: float) -> float:
    return tilt + offset


def analyse_two_states(series: Sequence[Tuple[float, object]], windows: TwoStateWindows, lam_nm: float, rinf: float,
                       angle_offset: float = 0.0, max_angle: float = DEFAULT_MAX_ANGLE, per_second: bool = True,
                       transmission: bool = True, align: bool = True, reference=None) -> ARXPSResult:
    """
    Mode 1: substrate and overlayer states in one core-level region (series = [(tilt, Block), ...], lowest angle first).
    rinf = (RSF N lambda)_overlayer / (RSF N lambda)_substrate, built exactly as the Calculator tab does.
    `reference` is the (be, intensity) the windows were set on (default: the first spectrum of the series).
    """
    if not (math.isfinite(rinf) and rinf > 0):
        raise ValueError("Rinf must be a positive number.")
    ref_be = series[0][1].be_ev
    if reference is None:
        reference = (ref_be, intensity(series[0][1], per_second, transmission))
    rows = []
    for tilt, blk in series:
        y = intensity(blk, per_second, transmission)
        ex = two_state_areas(blk.be_ev, y, windows, reference=reference, grid=ref_be, align=align)
        reference = (ex.grid, ex.y)                            # the next angle is matched to this one
        ratio = ex.area_ov / ex.area_sub if ex.area_sub > 0 else None
        rows.append(AngleRow(_emission(tilt, angle_offset), tilt, ex.area_sub, ex.area_ov, ratio,
                             shift_ev=ex.shift_ev, sub=ex, ov=ex))
    res = fit_series(rows, lam_nm, lambda r: math.log1p(r / rinf), max_angle, mode=1, rinf=rinf)
    probe = next((r.sub for r in rows if r.used), rows[0].sub)
    if probe.peak_sub is not None and probe.peak_ov is not None:
        if probe.peak_ov - probe.peak_sub < MIN_STATE_SEPARATION_EV or probe.peak_ov - windows.split < 0.15:
            res.warnings.append(f"The two states are less than {MIN_STATE_SEPARATION_EV:g} eV apart (or merge at the "
                                "split): a window split cannot separate overlapping peaks, so treat the thickness as a "
                                "rough estimate (peak fitting would be needed).")
    return res


def analyse_two_peaks(sub_series: Sequence[Tuple[float, object]], ov_series: Sequence[Tuple[float, object]],
                      sub_win: Window, ov_win: Window, so: float, ss: float, be_o: float, be_s: float,
                      source_ev: float, lam_nm: float, angle_offset: float = 0.0,
                      max_angle: float = DEFAULT_MAX_ANGLE, per_second: bool = True, transmission: bool = True,
                      align: bool = True, sub_reference=None, ov_reference=None) -> ARXPSResult:
    """
    Mode 2: overlayer and substrate are different peaks (different regions of the file, matched by tilt). RSFs are
    So (overlayer) and Ss (substrate); be_o/be_s set B = KE ratio. References as for analyse_two_states.
    """
    b = TG.energy_ratio_b(be_o, be_s, source_ev)
    ov_by_tilt = {round(t, 3): blk for t, blk in ov_series}
    if sub_reference is None:
        sub_reference = (sub_series[0][1].be_ev, intensity(sub_series[0][1], per_second, transmission))
    if ov_reference is None and ov_series:
        ov_reference = (ov_series[0][1].be_ev, intensity(ov_series[0][1], per_second, transmission))
    rows = []
    for tilt, blk_s in sub_series:
        blk_o = ov_by_tilt.get(round(tilt, 3))
        if blk_o is None:
            continue
        es = single_peak_area(blk_s.be_ev, intensity(blk_s, per_second, transmission), sub_win, sub_reference,
                              align=align)
        eo = single_peak_area(blk_o.be_ev, intensity(blk_o, per_second, transmission), ov_win, ov_reference,
                              align=align)
        sub_reference, ov_reference = (es.grid, es.y), (eo.grid, eo.y)       # match the next angle to these
        ratio = None
        if es.area_sub > 0 and eo.area_sub > 0:
            ratio = TG.ratio_a(eo.area_sub, es.area_sub, so, ss)
        rows.append(AngleRow(_emission(tilt, angle_offset), tilt, es.area_sub, eo.area_sub, ratio,
                             shift_ev=max(abs(es.shift_ev), abs(eo.shift_ev)), sub=es, ov=eo))
    if not rows:
        raise ValueError("The two regions share no tilt angle.")
    return fit_series(rows, lam_nm, lambda a: TG.solve_x(a, b), max_angle, mode=2, b=b)


# ---------------------------------------------------------------------------------------------------------------
# Areas from stored CasaXPS fits
# ---------------------------------------------------------------------------------------------------------------

SUBSTRATE, OVERLAYER, IGNORE = "substrate", "overlayer", "ignore"
MIN_SHARE_CHANGE = 0.02         # a state whose share of the fitted area moves less than this is not assigned


@dataclass
class FitAngle:
    """The chemical states of one fitted block (one angle of one region)."""
    tilt: float
    block: object
    states: Dict[str, casa.StateArea]


def fit_regions(series) -> List[str]:
    """Names of the CASA fit regions held in the first block of a series (e.g. 'Ga 3d', 'As 3d')."""
    for _, blk in series:
        fit = casa.fits_of(blk)
        if fit is not None:
            return [r.name for r in fit.regions]
    return []


def read_fits(series, region: Optional[str] = None) -> List[FitAngle]:
    """The fitted states of every block of a series, lowest angle first. Raises ValueError for a block without a fit."""
    out = []
    for tilt, blk in series:
        fit = casa.fits_of(blk)
        if fit is None:
            raise ValueError(f"The {tilt:g} degree block of '{blk.block_id}' has no CasaXPS fit: save the file from "
                             "CasaXPS with the fit in place, or use the windows method.")
        out.append(FitAngle(tilt, blk, casa.states_of(fit, blk.hv_ev, region)))
    return out


def state_names(fits: Sequence[FitAngle]) -> Dict[str, str]:
    """{state key: name} over all angles, in order of first appearance."""
    names: Dict[str, str] = {}
    for f in fits:
        for key, st in f.states.items():
            names.setdefault(key, st.name)
    return names


def state_shares(fits: Sequence[FitAngle]) -> Dict[str, Tuple[float, float]]:
    """{state key: (share of the total fitted area at the lowest angle, at the highest angle)}; missing counts as 0."""
    def shares(f):
        total = sum(st.area for st in f.states.values())
        return {k: (st.area / total if total > 0 else 0.0) for k, st in f.states.items()}
    lo, hi = shares(fits[0]), shares(fits[-1])
    return {k: (lo.get(k, 0.0), hi.get(k, 0.0)) for k in state_names(fits)}


def state_trends(fits: Sequence[FitAngle]) -> Dict[str, float]:
    """
    How each state's share of the total fitted area changes from the lowest to the highest angle (positive = grows
    with angle, as an overlayer does).
    """
    return {k: hi - lo for k, (lo, hi) in state_shares(fits).items()}


def auto_roles(fits: Sequence[FitAngle], min_change: float = MIN_SHARE_CHANGE) -> Dict[str, str]:
    """Substrate where a state's share falls with angle, overlayer where it rises, ignore where it hardly moves."""
    return {k: (SUBSTRATE if d < -min_change else OVERLAYER if d > min_change else IGNORE)
            for k, d in state_trends(fits).items()}


def auto_use(fits: Sequence[FitAngle], role: str, min_change: float = MIN_SHARE_CHANGE) -> Dict[str, str]:
    """
    For the two-peak mode, which states of a region to use for `role` (substrate or overlayer): the only state of a
    region is always used; with several, those whose share moves the right way with angle.
    """
    names = state_names(fits)
    if len(names) == 1:
        return {k: role for k in names}
    roles = auto_roles(fits, min_change)
    return {k: (role if roles[k] == role else IGNORE) for k in names}


def _transmission_factor(block, ke_raw: float) -> float:
    """Transmission function of the block at a kinetic energy, in units of TRANSMISSION_SCALE (1 if there is none)."""
    t = getattr(block, "transmission", None)
    if t is None or not np.all(np.asarray(t) > 0):
        return 1.0
    ke = block.hv_ev - np.asarray(block.be_ev, dtype=float)               # descending with ascending BE
    return float(np.interp(ke_raw, ke[::-1], np.asarray(t, dtype=float)[::-1])) / TRANSMISSION_SCALE


def _role_sum(f: FitAngle, roles: Dict[str, str], role: str, transmission: bool):
    """(area, raw area, area-weighted BE, [(name, BE, role)]) of the states of `f` assigned to `role`."""
    area = raw = wbe = 0.0
    marks = []
    for key, st in f.states.items():
        if roles.get(key) != role:
            continue
        for c in st.comps:
            corr = c.area / _transmission_factor(f.block, c.ke_raw) if transmission else c.area
            area += corr
            raw += c.area
            wbe += c.area * c.be
        marks.append((st.name, st.be, role))
    return area, raw, (wbe / raw if raw > 0 else None), marks


def _check_roles(fits, roles, want: Sequence[str]) -> None:
    present = {roles.get(k) for k in state_names(fits)}
    missing = [w for w in want if w not in present]
    if missing:
        raise ValueError("Assign at least one " + " and one ".join(missing) + " state (Components page).")


def analyse_fitted_states(series, roles: Dict[str, str], lam_nm: float, rinf: float, angle_offset: float = 0.0,
                          max_angle: float = DEFAULT_MAX_ANGLE, transmission: bool = True,
                          region: Optional[str] = None, fits: Optional[List[FitAngle]] = None) -> ARXPSResult:
    """
    Mode 1 from stored fits: the substrate and overlayer states are sums of the assigned states of one fitted region
    (series = [(tilt, Block), ...], lowest angle first); everything after the areas is shared with the window method.
    """
    if not (math.isfinite(rinf) and rinf > 0):
        raise ValueError("Rinf must be a positive number.")
    fits = fits or read_fits(series, region)
    _check_roles(fits, roles, (SUBSTRATE, OVERLAYER))
    rows, ref_be = [], None
    for f in fits:
        s_area, s_raw, s_be, s_marks = _role_sum(f, roles, SUBSTRATE, transmission)
        o_area, _, _, o_marks = _role_sum(f, roles, OVERLAYER, transmission)
        if s_be is not None and ref_be is None:
            ref_be = s_be
        shift = (s_be - ref_be) if (s_be is not None and ref_be is not None) else 0.0
        view = FitView(f.block, shift, s_raw, [(n, b - shift, r) for n, b, r in s_marks + o_marks])
        ratio = o_area / s_area if s_area > 0 else None
        rows.append(AngleRow(_emission(f.tilt, angle_offset), f.tilt, s_area, o_area, ratio, shift_ev=shift,
                             fit_sub=view, fit_ov=view))
    return fit_series(rows, lam_nm, lambda r: math.log1p(r / rinf), max_angle, mode=1, rinf=rinf)


def analyse_fitted_peaks(sub_series, ov_series, sub_roles: Dict[str, str], ov_roles: Dict[str, str], so: float,
                         ss: float, be_o: float, be_s: float, source_ev: float, lam_nm: float,
                         angle_offset: float = 0.0, max_angle: float = DEFAULT_MAX_ANGLE, transmission: bool = True,
                         sub_region: Optional[str] = None, ov_region: Optional[str] = None) -> ARXPSResult:
    """
    Mode 2 from stored fits: the substrate states of one fitted region against the overlayer states of another,
    matched by tilt, through the Thickogram equation. With `transmission` each stored area is divided by the block's
    transmission function at the component's kinetic energy, which matters here because the two peaks differ in KE.
    """
    b = TG.energy_ratio_b(be_o, be_s, source_ev)
    fs, fo = read_fits(sub_series, sub_region), read_fits(ov_series, ov_region)
    _check_roles(fs, sub_roles, (SUBSTRATE,))
    _check_roles(fo, ov_roles, (OVERLAYER,))
    ov_by_tilt = {round(f.tilt, 3): f for f in fo}
    rows, ref_s, ref_o = [], None, None
    for f in fs:
        g = ov_by_tilt.get(round(f.tilt, 3))
        if g is None:
            continue
        s_area, s_raw, s_be, s_marks = _role_sum(f, sub_roles, SUBSTRATE, transmission)
        o_area, o_raw, o_be, o_marks = _role_sum(g, ov_roles, OVERLAYER, transmission)
        ref_s = s_be if (ref_s is None and s_be is not None) else ref_s
        ref_o = o_be if (ref_o is None and o_be is not None) else ref_o
        sh_s = (s_be - ref_s) if (s_be is not None and ref_s is not None) else 0.0
        sh_o = (o_be - ref_o) if (o_be is not None and ref_o is not None) else 0.0
        ratio = TG.ratio_a(o_area, s_area, so, ss) if (s_area > 0 and o_area > 0) else None
        rows.append(AngleRow(_emission(f.tilt, angle_offset), f.tilt, s_area, o_area, ratio,
                             shift_ev=max(abs(sh_s), abs(sh_o)),
                             fit_sub=FitView(f.block, sh_s, s_raw, [(n, e - sh_s, r) for n, e, r in s_marks]),
                             fit_ov=FitView(g.block, sh_o, o_raw, [(n, e - sh_o, r) for n, e, r in o_marks])))
    if not rows:
        raise ValueError("The two regions share no tilt angle.")
    return fit_series(rows, lam_nm, lambda a: TG.solve_x(a, b), max_angle, mode=2, b=b)


def summarise_azimuths(results: Dict[Optional[float], ARXPSResult]) -> Tuple[float, Optional[float]]:
    """Mean thickness over azimuths and the standard deviation between azimuths (None for a single azimuth)."""
    d = [r.d_nm for r in results.values() if r.d_nm is not None]
    if not d:
        raise ValueError("No azimuth gave a thickness.")
    return float(np.mean(d)), (float(np.std(d, ddof=1)) if len(d) > 1 else None)
