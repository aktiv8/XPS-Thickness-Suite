"""
Tougaard inelastic-scattering background parameter B1 of a photoelectron peak (numpy only, no Tk). The background-based
DAXPES method of B. J. Murdoch, Surf. Interface Anal. 58, 125 (2026), Eqs. 1-2, uses B1 measured with a soft and a hard
X-ray source.

With the measured spectrum J(E) and the primary spectrum F(E),

    F(E) = J(E) - B1 integral_E^Emax J(E') K(E' - E) dE'                                           (1)
    K(T) = T / ((C - T^2)^2 + D T^2)                                                              (2)

(E kinetic energy, T the energy lost; Murdoch's oxide parameters are C = -1643 eV^2, D = 0, silicon C = 325 eV^2, D = 96 eV^2).
B1 (eV^2) is the single number that makes F(E) zero at a point DeltaE (15-30 eV) below the centroid of the peak in kinetic
energy, that is DeltaE ABOVE it in binding energy. It measures how much of the peak has lost about one inelastic collision, so
it grows with the thickness of a film that carries the peak and falls with the electron energy.

Everything here works on a binding-energy axis (ascending), where "lower kinetic energy" means a higher binding energy.
"""
import math
from dataclasses import dataclass, field
from typing import Dict, Optional, Sequence, Tuple

import numpy as np

_trapezoid = getattr(np, "trapezoid", None) or np.trapz

# (C, D) in eV^2 for K(T) of Eq. 2. The oxide and silicon values are Murdoch's; "universal" is Tougaard's usual C = 1643,
# D = 1 written in the form of Eq. 2 (C -> -C).
CROSS_SECTIONS: Dict[str, Tuple[float, float]] = {
    "Metal oxides (C = -1643, D = 0)": (-1643.0, 0.0),
    "Silicon (C = 325, D = 96)": (325.0, 96.0),
    "Universal (C = -1643, D = 1)": (-1643.0, 1.0),
}
DEFAULT_DELTA_E = 30.0          # eV; Murdoch used 30 eV (the allowed range is 15-30 eV)
DEFAULT_BAND = (15.0, 45.0)     # eV above the centroid for the optional band fit (a low-noise variant, not Murdoch's)
SMOOTH_POINTS = 5               # five-point adjacent averaging for the fit (Murdoch, Fig. 9)
BASELINE_POINTS = 10            # points at the low-binding-energy end taken as the flat zero level


def cross_section(t_ev, c: float, d: float):
    """K(T) of Eq. 2 (eV^-3); zero at T = 0."""
    t = np.asarray(t_ev, dtype=float)
    return t / ((c - t ** 2) ** 2 + d * t ** 2)


def kernel_matrix(be: np.ndarray, c: float, d: float) -> np.ndarray:
    """M[i, k] = K(be_i - be_k) w_k for k <= i (trapezoid weights), so that  integral = M @ y  for the electrons that lose
    energy to arrive at binding energy be_i from lower binding energy be_k."""
    x = np.asarray(be, dtype=float)
    n = len(x)
    w = np.zeros(n)
    dx = np.diff(x)
    w[:-1] += 0.5 * dx
    w[1:] += 0.5 * dx
    t = x[:, None] - x[None, :]
    m = np.where(t > 0, cross_section(np.where(t > 0, t, 1.0), c, d), 0.0) * w[None, :]
    m[np.arange(n), np.arange(n)] = 0.0                 # K(0) = 0 anyway
    return m


def smooth(y: np.ndarray, n: int = SMOOTH_POINTS) -> np.ndarray:
    """Moving average over n points (edges use the points that exist)."""
    y = np.asarray(y, dtype=float)
    if n <= 1:
        return y.copy()
    k = np.ones(n)
    return np.convolve(y, k, "same") / np.convolve(np.ones_like(y), k, "same")


@dataclass
class B1Fit:
    b1: float                              # eV^2
    centroid_be: float                     # centroid of the primary spectrum F inside the window
    eval_be: float                         # the point where F is forced to zero
    area: float                            # integral of F over the window (input units x eV)
    baseline: float
    be: np.ndarray
    j: np.ndarray                          # baseline-subtracted spectrum used
    background: np.ndarray                 # b1 * integral J K
    f: np.ndarray                          # J - background
    notes: list = field(default_factory=list)
    band_be: Optional[Tuple[float, float]] = None      # binding-energy range used when B1 was fitted over a band


def fit_b1(be: Sequence[float], y: Sequence[float], window: Tuple[float, float], c: float, d: float,
           delta_e: float = DEFAULT_DELTA_E, eval_be: Optional[float] = None, smooth_n: int = SMOOTH_POINTS,
           n_base: int = BASELINE_POINTS, matrix: Optional[np.ndarray] = None, iterations: int = 12,
           band: Optional[Tuple[float, float]] = None) -> B1Fit:
    """
    B1 of one spectrum. `window` = (lo, hi) binding energies around the peak (for the centroid and the area); the spectrum
    must run from a flat, peak-free zero level at the low-BE end to at least `delta_e` above the centroid (or `eval_be`, an
    absolute binding energy, if given). The zero level is the mean of the first `n_base` points; the fit uses the spectrum
    smoothed over `smooth_n` points. The centroid is found self-consistently from F. `matrix` can carry a precomputed
    `kernel_matrix` for the same axis. With `band` = (e1, e2) (eV above the centroid) B1 is instead the least-squares value
    over that range, B1 = sum(J I) / sum(I^2) with I the integral of Eq. 1, which is the same number for a spectrum that
    obeys Eq. 1 but uses many points rather than one, so it is far less noisy; it is a variant of Murdoch's definition.
    """
    x = np.asarray(be, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 20 or np.any(np.diff(x) <= 0):
        raise ValueError("The spectrum needs at least 20 points in ascending binding energy.")
    lo, hi = window
    if not x[0] < lo < hi < x[-1]:
        raise ValueError("The peak window must lie inside the spectrum.")
    base = float(np.mean(y[:n_base]))
    yb = y - base
    ys = smooth(yb, smooth_n)
    m = matrix if matrix is not None else kernel_matrix(x, c, d)
    integral = m @ ys
    win = (x >= lo) & (x <= hi)
    if win.sum() < 5:
        raise ValueError("The peak window holds too few points.")
    centroid = float(x[win][int(np.argmax(ys[win]))])
    b1 = 0.0
    xf = centroid
    band_be = None
    for _ in range(iterations):
        if band is not None and eval_be is None:
            lo_b, hi_b = centroid + band[0], centroid + band[1]
            if hi_b > x[-1]:
                raise ValueError(f"The spectrum ends at {x[-1]:.1f} eV but the background band reaches {hi_b:.1f} eV "
                                 f"({band[1]:g} eV above the peak): record a wider range.")
            sel = (x >= lo_b) & (x <= hi_b)
            if sel.sum() < 5:
                raise ValueError("The background band holds too few points.")
            den = float(integral[sel] @ integral[sel])
            b1 = float(ys[sel] @ integral[sel]) / den if den > 0 else 0.0
            if b1 <= 0:
                raise ValueError("There is no inelastic tail to fit in the band: the signal there is at or below the "
                                 "zero level (too noisy, or the baseline end of the spectrum is not flat).")
            xf, band_be = 0.5 * (lo_b + hi_b), (lo_b, hi_b)
        else:
            xf = eval_be if eval_be is not None else centroid + delta_e
            if xf > x[-1]:
                raise ValueError(f"The spectrum ends at {x[-1]:.1f} eV but the background must be fitted at {xf:.1f} eV "
                                 f"({delta_e:g} eV above the peak): record a wider range.")
            jf, i_f = float(np.interp(xf, x, ys)), float(np.interp(xf, x, integral))
            if i_f <= 0 or jf <= 0:
                raise ValueError("There is no inelastic tail to fit at the evaluation point: the signal there is at or "
                                 "below the zero level (too noisy, or the baseline end of the spectrum is not flat).")
            b1 = jf / i_f
        f = ys - b1 * integral
        w = np.where(win, f, 0.0)
        tot = float(_trapezoid(w, x))
        if tot <= 0:
            raise ValueError("The primary spectrum has no positive area in the window.")
        new = float(_trapezoid(w * x, x)) / tot
        if abs(new - centroid) < 1e-3:
            centroid = new
            break
        centroid = new
    f_raw = yb - b1 * (m @ yb)
    area = float(_trapezoid(np.where(win, f_raw, 0.0), x))
    return B1Fit(b1, centroid, xf, area, base, x, yb, b1 * integral, ys - b1 * integral, band_be=band_be)


def b1_sigma(be: Sequence[float], counts: Sequence[float], scale: float, window: Tuple[float, float], c: float, d: float,
             draws: int = 60, seed: int = 1, **kwargs) -> Optional[float]:
    """
    1-sigma of B1 from counting statistics by Monte Carlo: `counts` are raw counts, `scale` converts counts to the
    intensity unit used for the fit (for example 1 / (scans x dwell x transmission)); each draw replaces every count by a
    Poisson sample. Returns None when too few draws can be fitted.
    """
    x = np.asarray(be, dtype=float)
    counts = np.asarray(counts, dtype=float)
    m = kernel_matrix(x, c, d)
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(draws):
        try:
            vals.append(fit_b1(x, rng.poisson(np.maximum(counts, 0.0)) * scale, window, c, d, matrix=m, **kwargs).b1)
        except ValueError:
            continue
    return float(np.std(vals, ddof=1)) if len(vals) >= max(10, draws // 3) else None


def sensitivity_table(be_soft, y_soft, be_hard, y_hard, window: Tuple[float, float],
                      classes: Optional[Dict[str, Tuple[float, float]]] = None,
                      deltas: Sequence[float] = (15.0, 20.0, 30.0), band: Optional[Tuple[float, float]] = DEFAULT_BAND):
    """
    B1 of the soft and hard spectrum for each cross-section class and each DeltaE (and the band fit), as
    [(label, b1_soft, b1_hard)] with fits that fail left out. It shows how much of the B1 ratio, and so of the thickness,
    rests on the choices of C, D and DeltaE rather than on the data.
    """
    classes = classes or {k: v for k, v in CROSS_SECTIONS.items() if "Universal" not in k}
    rows = []
    for name, (c, d) in classes.items():
        mats = (kernel_matrix(np.asarray(be_soft, dtype=float), c, d), kernel_matrix(np.asarray(be_hard, dtype=float), c, d))
        variants = [(f"DeltaE {de:g} eV", dict(delta_e=de)) for de in deltas]
        if band is not None:
            variants.append((f"band {band[0]:g}-{band[1]:g} eV", dict(band=band)))
        for label, kw in variants:
            try:
                bs = fit_b1(be_soft, y_soft, window, c, d, matrix=mats[0], **kw).b1
                bh = fit_b1(be_hard, y_hard, window, c, d, matrix=mats[1], **kw).b1
            except ValueError:
                continue
            rows.append((f"{name.split(' (')[0]}, {label}", bs, bh))
    return rows
