"""
Layers (the Gries signal-ratio idea): thicknesses of a layer stack from the angle dependence of each peak on its own
(numpy only, no Tk). Part of the Depth profile tab; independent implementation, see depthprofile.py.

Layers 1..K from the top (thicknesses d_k) on a semi-infinite substrate. A species in layer j gives, for unit
concentration and including the attenuation by every layer above it,

    I(theta) = exp(-sum_{m<j} d_m / (lambda_m cos theta)) (1 - exp(-d_j / (lambda_j cos theta)))

and a substrate species exp(-sum_{m<=K} d_m / (lambda_m cos theta)). Every measured curve and every prediction is divided
by its own value at the first (smallest) angle, so no sensitivity factors and no common geometry factor are needed: only
the SHAPE of each peak's intensity against angle is used (W. H. Gries, J. Vac. Sci. Technol. A 7, 1655 (1989)).

Gries adjusts the thicknesses by eye until the predicted curves match; here that is also done automatically by least
squares, and the range of each thickness that fits within the noise is scanned, which shows how information fades with
depth.
"""
from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

import numpy as np

from depthprofile import DEFAULT_MAX_ANGLE, Species, stratification

SUBSTRATE_LAYER = -1


def layer_intensity(theta_deg, d: Sequence[float], layer: int, lam: Sequence[float]) -> np.ndarray:
    """
    Unnormalised intensity against angle of a species in `layer` (0-based from the top; SUBSTRATE_LAYER for the
    substrate) for layer thicknesses `d` (nm) and the species' attenuation length `lam[m]` (nm) in each layer m.
    """
    c = np.cos(np.radians(np.asarray(theta_deg, dtype=float)))
    above = layer if layer != SUBSTRATE_LAYER else len(d)
    atten = np.exp(-sum(d[m] / lam[m] for m in range(above)) / c)
    if layer == SUBSTRATE_LAYER:
        return atten
    return atten * (1.0 - np.exp(-d[layer] / (lam[layer] * c)))


def model_curves(species: Sequence[Species], assignment: Sequence[int], d: Sequence[float],
                 lam_layers: Optional[Sequence[Sequence[float]]] = None) -> List[np.ndarray]:
    """Unnormalised model intensities at each species' own angles. `lam_layers[i][m]` defaults to species i's lambda."""
    out = []
    for i, (sp, layer) in enumerate(zip(species, assignment)):
        lam = lam_layers[i] if lam_layers is not None else [sp.lam_nm] * len(d)
        out.append(layer_intensity(sp.angles, d, layer, lam))
    return out


def measured_normalised(sp: Species, max_angle: float = DEFAULT_MAX_ANGLE, floor: float = 0.03):
    """
    (angles, intensity / its value at the first usable angle, sigma, indices) over the positive intensities up to
    `max_angle`, or None with fewer than two. Sigma comes from the counting statistics when known (at least 0.1%), else
    `floor` (3%) of the value.
    """
    th = np.asarray(sp.angles, dtype=float)
    ok = np.flatnonzero((th <= max_angle) & (sp.intensity > 0))
    if len(ok) < 2:
        return None
    i0 = ok[0]
    y = sp.intensity[ok] / sp.intensity[i0]
    if sp.sigma is not None:
        rel = np.sqrt((sp.sigma[ok] / sp.intensity[ok]) ** 2 + (sp.sigma[i0] / sp.intensity[i0]) ** 2)
        sig = y * np.maximum(rel, 1e-3)
    else:
        sig = y * floor
    return th[ok], y, np.maximum(sig, 1e-6), ok


def gries_chi2(d: Sequence[float], species: Sequence[Species], assignment: Sequence[int],
               max_angle: float = DEFAULT_MAX_ANGLE, lam_layers=None, floor: float = 0.03) -> Tuple[float, int]:
    """Weighted sum of squared differences between measured and predicted normalised curves, and the points used."""
    chi, n = 0.0, 0
    for i, (sp, layer) in enumerate(zip(species, assignment)):
        meas = measured_normalised(sp, max_angle, floor)
        if meas is None:
            continue
        th, y, sig, _ = meas
        lam = lam_layers[i] if lam_layers is not None else [sp.lam_nm] * len(d)
        pred = layer_intensity(th, d, layer, lam)
        if pred[0] <= 1e-12:
            chi += 1e6
            continue
        pred = pred / pred[0]
        chi += float(np.sum(((y[1:] - pred[1:]) / sig[1:]) ** 2))        # the first point is the normalisation
        n += len(th) - 1
    return chi, n


def nelder_mead(f, x0: Sequence[float], step: Sequence[float], lo: Sequence[float], hi: Sequence[float],
                max_iter: int = 600, tol: float = 1e-9, restarts: int = 4) -> Tuple[np.ndarray, float]:
    """
    A small bounded Nelder-Mead simplex minimiser (points are clipped into [lo, hi]); no scipy needed. Clipping can
    flatten the simplex against a bound, so the search is restarted from its best point with a fresh simplex until a
    restart no longer improves it.
    """
    x, v = _nelder_mead_once(f, x0, step, lo, hi, max_iter, tol)
    for _ in range(restarts):
        x2, v2 = _nelder_mead_once(f, x, step, lo, hi, max_iter, tol)
        if v2 >= v - 1e-12 * (abs(v) + 1e-12):
            break
        x, v = x2, v2
    return x, v


def _nelder_mead_once(f, x0, step, lo, hi, max_iter, tol):
    n = len(x0)
    lo, hi = np.asarray(lo, dtype=float), np.asarray(hi, dtype=float)

    def clip(p):
        return np.minimum(np.maximum(p, lo), hi)

    pts = [clip(np.asarray(x0, dtype=float))]
    for k in range(n):
        p = pts[0].copy()
        p[k] = p[k] + step[k] if p[k] + step[k] <= hi[k] else p[k] - step[k]
        pts.append(clip(p))
    vals = [f(p) for p in pts]
    for _ in range(max_iter):
        order = np.argsort(vals)
        pts, vals = [pts[i] for i in order], [vals[i] for i in order]
        spread = float(np.max(np.abs(pts[-1] - pts[0])))
        if spread < 1e-8 or abs(vals[-1] - vals[0]) <= 1e-15 + tol * abs(vals[0]) and spread < 1e-6:
            break
        centre = np.mean(pts[:-1], axis=0)
        xr = clip(centre + (centre - pts[-1]))
        fr = f(xr)
        if fr < vals[0]:
            xe = clip(centre + 2.0 * (centre - pts[-1]))
            fe = f(xe)
            pts[-1], vals[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < vals[-2]:
            pts[-1], vals[-1] = xr, fr
        else:
            xc = clip(centre + 0.5 * ((xr if fr < vals[-1] else pts[-1]) - centre))
            fc = f(xc)
            if fc < min(fr, vals[-1]):
                pts[-1], vals[-1] = xc, fc
            else:
                for j in range(1, n + 1):
                    pts[j] = clip(pts[0] + 0.5 * (pts[j] - pts[0]))
                    vals[j] = f(pts[j])
    i = int(np.argmin(vals))
    return pts[i], float(vals[i])


@dataclass
class LayerFit:
    d_nm: np.ndarray                    # fitted thicknesses, top layer first
    chi2: float
    dof: int
    assignment: List[int]
    ranges: List[Optional[Tuple[float, float]]] = field(default_factory=list)    # per layer, from scan_ranges
    warnings: List[str] = field(default_factory=list)

    @property
    def chi2_red(self) -> Optional[float]:
        return self.chi2 / self.dof if self.dof > 0 else None


def fit_layers(species: Sequence[Species], assignment: Sequence[int], n_layers: int,
               d0: Optional[Sequence[float]] = None, d_max: float = 10.0, d_min: float = 0.02,
               max_angle: float = DEFAULT_MAX_ANGLE, lam_layers=None, floor: float = 0.03) -> LayerFit:
    """
    Least-squares layer thicknesses (nm) from the angle dependence of each species, by multi-start Nelder-Mead inside
    [d_min, d_max]. Raises ValueError when there is nothing to fit (no layer holds a species, or too few angles).
    """
    used = [k for k in range(n_layers) if k in assignment]
    if not used:
        raise ValueError("Assign at least one species to a layer.")
    if all(measured_normalised(sp, max_angle, floor) is None for sp in species):
        raise ValueError("Need at least two usable angles.")

    def to_d(x):
        d = np.full(n_layers, d_min, dtype=float)
        for j, k in enumerate(used):
            d[k] = x[j]
        return d

    def f(x):
        return gries_chi2(to_d(x), species, assignment, max_angle, lam_layers, floor)[0]

    starts = [np.full(len(used), s) for s in (0.5, 1.5, 4.0)]
    if d0 is not None:
        starts.insert(0, np.array([max(d0[k], d_min) for k in used], dtype=float))
    best = None
    for s in starts:
        x, v = nelder_mead(f, s, [0.5] * len(used), [d_min] * len(used), [d_max] * len(used))
        if best is None or v < best[1]:
            best = (x, v)
    d = to_d(best[0])
    chi2, npts = gries_chi2(d, species, assignment, max_angle, lam_layers, floor)
    fit = LayerFit(d, chi2, max(npts - len(used), 0), list(assignment))
    for j, k in enumerate(used):
        if best[0][j] >= d_max - 1e-6:
            fit.warnings.append(f"Layer {k + 1} ran to the upper limit of {d_max:g} nm: the data cannot pin its thickness.")
        if best[0][j] <= d_min + 1e-6:
            fit.warnings.append(f"Layer {k + 1} is at the lower limit: it may not be there.")
    unused = [k + 1 for k in range(n_layers) if k not in assignment]
    if unused:
        fit.warnings.append("Layer(s) " + ", ".join(map(str, unused)) + " hold no species, so their thickness is not fitted.")
    return fit


def scan_ranges(fit: LayerFit, species: Sequence[Species], n_layers: int, d_max: float = 10.0, d_min: float = 0.02,
                max_angle: float = DEFAULT_MAX_ANGLE, lam_layers=None, floor: float = 0.03, points: int = 41,
                delta: Optional[float] = None) -> LayerFit:
    """
    The range of each layer thickness that fits within the noise: the thickness is held at each value on a grid while the
    others are re-fitted, and the range is where chi-squared stays within `delta` of its minimum (default 1, scaled by the
    reduced chi-squared when that exceeds 1, so unrealistically small errors do not give unrealistically narrow ranges).
    Stores the ranges in `fit.ranges` (None for a layer that holds no species) and returns the fit.
    """
    limit = (1.0 if delta is None else delta) * max(1.0, fit.chi2_red or 1.0)
    used = [k for k in range(n_layers) if k in fit.assignment]
    grid = np.linspace(d_min, d_max, points)
    fit.ranges = [None] * n_layers
    for k in used:
        others = [j for j in used if j != k]
        prof = []
        for g in grid:
            def f(x, g=g):
                d = fit.d_nm.copy()
                d[k] = g
                for jj, j in enumerate(others):
                    d[j] = x[jj]
                return gries_chi2(d, species, fit.assignment, max_angle, lam_layers, floor)[0]
            if others:
                _, v = nelder_mead(f, [fit.d_nm[j] for j in others], [0.3] * len(others), [d_min] * len(others),
                                   [d_max] * len(others), max_iter=150)
            else:
                v = f([])
            prof.append(v)
        ok = grid[np.array(prof) <= fit.chi2 + limit]
        fit.ranges[k] = (float(ok.min()), float(ok.max())) if len(ok) else (float(fit.d_nm[k]), float(fit.d_nm[k]))
    return fit


def default_assignment(species: Sequence[Species], n_layers: int, max_angle: float = DEFAULT_MAX_ANGLE) -> List[int]:
    """
    A first guess for the layer of each species: rank them by average depth (stratification), shallowest in layer 1 and
    the deepest in the substrate; with more species than layers the extra ones share the deepest finite layer.
    """
    rows = stratification(species, max_angle=max_angle)
    order = sorted(range(len(species)), key=lambda i: rows[i].z_nm if rows[i].z_nm is not None else 1e9)
    out = [0] * len(species)
    for rank, i in enumerate(order):
        out[i] = SUBSTRATE_LAYER if rank == len(order) - 1 else min(rank, n_layers - 1)
    return out
