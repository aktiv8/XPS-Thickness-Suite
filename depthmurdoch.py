"""
Concentration depth profiling from apparent atomic concentrations by Tikhonov regularisation with the squared concentration
gradient (numpy only, no Tk). Part of the Depth profile tab.

The method is that of B. J. Murdoch and D. G. McCulloch, Surf. Interface Anal. 55, 658-664 (2023),
DOI 10.1002/sia.7240, written here independently from the equations of the paper:

* the data are the apparent atomic concentrations X[j, k] of peak j at emission angle k (already corrected for sensitivity
  factors and transmission, so any factor common to all peaks, such as the analysed area or the tilt geometry, cancels);
* the model has slabs of thickness dz with concentrations c[i, j] (summing to 100 at% over the peaks in every slab) and gives
  the intensity of peak j at angle theta as  I_j = sum_i c[i, j] exp(-z_i / (lambda_j cos theta)),  compared with the same
  sum for a pure infinitely thick sample  Iinf_j = sum_i exp(-z_i / (lambda_j cos theta)),  so that

      X_calc[j] = (I_j / Iinf_j) / sum_l (I_l / Iinf_l);

* Q = chi2 + alpha R is minimised, with chi2 = sum_{j,k} ((X_calc - X_meas) / sigma)^2 and R the squared concentration
  GRADIENT summed over the peaks and the slabs (c in at%, depths in angstrom, as in the authors' spreadsheet, so a value of alpha
  from the paper can be used as it is). The gradient, unlike the squared concentration, does not push a substrate to zero.

The concentrations are parametrised by a softmax over the peaks in each slab, so they sum to 100 at% exactly. The authors'
spreadsheet instead optimises free amplitudes per slab and renormalises to 100 at% only afterwards (the stored variables of
its example sum to between 32 and 148 at%). That model is not well posed as an optimisation: the apparent concentrations do not
depend on a common scale of all the amplitudes, so the regulariser, which falls with the scale, drives them all towards zero
(tried and rejected), and the spreadsheet's result depends on where Solver stops. With closure the same data give the same layer
structure (see the tests) at a higher chi2, because the model has fewer degrees of freedom. The gradient of Q is analytic and a
small L-BFGS minimiser replaces Solver. The slab depths are the slab centres by default; offset=1 gives the spreadsheet's
convention (the slabs lie at depths dz, 2 dz, ...).
"""
import math
from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

import numpy as np

ANGSTROM_PER_NM = 10.0


def eal_murdoch_nm(ke_ev: float) -> float:
    """
    The kinetic-energy-only effective attenuation length used in the spreadsheet that accompanies the paper (attributed there to
    Seah and Dench, Surf. Interface Anal. 1, 2 (1979)): 143 / E^2 + 0.054 sqrt(E) in nm, E the kinetic energy in eV. Semi-
    quantitative: the authors recommend EALs derived from material properties where they are available.
    """
    if ke_ev <= 0:
        raise ValueError("The kinetic energy must be positive.")
    return 143.0 / ke_ev ** 2 + 0.054 * math.sqrt(ke_ev)


def apparent_at_pct(intensities: Sequence[Sequence[float]], rsf: Optional[Sequence[float]] = None) -> np.ndarray:
    """
    Apparent atomic fractions X[j, k] (each angle sums to 1) from corrected intensities[j][k] divided by the sensitivity
    factors (the transmission correction and per-second scaling are the caller's business).
    """
    I = np.asarray(intensities, dtype=float)
    s = np.ones(I.shape[0]) if rsf is None else np.asarray(rsf, dtype=float)
    a = I / s[:, None]
    tot = a.sum(axis=0)
    if np.any(tot <= 0):
        raise ValueError("Every angle needs a positive total intensity.")
    return a / tot


def load_at_pct_csv(path: str) -> Tuple[List[str], np.ndarray, np.ndarray]:
    """
    Read measured apparent concentrations: a header row of peak names after an angle column, then one row per angle with the
    concentration of each peak in at% (or fractions). Returns (names, angles, X[j, k] as fractions summing to 1 per angle).
    """
    rows = []
    with open(path, "r", encoding="utf-8-sig", errors="replace") as fh:
        for line in fh:
            parts = [p.strip() for p in line.replace(";", ",").replace("\t", ",").split(",")]
            if any(parts):
                rows.append(parts)
    if len(rows) < 3:
        raise ValueError("The file needs a header row and at least two angles.")
    names = [n for n in rows[0][1:] if n]
    data = []
    for r in rows[1:]:
        try:
            data.append([float(x) for x in r[:1 + len(names)]])
        except ValueError:
            raise ValueError(f"Cannot read the numbers in the row {r[:1 + len(names)]}.")
    arr = np.array(data)
    if arr.shape[1] != 1 + len(names):
        raise ValueError("Every row needs an angle and one concentration per peak named in the header.")
    angles, x = arr[:, 0], arr[:, 1:].T
    tot = x.sum(axis=0)
    if np.any(tot <= 0):
        raise ValueError("Every angle needs a positive total concentration.")
    return names, angles, x / tot


# ---------------------------------------------------------------------------------------------------------------
# Forward model, objective and gradient
# ---------------------------------------------------------------------------------------------------------------

def depths(n: int, dz_nm: float, offset: float = 0.5) -> np.ndarray:
    """Depths of the n slabs in nm: (i + offset) * dz for i = 0..n-1 (centres by default)."""
    return (np.arange(n) + offset) * dz_nm


def attenuation_kernel(z_nm: np.ndarray, lam_nm: Sequence[float], angles_deg: Sequence[float]) -> np.ndarray:
    """K[j, i, k] = exp(-z_i / (lambda_j cos theta_k))."""
    lam = np.asarray(lam_nm, dtype=float)
    c = np.cos(np.radians(np.asarray(angles_deg, dtype=float)))
    return np.exp(-z_nm[None, :, None] / (lam[:, None, None] * c[None, None, :]))


def apparent_from_profile(c: np.ndarray, K: np.ndarray) -> np.ndarray:
    """X_calc[j, k] for concentrations c[i, j] (any scale per slab is fine) and the kernel."""
    inf = K.sum(axis=1)
    r = np.einsum("ij,jik->jk", c, K) / inf
    return r / r.sum(axis=0)


def _softmax(v: np.ndarray) -> np.ndarray:
    e = np.exp(v - v.max(axis=1, keepdims=True))
    return e / e.sum(axis=1, keepdims=True)


class _Problem:
    """Q(v) and its gradient for the softmax variables v[i, j] (c = softmax(v) over the peaks of each slab)."""

    def __init__(self, K, x_meas, sigma, alpha, dz_nm):
        self.K, self.xm = K, x_meas
        self.sigma = np.broadcast_to(np.asarray(sigma, dtype=float), x_meas.shape)
        self.alpha, self.dz_ang = alpha, dz_nm * ANGSTROM_PER_NM
        self.inf = K.sum(axis=1)
        self.n, self.s = K.shape[1], K.shape[0]

    def parts(self, v):
        c = _softmax(v.reshape(self.n, self.s))
        r = np.einsum("ij,jik->jk", c, self.K) / self.inf
        tot = r.sum(axis=0)
        x = r / tot
        e = (x - self.xm) / self.sigma
        chi2 = float(np.sum(e * e))
        d = 100.0 * (c[1:] - c[:-1])
        grad_term = float(np.sum(d * d) / self.dz_ang)
        return c, r, tot, x, e, chi2, d, grad_term

    def value(self, v):
        _, _, _, _, _, chi2, _, gt = self.parts(v)
        return chi2 + self.alpha * gt

    def value_and_grad(self, v):
        c, r, tot, x, e, chi2, d, gt = self.parts(v)
        g = 2.0 * e / self.sigma                                  # dchi2/dX
        gr = (g * tot - np.sum(g * r, axis=0)) / tot ** 2           # dchi2/dr  (the quotient rule through the closure)
        dc = np.einsum("jk,jik->ij", gr / self.inf, self.K)         # dchi2/dc
        dr = np.zeros_like(c)
        coef = 2.0 * 1e4 / self.dz_ang
        dr[1:] += coef * (c[1:] - c[:-1])
        dr[:-1] -= coef * (c[1:] - c[:-1])
        gc = dc + self.alpha * dr
        gv = c * (gc - np.sum(c * gc, axis=1, keepdims=True))       # through the softmax
        return chi2 + self.alpha * gt, gv.ravel()


def lbfgs(fg, x0: np.ndarray, max_iter: int = 1500, m: int = 10, gtol: float = 1e-8, ftol: float = 1e-10):
    """A compact L-BFGS minimiser with Armijo backtracking (no scipy). Returns (x, f, iterations, converged)."""
    x = np.asarray(x0, dtype=float).copy()
    f, g = fg(x)
    s_list, y_list = [], []
    history = [f]
    for it in range(1, max_iter + 1):
        if np.max(np.abs(g)) < gtol:
            return x, f, it - 1, True
        q = g.copy()
        alphas = []
        for s, y in reversed(list(zip(s_list, y_list))):
            rho = 1.0 / float(y @ s)
            a = rho * float(s @ q)
            alphas.append((a, rho, s, y))
            q -= a * y
        if s_list:
            q *= float(s_list[-1] @ y_list[-1]) / float(y_list[-1] @ y_list[-1])
        for a, rho, s, y in reversed(alphas):
            b = rho * float(y @ q)
            q += s * (a - b)
        direction = -q
        slope = float(g @ direction)
        if slope >= 0:                                           # not a descent direction: reset to steepest descent
            direction, slope = -g, -float(g @ g)
            s_list, y_list = [], []
        step, ok = 1.0 if s_list else min(1.0, 1.0 / max(float(np.linalg.norm(g)), 1e-12)), False
        for _ in range(40):
            xn = x + step * direction
            fn, gn = fg(xn)
            if fn <= f + 1e-4 * step * slope:
                ok = True
                break
            step *= 0.5
        if not ok:
            return x, f, it, False
        s_vec, y_vec = xn - x, gn - g
        if float(y_vec @ s_vec) > 1e-12:
            s_list.append(s_vec)
            y_list.append(y_vec)
            if len(s_list) > m:
                s_list.pop(0)
                y_list.pop(0)
        converged = abs(f - fn) <= ftol * (abs(f) + 1e-30)
        x, f, g = xn, fn, gn
        history.append(f)
        # a long flat valley: no worthwhile gain over the last 40 iterations counts as converged
        stalled = len(history) > 40 and history[-41] - f <= 1e-8 * (abs(f) + 1e-30)
        if converged or stalled:
            return x, f, it, True
    return x, f, max_iter, False


# ---------------------------------------------------------------------------------------------------------------
# The initial-guess layer model, the fit and its summaries
# ---------------------------------------------------------------------------------------------------------------

def layer_guess(order: Sequence[int], thickness_nm: Sequence[float], n_species: int, z_nm: np.ndarray,
                floor: float = 0.01) -> np.ndarray:
    """
    c0[i, j] (fractions summing to 1 per slab) for a stack: the peaks in `order` (indices) from the top, each occupying its
    thickness (nm); the last one in `order` is the substrate and fills the rest. A small `floor` keeps every peak
    present everywhere so the optimiser can move it.
    """
    c0 = np.full((len(z_nm), n_species), floor)
    top = 0.0
    for pos, j in enumerate(order):
        bottom = np.inf if pos == len(order) - 1 else top + thickness_nm[pos]
        c0[(z_nm >= top) & (z_nm < bottom), j] = 1.0
        top = bottom
    return c0 / c0.sum(axis=1, keepdims=True)


@dataclass
class MurdochFit:
    names: List[str]
    z_nm: np.ndarray
    c_pct: np.ndarray                   # c[i, j] in at%, each slab summing to 100
    c0_pct: np.ndarray                  # the initial guess
    x_meas: np.ndarray                  # X[j, k] fractions
    x_calc: np.ndarray
    angles: np.ndarray
    alpha: float
    chi2: float
    r_term: float                       # the squared-gradient term R (before multiplying by alpha)
    q: float
    iterations: int
    converged: bool
    lam_nm: np.ndarray
    notes: List[str] = field(default_factory=list)

    def peak_position_nm(self, j: int) -> float:
        return float(self.z_nm[int(np.argmax(self.c_pct[:, j]))])

    def fwhm_nm(self, j: int) -> Optional[float]:
        """Width at half maximum of species j's profile around its maximum (None if it never falls below half)."""
        c, z = self.c_pct[:, j], self.z_nm
        peak = int(np.argmax(c))
        half = 0.5 * c[peak]
        lo = peak
        while lo > 0 and c[lo - 1] >= half:
            lo -= 1
        hi = peak
        while hi < len(c) - 1 and c[hi + 1] >= half:
            hi += 1

        def edge(i_in, i_out):
            if i_out < 0 or i_out >= len(c):
                return None
            a, b = c[i_in], c[i_out]
            return float(z[i_in] + (half - a) / (b - a) * (z[i_out] - z[i_in])) if b != a else float(z[i_in])
        left, right = edge(lo, lo - 1), edge(hi, hi + 1)
        left = z[lo] - 0.5 * (z[1] - z[0]) if left is None else left          # the profile reaches the surface
        if right is None:
            return None                                                        # still above half at the bottom
        return right - left


def fit(names: Sequence[str], angles_deg: Sequence[float], x_meas: np.ndarray, lam_nm: Sequence[float],
        dz_nm: float = 0.1, depth_nm: float = 6.0, alpha: float = 1e-6, sigma=1.0, c0: Optional[np.ndarray] = None,
        offset: float = 0.5, max_iter: int = 1500, v_start: Optional[np.ndarray] = None) -> MurdochFit:
    """
    Minimise Q = chi2 + alpha R for the profile c(z) of every peak. `c0` (fractions, slabs x peaks) is the initial guess; None
    starts flat with equal shares (no prior knowledge: use a larger alpha). `sigma` is the 1-sigma error of the apparent
    concentrations in FRACTIONS (1.0, as in the spreadsheet, weights all points alike). Returns a MurdochFit.
    """
    names = list(names)
    x_meas = np.asarray(x_meas, dtype=float)
    angles = np.asarray(angles_deg, dtype=float)
    if x_meas.shape != (len(names), len(angles)):
        raise ValueError("x_meas must have one row per peak and one column per angle.")
    if len(names) < 2:
        raise ValueError("At least two peaks are needed (concentrations are relative).")
    if alpha < 0 or dz_nm <= 0 or depth_nm <= dz_nm:
        raise ValueError("alpha must be non-negative and the depth larger than one slab.")
    n = int(round(depth_nm / dz_nm))
    z = depths(n, dz_nm, offset)
    K = attenuation_kernel(z, lam_nm, angles)
    prob = _Problem(K, x_meas, sigma, alpha, dz_nm)
    if v_start is None:
        base = np.full((n, len(names)), 1.0 / len(names)) if c0 is None else np.asarray(c0, dtype=float)
        v_start = np.log(np.maximum(base, 1e-6)).ravel()
    v, q, iters, ok = lbfgs(prob.value_and_grad, v_start, max_iter)
    c, _, _, x, _, chi2, _, gt = prob.parts(v)
    c0_pct = 100.0 * (np.full((n, len(names)), 1.0 / len(names)) if c0 is None else np.asarray(c0, dtype=float))
    notes = []
    if not ok:
        notes.append("The minimiser stopped before converging: the result may not be a minimum.")
    if len(angles) * len(names) < n:
        notes.append("There are far more depth slabs than data points: the profile is determined by the regularisation "
                     "and the initial guess as much as by the data.")
    return MurdochFit(names, z, 100.0 * c, c0_pct, x_meas, x, angles, float(alpha), chi2, gt, q, iters, ok,
                      np.asarray(lam_nm, dtype=float), notes)


def s_curve(names, angles_deg, x_meas, lam_nm, log10_alphas: Sequence[float], dz_nm: float = 0.1,
            depth_nm: float = 6.0, sigma=1.0, c0=None, offset: float = 0.5) -> List[Tuple[float, float, float]]:
    """
    chi2 against log10(alpha), the paper's S-curve: [(log10 alpha, chi2, R)], each fit warm-started from the previous one
    (largest alpha first, so the smooth profile seeds the sharper ones).
    """
    out, start = [], None
    for la in sorted(log10_alphas, reverse=True):
        f = fit(names, angles_deg, x_meas, lam_nm, dz_nm, depth_nm, 10.0 ** la, sigma, c0, offset,
                v_start=start)
        start = np.log(np.maximum(f.c_pct / 100.0, 1e-6)).ravel()
        out.append((float(la), f.chi2, f.r_term))
    return sorted(out)


def suggest_alpha(curve: Sequence[Tuple[float, float, float]], tolerance: float = 0.25) -> float:
    """
    A log10(alpha) from an S-curve: the LARGEST alpha (smoothest profile) whose chi2 is within `tolerance` (25%, plus a
    tiny absolute floor) of the minimum chi2 on the curve. This automates the reading of the paper's S-curve (the
    minimum chi2 on the plateau of small alpha, without chasing the noise); it is a rule of thumb, not part of the paper.
    """
    best = min(c[1] for c in curve)
    ok = [c[0] for c in curve if c[1] <= best * (1.0 + tolerance) + 1e-12]
    return float(max(ok))


def ensemble(names, angles_deg, x_meas, lam_nm, order: Sequence[int], thickness_nm: Sequence[float],
             spread: float = 0.5, n: int = 12, seed: int = 0, dz_nm: float = 0.1, depth_nm: float = 6.0,
             alpha: float = 1e-6, sigma=1.0, offset: float = 0.5):
    """
    Profiles from `n` initial guesses whose layer thicknesses are multiplied by random factors between 1 - spread and
    1 + spread (the paper's test of the dependence on the initial guess). Returns (z, mean at%, sd at%, fits).
    """
    rng = np.random.default_rng(seed)
    nz = int(round(depth_nm / dz_nm))
    z = depths(nz, dz_nm, offset)
    fits = []
    for _ in range(n):
        th = [t * (1.0 + spread * (2.0 * rng.random() - 1.0)) for t in thickness_nm]
        c0 = layer_guess(order, th, len(names), z)
        fits.append(fit(names, angles_deg, x_meas, lam_nm, dz_nm, depth_nm, alpha, sigma, c0, offset))
    stack = np.array([f.c_pct for f in fits])
    return z, stack.mean(axis=0), stack.std(axis=0), fits


# ---------------------------------------------------------------------------------------------------------------
# Inputs from the Depth profile tab's species
# ---------------------------------------------------------------------------------------------------------------

def inputs_from_species(species, max_angle: float = 70.0):
    """
    The apparent concentrations of a set of species (depthprofile.Species) at the angles they all share, up to `max_angle`:
    (names, angles, X[j, k] as fractions summing to 1 per angle, sigma_counting[j, k] propagated from the counting
    statistics of the intensities when every species has them, else None). The intensities are divided by the species'
    sensitivity factors; the transmission correction was applied when the species were built.
    """
    if len(species) < 2:
        raise ValueError("Tick at least two species: concentrations are relative.")
    common = sorted({round(float(a), 6) for a in species[0].angles if a <= max_angle})
    for sp in species[1:]:
        have = {round(float(a), 6) for a in sp.angles}
        common = [a for a in common if a in have]
    if len(common) < 2:
        raise ValueError("The species share fewer than two usable angles.")
    idx = [[int(np.argmin(np.abs(np.asarray(sp.angles, dtype=float) - a))) for a in common] for sp in species]
    a = np.array([[sp.intensity[i] / sp.rsf for i in ids] for sp, ids in zip(species, idx)])
    if np.any(a <= 0):
        raise ValueError("Every species needs a positive intensity at every shared angle.")
    tot = a.sum(axis=0)
    x = a / tot
    sig = None
    if all(sp.sigma is not None for sp in species):
        s_a = np.array([[sp.sigma[i] / sp.rsf for i in ids] for sp, ids in zip(species, idx)])
        sig = np.zeros_like(x)
        for j in range(len(species)):
            d = -a[j][None, :] / tot ** 2 * np.ones((len(species), 1))          # dX_j / da_l = (delta_jl tot - a_j) / tot^2
            d[j] += 1.0 / tot
            sig[j] = np.sqrt(np.sum((d * s_a) ** 2, axis=0))
    return [sp.name for sp in species], np.array(common), x, sig
