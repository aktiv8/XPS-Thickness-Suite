"""
Regularised depth profiles from angle-resolved intensities (numpy only, no Tk). Part of the Depth profile tab; an
independent implementation of the idea in P. J. Cumpson, J. Electron Spectrosc. Relat. Phenom. 73, 25 (1995).

For one species the intensity at emission angle theta is a sum over depth slabs,

    I(theta) / S = sum_n c_n A_n(theta),      A_n(theta) = integral over slab n of exp(-z / (lambda cos theta)) dz

with c_n the concentration (in units of intensity per sensitivity factor per nm, so only ratios between species and
depths matter) and lambda the attenuation length. A few angles cannot determine many slabs, so the profile is found by
Tikhonov regularisation: minimise the weighted misfit plus a smoothness penalty on the second difference,

    || W (A c - b) ||^2 + alpha || L c ||^2 ,          W = diag(1 / sigma),  L = second-difference matrix,

and the strength alpha is chosen from the data (the smoothest profile that fits within the noise, or generalised
cross-validation) unless it is given. Optional prior knowledge sharpens the result: concentrations can be forced
non-negative, and a species can be held constant, or zero, below a chosen depth. The uncertainty of the profile from the
noise comes from the covariance of the linear estimator,

    cov(c) = G G^T,   G = (A^T W^T W A + alpha L^T L)^-1 A^T W^T W,

which excludes the bias that regularisation always adds, so the bands say how far noise alone can move the profile.

A second method, `fit_three_point`, follows the ARCtick spreadsheet (Cumpson 1995, Section 7.1.1, Paynter's model): the profile
is piecewise linear with just three concentrations c1, c2, c3 at three knots whose spacing is set by the depth resolution
the data can support, so the least-squares fit is stable without any smoothing parameter. See that function.
"""
import math
from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

import numpy as np

from depthprofile import DEFAULT_MAX_ANGLE, Species

BULK_FREE, BULK_CONSTANT, BULK_ZERO = "free", "constant", "zero"


def depth_grid(max_depth_nm: float, n: int) -> Tuple[np.ndarray, float]:
    """Slab centres (nm) and the slab thickness for n slabs down to max_depth_nm."""
    dz = max_depth_nm / n
    return (np.arange(n) + 0.5) * dz, dz


def kernel(angles_deg: Sequence[float], lam_nm: float, n: int, dz: float) -> np.ndarray:
    """A[m, n] = lam cos(theta_m) [exp(-z1 / (lam cos)) - exp(-z2 / (lam cos))] for slab n between z1 and z2."""
    lc = lam_nm * np.cos(np.radians(np.asarray(angles_deg, dtype=float)))
    z1 = np.arange(n) * dz
    z2 = z1 + dz
    return lc[:, None] * (np.exp(-z1[None, :] / lc[:, None]) - np.exp(-z2[None, :] / lc[:, None]))


def second_difference(n: int) -> np.ndarray:
    """(n-2) x n second-difference matrix, rows [1, -2, 1] (a straight line has zero penalty)."""
    L = np.zeros((max(n - 2, 0), n))
    for i in range(n - 2):
        L[i, i:i + 3] = (1.0, -2.0, 1.0)
    return L


def nnls(A: np.ndarray, b: np.ndarray, max_iter: Optional[int] = None) -> np.ndarray:
    """Non-negative least squares by the Lawson-Hanson active-set method (no scipy)."""
    m, n = A.shape
    x = np.zeros(n)
    passive = np.zeros(n, dtype=bool)
    w = A.T @ (b - A @ x)
    limit = max_iter or 5 * n + 20
    it = 0
    while not passive.all() and it < limit:
        cand = np.where(~passive, w, -np.inf)
        j = int(np.argmax(cand))
        if cand[j] <= 1e-12:
            break
        passive[j] = True
        it += 1
        while True:
            s = np.zeros(n)
            s[passive] = np.linalg.lstsq(A[:, passive], b, rcond=None)[0]
            if s[passive].min() > 0:
                break
            neg = passive & (s <= 0)
            ratio = x[neg] / np.maximum(x[neg] - s[neg], 1e-300)
            step = float(np.min(ratio)) if len(ratio) else 0.0
            x = x + step * (s - x)
            passive &= x > 1e-14
            if not passive.any():
                s = np.zeros(n)
                break
        x = s
        w = A.T @ (b - A @ x)
    return x


def _solve(Aw: np.ndarray, bw: np.ndarray, L: np.ndarray, alpha: float, nonneg: bool) -> np.ndarray:
    """Weighted Tikhonov solution; non-negative if asked (Lawson-Hanson on the stacked system)."""
    stack = np.vstack([Aw, math.sqrt(alpha) * L]) if len(L) else Aw
    rhs = np.concatenate([bw, np.zeros(len(L))]) if len(L) else bw
    if nonneg:
        return nnls(stack, rhs)
    return np.linalg.lstsq(stack, rhs, rcond=None)[0]


def choose_alpha(Aw: np.ndarray, bw: np.ndarray, L: np.ndarray, method: str = "discrepancy",
                 grid: Optional[np.ndarray] = None) -> float:
    """
    The regularisation strength, from the data. 'discrepancy' (default) takes the SMOOTHEST profile whose weighted misfit
    is no larger than the noise allows (||Aw c - bw||^2 <= number of points): the noise level is known from the counting
    statistics, and with only a handful of angles this is far more reliable than 'gcv', which minimises the generalised
    cross-validation score n ||Aw c - bw||^2 / (n - tr H)^2 and tends to under-smooth badly when there are few points.
    Both search a logarithmic grid scaled to the size of Aw.
    """
    n = len(bw)
    scale = float(np.linalg.norm(Aw, 2)) ** 2 / max(float(np.linalg.norm(L, 2)) ** 2, 1e-12) if len(L) else 1.0
    grid = grid if grid is not None else np.geomspace(1e-10, 1e3, 90) * scale
    scores, resids = [], []
    for a in grid:
        M = Aw.T @ Aw + a * (L.T @ L)
        try:
            Minv = np.linalg.inv(M)
        except np.linalg.LinAlgError:
            scores.append(math.inf)
            resids.append(math.inf)
            continue
        H = Aw @ Minv @ Aw.T
        resid = bw - H @ bw
        r2 = float(resid @ resid)
        dof = n - float(np.trace(H))
        resids.append(r2)
        scores.append(n * r2 / dof ** 2 if dof > 1e-6 else math.inf)
    if method == "gcv":
        return float(grid[int(np.argmin(scores))])
    ok = [i for i, r in enumerate(resids) if r <= n]
    return float(grid[max(ok)]) if ok else float(grid[0])


def _reduced(A: np.ndarray, z: np.ndarray, bulk: str, z_c: Optional[float]):
    """
    Merge or drop the slabs deeper than z_c: 'constant' -> one column holding the whole bulk, 'zero' -> removed. Returns
    the reduced matrix and, for each full slab, the index of its reduced parameter (-1 where it is held at zero).
    """
    n = len(z)
    if bulk == BULK_FREE or z_c is None or z_c >= z[-1]:
        return A, list(range(n))
    deep = z > z_c
    if not deep.any() or deep.all():
        return A, list(range(n))
    top = np.flatnonzero(~deep)
    pos = {int(i): j for j, i in enumerate(top)}
    if bulk == BULK_ZERO:
        return A[:, top], [pos.get(i, -1) for i in range(n)]
    merged = np.column_stack([A[:, top], A[:, deep].sum(axis=1)])
    return merged, [pos.get(i, len(top)) for i in range(n)]


@dataclass
class ProfileFit:
    name: str
    z_nm: np.ndarray
    c: np.ndarray                       # concentration (intensity / S per nm), per slab
    c_se: Optional[np.ndarray]          # 1-sigma from the noise (unconstrained estimator), per slab
    alpha: float
    chi2_red: Optional[float]
    predicted: np.ndarray               # I / S at the measured angles
    measured: np.ndarray                # I / S at the measured angles
    sigma: np.ndarray
    angles: np.ndarray
    negative: bool = False              # the unconstrained solution had negative parts
    notes: List[str] = field(default_factory=list)


def fit_profile(sp: Species, max_depth_nm: float = 12.0, n_depth: int = 20, alpha: Optional[float] = None,
                nonneg: bool = False, bulk: str = BULK_FREE, z_c_nm: Optional[float] = None,
                max_angle: float = DEFAULT_MAX_ANGLE, floor: float = 0.03, method: str = "discrepancy") -> ProfileFit:
    """
    Regularised concentration-depth profile of one species. `alpha` None chooses it from the data (`method`).
    `bulk` is 'free', 'constant' (one value below z_c_nm) or 'zero' (nothing below z_c_nm).
    """
    th = np.asarray(sp.angles, dtype=float)
    ok = np.flatnonzero((th <= max_angle) & (sp.intensity > 0))
    if len(ok) < 3:
        raise ValueError(f"{sp.name}: a depth profile needs at least three usable angles.")
    angles = th[ok]
    b = sp.intensity[ok] / sp.rsf
    sig = (sp.sigma[ok] / sp.rsf) if sp.sigma is not None else b * floor
    sig = np.maximum(sig, 1e-9 * float(np.max(b)))
    z, dz = depth_grid(max_depth_nm, n_depth)
    A_full = kernel(angles, sp.lam_nm, n_depth, dz)
    A, mapping = _reduced(A_full, z, bulk, z_c_nm)
    L = second_difference(A.shape[1])
    Aw, bw = A / sig[:, None], b / sig
    a = alpha if alpha is not None else choose_alpha(Aw, bw, L, method)
    c_par = _solve(Aw, bw, L, a, nonneg)
    free = _solve(Aw, bw, L, a, False)
    notes: List[str] = []
    negative = bool(np.min(free) < -1e-9 * max(float(np.max(np.abs(free))), 1e-30))
    if negative and not nonneg:
        notes.append("The profile has negative concentrations: an assumption or the kernel (lambda, geometry) may be "
                     "wrong, or the data have lost the information; try the non-negative option.")
    M = Aw.T @ Aw + a * (L.T @ L)
    try:
        G = np.linalg.solve(M, Aw.T)
        se_par = np.sqrt(np.sum(G * G, axis=1))
        dof = len(b) - float(np.trace(Aw @ G))
    except np.linalg.LinAlgError:
        se_par, dof = None, len(b) - 1.0

    def expand(arr):
        return np.array([arr[i] if i >= 0 else 0.0 for i in mapping])

    pred = A @ c_par
    chi2 = float(np.sum(((pred - b) / sig) ** 2))
    notes.append("depth resolution is roughly the attenuation length times cos θ at the largest angle "
                 f"({sp.lam_nm * math.cos(math.radians(float(angles.max()))):.2f} nm)")
    return ProfileFit(sp.name, z, expand(c_par), expand(se_par) if se_par is not None else None, float(a),
                      chi2 / dof if dof > 0.5 else None, pred, b, sig, angles, negative, notes)


def atomic_fraction(profiles: Sequence[ProfileFit]) -> np.ndarray:
    """Concentrations of several species on the same depth grid scaled so that they add up to 1 at every depth."""
    c = np.maximum(np.array([p.c for p in profiles]), 0.0)
    total = np.sum(c, axis=0)
    return c / np.where(total > 0, total, 1.0)


# ---------------------------------------------------------------------------------------------------------------
# Three-point (Paynter) fit, as in the ARCtick spreadsheet
# ---------------------------------------------------------------------------------------------------------------
#
# P. J. Cumpson, J. Electron Spectrosc. Relat. Phenom. 73, 25 (1995), Section 4.2 and 7.1.1 (Eqs. 11, 12, 23-26). The data
# hold about three degrees of freedom, so fit a profile with three parameters: c(z) = c1 down to z1, straight lines to c2
# at z2 and c3 at z3, then constant at c3 (the bulk) or, if the species is absent from the bulk, zero. With L = lambda cos
# theta and the knots' slopes s1 = (c2 - c1)/(z2 - z1), s2 = (c3 - c2)/(z3 - z2),
#
#     I(theta)/S = integral c(z) exp(-z/L) dz = L [ c1 + s1 L e^(-z1/L) + (s2 - s1) L e^(-z2/L) - s2 L e^(-z3/L) ]
#
# which is linear in c1, c2, c3. The knots are z1 = lambda cos(theta_max), z2 = z1 e^(pi/w), z3 = z2 e^(pi/w) with the highest
# resolvable generalised frequency w = acosh(pi (n - 1) (I/sigma)^2) / pi for n angles of relative uncertainty sigma/I.

DEFAULT_REL_SIGMA = 0.05       # the spreadsheet's "approx % uncertainty" cell


def omega_max(n_angles: int, rel_sigma: float) -> float:
    """Eq. 11: the highest generalised frequency the data can resolve (the spreadsheet's K27)."""
    return math.acosh(max(math.pi * (n_angles - 1) / rel_sigma ** 2, 1.0)) / math.pi


def three_point_knots(lam_nm: float, theta_max_deg: float, n_angles: int, rel_sigma: float = DEFAULT_REL_SIGMA
                      ) -> Tuple[Tuple[float, float, float], float]:
    """((z1, z2, z3) in nm, omega_max): Eqs. 24-26 with the geometric spacing of Eq. 12."""
    w = omega_max(n_angles, rel_sigma)
    z1 = lam_nm * math.cos(math.radians(theta_max_deg))
    ratio = math.exp(math.pi / w)
    return (z1, z1 * ratio, z1 * ratio * ratio), w


def three_point_matrix(angles_deg: Sequence[float], lam_nm: float, knots: Sequence[float], tail: bool = True) -> np.ndarray:
    """G[m, j] with I_m / S = sum_j G[m, j] c_j for the model above; `tail` False sets the concentration to zero below z3."""
    z1, z2, z3 = knots
    h1, h2 = z2 - z1, z3 - z2
    L = lam_nm * np.cos(np.radians(np.asarray(angles_deg, dtype=float)))
    e1, e2, e3 = np.exp(-z1 / L), np.exp(-z2 / L), np.exp(-z3 / L)
    g1 = L * (1.0 - L * e1 / h1 + L * e2 / h1)
    g2 = L * L * (e1 / h1 - e2 / h2 - e2 / h1 + e3 / h2)
    g3 = L * L * (e2 - e3) / h2
    if not tail:
        g3 = g3 - L * e3
    return np.column_stack([g1, g2, g3])


def three_point_profile(z_nm: np.ndarray, c: Sequence[float], knots: Sequence[float], tail: bool = True) -> np.ndarray:
    """The piecewise-linear concentration at depths z (constant c1 above z1, c3 below z3 or zero when `tail` is False)."""
    z = np.asarray(z_nm, dtype=float)
    out = np.interp(z, knots, c)
    if not tail:
        out = np.where(z > knots[2], 0.0, out)
    return out


def fit_three_point(sp: Species, max_depth_nm: float = 12.0, n_depth: int = 60, tail: bool = True,
                    max_angle: float = DEFAULT_MAX_ANGLE, rel_sigma: Optional[float] = None, floor: float = 0.03
                    ) -> ProfileFit:
    """
    ARCtick-style three-point depth profile of one species: weighted linear least squares for c1..c3 (no smoothing
    parameter). `rel_sigma` (the data's relative uncertainty, which sets the knots) defaults to the median relative
    sigma of the areas when known, else 5%. The profile is returned on a grid of `n_depth` slab centres down to
    `max_depth_nm` so it plots with the other method; `c_se` is the noise-only 1-sigma band from the covariance of
    c1..c3, (A_w^T A_w)^-1 as on the spreadsheet.
    """
    th = np.asarray(sp.angles, dtype=float)
    ok = np.flatnonzero((th <= max_angle) & (sp.intensity > 0))
    if len(ok) < 3:
        raise ValueError(f"{sp.name}: a three-point profile needs at least three usable angles.")
    angles = th[ok]
    b = sp.intensity[ok] / sp.rsf
    if sp.sigma is not None:
        sig = sp.sigma[ok] / sp.rsf
        rel = rel_sigma if rel_sigma else max(float(np.median(sig / b)), 0.01)
    else:
        rel = rel_sigma if rel_sigma else DEFAULT_REL_SIGMA
        sig = b * (rel_sigma if rel_sigma else floor)
    sig = np.maximum(sig, 1e-9 * float(np.max(b)))
    knots, w = three_point_knots(sp.lam_nm, float(angles.max()), len(ok), rel)
    G = three_point_matrix(angles, sp.lam_nm, knots, tail)
    Gw, bw = G / sig[:, None], b / sig
    c = np.linalg.lstsq(Gw, bw, rcond=None)[0]
    cov = np.linalg.pinv(Gw.T @ Gw)
    z, _ = depth_grid(max_depth_nm, n_depth)
    basis = np.column_stack([three_point_profile(z, e, knots, tail) for e in np.eye(3)])
    se = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", basis, cov, basis), 0.0))
    pred = G @ c
    chi2 = float(np.sum(((pred - b) / sig) ** 2))
    dof = len(b) - 3
    notes = [f"knots {knots[0]:.2f}, {knots[1]:.2f}, {knots[2]:.2f} nm (resolvable frequency {w:.2f}, "
             f"{100 * rel:.1f}% uncertainty, {len(ok)} angles)"]
    if dof <= 0:
        notes.append("three angles fix three concentrations exactly: there is no check on the fit")
    if float(np.min(c)) < -1e-9 * max(float(np.max(np.abs(c))), 1e-30):
        notes.append("a concentration is negative: an assumption (lambda, geometry, bulk) may be wrong")
    if knots[2] > max_depth_nm:
        notes.append(f"the third knot is deeper than the plotted {max_depth_nm:g} nm")
    return ProfileFit(sp.name, z, three_point_profile(z, c, knots, tail), se, float("nan"),
                      chi2 / dof if dof > 0 else None, pred, b, sig, angles, bool(np.min(c) < 0), notes)
