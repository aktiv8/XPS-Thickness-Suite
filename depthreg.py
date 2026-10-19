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
