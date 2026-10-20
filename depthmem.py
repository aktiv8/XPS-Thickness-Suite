"""
Maximum-entropy depth profiles from angle-resolved XPS (numpy only, no Tk). Part of the Depth profile tab, page "ARXPS MEMS
(Rossi)". After M. A. Scorciapino, G. Navarra, B. Elsener and A. Rossi, J. Phys. Chem. C 113, 21328 (2009),
DOI 10.1021/jp906326m (and its Supporting Information), written here from the equations of the paper.

Model (SI Eqs. 1-9). The sample is cut into layers of equal thickness t at depths z_i = i t (i = 0, 1, ...), layer i holding the
atomic fraction n[i, j] of species j with sum_j n[i, j] = 1. The peak of species j at emission angle theta (from the surface
normal) gives

    I_j = k_j sum_i n[i, j] exp(-z_i / (lambda[j, i] cos theta)),      Iinf_j = k_j sum_i exp(-z_i / (lambda[j, i] cos theta))

for the layered sample and for a thick sample of the pure species, and the **apparent concentration** (what an analyst reads off
the sensitivity-corrected intensities, the "apparent concentration diagram", ACD) is

    X[j, theta] = (I_j / Iinf_j) / sum_l (I_l / Iinf_l),      the k_j cancel.

**Entropy.** With m[i, j] the initial estimate of the profile,

    S = sum_ij [ n - m - n ln(n / m) ],     C = chi2 = sum_jk ((X_calc - X_obs) / sigma)^2,     Q = alpha S - C / 2,

and Q is maximised. A large alpha gives a smooth profile close to m that fits the data badly; a small alpha fits the data, noise
included. The paper lets alpha fall from 100 towards 1 over successive calculation cycles. Maximising Q is minimising
F = C / 2 - alpha S, done here over softmax variables (so every layer sums to 1 and no fraction is negative) with an analytic
gradient and the L-BFGS minimiser of depthmurdoch.py.

**Composition-dependent attenuation length** (paper Eqs. 4-7). The IMFP of the electrons of species j in layer i is the mean of
its IMFPs in the materials traversed, weighted by thickness, so it depends on the profile and is recomputed as the profile
changes (an inner iteration here). The paper takes the IMFP of an electron in each material from the G-1 formula of W. H. Gries,
Surf. Interface Anal. 24, 38 (1996) (through NIST SRD 71), implemented here as `g1_imfp_nm`: lambda[nm] = k1 (Va / Z*) E /
(lg E - k2), with Va = mean atomic mass / density (cm3 mol-1), Z* the mean of sqrt(Z), and (k1, k2) from the material class. Each
species is given a material (formula and density) and lambda_mat[j, m] = G-1 IMFP of species j's electrons in material m. As an
alternative, when no materials are given, the attenuation length of species m (typed or estimated elsewhere in the Suite)
stands for material m and is scaled to the energy of species j's electrons as (KE_j / KE_m)^0.75. In layer l a mixed material
has the composition-weighted mean of those lengths and lambda[j, i] is the mean over layers 0..i. NIST warns that G-1 can be off
by up to about 50 % for some materials and that TPP-2M is generally more reliable; G-1 is used because the paper uses it.

The paper's protocol is implemented in the last section of this file: `simulate_start` (the "simulator": a search for the layered
start model against the ACD), `fit_components` ("version 2": four-parameter flat-topped components with Gaussian edges, depths
free within +-0.6 nm), `run_protocol` (version 2, then the free-layer version above, repeated until the profile changes by less
than 1 %) and `minimum_depth_nm` (the minimum-depth test). The start model of `fit_mem` is the entropy prior m.
"""
import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

import depthlayers as DL
import depthmurdoch as DM

DEFAULT_ALPHAS = (100.0, 30.0, 10.0, 3.0, 1.0)
ENERGY_EXPONENT = 0.75               # lambda ~ KE^0.75 when moving an attenuation length to another electron energy
PRIOR_FLOOR = 1e-4                   # the prior never goes below this, so the logarithm of the entropy stays finite
LAMBDA_TOL = 0.01                    # the IMFP iteration stops when no length changes by more than 1 %
MAX_IMFP_ITERATIONS = 8
DEFAULT_REL_SIGMA = 0.06             # about the standard deviation of the paper's +-10 % uniform error on the ACD


# ---------------------------------------------------------------------------------------------------------------
# The model
# ---------------------------------------------------------------------------------------------------------------

def depths(n: int, dz_nm: float) -> np.ndarray:
    """Depths z_i = i dz of the n layers (the first layer is at the surface, as in the paper)."""
    return np.arange(n) * dz_nm


def lambda_matrix(lam_nm: Sequence[float], ke_ev: Optional[Sequence[Optional[float]]] = None) -> np.ndarray:
    """Lam[j, m] = lambda_m (KE_j / KE_m)^0.75: species m's attenuation length for the electrons of species j (diagonal = own)."""
    lam = np.asarray(lam_nm, dtype=float)
    if ke_ev is None or any(k is None or k <= 0 for k in ke_ev):
        ratio = np.ones((len(lam), len(lam)))
    else:
        ke = np.asarray(ke_ev, dtype=float)
        ratio = (ke[:, None] / ke[None, :]) ** ENERGY_EXPONENT
    return lam[None, :] * ratio


def layer_lambda(c: np.ndarray, lam_mat: np.ndarray) -> np.ndarray:
    """
    lambda_ji[j, i]: the attenuation length of species j's electrons from layer i. Layer l of composition c[l, :] is made of
    the species' materials, so its length for species j is sum_m c[l, m] Lam[j, m]; the electrons cross layers 0..i and the
    length is their (equal-thickness) mean (paper Eqs. 4-7).
    """
    per_layer = c @ lam_mat.T                                  # [layer, j]
    return (np.cumsum(per_layer, axis=0) / np.arange(1, len(c) + 1)[:, None]).T


def kernel(z_nm: np.ndarray, lam_ji: np.ndarray, angles_deg: Sequence[float]) -> np.ndarray:
    """K[j, i, k] = exp(-z_i / (lambda[j, i] cos theta_k))."""
    cos = np.cos(np.radians(np.asarray(angles_deg, dtype=float)))
    return np.exp(-z_nm[None, :, None] / (lam_ji[:, :, None] * cos[None, None, :]))


def apparent(c: np.ndarray, k: np.ndarray) -> np.ndarray:
    """X[j, k] for the profile c[i, j] and the kernel K[j, i, k] (the apparent concentration diagram)."""
    return DM.apparent_from_profile(c, k)


def entropy(c: np.ndarray, m: np.ndarray) -> float:
    """S = sum [n - m - n ln(n / m)] (paper Eq. 1): zero at n = m, negative elsewhere."""
    c = np.maximum(c, 1e-300)
    return float(np.sum(c - m - c * np.log(c / m)))


def layered_model(thickness_nm: Sequence[float], compositions: Sequence[Sequence[float]], z_nm: np.ndarray,
                  floor: float = 0.0) -> np.ndarray:
    """
    n0[i, j] for a stack: layer 1 of thickness thickness_nm[0] on top of layer 2, ... and the last layer (no thickness) below
    them all, each with its composition (any scale; normalised). `floor` adds a small share of every species everywhere.
    """
    comp = np.asarray(compositions, dtype=float)
    if comp.ndim != 2 or len(comp) != len(thickness_nm) + 1:
        raise ValueError("Give one composition per layer, with the last one for the bulk: one more than the thicknesses.")
    if np.any(comp < 0) or np.any(comp.sum(axis=1) <= 0):
        raise ValueError("Every layer needs a non-negative composition that is not all zero.")
    comp = comp / comp.sum(axis=1, keepdims=True)
    out = np.zeros((len(z_nm), comp.shape[1]))
    top = 0.0
    for idx, row in enumerate(comp):
        bottom = np.inf if idx == len(comp) - 1 else top + thickness_nm[idx]
        out[(z_nm >= top) & (z_nm < bottom)] = row
        top = bottom
    if floor:
        out = out + floor
        out = out / out.sum(axis=1, keepdims=True)
    return out


def acd_from_layers(thickness_nm: Sequence[float], compositions: Sequence[Sequence[float]], lam_nm: Sequence[float],
                    angles_deg: Sequence[float], dz_nm: float = 0.1, depth_nm: float = 15.0,
                    ke_ev: Optional[Sequence[float]] = None, iterate_imfp: bool = False) -> np.ndarray:
    """
    The apparent concentration diagram X[j, k] of a layered structure (the paper's "simulator"). With fixed lengths (the paper's
    numerical experiments) each species keeps its own lambda; with `iterate_imfp` the lengths follow the composition.
    """
    z = depths(int(round(depth_nm / dz_nm)), dz_nm)
    c = layered_model(thickness_nm, compositions, z)
    lam = np.asarray(lam_nm, dtype=float)
    if iterate_imfp:
        lam_ji = layer_lambda(c, lambda_matrix(lam, ke_ev))
    else:
        lam_ji = np.repeat(lam[:, None], len(z), axis=1)
    return apparent(c, kernel(z, lam_ji, angles_deg))


def add_acd_error(x: np.ndarray, rel: float = 0.10, seed: int = 0) -> np.ndarray:
    """The paper's error model: every ACD value times (1 + u), u uniform in [-rel, +rel]."""
    rng = np.random.default_rng(seed)
    return x * (1.0 + rel * (2.0 * rng.random(x.shape) - 1.0))


# ---------------------------------------------------------------------------------------------------------------
# The fit
# ---------------------------------------------------------------------------------------------------------------

class MemProblem:
    """F(v) = C/2 - alpha S and its gradient for the softmax variables v of the free layers (the rest are held fixed)."""

    def __init__(self, k: np.ndarray, x_meas: np.ndarray, sigma, alpha: float, prior: np.ndarray,
                 fixed: Optional[np.ndarray] = None):
        self.k, self.xm, self.alpha, self.m = k, x_meas, alpha, prior
        self.sigma = np.broadcast_to(np.asarray(sigma, dtype=float), x_meas.shape)
        self.inf = k.sum(axis=1)
        self.n, self.s = k.shape[1], k.shape[0]
        self.fixed = fixed                                   # n x s array of the fixed layers (rows past n_free), or None
        self.n_free = self.n if fixed is None else self.n - len(fixed)

    def profile(self, v: np.ndarray) -> np.ndarray:
        c = DM._softmax(v.reshape(self.n_free, self.s))
        return c if self.fixed is None else np.vstack([c, self.fixed])

    def parts(self, v):
        c = self.profile(v)
        r = np.einsum("ij,jik->jk", c, self.k) / self.inf
        tot = r.sum(axis=0)
        x = r / tot
        e = (x - self.xm) / self.sigma
        chi2 = float(np.sum(e * e))
        return c, r, tot, x, e, chi2

    def value(self, v):
        c, _, _, _, _, chi2 = self.parts(v)
        return 0.5 * chi2 - self.alpha * entropy(c, self.m)

    def value_and_grad(self, v):
        c, r, tot, x, e, chi2 = self.parts(v)
        g = 2.0 * e / self.sigma                                       # dchi2/dX
        gr = (g * tot - np.sum(g * r, axis=0)) / tot ** 2               # through the closure of X
        dc = np.einsum("jk,jik->ij", gr / self.inf, self.k)             # dchi2/dc
        gc = 0.5 * dc + self.alpha * np.log(np.maximum(c, 1e-300) / self.m)       # dF/dc: the entropy gives alpha ln(c/m)
        gc = gc[:self.n_free]
        cf = c[:self.n_free]
        gv = cf * (gc - np.sum(cf * gc, axis=1, keepdims=True))         # through the softmax
        return 0.5 * chi2 - self.alpha * entropy(c, self.m), gv.ravel()


@dataclass
class MemFit:
    names: List[str]
    z_nm: np.ndarray
    n: np.ndarray                          # n[i, j], fractions summing to 1 per layer
    prior: np.ndarray                      # the starting profile used as the entropy prior m
    x_meas: np.ndarray
    x_calc: np.ndarray
    angles: np.ndarray
    lam_ji: np.ndarray                     # [j, i]: the attenuation lengths at the solution
    lam_nm: np.ndarray                     # the species' own lengths
    alpha: float
    chi2: float
    entropy: float
    converged: bool
    cycles: List[Dict[str, float]] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)
    components: Optional[list] = None      # the components of a version-2 fit (None for a free-layer fit)

    @property
    def dz_nm(self) -> float:
        return float(self.z_nm[1] - self.z_nm[0]) if len(self.z_nm) > 1 else 0.1

    @property
    def n_points(self) -> int:
        return int(self.x_meas.size)

    def curve_area_nm(self, j: int) -> float:
        """The area under species j's profile (nm of the species): sum_i n[i, j] dz (the paper's "curve area")."""
        return float(np.sum(self.n[:, j]) * self.dz_nm)

    def peak_position_nm(self, j: int) -> float:
        return float(self.z_nm[int(np.argmax(self.n[:, j]))])

    def crossing_depth_nm(self, species: Sequence[int], level: float = 0.5) -> Optional[float]:
        """Depth where the summed fraction of `species` first falls below `level` (interpolated); None if it never does."""
        total = self.n[:, list(species)].sum(axis=1)
        below = np.flatnonzero(total < level)
        if len(below) == 0:
            return None
        i = int(below[0])
        if i == 0:
            return float(self.z_nm[0])
        a, b = total[i - 1], total[i]
        return float(self.z_nm[i - 1] + (a - level) / (a - b) * (self.z_nm[i] - self.z_nm[i - 1]))

    def interfaces(self) -> List[Tuple[float, str, str]]:
        """[(depth nm, dominant species above, below)] wherever the dominant species of the layers changes."""
        dom = np.argmax(self.n, axis=1)
        out = []
        for i in range(1, len(dom)):
            if dom[i] != dom[i - 1]:
                out.append((float(0.5 * (self.z_nm[i] + self.z_nm[i - 1])), self.names[dom[i - 1]], self.names[dom[i]]))
        return out

    def rdp(self) -> np.ndarray:
        """Relative depth of each species as defined in the paper: X at the largest angle over X at the smallest (model values)."""
        lo, hi = int(np.argmin(self.angles)), int(np.argmax(self.angles))
        return self.x_calc[:, hi] / np.maximum(self.x_calc[:, lo], 1e-12)


def _solve_fixed_lambda(prob: MemProblem, v0: np.ndarray, max_iter: int):
    v, f, iters, ok = DM.lbfgs(prob.value_and_grad, v0, max_iter)
    return v, iters, ok


def fit_mem(names: Sequence[str], angles_deg: Sequence[float], x_meas: np.ndarray, lam_nm: Sequence[float],
            start: np.ndarray, sigma=DEFAULT_REL_SIGMA, ke_ev: Optional[Sequence[Optional[float]]] = None,
            alphas: Sequence[float] = DEFAULT_ALPHAS, iterate_imfp: bool = True, dz_nm: float = 0.1,
            free_depth_nm: Optional[float] = None, target_chi2_per_point: Optional[float] = None,
            relative_sigma: bool = True, max_iter: int = 800, lam_mat: Optional[np.ndarray] = None) -> MemFit:
    """
    Maximise Q = alpha S - C/2 for the profile n[i, j]. `start` (layers x species, rows summing to 1) is the entropy prior m and
    the starting point; its number of rows sets the depth range (rows x dz_nm). `x_meas[j, k]` are the apparent concentrations
    (fractions summing to 1 per angle). `sigma` is the 1-sigma error of X: a relative error if `relative_sigma` (the default,
    sigma = rel X, never below 0.002) else the array of absolute errors. `alphas` is the schedule of calculation cycles, each
    started from the previous solution; with `target_chi2_per_point` the schedule stops at the first (largest) alpha whose
    chi2 per data point is at or below it. `iterate_imfp` recomputes the attenuation lengths from the profile of each cycle
    until none changes by more than 1 %. `free_depth_nm` holds every layer below that depth at the composition of `start` (a
    known bulk. `lam_mat[j, m]` (e.g. from `g1_lambda_matrix`) replaces the energy-scaled own lengths as the length of species
    j's electrons in material m; without the IMFP iteration each species then keeps its own-material length lam_mat[j, j].
    Returns a MemFit.
    """
    names = list(names)
    x_meas = np.asarray(x_meas, dtype=float)
    angles = np.asarray(angles_deg, dtype=float)
    start = np.asarray(start, dtype=float)
    s = len(names)
    if x_meas.shape != (s, len(angles)) or start.ndim != 2 or start.shape[1] != s:
        raise ValueError("x_meas needs one row per species and one column per angle, start one column per species.")
    if s < 2 or len(angles) < 2:
        raise ValueError("At least two species and two angles are needed.")
    if np.any(np.asarray(alphas) <= 0) or dz_nm <= 0:
        raise ValueError("The alphas and the layer thickness must be positive.")
    n = len(start)
    z = depths(n, dz_nm)
    prior = np.maximum(start, PRIOR_FLOOR)
    prior = prior / prior.sum(axis=1, keepdims=True)
    n_free = n if free_depth_nm is None else int(min(max(round(free_depth_nm / dz_nm), 2), n))
    fixed = None if n_free == n else prior[n_free:]
    sig = np.maximum(sigma * x_meas, 0.002) if relative_sigma else np.asarray(sigma, dtype=float)
    lam_own = np.asarray(lam_nm, dtype=float)
    if lam_mat is None:
        lam_mat = lambda_matrix(lam_own, ke_ev)
    else:
        lam_mat = np.asarray(lam_mat, dtype=float)
        if lam_mat.shape != (s, s):
            raise ValueError("lam_mat needs one row and one column per species.")
        lam_own = np.diag(lam_mat).copy()
    v = np.log(prior[:n_free]).ravel()
    cycles, notes = [], []
    c = prior.copy()
    converged = True
    chosen = None
    for alpha in alphas:
        inner = 0
        lam_prev = None
        while True:
            lam_ji = (layer_lambda(c, lam_mat) if iterate_imfp else np.repeat(lam_own[:, None], n, axis=1))
            prob = MemProblem(kernel(z, lam_ji, angles), x_meas, sig, alpha, prior, fixed)
            v, iters, ok = _solve_fixed_lambda(prob, v, max_iter)
            c = prob.profile(v)
            converged = converged and ok
            inner += 1
            if not iterate_imfp or inner >= MAX_IMFP_ITERATIONS:
                break
            new = layer_lambda(c, lam_mat)
            if float(np.max(np.abs(new - lam_ji) / lam_ji)) < LAMBDA_TOL:
                break
        _, _, _, x_calc, _, chi2 = prob.parts(v)
        cycles.append(dict(alpha=float(alpha), chi2=chi2, chi2_per_point=chi2 / x_meas.size, s=entropy(c, prior),
                           iterations=float(iters), imfp_iterations=float(inner), converged=float(ok)))
        chosen = (c.copy(), x_calc, lam_ji, float(alpha), chi2)
        if target_chi2_per_point is not None and chi2 / x_meas.size <= target_chi2_per_point:
            break
    c, x_calc, lam_ji, alpha, chi2 = chosen
    if not converged:
        notes.append("The minimiser stopped before converging in at least one cycle: the result may not be a maximum of Q.")
    if x_meas.size < 0.1 * n:
        notes.append("There are far fewer data points than layers: the profile is determined by the entropy prior (the "
                     "starting model) as much as by the data.")
    return MemFit(names, z, c, prior, x_meas, x_calc, angles, lam_ji, lam_own, alpha, chi2, entropy(c, prior), converged,
                  cycles, notes)


# ---------------------------------------------------------------------------------------------------------------
# The G-1 inelastic mean free path of Gries (1996), as used by the paper through NIST SRD 71
# ---------------------------------------------------------------------------------------------------------------

# symbol -> (atomic number, atomic weight in g/mol, IUPAC values to the precision needed here)
ELEMENTS: Dict[str, Tuple[int, float]] = {
    "H": (1, 1.008), "He": (2, 4.0026), "Li": (3, 6.94), "Be": (4, 9.0122), "B": (5, 10.81), "C": (6, 12.011),
    "N": (7, 14.007), "O": (8, 15.999), "F": (9, 18.998), "Ne": (10, 20.180), "Na": (11, 22.990), "Mg": (12, 24.305),
    "Al": (13, 26.982), "Si": (14, 28.086), "P": (15, 30.974), "S": (16, 32.06), "Cl": (17, 35.45), "Ar": (18, 39.948),
    "K": (19, 39.098), "Ca": (20, 40.078), "Sc": (21, 44.956), "Ti": (22, 47.867), "V": (23, 50.942), "Cr": (24, 51.996),
    "Mn": (25, 54.938), "Fe": (26, 55.845), "Co": (27, 58.933), "Ni": (28, 58.693), "Cu": (29, 63.546), "Zn": (30, 65.38),
    "Ga": (31, 69.723), "Ge": (32, 72.630), "As": (33, 74.922), "Se": (34, 78.971), "Br": (35, 79.904), "Kr": (36, 83.798),
    "Rb": (37, 85.468), "Sr": (38, 87.62), "Y": (39, 88.906), "Zr": (40, 91.224), "Nb": (41, 92.906), "Mo": (42, 95.95),
    "Tc": (43, 98.0), "Ru": (44, 101.07), "Rh": (45, 102.91), "Pd": (46, 106.42), "Ag": (47, 107.87), "Cd": (48, 112.41),
    "In": (49, 114.82), "Sn": (50, 118.71), "Sb": (51, 121.76), "Te": (52, 127.60), "I": (53, 126.90), "Xe": (54, 131.29),
    "Cs": (55, 132.91), "Ba": (56, 137.33), "La": (57, 138.91), "Ce": (58, 140.12), "Pr": (59, 140.91), "Nd": (60, 144.24),
    "Pm": (61, 145.0), "Sm": (62, 150.36), "Eu": (63, 151.96), "Gd": (64, 157.25), "Tb": (65, 158.93), "Dy": (66, 162.50),
    "Ho": (67, 164.93), "Er": (68, 167.26), "Tm": (69, 168.93), "Yb": (70, 173.05), "Lu": (71, 174.97), "Hf": (72, 178.49),
    "Ta": (73, 180.95), "W": (74, 183.84), "Re": (75, 186.21), "Os": (76, 190.23), "Ir": (77, 192.22), "Pt": (78, 195.08),
    "Au": (79, 196.97), "Hg": (80, 200.59), "Tl": (81, 204.38), "Pb": (82, 207.2), "Bi": (83, 208.98), "Th": (90, 232.04),
    "U": (92, 238.03),
}
# (k1, k2) of Gries' Eq. 2 by material class
G1_PARAMS = {"3d": (0.0020, 1.30), "4d": (0.0019, 1.35), "5d": (0.0019, 1.45), "main": (0.0014, 1.10),
             "organic": (0.0018, 1.00), "inorganic": (0.0019, 1.30)}
G1_RANGE_EV = (200.0, 2000.0)         # the energy range over which the parameters were fitted


@dataclass
class Material:
    """A material for the G-1 equation: a formula such as 'SiO2' or 'C6H8O2' and a density in g cm-3."""
    formula: str
    density: float
    category: str = "auto"             # 'auto', 'inorganic' or 'organic' (single elements use their own class)


def parse_formula(formula: str) -> List[Tuple[str, float]]:
    """[(symbol, count)] from a formula without parentheses, e.g. 'Al2O3' -> [('Al', 2), ('O', 3)]; counts may be decimals."""
    out: List[Tuple[str, float]] = []
    i, text = 0, formula.strip()
    if not text:
        raise ValueError("Give a formula such as SiO2.")
    while i < len(text):
        if not text[i].isupper():
            raise ValueError(f"Cannot read the formula '{formula}' (expected an element symbol at '{text[i:]}').")
        j = i + 1
        while j < len(text) and text[j].islower():
            j += 1
        symbol = text[i:j]
        if symbol not in ELEMENTS:
            raise ValueError(f"Unknown element '{symbol}' in '{formula}'.")
        k = j
        while k < len(text) and (text[k].isdigit() or text[k] == "."):
            k += 1
        count = float(text[j:k]) if k > j else 1.0
        if count <= 0:
            raise ValueError("Element counts must be positive.")
        out.append((symbol, count))
        i = k
    return out


def g1_class(composition: Sequence[Tuple[str, float]], category: str = "auto") -> str:
    """The (k1, k2) class: a single element by its place in the periodic table, a compound 'inorganic' or 'organic'."""
    if len(composition) == 1:
        z = ELEMENTS[composition[0][0]][0]
        return "3d" if 22 <= z <= 29 else "4d" if 40 <= z <= 47 else "5d" if 72 <= z <= 79 else "main"
    if category in ("inorganic", "organic"):
        return category
    symbols = {sym for sym, _ in composition}
    return "organic" if {"C", "H"} <= symbols else "inorganic"


def g1_imfp_nm(ke_ev: float, material: Material) -> float:
    """
    The G-1 inelastic mean free path in nm (Gries 1996, Eqs. 2-4; NIST SRD 71 Eqs. A.7-A.9):
    lambda = k1 (Va / Z*) E / (lg E - k2) with Va = sum(q M) / (rho sum q), Z* = sum(q sqrt(Z)) / sum(q). Meant for 200-2000 eV.
    """
    if ke_ev <= 10 ** G1_PARAMS["inorganic"][1] * 1.01 or material.density <= 0:
        raise ValueError("G-1 needs a kinetic energy above about 30 eV (it is fitted for 200-2000 eV) and a positive density.")
    comp = parse_formula(material.formula)
    total = sum(q for _, q in comp)
    z_star = sum(q * math.sqrt(ELEMENTS[sym][0]) for sym, q in comp) / total
    v_a = sum(q * ELEMENTS[sym][1] for sym, q in comp) / (material.density * total)
    k1, k2 = G1_PARAMS[g1_class(comp, material.category)]
    if math.log10(ke_ev) <= k2:
        raise ValueError("The kinetic energy is too low for the G-1 equation of this material class.")
    return k1 * (v_a / z_star) * ke_ev / (math.log10(ke_ev) - k2)


# Common materials for the G-1 input (densities in g cm-3 from standard handbooks, not from the paper; thin films differ).
MATERIAL_PRESETS: Dict[str, Material] = {
    "Si": Material("Si", 2.33), "SiO2": Material("SiO2", 2.196), "Si3N4": Material("Si3N4", 3.17),
    "Al": Material("Al", 2.70), "Al2O3": Material("Al2O3", 3.97), "Ti": Material("Ti", 4.51), "TiO2": Material("TiO2", 4.26),
    "Cr": Material("Cr", 7.19), "Cr2O3": Material("Cr2O3", 5.22), "Fe": Material("Fe", 7.87), "Fe2O3": Material("Fe2O3", 5.24),
    "Fe3O4": Material("Fe3O4", 5.17), "FeOOH": Material("FeOOH", 4.3), "Ni": Material("Ni", 8.91), "NiO": Material("NiO", 6.67),
    "Cu": Material("Cu", 8.96), "Cu2O": Material("Cu2O", 6.0), "CuO": Material("CuO", 6.31), "Zn": Material("Zn", 7.14),
    "ZnO": Material("ZnO", 5.61), "Ag": Material("Ag", 10.5), "W": Material("W", 19.3), "WO3": Material("WO3", 7.16),
    "Au": Material("Au", 19.3), "GaAs": Material("GaAs", 5.32), "Ga2O3": Material("Ga2O3", 5.88),
    "HfO2": Material("HfO2", 9.68), "Carbon (graphite)": Material("C", 2.2),
    "Hydrocarbon (CH2)": Material("CH2", 0.93, "organic"), "PMMA (C5H8O2)": Material("C5H8O2", 1.19, "organic"),
    "Water": Material("H2O", 1.0), "Phosphate (AlPO4)": Material("AlPO4", 2.56),
}


def guess_preset(name: str) -> str:
    """A material for a species from its name: 'Si El' -> Si, 'Si Ox' -> SiO2, 'Al Ox' -> Al2O3, 'C 1s' -> hydrocarbon, 'O 1s' -> SiO2."""
    token = name.replace("_", " ").split()
    sym = token[0] if token else ""
    sym = sym[:1].upper() + sym[1:]
    oxide = len(token) > 1 and token[1].lower().startswith("ox")
    if not oxide and len(sym) > 2 and sym.lower().endswith("ox") and sym[:-2] in ELEMENTS:      # "Alox", "Siox"
        sym, oxide = sym[:-2], True
    if sym[:2] in ("Al", "Ti", "Cr", "Fe", "Ni", "Cu", "Zn", "Ga", "Hf", "W", "Si") and sym not in MATERIAL_PRESETS:
        sym = sym[:2] if sym[:2] in ELEMENTS else sym[:1]
    if sym == "C":
        return "Hydrocarbon (CH2)"
    if sym in ("O", "N"):
        return "SiO2" if sym == "O" else "Si3N4"
    if oxide:
        return {"Si": "SiO2", "Al": "Al2O3", "Ti": "TiO2", "Cr": "Cr2O3", "Fe": "Fe2O3", "Ni": "NiO", "Cu": "CuO", "Zn": "ZnO",
                "W": "WO3", "Ga": "Ga2O3", "Hf": "HfO2"}.get(sym, "SiO2")
    return sym if sym in MATERIAL_PRESETS else "SiO2"


def g1_lambda_matrix(ke_ev: Sequence[float], materials: Sequence[Material]) -> np.ndarray:
    """Lam[j, m] = G-1 IMFP of species j's electrons (kinetic energy ke_ev[j]) in the material of species m."""
    if len(ke_ev) != len(materials) or any(k is None or k <= 0 for k in ke_ev):
        raise ValueError("G-1 needs the kinetic energy of every species and one material per species.")
    return np.array([[g1_imfp_nm(float(k), m) for m in materials] for k in ke_ev])


def g1_range_warning(ke_ev: Sequence[float]) -> Optional[str]:
    """A note when some electron energies lie outside the 200-2000 eV range over which the G-1 parameters were fitted."""
    out = [float(k) for k in ke_ev if k is not None and not G1_RANGE_EV[0] <= k <= G1_RANGE_EV[1]]
    if not out:
        return None
    return (f"Some electron energies ({', '.join(f'{k:.0f}' for k in out)} eV) lie outside 200-2000 eV, where the G-1 "
            "parameters were fitted; the IMFPs there are extrapolations.")


def ensemble_mem(names, angles_deg, x_meas, lam_nm, thickness_nm: Sequence[float], compositions, spread: float = 0.3,
                 n: int = 8, seed: int = 0, depth_nm: float = 15.0, dz_nm: float = 0.1, **kwargs):
    """
    Profiles from `n` starting models whose layer thicknesses are multiplied by random factors between 1 - spread and
    1 + spread (the dependence of the result on the starting profile). Returns (z, mean n, sd n, fits).
    """
    rng = np.random.default_rng(seed)
    z = depths(int(round(depth_nm / dz_nm)), dz_nm)
    fits = []
    for _ in range(n):
        th = [t * (1.0 + spread * (2.0 * rng.random() - 1.0)) for t in thickness_nm]
        fits.append(fit_mem(names, angles_deg, x_meas, lam_nm, layered_model(th, compositions, z, floor=0.01), dz_nm=dz_nm,
                            **kwargs))
    stack = np.array([f.n for f in fits])
    return z, stack.mean(axis=0), stack.std(axis=0), fits


# ---------------------------------------------------------------------------------------------------------------
# The paper's synthetic structures (Table 1 and the text), lengths in nm (the paper gives angstrom)
# ---------------------------------------------------------------------------------------------------------------

PAPER_ANGLES = (0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0)       # the paper does not list the angles; up to 60 deg is its limit
PAPER_STRUCTURES: Dict[str, dict] = {
    "4_2+2": dict(species="ABCD", lam=(4.0, 3.0, 2.0, 1.0), thickness=(1.6,),
                  compositions=((0.2, 0.8, 0, 0), (0, 0, 0.4, 0.6))),
    "7_b": dict(species="ABCDEFG", lam=(4.0, 3.0, 2.0, 1.0, 0.9, 3.5, 2.5), thickness=(0.7, 1.2, 0.6, 1.6),
                compositions=((0, 0.8, 0, 0, 0, 0, 0.2), (0, 0, 0, 0.23, 0, 0.15, 0.62), (1, 0, 0, 0, 0, 0, 0),
                              (0, 0, 0.65, 0, 0.35, 0, 0), (0, 0, 0.82, 0, 0.18, 0, 0))),
    "8": dict(species="ABCDEFGH", lam=(4.0, 3.0, 2.0, 1.0, 4.5, 3.5, 2.5, 2.7), thickness=(0.7, 1.2, 1.2),
              compositions=((0, 0.8, 0, 0, 0, 0, 0, 0.2), (0, 0, 0, 0.23, 0, 0.15, 0.62, 0),
                            (0.4, 0, 0.36, 0, 0.24, 0, 0, 0), (0, 0, 0.82, 0, 0.18, 0, 0, 0))),
}


# ---------------------------------------------------------------------------------------------------------------
# Phase 2: the paper's protocol (simulator, pseudo-Gaussian "version 2", cycles, minimum depth)
# ---------------------------------------------------------------------------------------------------------------

DEPTH_WINDOW_NM = 0.6                # the paper lets the component depths move by +-6 angstrom around the start model
EDGE_RANGE_NM = (0.05, 1.0)          # allowed edge widths (Gaussian sigma) of a component
CYCLE_TOL = 0.01                     # the protocol stops when the profile changes by less than 1 % between cycles
MAX_CYCLES = 6
_BOTTOM = 1.0e3                      # "end" of a component that runs to the bulk


def _erf(x: np.ndarray) -> np.ndarray:
    """erf by Abramowitz and Stegun 7.1.26 (absolute error below 1.5e-7), so that no scipy is needed."""
    sign = np.sign(x)
    ax = np.abs(x)
    t = 1.0 / (1.0 + 0.3275911 * ax)
    poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))))
    return sign * (1.0 - poly * np.exp(-ax * ax))


@dataclass
class Component:
    """One four-parameter component of species `species`: a flat top of height `amp` between `start` and `end` (nm) whose two
    edges are Gaussian error functions of width `width` (the paper's pseudo-Gaussian, defined here since the paper gives no formula)."""
    species: int
    start: float
    end: float
    amp: float
    width: float

    def values(self, z_nm: np.ndarray) -> np.ndarray:
        r = 1.0 / (max(self.width, 1e-6) * math.sqrt(2.0))
        return self.amp * 0.5 * (_erf((z_nm - self.start) * r) - _erf((z_nm - self.end) * r))

    @property
    def to_bottom(self) -> bool:
        return self.end >= 0.5 * _BOTTOM


def components_profile(comps: Sequence[Component], z_nm: np.ndarray, n_species: int, floor: float = 1e-4) -> np.ndarray:
    """n[i, j]: the components of each species added, then normalised so that every layer sums to 1 (floor keeps it positive)."""
    f = np.full((len(z_nm), n_species), floor)
    for c in comps:
        f[:, c.species] += c.values(z_nm)
    return f / f.sum(axis=1, keepdims=True)


def initial_components(start: np.ndarray, z_nm: np.ndarray, width_nm: float = 0.15, min_share: float = 0.05) -> List[Component]:
    """Components from a start profile: for each species every run of layers where it exceeds half its maximum (species that
    never reach `min_share` are left out: version 2 can only reshape what the start model contains)."""
    dz = float(z_nm[1] - z_nm[0]) if len(z_nm) > 1 else 0.1
    out: List[Component] = []
    for j in range(start.shape[1]):
        col = start[:, j]
        top = float(col.max())
        if top < min_share:
            continue
        on = col >= 0.5 * top
        i = 0
        while i < len(on):
            if not on[i]:
                i += 1
                continue
            k = i
            while k + 1 < len(on) and on[k + 1]:
                k += 1
            lo = 0.0 if i == 0 else float(z_nm[i]) - 0.5 * dz
            hi = _BOTTOM if k == len(on) - 1 else float(z_nm[k]) + 0.5 * dz
            out.append(Component(j, lo, hi, float(col[i:k + 1].max()), width_nm))
            i = k + 1
    return out


def _slots(comps: Sequence[Component]) -> List[Tuple[int, Optional[int], int, int]]:
    """Parameter positions in u for each component: (start slot, end slot or None, amplitude slot, width slot)."""
    out, p = [], 0
    for c in comps:
        a = p
        p += 1
        b = None
        if not c.to_bottom:
            b = p
            p += 1
        out.append((a, b, p, p + 1))
        p += 2
    return out


def _decode(u: np.ndarray, base: Sequence[Component], window: float) -> List[Component]:
    out = []
    for c, (a, b, m, w) in zip(base, _slots(base)):
        start = c.start + window * math.tanh(u[a])
        end = _BOTTOM if b is None else c.end + window * math.tanh(u[b])
        amp = c.amp * math.exp(min(max(u[m], -8.0), 8.0))
        width = EDGE_RANGE_NM[0] + (EDGE_RANGE_NM[1] - EDGE_RANGE_NM[0]) / (1.0 + math.exp(-u[w]))
        out.append(Component(c.species, start, end, amp, width))
    return out


def _encode_width(w: float) -> float:
    lo, hi = EDGE_RANGE_NM
    q = min(max((w - lo) / (hi - lo), 1e-3), 1 - 1e-3)
    return math.log(q / (1.0 - q))


def fit_components(names: Sequence[str], angles_deg: Sequence[float], x_meas: np.ndarray, lam_nm: Sequence[float],
                   start: np.ndarray, sigma=DEFAULT_REL_SIGMA, ke_ev: Optional[Sequence[Optional[float]]] = None,
                   iterate_imfp: bool = True, dz_nm: float = 0.1, window_nm: float = DEPTH_WINDOW_NM,
                   relative_sigma: bool = True, max_iter: int = 300, lam_mat: Optional[np.ndarray] = None,
                   components: Optional[Sequence[Component]] = None, rounds: int = 4) -> MemFit:
    """
    The paper's "version 2": the profile is a sum of four-parameter components (start depth, end depth, height, edge width), one
    per run of each species in `start` (or the `components` given), normalised to 1 in every layer, and chi2 is minimised over
    those parameters only (L-BFGS with numerical gradients). Each depth may move by `window_nm` (paper: 6 angstrom) from its
    start value; a depth that ends at the limit of its window is re-centred there and the fit repeated, up to `rounds` times,
    so a wrong start can still move further. The attenuation lengths follow the profile when `iterate_imfp` is set, as in
    `fit_mem`. Returns a MemFit whose `components` hold the fitted components and whose `prior` is `start`.
    """
    names = list(names)
    x_meas = np.asarray(x_meas, dtype=float)
    angles = np.asarray(angles_deg, dtype=float)
    start = np.asarray(start, dtype=float)
    s = len(names)
    if x_meas.shape != (s, len(angles)) or start.ndim != 2 or start.shape[1] != s:
        raise ValueError("x_meas needs one row per species and one column per angle, start one column per species.")
    if dz_nm <= 0 or window_nm <= 0:
        raise ValueError("The layer thickness and the depth window must be positive.")
    z = depths(len(start), dz_nm)
    sig = np.maximum(sigma * x_meas, 0.002) if relative_sigma else np.asarray(sigma, dtype=float)
    lam_own = np.asarray(lam_nm, dtype=float)
    if lam_mat is None:
        lam_mat = lambda_matrix(lam_own, ke_ev)
    else:
        lam_mat = np.asarray(lam_mat, dtype=float)
        if lam_mat.shape != (s, s):
            raise ValueError("lam_mat needs one row and one column per species.")
        lam_own = np.diag(lam_mat).copy()
    prior = np.maximum(start, PRIOR_FLOOR)
    prior = prior / prior.sum(axis=1, keepdims=True)
    comps = list(components) if components is not None else initial_components(prior, z)
    if len({c.species for c in comps}) < 2:
        raise ValueError("Version 2 needs a start model with at least two species present.")
    notes: List[str] = []
    c = components_profile(comps, z, s)
    lam_ji = layer_lambda(c, lam_mat) if iterate_imfp else np.repeat(lam_own[:, None], len(z), axis=1)
    ok_all, inner, total_iters = True, 0, 0
    while True:
        k = kernel(z, lam_ji, angles)
        for _ in range(max(rounds, 1)):
            base = comps

            def chi2_of(u, base=base, kk=k):
                x = apparent(components_profile(_decode(u, base, window_nm), z, s), kk)
                e = (x - x_meas) / sig
                return float(np.sum(e * e))

            def fg(u):
                g = np.empty_like(u)
                for q in range(len(u)):
                    d = np.zeros_like(u)
                    d[q] = 1e-5
                    g[q] = (chi2_of(u + d) - chi2_of(u - d)) / 2e-5
                return chi2_of(u), g

            slots = _slots(base)
            u0 = np.zeros(slots[-1][3] + 1)
            for cc, (_, _, _, w) in zip(base, slots):
                u0[w] = _encode_width(cc.width)
            u, f, iters, ok = DM.lbfgs(fg, u0, max_iter)
            total_iters += iters
            ok_all = ok_all and ok
            comps = _decode(u, base, window_nm)
            at_limit = any(abs(math.tanh(u[i])) > 0.97 for a, b, _, _ in slots for i in (a, b) if i is not None)
            if not at_limit:
                break
        else:
            notes.append("A component depth is still at the edge of its +-%.1f nm window: the start model may be too far off." %
                         window_nm)
        c = components_profile(comps, z, s)
        inner += 1
        if not iterate_imfp or inner >= MAX_IMFP_ITERATIONS:
            break
        new_lam = layer_lambda(c, lam_mat)
        done = float(np.max(np.abs(new_lam - lam_ji) / lam_ji)) < LAMBDA_TOL
        lam_ji = new_lam
        if done:
            break
    x_calc = apparent(c, kernel(z, lam_ji, angles))
    e = (x_calc - x_meas) / sig
    chi2 = float(np.sum(e * e))
    cycles = [dict(alpha=0.0, chi2=chi2, chi2_per_point=chi2 / x_meas.size, s=entropy(c, prior), iterations=float(total_iters),
                   imfp_iterations=float(inner), converged=float(ok_all))]
    if not ok_all:
        notes.append("The minimiser stopped before converging: the components may not be at the best fit.")
    return MemFit(names, z, c, prior, x_meas, x_calc, angles, lam_ji, lam_own, 0.0, chi2, entropy(c, prior), ok_all, cycles, notes,
                  comps)


def _layer_weights(thickness_nm: Sequence[float], lam_nm: np.ndarray, angles_deg: np.ndarray) -> np.ndarray:
    """W[j, l, k] = exp(-a_l / (lam_j cos theta_k)) - exp(-b_l / (lam_j cos theta_k)): the share of a thick pure sample's
    intensity that layer l (from depth a_l to b_l, the last one to infinity) supplies, for fixed lengths."""
    edges = np.concatenate([[0.0], np.cumsum(thickness_nm)])
    cos = np.cos(np.radians(angles_deg))
    arg = lam_nm[:, None] * cos[None, :]
    top = np.exp(-edges[None, :, None] / arg[:, None, :])
    bottom = np.concatenate([top[:, 1:, :], np.zeros_like(top[:, :1, :])], axis=1)
    return top - bottom


def _best_compositions(w: np.ndarray, x_meas: np.ndarray, sig: np.ndarray, start: Optional[np.ndarray], max_iter: int = 200):
    """Compositions c[l, j] (rows sum to 1) minimising chi2 for fixed layer weights w[j, l, k]; reuses MemProblem with alpha 0."""
    s, n_layers, n_ang = w.shape
    prob = MemProblem(w, x_meas, sig, 0.0, np.full((n_layers, s), 1.0 / s))
    prob.inf = np.ones((s, n_ang))
    v0 = np.log(np.full((n_layers, s), 1.0 / s) if start is None else np.maximum(start, 1e-3)).ravel()
    v, f, _, _ = DM.lbfgs(prob.value_and_grad, v0, max_iter)
    return prob.profile(v), 2.0 * f


def simulate_start(names: Sequence[str], angles_deg: Sequence[float], x_meas: np.ndarray, lam_nm: Sequence[float],
                   n_layers: int = 3, sigma=DEFAULT_REL_SIGMA, relative_sigma: bool = True, depth_nm: float = 15.0,
                   dz_nm: float = 0.1, initial_thickness_nm: Optional[Sequence[float]] = None,
                   sizes_nm: Sequence[float] = (0.7, 1.5), max_iter: int = 300):
    """
    The paper's "simulator", automated: the layered structure (n_layers - 1 layers of unknown thickness on a bulk, each with a
    composition) whose ACD is closest to `x_meas`, each species keeping its own attenuation length. The compositions of a stack
    are linear in the intensities, so for each set of thicknesses they are solved exactly (softmax L-BFGS), and a bounded
    Nelder-Mead simplex searches the thicknesses, started from `initial_thickness_nm` and from equal layers of each of
    `sizes_nm`. Returns (thickness_nm, compositions[n_layers, species], start profile on the grid, chi2).
    """
    x_meas = np.asarray(x_meas, dtype=float)
    angles = np.asarray(angles_deg, dtype=float)
    lam = np.asarray(lam_nm, dtype=float)
    if n_layers < 2:
        raise ValueError("The simulator needs at least two layers (one on a bulk).")
    if x_meas.shape != (len(names), len(angles)):
        raise ValueError("x_meas needs one row per species and one column per angle.")
    sig = np.maximum(sigma * x_meas, 0.002) if relative_sigma else np.asarray(sigma, dtype=float)
    state = {"comp": None}

    def cost(th):
        return _best_compositions(_layer_weights(th, lam, angles), x_meas, sig, state["comp"], max_iter)

    def objective(th):
        comp, chi2 = cost(th)
        state["comp"] = comp
        return chi2

    n_th = n_layers - 1
    lo, hi = np.full(n_th, 0.05), np.full(n_th, 0.5 * depth_nm)
    starts = [np.asarray(initial_thickness_nm, dtype=float)] if initial_thickness_nm is not None else []
    starts += [np.full(n_th, float(sz)) for sz in sizes_nm]
    best = None
    for x0 in starts:
        state["comp"] = None
        x, v = DL.nelder_mead(objective, np.clip(x0, lo, hi), np.maximum(0.3 * x0, 0.1), lo, hi, max_iter=150, restarts=1)
        if best is None or v < best[1]:
            best = (x, v)
    th = np.asarray(best[0])
    state["comp"] = None
    comp, chi2 = cost(th)
    z = depths(int(round(depth_nm / dz_nm)), dz_nm)
    return th, comp, layered_model(th, comp, z), float(chi2)


def profile_change(a: np.ndarray, b: np.ndarray) -> float:
    """The mean absolute difference of two profiles' fractions (0.01 = 1 % of composition on average)."""
    return float(np.mean(np.abs(a - b)))


@dataclass
class ProtocolResult:
    fit: MemFit                            # the last free-layer (version 1) fit
    version2: Optional[MemFit]             # the last pseudo-Gaussian (version 2) fit
    history: List[Dict[str, float]]        # per cycle: chi2/point of version 2 and 1 and the change from the cycle before
    converged: bool                        # the profile changed by less than CYCLE_TOL in the last cycle


def run_protocol(names, angles_deg, x_meas, lam_nm, start: np.ndarray, max_cycles: int = MAX_CYCLES, tol: float = CYCLE_TOL,
                 alphas: Sequence[float] = DEFAULT_ALPHAS, dz_nm: float = 0.1, window_nm: float = DEPTH_WINDOW_NM,
                 use_version2: bool = True, fit_kwargs: Optional[dict] = None, **kwargs) -> ProtocolResult:
    """
    The paper's cycle: version 2 (`fit_components`) reshapes the start model, version 1 (`fit_mem` with that result as start and
    prior) lets every layer move, and the pair is repeated, each cycle starting from the last version-1 profile, until the
    profile changes by less than `tol` (mean absolute fraction) between cycles. `kwargs` (sigma, ke_ev, iterate_imfp, lam_mat,
    relative_sigma) go to both fits; `fit_kwargs` (free_depth_nm, target_chi2_per_point) only to the free-layer fit.
    """
    cur = np.asarray(start, dtype=float)
    history: List[Dict[str, float]] = []
    v2 = fit = None
    converged = False
    for cycle in range(1, max_cycles + 1):
        first = cur
        if use_version2:
            v2 = fit_components(names, angles_deg, x_meas, lam_nm, cur, dz_nm=dz_nm, window_nm=window_nm, **kwargs)
            first = v2.n
        fit = fit_mem(names, angles_deg, x_meas, lam_nm, first, alphas=alphas, dz_nm=dz_nm, **kwargs, **(fit_kwargs or {}))
        change = profile_change(fit.n, cur)
        history.append(dict(cycle=float(cycle), chi2_v2=(v2.cycles[0]["chi2_per_point"] if v2 else float("nan")),
                            chi2_v1=fit.cycles[-1]["chi2_per_point"], change=change))
        cur = fit.n
        if change < tol:
            converged = True
            break
    if not converged:
        fit.notes.append("The profile was still changing by %.1f %% after %d cycles." % (100 * history[-1]["change"], max_cycles))
    return ProtocolResult(fit, v2, history, converged)


def minimum_depth_nm(fit: MemFit, tol: float = CYCLE_TOL) -> float:
    """
    The paper's minimum-depth test: the shallowest depth D such that replacing everything below D by the composition of the
    deepest layer (a constrained bulk) changes no apparent concentration by more than `tol` (relative). Below it the profile
    only has to be "bulk" as far as these data can tell.
    """
    n = fit.n
    k = kernel(fit.z_nm, fit.lam_ji, fit.angles)
    full = apparent(n, k)
    for i in range(2, len(n)):
        cut = n.copy()
        cut[i:] = n[-1]
        if float(np.max(np.abs(apparent(cut, k) - full) / np.maximum(full, 1e-9))) < tol:
            return float(fit.z_nm[i - 1] + fit.dz_nm)
    return float(fit.z_nm[-1] + fit.dz_nm)
