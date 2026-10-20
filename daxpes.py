"""
Dual-anode XPS (DAXPES) thickness and coverage of a thin film from the peak area of ONE core level measured with a soft
and a hard X-ray source (standard library only, no Tk). B. J. Murdoch, Surf. Interface Anal. 58, 125 (2026), Eqs. 3, 6, 12.

A film of thickness d covering a fraction `cov` of the analysed area gives, for each source (attenuation length lam,
take-off angle theta from the surface normal), the peak area relative to a thick reference of the same material

    r = I / I_inf = cov (1 - exp(-d / (lam cos theta)))                                           (3), (12)

A single source therefore cannot separate d from the coverage. But the coverage cancels in the ratio of the two sources,

    (I_H / I_inf,H) / (I_X / I_inf,X) = (1 - exp(-d / (lam_H cos theta))) / (1 - exp(-d / (lam_X cos theta)))      (6)

(X = soft, H = hard source), which gives d on its own; the coverage is then what is left over. Equation 6 needs the two
sources to see the SAME analysed area (the analysis area must not depend on the source), and the reference ratio
rho = I_inf,H / I_inf,X, which can be measured on a thick reference of the same material with both sources or typed in.
The ratio r_H / r_X lies between lam_X / lam_H (a vanishingly thin film) and 1 (a film thicker than the soft source's
information depth), so the method is sensitive only in between, and d cannot be trusted much beyond the elastic limit
3 lam_H cos theta of the hard electrons (Murdoch, Section 3).
"""
import math
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

# Photon energies (eV) of common laboratory sources (Murdoch used Al Ka, Ag La and Ga Ka).
SOURCES: Dict[str, float] = {"Al Kα": 1486.7, "Mg Kα": 1253.6, "Ag Lα": 2984.3, "Ga Kα": 9252.0}

ELASTIC_LIMIT_FACTOR = 3.0      # the thickest film the hard-source electrons can report is about 3 lam_H cos theta
NEAR_LIMIT = 0.9                # warn above this fraction of that limit
SOFT_SATURATED = 0.95           # warn when the soft source has already reached this fraction of its thick-film signal

# Film materials for the TPP-2M starting attenuation lengths: (density g/cm3, molar mass g/mol, valence electrons per
# formula unit, band gap eV). The band gaps of TiO2, Ga2O3, SnO2 and Si are those Murdoch used (Table 1); densities and
# valence-electron counts are handbook values, not from his paper. TPP-2M gives IMFPs; attenuation lengths are shorter.
FILM_MATERIALS: Dict[str, Tuple[float, float, float, float]] = {
    "TiO2": (4.26, 79.87, 16, 3.2),
    "Ga2O3": (5.88, 187.44, 24, 4.8),
    "SnO2": (6.95, 150.71, 16, 3.8),
    "SiO2": (2.196, 60.08, 16, 8.9),
    "Si": (2.33, 28.09, 4, 1.1),
}


def tpp2m_start_nm(tpp2m: Callable[..., float], material: str, ke_ev: float, factor: float = 1.0) -> float:
    """Starting attenuation length: the TPP-2M IMFP of `material` at `ke_ev`, times `factor` (AL/IMFP)."""
    rho, mass, nv, gap = FILM_MATERIALS[material]
    return factor * float(tpp2m(ke_ev, rho, mass, nv, gap))


# Seah's universal attenuation length (M. P. Seah, Surf. Interface Anal. 44, 1353 (2012), Eq. 11, "S3"): only the average
# atomic number Z, the average atomic size a (Eq. 7) and, for strongly bonded compounds, a small term W are needed, and it
# agrees with the attenuation lengths computed from optical data to an RMS 8 % (100 eV to 30 keV, emission up to 65 deg).
#     L = (5.8 + 0.0041 Z^1.7 + 0.088 E^0.93) a^1.82 / (Z^0.38 (1 - W))   nm,    a^3 = 1e21 M / (rho N_A (g + h)) nm^3
# with E the kinetic energy in eV, W = 0 for an element and 0.06 H (H the heat of formation in eV per atom) or 0.02 E_g for a
# compound. Z is the mean atomic number over the atoms of the formula unit. Composition (Z, atoms) and the heat of
# formation per atom (kJ/mol from standard tables, divided by the atoms per formula unit and by 96.485) per film material;
# the heats are handbook values, not from Seah's paper.
AVOGADRO = 6.02214076e23
SEAH_INPUTS: Dict[str, Tuple[Tuple[Tuple[int, int], ...], float]] = {
    "TiO2": (((22, 1), (8, 2)), 944.0 / 96.485 / 3),
    "Ga2O3": (((31, 2), (8, 3)), 1089.0 / 96.485 / 5),
    "SnO2": (((50, 1), (8, 2)), 577.6 / 96.485 / 3),
    "SiO2": (((14, 1), (8, 2)), 910.7 / 96.485 / 3),
    "Si": (((14, 1),), 0.0),
}


def seah_a_nm(molar_mass: float, density_g_cm3: float, atoms_per_formula: float) -> float:
    """Seah's average atomic size a (nm) from Eq. 7: a^3 = 1e21 M / (rho N_A (g + h))."""
    return (1e21 * molar_mass / (density_g_cm3 * AVOGADRO * atoms_per_formula)) ** (1.0 / 3.0)


def seah_al_nm(ke_ev: float, z_avg: float, a_nm: float, w: float = 0.0) -> float:
    """Seah (2012) Eq. 11 ("S3"): the attenuation length in nm at kinetic energy `ke_ev`; `w` is 0 for an element."""
    if ke_ev <= 0 or z_avg <= 0 or a_nm <= 0 or not 0.0 <= w < 1.0:
        raise ValueError("Seah's equation needs a positive energy, atomic number and atomic size, and 0 <= W < 1.")
    return (5.8 + 0.0041 * z_avg ** 1.7 + 0.088 * ke_ev ** 0.93) * a_nm ** 1.82 / (z_avg ** 0.38 * (1.0 - w))


def seah_start_nm(material: str, ke_ev: float) -> float:
    """Seah's S3 attenuation length of a film material in FILM_MATERIALS at `ke_ev` (W = 0.06 H; 0 for an element)."""
    comp, heat = SEAH_INPUTS[material]
    rho, mass, _, _ = FILM_MATERIALS[material]
    atoms = sum(n for _, n in comp)
    z_avg = sum(z * n for z, n in comp) / atoms
    return seah_al_nm(ke_ev, z_avg, seah_a_nm(mass, rho, atoms), 0.06 * heat)


def film_factor(d_nm: float, lam_nm: float, theta_deg: float) -> float:
    """1 - exp(-d / (lam cos theta)): the signal of a complete film of thickness d relative to a thick one."""
    return -math.expm1(-d_nm / (lam_nm * math.cos(math.radians(theta_deg))))


def forward_ratios(d_nm: float, cov: float, lam_soft: float, lam_hard: float, theta_deg: float = 0.0
                   ) -> Tuple[float, float]:
    """(r_soft, r_hard) = cov (1 - e^(-d/(lam cos theta))) for a film of thickness d covering a fraction `cov`."""
    return (cov * film_factor(d_nm, lam_soft, theta_deg), cov * film_factor(d_nm, lam_hard, theta_deg))


def single_source_thickness(r: float, lam_nm: float, theta_deg: float = 0.0) -> float:
    """Eq. 3 inverted: the thickness a single source reports, t = -lam cos(theta) ln(1 - I/I_inf), assuming full coverage."""
    if not (math.isfinite(r) and 0.0 < r < 1.0):
        raise ValueError("I/I∞ must be between 0 and 1 for a thickness (a film cannot give more than the thick reference).")
    return -lam_nm * math.cos(math.radians(theta_deg)) * math.log1p(-r)


def ratio_limits(lam_soft: float, lam_hard: float) -> Tuple[float, float]:
    """The range (lam_X / lam_H, 1) that r_H / r_X can take; the hard source must be the more penetrating one."""
    if not (lam_soft > 0 and lam_hard > lam_soft):
        raise ValueError("The hard source needs the longer attenuation length (check which source is which).")
    return lam_soft / lam_hard, 1.0


def daxpes_thickness(q: float, lam_soft: float, lam_hard: float, theta_deg: float = 0.0,
                     label: str = "area ratio") -> float:
    """
    Eq. 6 solved for d by bisection, where q = (I_H/I_inf,H) / (I_X/I_inf,X). q rises monotonically from lam_X/lam_H
    (d -> 0) to 1 (d -> infinity); outside that range there is no solution and a ValueError says why.
    """
    lo, hi = ratio_limits(lam_soft, lam_hard)
    if not math.isfinite(q):
        raise ValueError("The area ratio is not a number.")
    if q <= lo:
        raise ValueError(f"The {label} {q:.3g} is at or below the thin-film limit λ_soft/λ_hard = {lo:.3g}: the film is "
                         "thinner than the two sources can resolve, or the reference ratio is wrong.")
    if q >= hi:
        raise ValueError(f"The {label} {q:.3g} is at or above 1: both sources already see a thick film (or the "
                         "reference ratio is wrong), so the thickness is not determined.")
    c = math.cos(math.radians(theta_deg))

    def ratio(d):
        return film_factor(d, lam_hard, theta_deg) / film_factor(d, lam_soft, theta_deg)

    a, b = 1e-9 * lam_soft * c, lam_soft * c
    while ratio(b) < q:                       # grow the bracket until it holds the answer
        b *= 2.0
        if b > 1e6 * lam_hard:
            raise ValueError("The thickness is beyond anything the two sources can resolve.")
    for _ in range(200):
        m = 0.5 * (a + b)
        if ratio(m) < q:
            a = m
        else:
            b = m
        if b - a < 1e-13 * b:
            break
    return 0.5 * (a + b)


@dataclass
class DaxpesResult:
    d_nm: float                              # DAXPES thickness, independent of coverage
    q: float                                 # (I_H/I_inf,H) / (I_X/I_inf,X)
    d_se_nm: Optional[float] = None
    t_soft_nm: Optional[float] = None        # single-source (full coverage assumed) thicknesses
    t_hard_nm: Optional[float] = None
    coverage: Optional[float] = None         # exact within the model: r / (1 - e^(-d/(lam cos theta)))
    coverage_se: Optional[float] = None
    coverage_paper: Optional[float] = None   # hard-source single-source thickness / DAXPES thickness (Murdoch, Fig. 11)
    limit_nm: float = 0.0                    # elastic limit 3 lam_H cos theta
    warnings: List[str] = field(default_factory=list)
    method: str = "area"                     # "area" (Eq. 6) or "background" (Eqs. 7-8)
    d_star_nm: Optional[float] = None        # background method: the thickness measure d* of Eq. 7
    d_simple_nm: Optional[float] = None      # background method: 2 d* cos(theta), the thin-film reading without the 3-lambda cap
    b1_soft: Optional[float] = None          # background method: B1 of each source (eV^2)
    b1_hard: Optional[float] = None
    lam_band: Dict[str, Optional[Tuple[float, float, Optional[float], Optional[float]]]] = field(default_factory=dict)
    # lam_band["ratio"] / ["common"] = (d at -band, d at +band, coverage at -band, coverage at +band) when the hard-source
    # attenuation length only / both lengths together are changed by +-`band` (None if a variant has no solution)


def _propagate(f: Callable[[Sequence[float]], float], x: Sequence[float], sig: Sequence[Optional[float]]) -> Optional[float]:
    """1-sigma of f(x) from independent 1-sigma errors `sig` (None = exact), by central differences."""
    var = 0.0
    for i, s in enumerate(sig):
        if not s:
            continue
        h = 1e-6 * abs(x[i]) if x[i] else 1e-9
        up, dn = list(x), list(x)
        up[i] += h
        dn[i] -= h
        var += ((f(up) - f(dn)) / (2.0 * h) * s) ** 2
    return math.sqrt(var) if var > 0 else None


def analyse(i_soft: float, i_hard: float, lam_soft: float, lam_hard: float, theta_deg: float = 0.0,
            iinf_soft: Optional[float] = None, iinf_hard: Optional[float] = None, rho: Optional[float] = None,
            sig_soft: Optional[float] = None, sig_hard: Optional[float] = None,
            sig_inf_soft: Optional[float] = None, sig_inf_hard: Optional[float] = None,
            band: float = 0.10) -> DaxpesResult:
    """
    DAXPES thickness and, when the thick-reference areas are known, the single-source thicknesses and the coverage.
    Give either both reference areas (`iinf_soft`, `iinf_hard`) or just their ratio `rho` = I_inf,H / I_inf,X (then only d
    is returned: without a reference for each source there is no coverage). The `sig_*` arguments are 1-sigma
    uncertainties of the same quantities (None = ignore); the thickness error is the counting-statistics part only.
    """
    if i_soft <= 0 or i_hard <= 0:
        raise ValueError("Both peak areas must be positive.")
    have_inf = iinf_soft is not None and iinf_hard is not None
    if have_inf:
        if iinf_soft <= 0 or iinf_hard <= 0:
            raise ValueError("The reference areas must be positive.")
        rho = iinf_hard / iinf_soft
    elif rho is None or not rho > 0:
        raise ValueError("Give the reference areas of both sources, or their ratio.")
    q = (i_hard / i_soft) / rho
    d = daxpes_thickness(q, lam_soft, lam_hard, theta_deg)
    res = DaxpesResult(d, q, limit_nm=ELASTIC_LIMIT_FACTOR * lam_hard * math.cos(math.radians(theta_deg)))

    def d_of(v):                                          # v = [I_X, I_H, rho]
        return daxpes_thickness((v[1] / v[0]) / v[2], lam_soft, lam_hard, theta_deg)

    rel_rho = None
    if have_inf and (sig_inf_soft or sig_inf_hard):
        rel_rho = math.hypot((sig_inf_soft or 0.0) / iinf_soft, (sig_inf_hard or 0.0) / iinf_hard)
    elif not have_inf:
        rel_rho = None
    sigs = [sig_soft, sig_hard, rho * rel_rho if rel_rho else None]
    res.d_se_nm = _propagate(d_of, [i_soft, i_hard, rho], sigs)

    c = math.cos(math.radians(theta_deg))
    if have_inf:
        r_x, r_h = i_soft / iinf_soft, i_hard / iinf_hard
        if r_x < 1.0:
            res.t_soft_nm = single_source_thickness(r_x, lam_soft, theta_deg)
        if r_h < 1.0:
            res.t_hard_nm = single_source_thickness(r_h, lam_hard, theta_deg)
            res.coverage_paper = res.t_hard_nm / d
        res.coverage = r_x / film_factor(d, lam_soft, theta_deg)

        def cov_of(v):                                    # v = [I_X, I_H, Iinf_X, Iinf_H]
            dd = daxpes_thickness(((v[1] / v[0]) / (v[3] / v[2])), lam_soft, lam_hard, theta_deg)
            return (v[0] / v[2]) / film_factor(dd, lam_soft, theta_deg)

        res.coverage_se = _propagate(cov_of, [i_soft, i_hard, iinf_soft, iinf_hard],
                                     [sig_soft, sig_hard, sig_inf_soft, sig_inf_hard])
        if res.coverage > 1.05:
            res.warnings.append(f"The coverage comes out at {100 * res.coverage:.0f} %: above 100 % means the areas, the "
                                "reference or the attenuation lengths are inconsistent.")
    if band:
        def at(ls, lh):
            """(d, coverage) for other attenuation lengths, or None where Eq. 6 has no solution for them."""
            try:
                dd = daxpes_thickness(q, ls, lh, theta_deg)
            except ValueError:
                return None
            cov = (i_soft / iinf_soft) / film_factor(dd, ls, theta_deg) if have_inf else None
            return dd, cov

        for name, variants in (("ratio", ((lam_soft, lam_hard * (1 - band)), (lam_soft, lam_hard * (1 + band)))),
                               ("common", ((lam_soft * (1 - band), lam_hard * (1 - band)),
                                           (lam_soft * (1 + band), lam_hard * (1 + band))))):
            lo, hi = at(*variants[0]), at(*variants[1])
            res.lam_band[name] = (None if lo is None or hi is None else
                                  (lo[0], hi[0], lo[1], hi[1]))
    if d > NEAR_LIMIT * res.limit_nm:
        res.warnings.append(f"The thickness is close to the elastic limit 3 λ_hard cos θ = {res.limit_nm:.1f} nm of the "
                            "hard-source electrons; a thicker film would read about the same.")
    if film_factor(d, lam_soft, theta_deg) > SOFT_SATURATED:
        res.warnings.append("The soft source is almost saturated (the film is thicker than about 3 λ_soft cos θ): the "
                            "thickness then depends on a small difference and is poorly determined.")
    if d < 0.3 * lam_soft * c:
        res.warnings.append("Very thin film: the two sources differ little, so the thickness is sensitive to errors in "
                            "the reference ratio and the attenuation lengths.")
    return res


# ---------------------------------------------------------------------------------------------------------------
# Background method (Murdoch 2026, Eqs. 5, 7, 8): the ratio of the inelastic-background parameters B1 of the two sources
# ---------------------------------------------------------------------------------------------------------------
#
# B1 / B0 = 1 - exp(-d / (2 lam cos theta))                                                          (5)
# B1_H / B1_X = (1 - e^(-d*/lam_H)) / (1 - e^(-d*/lam_X)),    d* = d / (2 cos theta) if elastic scattering is ignored  (7)
# d = (1 - exp(-d* / (lam_H cos theta))) 3 lam_H cos theta                                           (8)
#
# B0 (the thick-sample value) cancels in the ratio, which needs no reference sample and no knowledge of the flux or
# cross-section of either source, only that the two sources analyse the same material (B0 equal for both). Eq. 7 has the
# same form as Eq. 6 with theta = 0, so d* is found by the same solver. Eq. 8 then converts d* to the film thickness and
# caps it at 3 lam_H cos theta, the elastic limit of the hard-source electrons. For a thin film Eq. 8 gives about 3 d*
# where the simple relation of Eq. 5 gives 2 d* cos theta, so both are returned.


def background_thickness(b1_soft: float, b1_hard: float, lam_soft: float, lam_hard: float, theta_deg: float = 0.0
                         ) -> Tuple[float, float, float, float]:
    """(d, d*, 2 d* cos(theta), B1 ratio) from the B1 values (eV^2) of the soft and hard source; d is Eq. 8."""
    if not (b1_soft > 0 and b1_hard > 0):
        raise ValueError("Both B1 values must be positive.")
    q = b1_hard / b1_soft
    d_star = daxpes_thickness(q, lam_soft, lam_hard, 0.0, label="B1 ratio")
    c = math.cos(math.radians(theta_deg))
    d = 3.0 * lam_hard * c * (-math.expm1(-d_star / (lam_hard * c)))
    return d, d_star, 2.0 * d_star * c, q


def analyse_background(b1_soft: float, b1_hard: float, lam_soft: float, lam_hard: float, theta_deg: float = 0.0,
                       sig_b1_soft: Optional[float] = None, sig_b1_hard: Optional[float] = None,
                       i_soft: Optional[float] = None, iinf_soft: Optional[float] = None,
                       i_hard: Optional[float] = None, iinf_hard: Optional[float] = None,
                       band: float = 0.10) -> DaxpesResult:
    """
    Thickness from the B1 values of the two sources and, if a peak area and a thick-reference area are given for the soft
    source (and the hard one), the coverage c = (I_X/I_inf,X) / (1 - exp(-d/(lam_X cos theta))) of the film.
    """
    d, d_star, d_simple, q = background_thickness(b1_soft, b1_hard, lam_soft, lam_hard, theta_deg)
    c = math.cos(math.radians(theta_deg))
    res = DaxpesResult(d, q, method="background", d_star_nm=d_star, d_simple_nm=d_simple, b1_soft=b1_soft,
                       b1_hard=b1_hard, limit_nm=ELASTIC_LIMIT_FACTOR * lam_hard * c)
    res.d_se_nm = _propagate(lambda v: background_thickness(v[0], v[1], lam_soft, lam_hard, theta_deg)[0],
                             [b1_soft, b1_hard], [sig_b1_soft, sig_b1_hard])
    if i_soft is not None and iinf_soft is not None and i_soft > 0 and iinf_soft > 0:
        r_x = i_soft / iinf_soft
        res.coverage = r_x / film_factor(d, lam_soft, theta_deg)
        if r_x < 1.0:
            res.t_soft_nm = single_source_thickness(r_x, lam_soft, theta_deg)
        if i_hard is not None and iinf_hard is not None and i_hard > 0 and 0 < i_hard / iinf_hard < 1.0:
            res.t_hard_nm = single_source_thickness(i_hard / iinf_hard, lam_hard, theta_deg)
            res.coverage_paper = res.t_hard_nm / d
        if res.coverage > 1.05:
            res.warnings.append(f"The coverage comes out at {100 * res.coverage:.0f} %: above 100 % means the B1 values, "
                                "the areas or the attenuation lengths are inconsistent.")
    if band:
        for name, variants in (("ratio", ((lam_soft, lam_hard * (1 - band)), (lam_soft, lam_hard * (1 + band)))),
                               ("common", ((lam_soft * (1 - band), lam_hard * (1 - band)),
                                           (lam_soft * (1 + band), lam_hard * (1 + band))))):
            got = []
            for ls, lh in variants:
                try:
                    dd = background_thickness(b1_soft, b1_hard, ls, lh, theta_deg)[0]
                except ValueError:
                    got.append(None)
                    continue
                cov = ((i_soft / iinf_soft) / film_factor(dd, ls, theta_deg)
                       if i_soft and iinf_soft else None)
                got.append((dd, cov))
            res.lam_band[name] = None if None in got else (got[0][0], got[1][0], got[0][1], got[1][1])
    if d > NEAR_LIMIT * res.limit_nm:
        res.warnings.append(f"The thickness is close to the elastic limit 3 λ_hard cos θ = {res.limit_nm:.1f} nm; a "
                            "thicker film would read about the same (Eq. 8 saturates there).")
    if d_star < 0.2 * lam_soft:
        res.warnings.append("Very thin film: the two B1 values differ little, and the paper found the background method "
                            "underestimates the thickness of ultra-thin films (below about 3 nm).")
    return res
