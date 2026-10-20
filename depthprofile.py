"""
Depth-profile analysis of angle-resolved XPS (numpy only, no Tk): the species behind the depth-profile panels, and the
Surrey Trend Plot statistics.

Inspired by the methods of the NPL ARCtick package (P. J. Cumpson, J. Electron Spectrosc. Relat. Phenom. 73, 25
(1995)); everything here is an independent implementation from the published physics, not a copy of that software.

A *species* is anything whose intensity can be followed against emission angle: one chemical state of a CasaXPS fit
(casa.py), or a whole region integrated with a Shirley background (arxps.py). Angles are from the surface NORMAL
(0 = normal emission), the same convention as the rest of the Suite.

For a species that sits in a thin layer at depth z, I(theta) = I(0) exp(-z sec(theta) / lambda), so ln I is a straight
line in sec(theta) with slope -z/lambda. That is the idea behind the Surrey Trend Plot: species with parallel lines
(and similar lambda) may be distributed together, and the scatter about the line shows the experimental precision.
"""
import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

import arxps as AR
import casa

DEFAULT_MAX_ANGLE = AR.DEFAULT_MAX_ANGLE
SIMILAR_KE_RATIO = (0.7, 1.4)       # kinetic-energy ratios outside this give quite different attenuation lengths


@dataclass
class Species:
    name: str
    region: str                         # the block (region) it comes from
    angles: np.ndarray                  # emission angles from the surface normal, ascending
    intensity: np.ndarray               # corrected areas (counts/s x eV, transmission-corrected if asked)
    sigma: Optional[np.ndarray] = None  # 1-sigma counting-statistics uncertainty of each intensity, if known
    lam_nm: float = 2.0                 # attenuation length of this species' photoelectrons (editable)
    rsf: float = 1.0                    # relative sensitivity factor
    ke_ev: Optional[float] = None
    source: str = "fit"                 # "fit" or "region"
    key: str = ""

    @property
    def sec(self) -> np.ndarray:
        return 1.0 / np.cos(np.radians(self.angles))


def area_sigma(block, ke_lo: Optional[float] = None, ke_hi: Optional[float] = None, transmission_factor: float = 1.0
               ) -> float:
    """
    1-sigma counting-statistics uncertainty of an integrated area in counts/s x eV: area = step x sum(y_i) with
    y_i = counts_i / (scans x dwell), so var(area) = step^2 x sum(counts_i) / (scans x dwell)^2 over the integrated
    points (all points of the block when no kinetic-energy limits are given or none fall inside them).
    """
    counts = np.asarray(block.counts, dtype=float)
    be = np.asarray(block.be_ev, dtype=float)
    if len(be) < 2:
        return float("nan")
    step = float(np.mean(np.abs(np.diff(be))))
    mask = np.ones(len(be), dtype=bool)
    if ke_lo is not None and ke_hi is not None:
        ke = block.hv_ev - be
        inside = (ke >= ke_lo) & (ke <= ke_hi)
        if inside.sum() >= 3:
            mask = inside
    t = (block.scans or 1) * (block.dwell_s or 1.0)
    return step * math.sqrt(max(float(counts[mask].sum()), 0.0)) / t / max(transmission_factor, 1e-12)


def species_from_fits(series: Sequence[Tuple[float, object]], angle_offset: float = 0.0, transmission: bool = True,
                      fit_region: Optional[str] = None) -> List[Species]:
    """
    One species per chemical state of the CasaXPS fit stored in a region's blocks (series = [(tilt, Block), ...], lowest
    angle first). Areas are used as stored; with `transmission` each component area is divided by the block's
    transmission function at the component's kinetic energy. Raises ValueError for a block without a fit.
    """
    fits = AR.read_fits(series, fit_region)
    names = AR.state_names(fits)
    angles = np.array([f.tilt + angle_offset for f in fits], dtype=float)
    out = []
    for key, name in names.items():
        area, sig, ke_w, ke_a = [], [], 0.0, 0.0
        for f in fits:
            st = f.states.get(key)
            a = s = 0.0
            if st is not None:
                for c in st.comps:
                    tf = AR._transmission_factor(f.block, c.ke_raw) if transmission else 1.0
                    a += c.area / tf
                    ke_w += c.area * c.ke_raw
                    ke_a += c.area
                s = area_sigma(f.block, transmission_factor=AR._transmission_factor(f.block, st.comps[0].ke_raw)
                               if transmission else 1.0)
            area.append(a)
            sig.append(s)
        out.append(Species(name=name, region=fits[0].block.block_id, angles=angles, intensity=np.array(area),
                           sigma=np.array(sig), ke_ev=(ke_w / ke_a if ke_a > 0 else None), source="fit", key=key))
    return out


def auto_envelope(be: np.ndarray, y: np.ndarray, edge_ev: float = 0.5, level: float = 0.10, pad_ev: float = 1.5
                  ) -> "AR.Window":
    """
    The part of a region to integrate: from where the smoothed spectrum first rises above `level` of its range to where it
    last falls below it, plus `pad_ev` each side (never beyond the data less `edge_ev`). A Shirley background over a
    whole wide region, with a long loss tail at the high-BE end, can fail to converge and give nonsense areas.
    """
    sm = AR.smooth(np.asarray(y, dtype=float), 7)
    cut = float(sm.min()) + level * float(np.ptp(sm))
    above = np.flatnonzero(sm > cut)
    lo = float(be[above[0]]) - pad_ev if len(above) else float(be[0])
    hi = float(be[above[-1]]) + pad_ev if len(above) else float(be[-1])
    return AR.Window(max(lo, float(be[0]) + edge_ev), min(hi, float(be[-1]) - edge_ev))


def species_from_region(series: Sequence[Tuple[float, object]], angle_offset: float = 0.0, per_second: bool = True,
                        transmission: bool = True, edge_ev: float = 0.5) -> Species:
    """
    A whole region as one species: Shirley-corrected area over the peaks of the region (`auto_envelope`), each spectrum
    aligned to the previous one so a charging drift does not matter.
    """
    first = series[0][1]
    be0 = np.asarray(first.be_ev, dtype=float)
    reference = (be0, AR.intensity(first, per_second, transmission))
    win = auto_envelope(be0, reference[1], edge_ev)
    area, sig, angles = [], [], []
    for tilt, blk in series:
        ex = AR.single_peak_area(blk.be_ev, AR.intensity(blk, per_second, transmission), win, reference, grid=be0)
        reference = (ex.grid, ex.y)
        tf = 1.0
        if transmission and blk.transmission is not None and np.all(np.asarray(blk.transmission) > 0):
            tf = float(np.mean(blk.transmission)) / AR.TRANSMISSION_SCALE
        s = area_sigma(blk, transmission_factor=tf)
        area.append(ex.area_sub)
        sig.append(s if not per_second else s)
        angles.append(tilt + angle_offset)
    ke = first.hv_ev - float(np.mean(be0))
    return Species(name=first.species or first.block_id, region=first.block_id, angles=np.array(angles, dtype=float),
                   intensity=np.array(area), sigma=np.array(sig), ke_ev=ke, source="region", key=first.block_id)


# ---------------------------------------------------------------------------------------------------------------
# Surrey Trend Plot statistics
# ---------------------------------------------------------------------------------------------------------------

def apply_geometry(species: Sequence[Species], power: float) -> List[Species]:
    """
    Copies of the species with every intensity (and its uncertainty) divided by cos(theta) ** power. Absolute intensities
    carry any change of the analysed area with tilt, which is common to all species. Power 0 when the X-ray beam defines
    the area (it does not change with tilt). When the analyser defines it the area is A0 sec(theta) (Cumpson, J. Electron
    Spectrosc. 73, 25 (1995), Section 2) so the raw intensity rises as sec(theta) and the reduced intensity I = F / (T A
    sigma L) needs it removed: power -1, i.e. multiply by cos(theta). That is also what the ARCtick spreadsheet does.
    Ratios of species at one angle, as on the ARXPS tab, do not need this.
    """
    out = []
    for sp in species:
        if not power:
            out.append(sp)
            continue
        g = np.cos(np.radians(sp.angles)) ** power
        out.append(Species(sp.name, sp.region, sp.angles, sp.intensity / g,
                           None if sp.sigma is None else sp.sigma / g, sp.lam_nm, sp.rsf, sp.ke_ev, sp.source, sp.key))
    return out


@dataclass
class TrendRow:
    name: str
    slope: Optional[float]              # d ln I / d sec(theta) = -z/lambda
    slope_se: Optional[float]
    intercept: Optional[float]
    rms: Optional[float]                # rms residual of ln I about the line (the experimental scatter)
    depth_nm: Optional[float]           # apparent depth -lambda x slope
    n: int = 0
    ke_ev: Optional[float] = None
    note: str = ""
    similar: List[str] = field(default_factory=list)       # names of species whose trend matches this one


def trend_fit(angles: Sequence[float], intensity: Sequence[float], max_angle: float = DEFAULT_MAX_ANGLE
              ) -> Tuple[Optional[float], Optional[float], Optional[float], Optional[float], int]:
    """
    Straight-line fit of ln I against sec(theta) over the positive intensities at angles up to `max_angle`:
    (slope, standard error of the slope, intercept, rms residual of ln I, points used). Fewer than two points give
    all None; two points give no error or scatter.
    """
    th, y = np.asarray(angles, dtype=float), np.asarray(intensity, dtype=float)
    ok = (y > 0) & (th <= max_angle) & (th >= 0) & (th < 90)
    x, ly = 1.0 / np.cos(np.radians(th[ok])), np.log(y[ok])
    n = int(ok.sum())
    if n < 2 or np.ptp(x) <= 0:
        return None, None, None, None, n
    sxx = float(np.sum((x - x.mean()) ** 2))
    slope = float(np.sum((x - x.mean()) * (ly - ly.mean())) / sxx)
    icpt = float(ly.mean() - slope * x.mean())
    if n < 3:
        return slope, None, icpt, None, n
    res = ly - icpt - slope * x
    s = math.sqrt(float(np.sum(res ** 2)) / (n - 2))
    return slope, s / math.sqrt(sxx), icpt, math.sqrt(float(np.mean(res ** 2))), n


def trend_table(species: Sequence[Species], max_angle: float = DEFAULT_MAX_ANGLE) -> List[TrendRow]:
    """
    Trend statistics for every species, with the co-location hint: two species look alike when their slopes differ
    by less than twice the combined standard error, their kinetic energies are within a factor of 0.7-1.4 (so their
    attenuation lengths are similar) and neither is a single-point fit. A hint, never a proof.
    """
    rows = []
    for sp in species:
        slope, se, icpt, rms, n = trend_fit(sp.angles, sp.intensity, max_angle)
        depth = -sp.lam_nm * slope if slope is not None else None
        note = "" if n >= 3 else "needs three or more angles for a scatter estimate"
        rows.append(TrendRow(sp.name, slope, se, icpt, rms, depth, n, sp.ke_ev, note))
    for i, a in enumerate(rows):
        for j in range(i + 1, len(rows)):
            b = rows[j]
            if None in (a.slope, b.slope, a.slope_se, b.slope_se):
                continue
            tol = 2.0 * math.hypot(a.slope_se, b.slope_se)
            ke_ok = (a.ke_ev is None or b.ke_ev is None
                     or SIMILAR_KE_RATIO[0] <= a.ke_ev / b.ke_ev <= SIMILAR_KE_RATIO[1])
            if abs(a.slope - b.slope) <= tol and ke_ok:
                a.similar.append(b.name)
                b.similar.append(a.name)
    for r in rows:
        ke_far = [o.name for o in rows if o is not r and r.ke_ev and o.ke_ev
                  and not SIMILAR_KE_RATIO[0] <= r.ke_ev / o.ke_ev <= SIMILAR_KE_RATIO[1]]
        if ke_far and not r.note:
            r.note = "kinetic energy differs a lot from " + ", ".join(ke_far[:2]) + ": compare slopes with care"
    return rows


def normalised(sp: Species) -> np.ndarray:
    """Intensities divided by the value at the first angle with a positive intensity (for parallel-line comparison)."""
    y = np.asarray(sp.intensity, dtype=float)
    pos = np.flatnonzero(y > 0)
    return y / y[pos[0]] if len(pos) else y


# ---------------------------------------------------------------------------------------------------------------
# Stratification (average depth and amount of each species)
# ---------------------------------------------------------------------------------------------------------------
#
# A thin layer of areal amount n at depth z gives I(theta) = S n exp(-z sec(theta) / lambda), S the sensitivity factor.
# Two angles a and b therefore give its depth and amount:
#
#     R = I_b / I_a = exp(-z (sec_b - sec_a) / lambda)     ->   z = -lambda ln R / (sec_b - sec_a)
#     n = I_a exp(z sec_a / lambda) / S                      (the intensity corrected for the overlying attenuation)
#
# A species spread over several depths is reported at one average depth, so a species in two distinct layers is
# misleading; curvature of ln I against sec(theta) is the warning.

@dataclass
class StratRow:
    name: str
    z_nm: Optional[float]
    z_se_nm: Optional[float]
    amount: Optional[float]             # areal amount, arbitrary units (intensity / S, corrected for attenuation)
    amount_pct: Optional[float] = None  # share of the total amount of the species used
    ratio: Optional[float] = None       # I_b / I_a (pair mode)
    curved_t: Optional[float] = None    # significance of the curvature of ln I against sec(theta)
    note: str = ""


def curvature_t(angles: Sequence[float], intensity: Sequence[float], max_angle: float = DEFAULT_MAX_ANGLE
                ) -> Optional[float]:
    """t statistic of the quadratic term when ln I is fitted against sec(theta) (None with fewer than 4 points)."""
    th, y = np.asarray(angles, dtype=float), np.asarray(intensity, dtype=float)
    ok = (y > 0) & (th <= max_angle) & (th >= 0) & (th < 90)
    if ok.sum() < 4:
        return None
    x, ly = 1.0 / np.cos(np.radians(th[ok])), np.log(y[ok])
    x0 = x - x.mean()
    design = np.column_stack([np.ones_like(x0), x0, x0 ** 2])
    coef, res, rank, _ = np.linalg.lstsq(design, ly, rcond=None)
    if rank < 3:
        return None
    dof = len(x) - 3
    resid = ly - design @ coef
    s2 = float(resid @ resid) / dof
    if s2 < 1e-24:                                     # an exact straight line: nothing to be curved about
        return 0.0
    cov = s2 * np.linalg.inv(design.T @ design)
    se = math.sqrt(max(cov[2, 2], 0.0))
    return float(coef[2] / se) if se > 0 else (0.0 if abs(coef[2]) < 1e-12 else math.inf)


def stratification(species: Sequence[Species], theta_a: Optional[float] = None, theta_b: Optional[float] = None,
                   use_all: bool = False, max_angle: float = DEFAULT_MAX_ANGLE, curved_limit: float = 3.0
                   ) -> List[StratRow]:
    """
    Average depth and relative amount of each species. Pair mode (default) uses two emission angles: `theta_a` (default the
    smallest angle) and `theta_b` (default the angle nearest 67.5 degrees that is within `max_angle`). `use_all` fits
    ln I against sec(theta) over every angle instead. The uncertainty of the depth comes from the counting statistics
    (pair mode, when the areas have one) or from the scatter about the line (all angles).
    """
    rows: List[StratRow] = []
    for sp in species:
        th = np.asarray(sp.angles, dtype=float)
        usable = np.flatnonzero((th <= max_angle) & (sp.intensity > 0))
        if len(usable) < 2:
            rows.append(StratRow(sp.name, None, None, None, note="needs positive intensities at two usable angles"))
            continue
        t_curve = curvature_t(sp.angles, sp.intensity, max_angle)
        warn = ""
        if t_curve is not None and abs(t_curve) > curved_limit:
            warn = "ln I is curved against sec θ: the species may be in more than one layer"
        if use_all:
            slope, se, icpt, _, n = trend_fit(sp.angles, sp.intensity, max_angle)
            if slope is None:
                rows.append(StratRow(sp.name, None, None, None, note="not enough angles"))
                continue
            z = -sp.lam_nm * slope
            amount = math.exp(icpt) / sp.rsf
            rows.append(StratRow(sp.name, z, (sp.lam_nm * se if se is not None else None), amount, None, None, t_curve,
                                 warn if z >= 0 else "intensity rises with angle: species lies at the surface "
                                 "(depth is negative)"))
            continue
        ia = int(usable[np.argmin(np.abs(th[usable] - (theta_a if theta_a is not None else th[usable].min())))])
        want_b = theta_b if theta_b is not None else 67.5
        cand = [i for i in usable if i != ia]
        ib = int(min(cand, key=lambda i: abs(th[i] - want_b)))
        sec_a, sec_b = 1.0 / math.cos(math.radians(th[ia])), 1.0 / math.cos(math.radians(th[ib]))
        if abs(sec_b - sec_a) < 1e-9:
            rows.append(StratRow(sp.name, None, None, None, note="the two angles are the same"))
            continue
        ratio = float(sp.intensity[ib] / sp.intensity[ia])
        z = -sp.lam_nm * math.log(ratio) / (sec_b - sec_a)
        z_se = None
        if sp.sigma is not None and sp.intensity[ia] > 0 and sp.intensity[ib] > 0:
            rel = math.hypot(sp.sigma[ia] / sp.intensity[ia], sp.sigma[ib] / sp.intensity[ib])
            z_se = sp.lam_nm * rel / abs(sec_b - sec_a)
        amount = float(sp.intensity[ia]) * math.exp(z * sec_a / sp.lam_nm) / sp.rsf
        note = warn
        if ratio > 1.0 and sec_b > sec_a:
            note = "intensity rises with angle: species lies at the surface (depth is negative)"
        rows.append(StratRow(sp.name, z, z_se, amount, None, ratio, t_curve, note))
    total = sum(r.amount for r in rows if r.amount is not None and r.amount > 0)
    for r in rows:
        if r.amount is not None and total > 0:
            r.amount_pct = 100.0 * max(r.amount, 0.0) / total
    return rows


# ---------------------------------------------------------------------------------------------------------------
# Relative depth plot (RDP): which species lie nearer the surface, with no depth scale
# ---------------------------------------------------------------------------------------------------------------
#
# C. R. Brundle, G. Conti, P. Mack, J. Electron Spectrosc. Relat. Phenom. 178-179, 433 (2010), Section 4.1 (after Seah et
# al.); R. N. S. Sodhi et al., Surf. Interface Anal. 49, 1309 (2017). For each species the ratio of the logarithm of its
# peak intensity at a high (near-grazing) angle to that at a low (near-normal) angle: the nearer the surface a species
# sits, the more it gains with angle, so the larger the value. It orders layers; it gives no depth, no thickness and no
# way of seeing the same species in two layers (that is averaged), and species with very different kinetic energies have
# different attenuation lengths, which the plain index ignores.
#
# Two forms are offered. 'published' is the papers' ln(I_high) / ln(I_low), which only has a meaning for dimensionless
# intensities, so each species is expressed as its atomic % (I/S as a share of all the species used) at each angle; that
# also removes any geometry factor common to every species. It is undefined where a share is at or below about 1 % (the
# logarithm is zero or negative). 'ratio' is the unit-free ln(I_high / I_low) = -(sec_high - sec_low) z / lambda, the
# Stratification depth in units of lambda with the sign reversed; a factor common to all species only shifts every
# value by the same amount, so the order never changes.

RDP_PUBLISHED, RDP_RATIO = "published", "ratio"
RDP_MIN_SHARE_PCT = 1.1       # below this the logarithm of the percentage is about zero or negative


@dataclass
class RdpRow:
    name: str
    value: Optional[float]              # None when undefined
    se: Optional[float]
    pct_lo: Optional[float]             # atomic % at the low and high angle (published form)
    pct_hi: Optional[float]
    ratio: Optional[float]              # I_high / I_low
    curved_t: Optional[float] = None
    note: str = ""


def _nearest(angles: np.ndarray, usable: np.ndarray, target: float) -> int:
    return int(usable[np.argmin(np.abs(angles[usable] - target))])


def relative_depth(species: Sequence[Species], theta_lo: Optional[float] = None, theta_hi: Optional[float] = None,
                   mode: str = RDP_PUBLISHED, max_angle: float = DEFAULT_MAX_ANGLE, curved_limit: float = 3.0
                   ) -> List[RdpRow]:
    """
    Relative depth index of every species between two emission angles: `theta_lo` (default the smallest usable angle)
    and `theta_hi` (default the angle nearest 67.5 degrees within `max_angle`); each species uses its own angle nearest
    to those. `mode` is RDP_PUBLISHED or RDP_RATIO (see above). The uncertainty comes from the counting statistics when
    every species has them. A larger value means nearer the surface.
    """
    usable = {id(sp): np.flatnonzero((np.asarray(sp.angles, dtype=float) <= max_angle) & (sp.intensity > 0))
              for sp in species}
    pool = [float(sp.angles[i]) for sp in species for i in usable[id(sp)]]
    if len(set(pool)) < 2:
        return [RdpRow(sp.name, None, None, None, None, None, note="needs positive intensities at two usable angles")
                for sp in species]
    lo = theta_lo if theta_lo is not None else min(pool)
    hi = theta_hi if theta_hi is not None else min(sorted(set(pool)), key=lambda a: abs(a - 67.5))
    if hi == lo:
        hi = max(pool) if lo == min(pool) else min(pool)
    idx = {}
    for sp in species:
        u = usable[id(sp)]
        if len(u) < 2:
            idx[id(sp)] = None
            continue
        a = np.asarray(sp.angles, dtype=float)
        idx[id(sp)] = (_nearest(a, u, lo), _nearest(a, u, hi))
    have = [sp for sp in species if idx[id(sp)] is not None and idx[id(sp)][0] != idx[id(sp)][1]]
    sig_ok = bool(have) and all(sp.sigma is not None for sp in have)

    def column(which: int):
        """I/S, its relative sigma and the share (fraction) of each species that has data."""
        vals = {id(sp): float(sp.intensity[idx[id(sp)][which]]) / sp.rsf for sp in have}
        total = sum(vals.values())
        rel = {id(sp): (float(sp.sigma[idx[id(sp)][which]] / sp.intensity[idx[id(sp)][which]]) if sig_ok else None)
               for sp in have}
        return vals, rel, {k: v / total for k, v in vals.items()}

    v_lo, rel_lo, f_lo = column(0) if have else ({}, {}, {})
    v_hi, rel_hi, f_hi = column(1) if have else ({}, {}, {})

    def share_se(sp, rel, frac):
        """1-sigma of ln(share): d ln p_i / d ln I_i = 1 - p_i and d ln p_i / d ln I_j = -p_j."""
        if not sig_ok:
            return None
        var = sum(((1.0 if o is sp else 0.0) - frac[id(o)]) ** 2 * rel[id(o)] ** 2 for o in have)
        return math.sqrt(var)

    rows: List[RdpRow] = []
    for sp in species:
        if idx[id(sp)] is None or idx[id(sp)][0] == idx[id(sp)][1]:
            rows.append(RdpRow(sp.name, None, None, None, None, None, note="needs positive intensities at two usable angles"))
            continue
        ia, ib = idx[id(sp)]
        ratio = float(sp.intensity[ib] / sp.intensity[ia])
        p_lo, p_hi = 100.0 * f_lo[id(sp)], 100.0 * f_hi[id(sp)]
        t = curvature_t(sp.angles, sp.intensity, max_angle)
        notes = []
        value = se = None
        if mode == RDP_RATIO:
            value = math.log(ratio)
            if sig_ok:
                se = math.hypot(rel_lo[id(sp)], rel_hi[id(sp)])
        elif p_lo < RDP_MIN_SHARE_PCT or p_hi < RDP_MIN_SHARE_PCT:
            notes.append(f"undefined: a share of {min(p_lo, p_hi):.2g} % is too small for the ratio of logarithms "
                         "(use the log-of-ratio form)")
        else:
            value = math.log(p_hi) / math.log(p_lo)
            if sig_ok:
                se = abs(value) * math.hypot(share_se(sp, rel_hi, f_hi) / math.log(p_hi),
                                             share_se(sp, rel_lo, f_lo) / math.log(p_lo))
        if t is not None and abs(t) > curved_limit:
            notes.append("ln I is curved against sec θ: the species may sit in more than one layer (averaged here)")
        ke_far = [o.name for o in species if o is not sp and sp.ke_ev and o.ke_ev
                  and not SIMILAR_KE_RATIO[0] <= sp.ke_ev / o.ke_ev <= SIMILAR_KE_RATIO[1]]
        if ke_far:
            notes.append("kinetic energy differs a lot from " + ", ".join(ke_far[:2])
                         + ": different attenuation length, so compare with care (Stratification scales by lambda)")
        rows.append(RdpRow(sp.name, value, se, p_lo, p_hi, ratio, t, "; ".join(notes)))
    return rows
