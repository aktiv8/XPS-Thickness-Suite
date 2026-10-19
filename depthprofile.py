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
