"""
XPS Topofactors: correcting a Thickogram thickness for particle / fibre morphology.

Source: A. G. Shard, J. Wang and S. J. Spencer, "XPS topofactors: determining overlayer thickness on particles and
fibres", Surf. Interface Anal. 41, 541-548 (2009), DOI 10.1002/sia.3044.

The Thickogram (thickogram.py) assumes a flat sample. For a sample of known geometry, analysed with the mean surface
normal directed at the analyser (theta = 0), the planar ("apparent") thickness t from the Thickogram is multiplied by
the Topofactor TF = T/t to give the true conformal overlayer thickness T = TF * t.

Implemented from the paper:
  * thin overlayers (t < 0.1 L): the limiting values 2/pi (cylinder), 0.50 (sphere, X-rays reach the whole surface)
    and 0.57 (sphere with X-ray shadowing, source-analyser angle 45-65 deg)
  * 0.1 < t < 3 (t in units of L): the empirical Eq. (6), TF = c0 + (c1 + c2 ln B) t^k + c3 exp(-c4 B), accurate to
    ~1 % against the paper's numerical calculations, with the shape parameters of its Table 1
  * the paper's constant "useful" values (0.79 cylinder, 0.67 sphere, 0.71 shadowed sphere; within ~10 %)

NOT implemented: the paper's nanoscale model (its Eqs. 7-9, for particles/fibres whose radius is comparable to the
electron attenuation length). Topofactors are for objects much larger than L.

Units and conventions (easy to get wrong):
  * t is the planar thickness in units of L, the attenuation length of the overlayer electrons in the overlayer:
    t_over_l = t_nm / L.
  * The paper's B is the ratio of attenuation lengths overlayer : substrate electrons, B = (Eo/Es)^0.75 (Cumpson).
    thickogram.py's `b` is the plain kinetic-energy ratio Eo/Es (it applies ^0.75 internally), so the paper's B is
    `b ** 0.75` - use b_paper().
  * Table 1's caption in the paper says "Eqn (7)"; the parameters belong to Eqn (6).
  * The paper quotes 0.67 for spheres as the value at t = 1.5 in the text but at T = 1.5 in its Table 2; Eq. (6) gives
    0.66 at t = 1.5, B = 1.

Pure Python (math only): no Tk, no numpy.
"""
import math
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

T_THIN = 0.1          # below this (units of L) the thin-film limit is used
T_MAX = 3.0           # Eq. (6) is stated for t up to 3
B_USUAL = (0.7, 1.4)  # "in general, B lies between 0.7 and 1.4"

# Table 1 of the paper (parameters of Eq. 6): c0, c1, c2, c3, c4, k
TABLE_1 = {
    "cylinder": (0.46, 0.25, 0.09, 0.16, 1.1, 0.22),
    "sphere": (0.40, 0.20, 0.09, 0.14, 1.5, 0.33),
    "sphere_shadowed": (0.46, 0.19, 0.08, 0.13, 1.4, 0.30),
}

MODES = ("eq6", "constant")
MODE_LABELS = {"eq6": "Eq. (6): depends on thickness and B", "constant": "Constant average (within ~10 %)"}

NANOSCALE_NOTE = ("Topofactors are for objects much larger than the attenuation length L (micro- and macroscopic "
                  "fibres and particles). They are not valid for nanoparticles or nanofibres.")


@dataclass(frozen=True)
class Shape:
    key: str
    label: str
    params: Optional[Tuple[float, ...]]    # Table 1 row, None for a flat sample
    thin_limit: float                      # t -> 0 value
    constant: float                        # the paper's average value
    hint: str


SHAPES = {
    "flat": Shape("flat", "Flat (planar) - no correction", None, 1.0, 1.0,
                  "Flat sample: Topofactor = 1."),
    "cylinder": Shape("cylinder", "Cylinder (fibre)", TABLE_1["cylinder"], 2.0 / math.pi, 0.79,
                      "Fibre axis perpendicular to the plane of source, sample and analyser."),
    "sphere": Shape("sphere", "Sphere (particle), X-rays reach all of it", TABLE_1["sphere"], 0.50, 0.67,
                    "Source collinear with the analyser, so no X-ray shadowing."),
    "sphere_shadowed": Shape("sphere_shadowed", "Sphere (particle), X-ray shadowed", TABLE_1["sphere_shadowed"], 0.57,
                             0.71, "Source-analyser angle about 45-65 deg: half of the sphere is in shadow."),
}
SHAPE_KEYS = tuple(SHAPES)


def b_paper(b_energy_ratio: float) -> float:
    """The paper's B = (Eo/Es)^0.75 from thickogram.py's plain kinetic-energy ratio `b` = Eo/Es."""
    if not math.isfinite(b_energy_ratio) or b_energy_ratio <= 0:
        raise ValueError(f"The kinetic-energy ratio must be positive (got {b_energy_ratio!r}).")
    return b_energy_ratio ** 0.75


def eq6(shape_key: str, t_over_l: float, b: float) -> float:
    """Eq. (6): TF = c0 + (c1 + c2 ln B) t^k + c3 exp(-c4 B), with the paper's B and t in units of L."""
    params = SHAPES[shape_key].params
    if params is None:
        return 1.0
    c0, c1, c2, c3, c4, k = params
    return c0 + (c1 + c2 * math.log(b)) * t_over_l ** k + c3 * math.exp(-c4 * b)


@dataclass
class TopofactorResult:
    tf: float
    regime: str                    # planar / constant / thin limit / Eq. 6 / Eq. 6 (held at t = 3)
    warnings: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)


def topofactor(shape_key: str, t_over_l: float, b: float, mode: str = "eq6",
               theta_normal_deg: float = 0.0) -> TopofactorResult:
    """
    The Topofactor TF = T/t for a shape, planar thickness t (units of L), the paper's B and the analysis angle.

    mode "eq6": the thin-film limit below t = 0.1, Eq. (6) for 0.1-3, and Eq. (6) held at its t = 3 value above that
    (with a warning). mode "constant": the paper's single average value for the shape.
    """
    if shape_key not in SHAPES:
        raise ValueError(f"Unknown shape {shape_key!r}; choose from {', '.join(SHAPES)}.")
    if mode not in MODES:
        raise ValueError(f"Unknown mode {mode!r}; choose from {', '.join(MODES)}.")
    if not math.isfinite(t_over_l) or t_over_l < 0:
        raise ValueError(f"The planar thickness in units of L must be >= 0 (got {t_over_l!r}).")
    if not math.isfinite(b) or b <= 0:
        raise ValueError(f"B must be positive (got {b!r}).")
    shape = SHAPES[shape_key]
    if shape_key == "flat":
        return TopofactorResult(1.0, "planar")

    warnings: List[str] = []
    notes = [NANOSCALE_NOTE]
    if abs(theta_normal_deg) > 1e-9:
        warnings.append(f"Topofactors assume the sample's mean normal faces the analyser (theta = 0); at theta = "
                        f"{theta_normal_deg:g} deg the value is only approximate. Re-measure at theta = 0 if you can.")
    if mode == "constant":
        return TopofactorResult(shape.constant, "constant", warnings, notes)

    if not (B_USUAL[0] <= b <= B_USUAL[1]):
        warnings.append(f"B = {b:.2f} is outside the usual range {B_USUAL[0]}-{B_USUAL[1]}; Eq. (6) is best when B is "
                        "close to 1.")
    if t_over_l < T_THIN:
        return TopofactorResult(shape.thin_limit, "thin limit", warnings, notes)
    if t_over_l <= T_MAX:
        return TopofactorResult(eq6(shape_key, t_over_l, b), "Eq. 6", warnings, notes)
    warnings.append(f"t/L = {t_over_l:.2f} is above 3, beyond the range of Eq. (6): the factor is held at its t = 3 "
                    "value, and a layer this thick is poorly determined by XPS anyway.")
    return TopofactorResult(eq6(shape_key, T_MAX, b), "Eq. 6 (held at t = 3)", warnings, notes)


@dataclass
class CorrectedResult:
    shape: str
    mode: str
    tf: float
    planar_t_nm: float              # the Thickogram (flat-sample) thickness t
    t_nm: float                     # the corrected, true conformal thickness T = TF * t
    t_over_l: float
    b_paper: float
    regime: str
    warnings: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)


def apply_topofactor(result, shape_key: str, mode: str = "eq6") -> CorrectedResult:
    """Apply the Topofactor to a thickogram.ThickogramResult: T = TF * t (t from the Thickogram, planar assumption)."""
    t_over_l = result.t_nm / result.lam_nm
    bp = b_paper(result.b)
    tfr = topofactor(shape_key, t_over_l, bp, mode, result.theta_normal_deg)
    return CorrectedResult(shape=shape_key, mode=mode, tf=tfr.tf, planar_t_nm=result.t_nm, t_nm=tfr.tf * result.t_nm,
                           t_over_l=t_over_l, b_paper=bp, regime=tfr.regime, warnings=tfr.warnings, notes=tfr.notes)
