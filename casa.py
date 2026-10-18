"""
CasaXPS peak fits stored in VAMAS block comments (stdlib only, no Tk, no numpy).

When CasaXPS saves a fitted spectrum as VAMAS it writes the fit into each block's comment lines:

    Casa Info Follows CasaXPS Version 2.3.27PR7.6
    ...
    CASA region (*Si 2p*) (*Shirley*) 1379.1765 1392.2142 0.817 2 0 0 0 -450 0 0 (*Si 2p*) 28.0858 0 0.817
    2
    CASA comp (*Si 2p El*) (*LA(50)*) Area 73449.445 1e-20 748475.55 -1 1 MFWHM 1.0834043 ... Position 1387.8626 ...
        RSF 0.817 MASS 28.0858 INDEX 0 (*Si 2p*) CONST (**) UNCORR

Areas are counts/s x eV as stored (raw: not transmission-corrected) and are used as they are: nothing is re-integrated
and no line shape is evaluated. Positions are kinetic energies in the frame stored in the file; an optional
`Calib M = <measured> A = <assigned> ...` line carries the charge correction CasaXPS applied.

This module is a port of the comment parser of the author's eXPoSe SpectraDeck project (casafit.parse, calibration,
component_be and quant.state_of; https://github.com/aktiv8/eXPoSe-Spectra-Deck), trimmed to what the ARXPS tab needs.
Its line-shape and background code is deliberately not included. The behaviour that matters:

* a **state** is a set of components sharing an INDEX >= 0 (named by the group tag, unless that is only the label of
  the region), otherwise a component on its own; doublet partners therefore add up;
* component RSFs are never used (linked doublet partners often store RSF 0).
"""
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

_TAG = re.compile(r"\(\*(.*?)\*\)")
_NUM = re.compile(r"^-?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?$")
_CALIB = re.compile(r"^Calib\s+M\s*=\s*([-\d.eE+]+)\s+A\s*=\s*([-\d.eE+]+)(.*)$")
_REGION = re.compile(r"^CASA region \(\*(.*?)\*\) \(\*(.*?)\*\)\s*(.*)$")
_COMP = re.compile(r"^CASA comp \(\*(.*?)\*\) \(\*(.*?)\*\)\s*(.*)$")
_INDEX = re.compile(r"\bINDEX\s+(-?\d+)\s*(?:\(\*(.*?)\*\))?")


@dataclass
class FitRegion:
    name: str
    background: str = "Shirley"
    start_ke: float = 0.0
    end_ke: float = 0.0
    rsf: float = 1.0
    line: str = ""


@dataclass
class FitComponent:
    name: str
    shape: str = ""
    area: float = 0.0                # counts/s x eV, as stored
    fwhm: float = 1.0
    pos_ke: float = 0.0              # KE as stored (see Fit.shift_of_comps)
    rsf: float = 1.0
    index: int = -1                  # chemical-state group (>= 0) shared by linked components
    group: str = ""                  # the group's tag
    region: str = ""
    line: str = ""


@dataclass
class Fit:
    regions: List[FitRegion] = field(default_factory=list)
    components: List[FitComponent] = field(default_factory=list)
    calib_shift: float = 0.0         # eV added to BE by the Calib lines
    region_shift: Optional[float] = None    # part of it the stored region limits / component positions
    comp_shift: Optional[float] = None      # do NOT already include (None: all of calib_shift)

    def shift_of_comps(self) -> float:
        """eV to add to a stored component position to get the raw kinetic energy."""
        return self.calib_shift if self.comp_shift is None else self.comp_shift

    def region_components(self, name: str) -> List[FitComponent]:
        """Components of a region (all of them when the file has a single region)."""
        mine = [c for c in self.components if c.region == name]
        if mine or len(self.regions) != 1:
            return mine
        return list(self.components)


@dataclass
class CompArea:
    name: str
    area: float
    be: float                        # eV, in the calibrated frame CasaXPS shows
    ke_raw: float                    # eV, on the spectrum's own kinetic-energy axis
    fwhm: float


@dataclass
class StateArea:
    """One chemical state of one block: its components (doublet partners add up)."""
    key: str
    name: str
    comps: List[CompArea]

    @property
    def area(self) -> float:
        return sum(c.area for c in self.comps)

    @property
    def be(self) -> float:
        """Area-weighted binding energy of the state."""
        total = self.area
        if total <= 0:
            return sum(c.be for c in self.comps) / len(self.comps)
        return sum(c.area * c.be for c in self.comps) / total


def _floats(text: str) -> List[float]:
    return [float(t) for t in text.split() if _NUM.match(t)]


def _region(line: str) -> Optional[FitRegion]:
    m = _REGION.match(line)
    if not m:
        return None
    body = m.group(3)
    nums = _floats(_TAG.sub(" ", body))
    if len(nums) < 2:
        return None
    return FitRegion(name=m.group(1), background=m.group(2), start_ke=min(nums[0], nums[1]),
                     end_ke=max(nums[0], nums[1]), rsf=nums[2] if len(nums) > 2 else 1.0, line=line)


def _comp(line: str, region_name: str) -> Optional[FitComponent]:
    m = _COMP.match(line)
    if not m:
        return None
    rest = m.group(3)

    def num(key, default=None):
        r = re.search(r"\b" + key + r"\s+(-?[\d.]+(?:[eE][-+]?\d+)?)", rest)
        return float(r.group(1)) if r else default
    area, fwhm, pos = num("Area"), num("MFWHM"), num("Position")
    if area is None or pos is None:
        return None
    idx = _INDEX.search(rest)
    return FitComponent(name=m.group(1), shape=m.group(2), area=area, fwhm=fwhm if fwhm else 1.0, pos_ke=pos,
                        rsf=num("RSF", 1.0), index=int(idx.group(1)) if idx else -1,
                        group=(idx.group(2) or "") if idx else "", region=region_name, line=line)


def _calibration(lines: Sequence[str]) -> Optional[Tuple[float, float, float]]:
    """(total shift, region shift, component shift) from the Calib lines of a block, or None."""
    total = reg = comp = 0.0
    found = False
    for raw in lines or ():
        cm = _CALIB.match(str(raw).strip())
        if not cm:
            continue
        found = True
        d = float(cm.group(2)) - float(cm.group(1))
        flags = {w.upper() for w in cm.group(3).split()}
        total += d
        reg += 0.0 if "REGIONS" in flags else d
        comp += 0.0 if "COMPS" in flags else d
    return (total, reg, comp) if found else None


def _fix_region_attribution(fit: Fit) -> None:
    """When a file lists several region lines before their comps, give each component to the one other region whose
    KE window contains it (only when exactly one qualifies)."""
    by_name = {r.name: r for r in fit.regions}
    for comp in fit.components:
        cur = by_name.get(comp.region)
        if cur is not None and cur.start_ke <= comp.pos_ke <= cur.end_ke:
            continue
        hits = [r for r in fit.regions if r.start_ke <= comp.pos_ke <= r.end_ke]
        if len(hits) == 1:
            comp.region = hits[0].name


def parse(comment_lines: Sequence[str]) -> Optional[Fit]:
    """The fit held in a block's comment lines, or None when it has none."""
    fit = Fit()
    cal = _calibration(comment_lines)
    if cal is not None:
        fit.calib_shift, fit.region_shift, fit.comp_shift = cal
    current = ""
    for raw in comment_lines or ():
        line = str(raw).strip()
        if line.startswith("CASA region"):
            reg = _region(line)
            if reg:
                fit.regions.append(reg)
                current = reg.name
        elif line.startswith("CASA comp"):
            comp = _comp(line, current)
            if comp:
                fit.components.append(comp)
    if len(fit.regions) > 1:
        _fix_region_attribution(fit)
    return None if not (fit.regions or fit.components) else fit


def fits_of(block) -> Optional[Fit]:
    """The CasaXPS fit of a vamas.Block, or None."""
    return parse(block.comments)


def component_be(comp: FitComponent, hv: float, fit: Fit) -> float:
    """Binding energy of a component in the calibrated frame CasaXPS shows."""
    return hv - comp.pos_ke - fit.shift_of_comps() + fit.calib_shift


def state_key(comp: FitComponent, region_name: str = "") -> Tuple[str, str]:
    """(key, name) of the chemical state a component belongs to (see the module docstring)."""
    if comp.index >= 0:
        grp = comp.group.strip()
        named = grp and grp.lower() != (region_name or "").strip().lower()
        return f"i{comp.index}", (grp if named else comp.name)
    return f"n{comp.name}", comp.name


def states_of(fit: Fit, hv: float, region: Optional[str] = None) -> Dict[str, StateArea]:
    """
    The chemical states of a fit, in the order they first appear: {key: StateArea}. `region` (a CASA region name)
    restricts them to one region of the block; by default every component is used.
    """
    comps = fit.region_components(region) if region else list(fit.components)
    out: Dict[str, StateArea] = {}
    for c in comps:
        key, name = state_key(c, c.region)
        st = out.setdefault(key, StateArea(key, name, []))
        st.comps.append(CompArea(c.name, c.area, component_be(c, hv, fit), c.pos_ke + fit.shift_of_comps(), c.fwhm))
    return out
