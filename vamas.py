"""
VAMAS (ISO 14976) reader and a minimal writer, numpy only.

Only what an angle-resolved XPS tilt series needs is supported: experiment mode NORM, scan mode REGULAR and technique
XPS. The layout was worked out from Thermo K-Alpha+ exports and is the standard one:

* file header (5 text lines, comment lines, mode lines, experimental variables, four counts, number of blocks);
* per block: ids, date, comments, technique, one value per experimental variable, source, analyser, species,
  abscissa (start, step), corresponding variables, acquisition, sample-normal angles, extra parameters, then the data.

Quirks that matter (all covered by tests):

* 1e37 means "not specified" and is returned as None.
* The number of ordinate values in a block is the TOTAL count: with two corresponding variables (Counts and
  Transmission, as the K-Alpha writes) the values are interleaved, so n = 2 x the number of energy points.
* The K-Alpha stores the stage tilt in the first experimental variable, which is only labelled "Exp Variable 1"
  (unit d). Nothing in the file says it is a tilt, so the caller chooses the angle source (see ANGLE_SOURCES).
* The energy axis is Kinetic Energy; binding energy is hv - KE (the K-Alpha already refers KE to the Fermi level, so
  no work function is subtracted).
"""
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

NOT_SPECIFIED = 1e37
FORMAT_ID = "VAMAS Surface Chemical Analysis Standard Data Transfer Format 1988 May 4"

# Where an angle can come from. exp1/exp2 are the experimental variables, the rest are fixed block fields.
ANGLE_SOURCES = {
    "exp1": "Experimental variable 1",
    "exp2": "Experimental variable 2",
    "normal_polar": "Sample-normal polar angle (block field)",
    "normal_azimuth": "Sample-normal azimuth (block field)",
    "rotation": "Sample rotation (block field)",
}


class VamasError(ValueError):
    """The file is not a VAMAS file we can read."""


def _opt(text: str) -> Optional[float]:
    """A VAMAS number, with 1e37 meaning 'not specified'."""
    try:
        v = float(text)
    except ValueError:
        raise VamasError(f"Expected a number, found {text!r}.")
    return None if abs(v) >= NOT_SPECIFIED * 0.99 else v


@dataclass
class Block:
    block_id: str
    sample_id: str = ""
    species: str = ""
    transition: str = ""
    technique: str = "XPS"
    hv_ev: float = 1486.68
    be_ev: np.ndarray = field(default_factory=lambda: np.zeros(0))       # ascending binding energy
    counts: np.ndarray = field(default_factory=lambda: np.zeros(0))
    transmission: Optional[np.ndarray] = None
    dwell_s: Optional[float] = None
    scans: Optional[int] = None
    pass_energy_ev: Optional[float] = None
    exp_values: List[Optional[float]] = field(default_factory=list)
    exp_labels: List[str] = field(default_factory=list)
    exp_units: List[str] = field(default_factory=list)
    normal_polar: Optional[float] = None
    normal_azimuth: Optional[float] = None
    rotation: Optional[float] = None
    comments: List[str] = field(default_factory=list)
    timestamp: str = ""

    def value(self, source: str) -> Optional[float]:
        if source in ("exp1", "exp2"):
            i = 0 if source == "exp1" else 1
            return self.exp_values[i] if i < len(self.exp_values) else None
        if source in ("normal_polar", "normal_azimuth", "rotation"):
            return getattr(self, source)
        raise KeyError(source)

    @property
    def region(self) -> str:
        return self.block_id


@dataclass
class VamasFile:
    path: str
    blocks: List[Block]
    exp_labels: List[str] = field(default_factory=list)
    exp_units: List[str] = field(default_factory=list)

    def regions(self) -> List[str]:
        """Region (block id) names in file order, without repeats."""
        seen: Dict[str, None] = {}
        for b in self.blocks:
            seen.setdefault(b.block_id, None)
        return list(seen)

    def blocks_of(self, region: str) -> List[Block]:
        return [b for b in self.blocks if b.block_id == region]

    def sources_present(self, region: Optional[str] = None) -> List[str]:
        """Angle sources that have a value in at least one block (of `region`, or of the file)."""
        blocks = self.blocks_of(region) if region else self.blocks
        return [s for s in ANGLE_SOURCES if any(b.value(s) is not None for b in blocks)]

    def default_sources(self, region: Optional[str] = None) -> Tuple[Optional[str], Optional[str]]:
        """(angle source, azimuth source): the first source that varies is the tilt, a second, different one the azimuth."""
        blocks = self.blocks_of(region) if region else self.blocks
        order = ("exp1", "normal_polar", "exp2", "rotation", "normal_azimuth")
        angle = next((s for s in order if _distinct(blocks, s) > 1), None)
        if angle is None:
            angle = next((s for s in order if any(b.value(s) is not None for b in blocks)), None)
        azimuth = None
        for s in ("exp2", "normal_azimuth", "rotation", "exp1"):
            if s != angle and any(b.value(s) is not None for b in blocks):
                azimuth = s
                break
        return angle, azimuth

    def series(self, region: str, angle_source: Optional[str] = "exp1",
               azimuth_source: Optional[str] = None) -> Dict[Optional[float], List[Tuple[float, Block]]]:
        """
        Blocks of one region as {azimuth: [(angle, block), ...]} with each list sorted by angle. Without an azimuth source
        (or when it is unspecified) the key is None. Blocks without an angle are skipped.
        """
        out: Dict[Optional[float], List[Tuple[float, Block]]] = {}
        for b in self.blocks_of(region):
            angle = b.value(angle_source) if angle_source else None
            if angle is None:
                continue
            az = b.value(azimuth_source) if azimuth_source else None
            out.setdefault(az, []).append((angle, b))
        for rows in out.values():
            rows.sort(key=lambda r: r[0])
        return out


def _distinct(blocks: Sequence[Block], source: str) -> int:
    return len({b.value(source) for b in blocks if b.value(source) is not None})


# ---------------------------------------------------------------------------------------------------------------
# Reader
# ---------------------------------------------------------------------------------------------------------------

class _Lines:
    def __init__(self, lines: List[str]):
        self.lines, self.pos = lines, 0

    def next(self) -> str:
        if self.pos >= len(self.lines):
            raise VamasError("The file ends early (truncated VAMAS file?).")
        v = self.lines[self.pos]
        self.pos += 1
        return v.strip()

    def int(self) -> int:
        t = self.next()
        try:
            return int(float(t))
        except ValueError:
            raise VamasError(f"Expected an integer, found {t!r}.")

    def skip(self, n: int):
        for _ in range(n):
            self.next()


def read(path: str) -> VamasFile:
    """Read a VAMAS file. Raises VamasError with a readable message for anything it cannot handle."""
    try:
        with open(path, "r", encoding="latin-1") as fh:
            text = fh.read()
    except OSError as exc:
        raise VamasError(f"Cannot open the file: {exc}")
    return parse(text, path)


def parse(text: str, path: str = "") -> VamasFile:
    lines = text.splitlines()
    if not lines or not lines[0].strip().upper().startswith("VAMAS"):
        raise VamasError("This is not a VAMAS file (the first line should start with 'VAMAS Surface Chemical Analysis').")
    ln = _Lines(lines)
    ln.skip(5)                                             # format, institution, instrument, operator, experiment id
    ln.skip(ln.int())                                      # file comments
    exp_mode, scan_mode = ln.next().upper(), ln.next().upper()
    if exp_mode != "NORM":
        raise VamasError(f"Experiment mode {exp_mode} is not supported (only NORM: a list of spectra, not maps or depth "
                         "profiles).")
    if scan_mode != "REGULAR":
        raise VamasError(f"Scan mode {scan_mode} is not supported (only REGULAR, an evenly spaced energy axis).")
    ln.int()                                               # number of spectral regions
    nev = ln.int()
    exp_labels, exp_units = [], []
    for _ in range(nev):
        exp_labels.append(ln.next())
        exp_units.append(ln.next())
    ln.skip(4)                                             # inclusion list, manual items, future upgrade entries (x2)
    nblocks = ln.int()

    blocks = []
    for _ in range(nblocks):
        blocks.append(_read_block(ln, nev, exp_labels, exp_units))
    return VamasFile(path=path, blocks=blocks, exp_labels=exp_labels, exp_units=exp_units)


def _read_block(ln: _Lines, nev: int, exp_labels: List[str], exp_units: List[str]) -> Block:
    b = Block(block_id=ln.next(), sample_id=ln.next())
    y, mo, d, h, mi, s = (ln.next() for _ in range(6))
    ln.next()                                              # hours from GMT
    b.timestamp = f"{y}-{mo.zfill(2)}-{d.zfill(2)} {h.zfill(2)}:{mi.zfill(2)}:{s.zfill(2)}"
    b.comments = [ln.next() for _ in range(ln.int())]
    b.technique = ln.next()
    if b.technique.upper() != "XPS":
        raise VamasError(f"Block '{b.block_id}' uses technique {b.technique}; only XPS blocks are supported.")
    b.exp_labels, b.exp_units = list(exp_labels), list(exp_units)
    b.exp_values = [_opt(ln.next()) for _ in range(nev)]
    ln.next()                                              # analysis source label
    b.hv_ev = _opt(ln.next()) or 0.0                       # source characteristic energy
    ln.skip(5)                                             # strength, beam widths, incidence polar/azimuth
    ln.next()                                              # analyser mode
    b.pass_energy_ev = _opt(ln.next())
    ln.skip(7)                                             # magnification, work function, bias, widths, take-off angles
    b.species, b.transition = ln.next(), ln.next()
    ln.next()                                              # charge of detected particle
    abscissa = ln.next().lower()
    ln.next()                                              # abscissa unit
    start, step = float(ln.next()), float(ln.next())
    ncv = ln.int()
    names = []
    for _ in range(ncv):
        names.append(ln.next().lower())
        ln.next()
    ln.next()                                              # signal mode
    b.dwell_s = _opt(ln.next())
    scans = _opt(ln.next())
    b.scans = int(scans) if scans is not None else None
    ln.next()                                              # signal time correction
    b.normal_polar, b.normal_azimuth, b.rotation = _opt(ln.next()), _opt(ln.next()), _opt(ln.next())
    for _ in range(ln.int()):                              # additional numerical parameters (label, unit, value)
        ln.skip(3)
    n = ln.int()
    ln.skip(2 * ncv)                                       # min and max of each corresponding variable
    if ncv < 1 or n % ncv:
        raise VamasError(f"Block '{b.block_id}' has {n} values for {ncv} variables.")
    data = np.array([float(ln.next()) for _ in range(n)]).reshape(-1, ncv)
    axis = start + step * np.arange(data.shape[0])
    if "binding" in abscissa:
        be = axis
    elif "kinetic" in abscissa:
        be = b.hv_ev - axis
    else:
        raise VamasError(f"Block '{b.block_id}': unknown energy axis {abscissa!r}.")
    order = np.argsort(be)
    b.be_ev = be[order]
    columns = {name: data[order, i] for i, name in enumerate(names)}
    b.counts = columns.get("counts", data[order, 0])
    b.transmission = columns.get("transmission")
    return b


# ---------------------------------------------------------------------------------------------------------------
# Writer (round-trip tests, the docs sample file, synthetic tilt series)
# ---------------------------------------------------------------------------------------------------------------

def _num(v: Optional[float]) -> str:
    return "1e+37" if v is None else f"{v:.10g}"


def write(path: str, blocks: Sequence[Block], exp_labels: Sequence[str] = ("Tilt", "Azimuth"),
          exp_units: Sequence[str] = ("d", "d"), comment: str = "") -> None:
    """Write blocks as a minimal NORM/REGULAR XPS VAMAS file with Counts (and Transmission, if present) columns."""
    out = [FORMAT_ID, "Not Specified", "XPS Thickness Suite", "Not Specified", "Not Specified"]
    out.append("1" if comment else "0")
    if comment:
        out.append(comment)
    out += ["NORM", "REGULAR", str(len({b.block_id for b in blocks})), str(len(exp_labels))]
    for label, unit in zip(exp_labels, exp_units):
        out += [label, unit]
    out += ["0", "0", "0", "0", str(len(blocks))]
    for b in blocks:
        be = np.asarray(b.be_ev, dtype=float)
        if len(be) < 2:
            raise ValueError("A block needs at least two points.")
        step = float(be[1] - be[0])
        ke = b.hv_ev - be[::-1]                            # kinetic energy ascending, as the instrument writes it
        cols = [np.asarray(b.counts, dtype=float)[::-1]]
        names = ["Counts"]
        if b.transmission is not None:
            cols.append(np.asarray(b.transmission, dtype=float)[::-1])
            names.append("Transmission")
        ev = list(b.exp_values) + [None] * (len(exp_labels) - len(b.exp_values))
        out += [b.block_id, b.sample_id, "2025", "1", "1", "0", "0", "0", "0", str(len(b.comments))] + list(b.comments)
        out += ["XPS"] + [_num(v) for v in ev[:len(exp_labels)]]
        out += ["KA1486", _num(b.hv_ev), "1e+37", "1e+37", "1e+37", "1e+37", "1e+37"]
        out += ["FAT", _num(b.pass_energy_ev)] + ["1e+37"] * 7
        out += [b.species, b.transition, "-1", "Kinetic Energy", "eV", _num(float(ke[0])), _num(step), str(len(cols))]
        for name in names:
            out += [name, "d"]
        out += ["pulse counting", _num(b.dwell_s), _num(b.scans), "1e+37",
                _num(b.normal_polar), _num(b.normal_azimuth), _num(b.rotation), "0", str(len(be) * len(cols))]
        for _ in cols:
            out += ["0", "1"]
        for row in zip(*cols):
            out += [f"{v:.10g}" for v in row]
    out.append("end of experiment")
    with open(path, "w", encoding="latin-1", newline="") as fh:
        fh.write("\r\n".join(out) + "\r\n")
