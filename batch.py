"""
Batch Strohmeier thickness from a CasaXPS-fitted VAMAS file (stdlib only, no Tk).

One block per sample of a core-level region (e.g. "Al2p Scan"), each carrying the fit CasaXPS stored in its comment
lines (see casa.py). A chemical state of the fit is given a role (substrate / overlayer / ignore) and the Strohmeier
equation of app.strohmeier_thickness is applied to every sample in turn.

Areas are used as stored (raw counts/s x eV): the states belong to one core level, so the transmission function cancels
in the ratio. A sample that cannot be calculated (a state missing, zero area) becomes a row with an error message rather
than stopping the batch.
"""
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import casa

SUBSTRATE, OVERLAYER, IGNORE = "substrate", "overlayer", "ignore"


@dataclass
class Sample:
    index: int                        # block position in the file
    sample_id: str
    states: Dict[str, Tuple[str, float, float]] = field(default_factory=dict)   # key -> (name, area, BE)


@dataclass
class StateInfo:
    key: str
    name: str
    be: float                         # mean binding energy over the samples
    count: int                        # samples in which the state appears


@dataclass
class Row:
    sample: Sample
    i_sub: float = 0.0
    i_ov: float = 0.0
    d_nm: Optional[float] = None
    validity_limit_nm: Optional[float] = None
    within_validity: Optional[bool] = None
    error: str = ""


def collect(vf, region: str) -> Tuple[List[Sample], int]:
    """(samples, number of blocks of the region without a CasaXPS fit), in file order."""
    samples, skipped = [], 0
    for i, b in enumerate(vf.blocks):
        if b.region != region:
            continue
        fit = casa.fits_of(b)
        if fit is None or not fit.components:
            skipped += 1
            continue
        st = casa.states_of(fit, b.hv_ev)
        samples.append(Sample(i, b.sample_id, {k: (s.name, s.area, s.be) for k, s in st.items()}))
    return samples, skipped


def fitted_regions(vf) -> List[str]:
    """Regions with at least one fitted block, in file order."""
    out = []
    for b in vf.blocks:
        if b.region not in out:
            fit = casa.fits_of(b)
            if fit is not None and fit.components:
                out.append(b.region)
    return out


def fit_states(samples: Sequence[Sample]) -> List[StateInfo]:
    """The states seen across the samples, lowest binding energy first."""
    acc: Dict[str, list] = {}
    for s in samples:
        for key, (name, _area, be) in s.states.items():
            a = acc.setdefault(key, [name, 0.0, 0])
            a[1] += be
            a[2] += 1
    infos = [StateInfo(k, n, tot / cnt, cnt) for k, (n, tot, cnt) in acc.items()]
    return sorted(infos, key=lambda s: s.be)


def guess_roles(states: Sequence[StateInfo]) -> Dict[str, str]:
    """Lowest-BE state is the substrate (the metal), every other one an overlayer. All editable in the GUI."""
    ordered = sorted(states, key=lambda s: s.be)
    return {s.key: (SUBSTRATE if i == 0 else OVERLAYER) for i, s in enumerate(ordered)}


def compute(samples: Sequence[Sample], roles: Dict[str, str], substrate, overlayers: Dict[str, object],
            theta_deg: float, volume_density) -> List[Row]:
    """
    One Row per sample. `substrate` is a Compound preset; `overlayers` maps each overlayer state key to its Compound
    preset; `volume_density` is app.volume_density. Several overlayer states are intensity-weighted, as in the
    Calculator tab (oxide + hydroxide).
    """
    from app import strohmeier_thickness, OverlayerInput

    sub_keys = [k for k, r in roles.items() if r == SUBSTRATE]
    ov_keys = [k for k, r in roles.items() if r == OVERLAYER]
    if len(sub_keys) != 1:
        raise ValueError("Give exactly one state the Substrate role.")
    if not ov_keys:
        raise ValueError("Give at least one state the Overlayer role.")
    for k in ov_keys:
        if overlayers.get(k) is None:
            raise ValueError("Choose an overlayer material for every overlayer state.")
    sub_key = sub_keys[0]
    n_sub = volume_density(substrate.density_g_cm3, substrate.molar_mass_g_mol, substrate.atoms_per_formula)

    rows = []
    for s in samples:
        row = Row(s)
        rows.append(row)
        if sub_key not in s.states:
            row.error = f"no '{sub_key}' state in the fit"
            continue
        row.i_sub = s.states[sub_key][1]
        comps = []
        for k in ov_keys:
            if k not in s.states:
                continue
            m = overlayers[k]
            n = volume_density(m.density_g_cm3, m.molar_mass_g_mol, m.atoms_per_formula)
            comps.append(OverlayerInput(m.label, n, m.imfp_nm, m.rsf, s.states[k][1]))
        row.i_ov = sum(c.area for c in comps)
        try:
            res = strohmeier_thickness(theta_deg, n_sub, substrate.imfp_nm, substrate.rsf, row.i_sub, comps)
        except ValueError as exc:
            row.error = str(exc)
            continue
        row.d_nm = res["d_nm"]
        row.validity_limit_nm = res["validity_limit_nm"]
        row.within_validity = res["within_validity"]
    return rows
