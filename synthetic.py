"""
Synthetic ARXPS tilt series (numpy only), written as vamas.Block objects. Used by the tests, to make the sample file in
docs/sample/ and to check that the analysis recovers a thickness it is given. Nothing here is a measurement.

The intensities follow the same equations the analysis inverts, so with no noise and well-separated peaks the thickness
comes back to within the effect of the peak tails that fall outside the integration windows:

    two states of one peak:   I_sub = I_inf e^-x,   I_ov = Rinf I_inf (1 - e^-x),         x = d / (lambda cos theta)
    two different peaks:      I_s = Is_inf e^(-B' x), I_o = Ro Is_inf (1 - e^-x),  B' = (KEo / KEs)^0.75, Ro = So / Ss
"""
import math
from typing import List, Optional, Sequence

import numpy as np

import spectrum_sim as sim
import thickogram as TG
import vamas as VM

DEFAULT_TILTS = (0.0, 15.0, 30.0, 45.0, 60.0, 75.0)
SPIN_ORBIT_EV = 0.6          # 2p doublet: 2p3/2 and 2p1/2 (area ratio 2:1) this far apart
DWELL_S, SCANS, TRANSMISSION = 0.05, 4, 8800.0


def _doublet(be, centre, fwhm, area, lorentz):
    """A 2p doublet with the given TOTAL area (2p3/2 at `centre`)."""
    return area * (2 / 3 * sim.pseudo_voigt(be, centre, fwhm, lorentz)
                   + 1 / 3 * sim.pseudo_voigt(be, centre + SPIN_ORBIT_EV, fwhm, lorentz))


def casa_fit_lines(region_name, ke_lo, ke_hi, region_rsf, comps, calib=None):
    """
    The comment lines CasaXPS writes for a fit, so the tests and the docs sample can carry one. `comps` are dicts with
    name, area (counts/s x eV), ke (as stored), fwhm, index, group and rsf. `calib` = (measured BE, assigned BE) adds
    the charge-correction line (positions then stay in the frame they were written in).
    """
    lines = ["Casa Info Follows CasaXPS Version synthetic", "0", "0", "1"]
    if calib is not None:
        lines.append(f"Calib M = {calib[0]:.6g} A = {calib[1]:.6g} BE ADD")
    lines.append(f"CASA region (*{region_name}*) (*Shirley*) {ke_lo:.6g} {ke_hi:.6g} {region_rsf:g} 2 0 0 0 -450 0 0 "
                 f"(*{region_name}*) 28.0858 0 {region_rsf:g}")
    lines.append(str(len(comps)))
    for c in comps:
        lines.append(f"CASA comp (*{c['name']}*) (*LA(50)*) Area {c['area']:.8g} 1e-20 10000000 -1 1 MFWHM {c['fwhm']:.6g} "
                     f"0.2 6 -1 1 Position {c['ke']:.8g} 0 0 -1 1 RSF {c['rsf']:g} MASS 28.0858 INDEX {c['index']} "
                     f"(*{c['group']}*) CONST (**) UNCORR")
    return lines


def _doublet_comps(name, group, index, hv, be, fwhm, area, rsf):
    """The two spin-orbit components of a state, linked by INDEX (the partner stores RSF 0, as CasaXPS does)."""
    return [dict(name=f"{name} 2p3/2", area=area * 2 / 3, ke=hv - be, fwhm=fwhm, index=index, group=group, rsf=rsf),
            dict(name=f"{name} 2p1/2", area=area / 3, ke=hv - (be + SPIN_ORBIT_EV), fwhm=fwhm, index=index, group=group,
                 rsf=0)]


def _block(region, species, hv, be, rate, tilt, azimuth, noise, rng) -> VM.Block:
    counts = np.maximum(rate, 0.0) * SCANS * DWELL_S
    if noise:
        counts = rng.poisson(counts).astype(float)
    return VM.Block(block_id=region, sample_id="synthetic", species=species, technique="XPS", hv_ev=hv, be_ev=be,
                    counts=counts, transmission=np.full_like(be, TRANSMISSION), dwell_s=DWELL_S, scans=SCANS,
                    pass_energy_ev=50.0, exp_values=[tilt, azimuth], exp_labels=["Tilt", "Azimuth"],
                    exp_units=["d", "d"], comments=["Synthetic data: not a measurement"])


def two_state_series(d_nm: float = 1.5, lam_nm: float = 3.89, rinf: float = 0.5547, tilts: Sequence[float] = DEFAULT_TILTS,
                     azimuths: Sequence[Optional[float]] = (None,), sub_be: float = 98.9, ov_be: float = 102.9,
                     hv: float = 1486.68, fwhm_sub: float = 0.75, fwhm_ov: float = 1.6, lorentz: float = 0.0,
                     i_inf: float = 60000.0, drift_ev: float = 0.3, noise: bool = False, seed: int = 1,
                     region: str = "Si2p Scan", species: str = "Si 2p", d_by_azimuth: Optional[dict] = None,
                     step: float = 0.1, be_range=(90.0, 115.0), fitted: bool = False,
                     calib: float = 0.0) -> List[VM.Block]:
    """
    Si 2p-like region: a substrate doublet and an overlayer doublet `ov_be - sub_be` higher, on a Shirley background.
    `drift_ev` moves the whole spectrum up to that much with the tilt (charging). `d_by_azimuth` gives each azimuth its
    own thickness (default: all `d_nm`). Areas are in counts per second x eV. With `fitted=True` each block also carries
    a CasaXPS-style fit of the two states (areas exact, doublet partners linked by INDEX); `calib` adds a charge
    correction line (the displayed BE is then the measured BE plus `calib`).
    """
    rng = np.random.default_rng(seed)
    be = np.arange(be_range[0], be_range[1] + step / 2, step)
    out = []
    for az in azimuths:
        d = (d_by_azimuth or {}).get(az, d_nm)
        for tilt in tilts:
            x = d / (lam_nm * math.cos(math.radians(tilt)))
            area_sub = i_inf * math.exp(-x)
            area_ov = rinf * i_inf * (1.0 - math.exp(-x))
            shift = drift_ev * tilt / max(max(tilts), 1.0)
            signal = (_doublet(be, sub_be + shift, fwhm_sub, area_sub, lorentz)
                      + _doublet(be, ov_be + shift, fwhm_ov, area_ov, lorentz))
            bg = sim.shirley_background(signal, step=0.06 * (area_sub + area_ov), base=0.02 * i_inf / 10.0)
            blk = _block(region, species, hv, be, signal + bg, tilt, az, noise, rng)
            if fitted:
                shown = calib
                comps = (_doublet_comps("Si El", "Si El", 0, hv, sub_be + shift + shown, fwhm_sub, area_sub, 0.817)
                         + _doublet_comps("Si Ox", "Si Ox", 1, hv, ov_be + shift + shown, fwhm_ov, area_ov, 0.817))
                blk.comments += casa_fit_lines("Si 2p", hv - be_range[1], hv - be_range[0], 0.817, comps,
                                               (sub_be + shift, sub_be + shift + calib) if calib else None)
            out.append(blk)
    return out


def two_peak_series(d_nm: float = 1.5, lam_nm: float = 3.0, so: float = 2.93, ss: float = 0.82,
                    be_o: float = 532.0, be_s: float = 99.3, hv: float = 1486.68, tilts: Sequence[float] = DEFAULT_TILTS,
                    fwhm_o: float = 1.7, fwhm_s: float = 0.8, is_inf: float = 40000.0, noise: bool = False,
                    seed: int = 1, azimuths: Sequence[Optional[float]] = (None,), step: float = 0.1,
                    fitted: bool = False) -> List[VM.Block]:
    """
    An overlayer peak (O 1s-like) and a substrate peak (Si 2p-like) in two regions. RSFs `so`, `ss`; B' = (KEo/KEs)^0.75.
    Windows for the analysis: the peaks sit at be_o and be_s (no drift here). With `fitted=True` each block carries a
    single-component CasaXPS-style fit with the exact area.
    """
    rng = np.random.default_rng(seed)
    b_prime = ((hv - be_o) / (hv - be_s)) ** 0.75
    ro = so / ss
    out = []
    for az in azimuths:
        for tilt in tilts:
            x = d_nm / (lam_nm * math.cos(math.radians(tilt)))
            area_s = is_inf * math.exp(-b_prime * x)
            area_o = ro * is_inf * (1.0 - math.exp(-x))
            for region, species, centre, fwhm, area in (("Si2p Scan", "Si 2p", be_s, fwhm_s, area_s),
                                                        ("O1s Scan", "O 1s", be_o, fwhm_o, area_o)):
                be = np.arange(centre - 12.0, centre + 12.0 + step / 2, step)
                sig = area * sim.pseudo_voigt(be, centre, fwhm, 0.0)
                bg = sim.shirley_background(sig, step=0.05 * area, base=0.02 * is_inf / 10.0)
                blk = _block(region, species, hv, be, sig + bg, tilt, az, noise, rng)
                if fitted:
                    comp = dict(name=species, area=area, ke=hv - centre, fwhm=fwhm, index=0, group=species,
                                rsf=so if species == "O 1s" else ss)
                    blk.comments += casa_fit_lines(species, hv - be[-1], hv - be[0], comp["rsf"], [comp])
                out.append(blk)
    return out


def layered_series(d_nm: Sequence[float], species: Sequence[dict], tilts: Sequence[float] = DEFAULT_TILTS,
                   hv: float = 1486.68, noise: bool = False, seed: int = 1, scale: float = 1000.0,
                   step: float = 0.1, fwhm: float = 1.0) -> List[VM.Block]:
    """
    A layer stack on a substrate, one region (a single Gaussian peak) per species, with exactly the intensity the Gries
    layer model predicts. `species` are dicts with name, be (eV), layer (0-based from the top, -1 for the substrate) and
    optionally lam (nm, default 3). Intensities are in counts per second x eV and `scale` sets the signal level.
    """
    import depthlayers as DL
    rng = np.random.default_rng(seed)
    out = []
    for tilt in tilts:
        for spec in species:
            lam = spec.get("lam", 3.0)
            area = scale * float(DL.layer_intensity([tilt], list(d_nm), spec["layer"], [lam] * len(d_nm))[0])
            be = np.arange(spec["be"] - 8.0, spec["be"] + 8.0 + step / 2, step)
            sig = area * sim.pseudo_voigt(be, spec["be"], fwhm, 0.0)
            bg = sim.shirley_background(sig, step=0.04 * area, base=0.01 * scale / 10.0)
            blk = _block(spec["name"] + " Scan", spec["name"], hv, be, sig + bg, tilt, None, noise, rng)
            out.append(blk)
    return out


def write_sample(path: str, fitted: bool = False, **kwargs) -> None:
    """The docs sample files: a noisy synthetic Si/SiO2 tilt series, optionally with a CasaXPS-style fit stored in it."""
    blocks = two_state_series(noise=True, lorentz=0.2, fitted=fitted, **kwargs)
    what = "with a CasaXPS-style fit" if fitted else "without a fit"
    VM.write(path, blocks, comment=f"Synthetic Si 2p tilt series (Si + SiO2) {what}, not a measurement. "
                                   "Made by synthetic.py.")
