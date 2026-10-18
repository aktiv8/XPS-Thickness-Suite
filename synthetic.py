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
                     step: float = 0.1, be_range=(90.0, 115.0)) -> List[VM.Block]:
    """
    Si 2p-like region: a substrate doublet and an overlayer doublet `ov_be - sub_be` higher, on a Shirley background.
    `drift_ev` moves the whole spectrum up to that much with the tilt (charging). `d_by_azimuth` gives each azimuth its
    own thickness (default: all `d_nm`). Areas are in counts per second x eV.
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
            out.append(_block(region, species, hv, be, signal + bg, tilt, az, noise, rng))
    return out


def two_peak_series(d_nm: float = 1.5, lam_nm: float = 3.0, so: float = 2.93, ss: float = 0.82,
                    be_o: float = 532.0, be_s: float = 99.3, hv: float = 1486.68, tilts: Sequence[float] = DEFAULT_TILTS,
                    fwhm_o: float = 1.7, fwhm_s: float = 0.8, is_inf: float = 40000.0, noise: bool = False,
                    seed: int = 1, azimuths: Sequence[Optional[float]] = (None,), step: float = 0.1) -> List[VM.Block]:
    """
    An overlayer peak (O 1s-like) and a substrate peak (Si 2p-like) in two regions. RSFs `so`, `ss`; B' = (KEo/KEs)^0.75.
    Windows for the analysis: the peaks sit at be_o and be_s (no drift here).
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
                out.append(_block(region, species, hv, be, sig + bg, tilt, az, noise, rng))
    return out


def write_sample(path: str, **kwargs) -> None:
    """The docs sample file: a noisy synthetic Si/SiO2 tilt series."""
    blocks = two_state_series(noise=True, lorentz=0.2, **kwargs)
    VM.write(path, blocks, comment="Synthetic Si 2p tilt series (Si + SiO2), not a measurement. Made by synthetic.py.")
