"""
Generate docs/presets.md from the preset data in app.py, so the tables can never drift from the code.

    python docs/make_presets.py            # rewrite docs/presets.md
    python docs/make_presets.py --check    # exit 1 if docs/presets.md is out of date (tests/test_docs.py does this)

Only plain-text tables are used (no maths in cells), so the page renders identically on GitHub and in the app.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

import app as A  # noqa: E402  (importing app opens no window)

OUT = os.path.join(HERE, "presets.md")

INTRO = """# Presets and data

The Calculator tab and the Thickogram tab's preset drop-downs are filled from a reference table in `app.py`. **Every
value is an editable default to check against your own material and instrument; none is an authority.** This page lists
all of them, generated from the code, and says where each comes from.

Contents: [Provenance](#provenance), [Substrates](#substrates), [Overlayers](#overlayers),
[Band gaps](#band-gaps), [Notes on particular sets](#notes-on-particular-sets).

## Provenance

| Label | Meaning |
|---|---|
| Published | Density, IMFP and related values taken from the paper named: Strohmeier [1], Biesinger [5] or Carlson and McGuire [7]. The tests pin these presets to the published numbers. |
| Avantage density, computed IMFP | Density (and band gap) from the Thermo Avantage `Densities.txt` [11]; the IMFP is calculated with TPP-2M [6] for Al K-alpha electrons. |
| Handbook density, computed IMFP | A standard handbook density that is **not** from the supplied files; IMFP calculated with TPP-2M. |
| User-supplied density, computed IMFP | A density you supplied; IMFP calculated with TPP-2M. |
| Biesinger density, computed IMFP | Biesinger's densities with an IMFP calculated for Cu 2p3/2 electrons, which he does not publish. |
| Placeholder | A rough stand-in with no source. |

Binding energies are informational and position the simulated peaks; they do not enter the thickness calculation.
Many are typical literature values rather than values from the supplied files.
"""

NOTES = """## Notes on particular sets

**Copper: two sets of IMFPs.** Biesinger's IMFPs [5] are for the Cu L3M4,5M4,5 Auger line (kinetic energy about 918
eV), so the presets marked LMM must be used with areas from the Cu LMM fit. The presets marked 2p3/2 are for Cu 2p3/2
electrons (about 554 eV with Al K-alpha); their IMFPs are computed here, not published. Thickness scales almost
linearly with the IMFPs, so mixing the sets changes the answer by about 30 percent. The same TPP-2M code reproduces
Biesinger's three oxide IMFPs at 918 eV to 0.2 percent, but his value for the metal is about 6 percent above
TPP-2M at that energy, so the computed 2p3/2 metal value may sit slightly low.

**Aluminium.** Strohmeier's values [1] are IMFPs of 2.8 nm (oxide) and 2.6 nm (metal) for Al K-alpha, 2.4 and 2.2 nm
for Mg K-alpha, and a ratio of metal atom densities of 1.5. The Al2O3 density of 3.40 g/cm3 is back-calculated to
give that ratio; it is not a measured density (native and anodic films are nearer 3.0 to 3.2, and sapphire is 3.95).
The Avantage file lists 2.6989 and 3.97 g/cm3 for Al and Al2O3, which would give a ratio of about 1.29, so those
values are not used for the Al presets.

**Silicon family.** The Si3N4 density of 3.44 g/cm3 is the Avantage value and is high for thin amorphous films, which
are often 2.9 to 3.1; a lower density raises the density ratio and changes the thickness. SiON is only a rough
placeholder.

**Tungsten.** The W / WO3 presets follow Carlson and McGuire [7], who took the anodic oxide as 6.8 g/cm3 (bulk is
7.16) and used a 45 degree take-off.

**Extra metal and oxide pairs.** Mg / MgO, Ge / GeO2, Zn / ZnO, Zr / ZrO2, Hf / HfO2, Sn / SnO2, SiC with SiO2, GaAs and
GaP with Ga2O3, and CdTe with TeO2 were added with densities from the Avantage file where it has them and handbook
values otherwise. As a sanity check, the density ratios for the 1:1 metal oxides agree with the known
Pilling-Bedworth ratios (Mg 0.79, Zr 1.55, Hf 1.62, Sn 1.34, Zn 1.58). Choices that are judgement calls include the
density of GeO2 (3.64, amorphous; hexagonal is 4.23) and of TeO2 (6.04, tetragonal; orthorhombic is 5.67).

**Weaker pairs.** For InP the In 3d5/2 shift to In2O3 is only 0.3 to 0.5 eV, so the oxide component is hard to fit, and
the native oxide also contains InPO4, which is not modelled. PbS uses the S 2p line, where the sulfate shift is about
7.5 eV, so PbS with PbSO4 is clean but PbO is not modelled. For Zn and Mg the Zn LMM and Mg KLL Auger lines separate
metal and oxide better than the core lines.

**Not included.** ZnS (no clean oxide route), diamond and graphite (no meaningful oxide overlayer), As2O3 on GaAs, CdO
on CdTe (Cd 3d cannot resolve it) and PbO.

**Valence electrons for TPP-2M.** Valence-electron counts count s and p electrons plus unfilled d electrons, with
filled d shells treated as core (copper is kept at 11 to match Biesinger). They are listed in `TPP2M_INPUTS` in
`app.py`, together with the kinetic energy used for each computed IMFP.
"""


def provenance(key, comp):
    note = comp.note
    if "Strohmeier" in note:
        return "Published: Strohmeier [1]"
    if note.startswith("Biesinger"):
        return "Published: Biesinger [5]"
    if "Carlson" in note:
        return "Published: Carlson and McGuire [7]"
    if "COMPUTED" in note:
        if "(2p3/2)" in key:
            return "Biesinger density, computed IMFP"
        if "your message" in note:
            return "User-supplied density, computed IMFP"
        if "Avantage" in note:
            return "Avantage density, computed IMFP"
        return "Handbook density, computed IMFP"
    return "Placeholder"


def _fmt(v):
    return "" if v in (None, 0, 0.0) else f"{v:g}"


def _gap(v):
    if v is None:
        return "-"
    return "metal" if v == 0 else f"{v:g}"


def _table(presets):
    rows = ["| Preset | Core level | BE (eV) | Density (g/cm3) | Molar mass | Atoms | IMFP (nm) | Band gap (eV) | Provenance |",
            "|---|---|---|---|---|---|---|---|---|"]
    for key, c in presets.items():
        if c is None or key.startswith("Custom"):
            continue
        rows.append(f"| {key} | {c.core_level or '-'} | {_fmt(c.be_ev) or '-'} | {c.density_g_cm3:g} | "
                    f"{c.molar_mass_g_mol:g} | {c.atoms_per_formula:g} | {c.imfp_nm:g} | {_gap(c.band_gap_ev)} | "
                    f"{provenance(key, c)} |")
    return "\n".join(rows)


def build() -> str:
    gaps = ["| Material | Band gap (eV) | Source |", "|---|---|---|"]
    for name in sorted(A.BAND_GAPS_EV):
        value, source = A.BAND_GAPS_EV[name]
        gaps.append(f"| {name} | {value:g} | {source} |")
    parts = [
        INTRO,
        "## Substrates\n\nThe element or compound the overlayer sits on. Atoms is the number of atoms of the measured "
        "element per formula unit.\n\n" + _table(A.SUBSTRATES),
        "## Overlayers\n\nOxides, hydroxides and nitrides. A blank binding energy marks an Auger preset (no core-level "
        "binding energy).\n\n" + _table(A.OVERLAYERS),
        "## Band gaps\n\nBand gaps are used only by the TPP-2M IMFP estimate. Existing entries are never overwritten by "
        "later lists: the Avantage values take precedence.\n\n" + "\n".join(gaps),
        NOTES,
    ]
    return "\n\n".join(p.rstrip("\n") for p in parts) + "\n"


if __name__ == "__main__":
    text = build()
    if "--check" in sys.argv:
        current = open(OUT, encoding="utf-8").read() if os.path.exists(OUT) else ""
        sys.exit(0 if current == text else 1)
    with open(OUT, "w", encoding="utf-8", newline="") as fh:
        fh.write(text)
    print(f"wrote {OUT} ({text.count(chr(10))} lines)")
