# Presets and data

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

## Substrates

The element or compound the overlayer sits on. Atoms is the number of atoms of the measured element per formula unit.

| Preset | Core level | BE (eV) | Density (g/cm3) | Molar mass | Atoms | IMFP (nm) | Band gap (eV) | Provenance |
|---|---|---|---|---|---|---|---|---|
| Cu (metal, LMM) | Cu LMM | - | 8.96 | 63.546 | 1 | 1.54 | metal | Published: Biesinger [5] |
| Cu (metal, 2p3/2) | Cu 2p3/2 | 932.63 | 8.96 | 63.546 | 1 | 1.01 | metal | Handbook density, computed IMFP |
| W (metal) | W 4f7/2 | 31.4 | 19.3 | 183.84 | 1 | 1.28 | metal | Published: Carlson and McGuire [7] |
| Al (metal) | Al 2p | 72.9 | 2.7 | 26.98 | 1 | 2.6 | metal | Published: Strohmeier [1] |
| Al (metal, Mg Ka) | Al 2p | 72.9 | 2.7 | 26.98 | 1 | 2.2 | metal | Published: Strohmeier [1] |
| Si (metal/elemental) | Si 2p3/2 | 99.3 | 2.33 | 28.09 | 1 | 3.09 | 1.12 | Avantage density, computed IMFP |
| Ti (metal) | Ti 2p3/2 | 453.9 | 4.54 | 47.867 | 1 | 2.19 | metal | Avantage density, computed IMFP |
| Mg (metal) | Mg 2p | 49.7 | 1.7 | 24.305 | 1 | 3.85 | metal | Avantage density, computed IMFP |
| Ge (elemental) | Ge 3d | 29.4 | 5.323 | 72.63 | 1 | 2.99 | 0.67 | Avantage density, computed IMFP |
| Zn (metal, LMM) | Zn LMM | - | 7.14 | 65.38 | 1 | 2.39 | metal | Handbook density, computed IMFP |
| Zr (metal) | Zr 3d5/2 | 178.8 | 6.52 | 91.224 | 1 | 2.72 | metal | Handbook density, computed IMFP |
| Hf (metal) | Hf 4f7/2 | 14.2 | 13.31 | 178.49 | 1 | 2.79 | metal | Handbook density, computed IMFP |
| Sn (metal) | Sn 3d5/2 | 484.9 | 7.31 | 118.71 | 1 | 2.31 | metal | Avantage density, computed IMFP |
| SiC | Si 2p | 100.4 | 3.21 | 40.097 | 1 | 2.9 | 2.86 | Avantage density, computed IMFP |
| GaAs (Ga 3d) | Ga 3d | 19.3 | 5.31 | 144.645 | 1 | 3.02 | 1.35 | Avantage density, computed IMFP |
| GaP (Ga 3d) | Ga 3d | 19.2 | 4.14 | 100.697 | 1 | 3.04 | 2.25 | Handbook density, computed IMFP |
| CdTe (Te 3d5/2) | Te 3d5/2 | 572.1 | 5.85 | 240.014 | 1 | 2.37 | 1.44 | Handbook density, computed IMFP |
| InP (In 3d5/2) | In 3d5/2 | 444.6 | 4.81 | 145.792 | 1 | 2.4 | 1.27 | User-supplied density, computed IMFP |
| PbS (S 2p) | S 2p3/2 | 161 | 7.6 | 239.26 | 1 | 2.63 | 0.37 | User-supplied density, computed IMFP |

## Overlayers

Oxides, hydroxides and nitrides. A blank binding energy marks an Auger preset (no core-level binding energy).

| Preset | Core level | BE (eV) | Density (g/cm3) | Molar mass | Atoms | IMFP (nm) | Band gap (eV) | Provenance |
|---|---|---|---|---|---|---|---|---|
| CuO (LMM) | Cu LMM | - | 6.31 | 79.55 | 1 | 1.76 | 1.25 | Published: Biesinger [5] |
| Cu2O (LMM) | Cu LMM | - | 6 | 143.09 | 2 | 1.77 | 2.1 | Published: Biesinger [5] |
| Cu(OH)2 (LMM) | Cu LMM | - | 3.37 | 97.56 | 1 | 2.19 | 1.85 | Published: Biesinger [5] |
| CuO (2p3/2) | Cu 2p3/2 | 933.76 | 6.31 | 79.55 | 1 | 1.22 | 1.25 | Biesinger density, computed IMFP |
| Cu2O (2p3/2) | Cu 2p3/2 | 932.18 | 6 | 143.09 | 2 | 1.23 | 2.1 | Biesinger density, computed IMFP |
| Cu(OH)2 (2p3/2) | Cu 2p3/2 | 934.67 | 3.37 | 97.56 | 1 | 1.51 | 1.85 | Biesinger density, computed IMFP |
| WO3 | W 4f7/2 | 35.8 | 6.8 | 231.84 | 1 | 2.63 | - | Published: Carlson and McGuire [7] |
| Al2O3 | Al 2p | 74.4 | 3.4 | 101.96 | 2 | 2.8 | 8.7 | Published: Strohmeier [1] |
| Al2O3 (Mg Ka) | Al 2p | 74.4 | 3.4 | 101.96 | 2 | 2.4 | 8.7 | Published: Strohmeier [1] |
| SiO2 | Si 2p3/2 | 103.4 | 2.196 | 60.08 | 1 | 3.89 | 8.9 | Avantage density, computed IMFP |
| Si3N4 | Si 2p3/2 | 101.8 | 3.44 | 140.28 | 3 | 3 | 4.65 | Avantage density, computed IMFP |
| TiO2 | Ti 2p3/2 | 458.8 | 4.26 | 79.866 | 1 | 2.17 | 3.05 | Avantage density, computed IMFP |
| MgO | Mg 2p | 50.8 | 3.58 | 40.304 | 1 | 3.22 | 7.8 | Handbook density, computed IMFP |
| GeO2 | Ge 3d | 32.9 | 3.64 | 104.628 | 1 | 3.03 | 5.6 | Handbook density, computed IMFP |
| ZnO (LMM) | Zn LMM | - | 5.61 | 81.379 | 1 | 1.95 | 3.2 | Avantage density, computed IMFP |
| ZrO2 | Zr 3d5/2 | 182.3 | 5.68 | 123.222 | 1 | 2.42 | 6 | Avantage density, computed IMFP |
| HfO2 | Hf 4f7/2 | 16.8 | 9.68 | 210.488 | 1 | 2.26 | 5.7 | Avantage density, computed IMFP |
| SnO2 | Sn 3d5/2 | 486.6 | 6.95 | 150.708 | 1 | 1.8 | 3.6 | Handbook density, computed IMFP |
| Ga2O3 | Ga 3d | 20.7 | 5.88 | 187.443 | 2 | 2.56 | 4.8 | Handbook density, computed IMFP |
| TeO2 | Te 3d5/2 | 576 | 6.04 | 159.598 | 1 | 1.76 | 3.7 | Handbook density, computed IMFP |
| In2O3 | In 3d5/2 | 444.9 | 7.18 | 277.633 | 2 | 1.88 | 2.9 | Handbook density, computed IMFP |
| PbSO4 | S 2p3/2 | 168.6 | 6.29 | 303.256 | 1 | 2.32 | 4.5 | Handbook density, computed IMFP |
| SiON (oxynitride, approx.) | Si 2p3/2 | 102.3 | 2.5 | 87.09 | 1 | 3.3 | - | Placeholder |

## Band gaps

Band gaps are used only by the TPP-2M IMFP estimate. Existing entries are never overwritten by later lists: the Avantage values take precedence.

| Material | Band gap (eV) | Source |
|---|---|---|
| Al2O3 | 8.7 | Thermo Avantage Densities.txt |
| CdTe | 1.44 | user-supplied list (unsourced) |
| Cu(OH)2 | 1.85 | Biesinger, Surf. Interface Anal. 2017, 49, 1325 |
| Cu2O | 2.1 | Biesinger, Surf. Interface Anal. 2017, 49, 1325 |
| CuO | 1.25 | Biesinger, Surf. Interface Anal. 2017, 49, 1325 |
| Diamond | 5.5 | Thermo Avantage Densities.txt |
| Ga2O3 | 4.8 | standard literature value, not from the user's files |
| GaAs | 1.35 | Thermo Avantage Densities.txt |
| GaP | 2.25 | user-supplied list (unsourced) |
| Ge | 0.67 | Thermo Avantage Densities.txt |
| GeO2 | 5.6 | standard literature value, not from the user's files |
| Graphite | 0 | Thermo Avantage Densities.txt |
| HfO2 | 5.7 | Thermo Avantage Densities.txt |
| In2O3 | 2.9 | rough estimate, not from the user's files |
| InP | 1.27 | user-supplied list (unsourced) |
| MgO | 7.8 | standard literature value, not from the user's files |
| PbS | 0.37 | user-supplied list (unsourced) |
| PbSO4 | 4.5 | rough estimate, not from the user's files |
| Si | 1.12 | Thermo Avantage Densities.txt |
| Si3N4 | 4.65 | Thermo Avantage Densities.txt |
| SiC | 2.86 | Thermo Avantage Densities.txt |
| SiO2 | 8.9 | Thermo Avantage Densities.txt |
| Sn (gray) | 0.08 | Thermo Avantage Densities.txt |
| SnO2 | 3.6 | standard literature value, not from the user's files |
| TeO2 | 3.7 | standard literature value, not from the user's files |
| TiO2 | 3.05 | Thermo Avantage Densities.txt |
| ZnO | 3.2 | Thermo Avantage Densities.txt |
| ZnS | 3.6 | Thermo Avantage Densities.txt |
| ZrO2 | 6 | Thermo Avantage Densities.txt |

## Notes on particular sets

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
