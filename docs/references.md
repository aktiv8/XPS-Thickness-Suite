# References

The numbers match the bracketed citations in [Methods](methods.md) and the [User guide](user-guide.md). DOIs are given
only where they were checked against the paper itself or a bibliographic record.

## References

1. B. R. Strohmeier, An ESCA method for determining the oxide thickness on aluminum alloys, *Surf. Interface Anal.* **15**,
   51-56 (1990). DOI [10.1002/sia.740150109](https://doi.org/10.1002/sia.740150109).
2. HarwellXPS, Technical Note #1: Oxide thickness determination using the Strohmeier equation (the mixed
   oxide and hydroxide form of the equation).
   https://subsite.harwellxps.uk/wp-content/uploads/2018/02/HarwellXPS_TechNote_01.pdf
3. P. J. Cumpson, The Thickogram: a method for easy film thickness measurement in XPS, *Surf. Interface Anal.* **29**,
   403-406 (2000). DOI `10.1002/1096-9918(200006)29:6<403::AID-SIA884>3.0.CO;2-8`.
4. A. G. Shard, J. Wang and S. J. Spencer, **XPS topofactors: determining overlayer thickness on particles and
   fibres**, *Surf. Interface Anal.* **41**, 541-548 (2009). DOI [10.1002/sia.3044](https://doi.org/10.1002/sia.3044).
   This is the source of the Topofactor correction, its Eq. 6 and Table 1.
5. M. C. Biesinger, Advanced analysis of copper X-ray photoelectron spectra, *Surf. Interface Anal.* **49**, 1325-1334
   (2017). DOI [10.1002/sia.6239](https://doi.org/10.1002/sia.6239). Source of the Cu LMM presets, the Cu(II) satellite
   ratios and the worked Cu example.
6. S. Tanuma, C. J. Powell and D. R. Penn, Calculations of electron inelastic mean free paths. V. Data for 14 organic
   compounds over the 50-2000 eV range, *Surf. Interface Anal.* **21**, 165-176 (1994). DOI
   [10.1002/sia.740210302](https://doi.org/10.1002/sia.740210302). This paper introduced the TPP-2M predictive equation
   used for the computed IMFPs.
7. T. A. Carlson and G. E. McGuire, Study of the x-ray photoelectron spectrum of tungsten-tungsten oxide as a
   function of thickness of the surface oxide layer, *J. Electron Spectrosc. Relat. Phenom.* **1**, 161-168 (1972/73).
   Source of the W / WO3 presets and the original uniform-overlayer ratio.

## Software and data

8. M. C. Biesinger, Aluminum oxide-hydroxide thickness calculator, spreadsheet, Surface Science Western (2008). It uses
   Strohmeier's values and is reproduced by the Al presets.
9. J. W. Lee, Hallym University, the GlowScript program from which the Thickogram code derives.
   https://www.glowscript.org/#/user/jwsslee/folder/Private/program/thickogram
10. The standalone Python/Tk conversion of that program, https://github.com/aktiv8/Thickogram, kept unmodified in the
    suite's `legacy` folder.
11. Thermo Scientific Avantage, the `Densities.txt` table of densities and band gaps used for the Si, SiO2, Si3N4, Ti
    and TiO2 presets and several others.

## Citing

If you use the Topofactor correction in published work, please cite the paper by Shard and co-workers [4] and the
Thickogram paper by Cumpson [3]. For the Strohmeier calculator cite [1], and for the Cu presets [5].
