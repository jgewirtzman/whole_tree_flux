# Data

Inputs to the analysis. Everything in `flux/out/`, `scaling/out/` and `figures/` is rebuilt from these files by
`Rscript run_analysis.R`.

| Folder | Contents |
|---|---|
| `analyzer/` | Raw concentration records, one folder per site (Harvard Forest has one folder per analyzer) |
| `metadata/` | One row per chamber closure for each site, plus leaf areas and black gum chamber geometry |
| `meteorology/` | Harvard Forest Fisher station (HF001) records for the three campaigns |
| `photos/` | Field photographs (Figure 1) and chamber photographs (Figure S1) |
| `soil_jevon2023/` | Download instructions for the soil CH₄ fluxes used for comparison (not redistributed) |

## Sites and campaigns

| Site | Dates | Analyzer(s) | Closures | Trees |
|---|---|---|---|---|
| Harvard Forest, Petersham, MA (canopy lift; sites Swamp Rd and EMS) | 18–19 Jul and 16–17 Aug 2023 | LGR1, LGR2, LGR3 | 136 | *Nyssa sylvatica* (2), *Tsuga canadensis* (2), *Acer rubrum* (2), *Quercus rubra* (1)¹ |
| Yale Myers Forest, Connecticut (tree climbing) | 4 Oct 2022 | LGR3 | 9 | *Quercus velutina* (1) |
| Black Gum Swamp, Harvard Forest (saturated peat) | 28–29 Aug 2024 | LGR2 | 105 | *N. sylvatica* (1) |

¹ The second *N. sylvatica* (tag 3) has three basal stem closures only and is not part of the six-tree analysis.

All analyzers are ABB/Los Gatos Research GLA131-GGA microportable greenhouse gas analyzers (CO₂, CH₄, H₂O).

## analyzer/

The analyzer's concentration files (`micro_<date>_f<nnnn>.txt`), as logged. Times are the analyzer's local clock. The
analyzer's other log files (b, l, p and lims files) are not used and not included. Logging intervals are 1 s
(Harvard Forest LGR1 and LGR3, Black Gum Swamp), 10 s (Harvard Forest LGR2) and 5 s (Yale Myers).

## metadata/

Times are local clock times (`YYYY-MM-DD HH:MM:SS`, no time zone). `Field_*` times come from the field record;
`Window_*` times are the measurement window on the analyzer clock, the period used for flux fitting. Missing values
are empty cells.

### harvard_forest_closures.csv, yale_myers_closures.csv, black_gum_swamp_closures.csv

| Column | Description |
|---|---|
| `UniqueID` | Closure identifier (Harvard Forest: tree tag, nominal height, tissue, replicate; join key only) |
| `Date` | Measurement date |
| `Site` | Site (Harvard Forest: Swamp Rd or EMS) |
| `Analyzer` | Analyzer unit |
| `Tree_Tag` | Tree tag (Harvard Forest) |
| `Species` | `bg` *N. sylvatica*, `hem` *T. canadensis*, `rm` *A. rubrum*, `ro` *Q. rubra*, `qv` *Q. velutina* |
| `Height_m` | Chamber height above ground (m) |
| `Type` | Tissue: `stem`, `branch`, `leaf` or `leaf (shaded)` |
| `Replicate` | Repeat closure at the same position (Harvard Forest) |
| `ChamberID` | Chamber used |
| `Chamber_volume_L` | Chamber volume (L) |
| `Tubing_volume_L` | Tubing volume (L). The system volume is chamber + tubing + 0.028 L analyzer cell |
| `Chamber_area_cm2` | Enclosed bark area (cm²). Foliage closures use leaf areas instead (see below) |
| `Field_start`, `Field_end` | Closure start and end in the field record (local time) |
| `Analyzer_start`, `Analyzer_end` | Field record converted to the analyzer clock (Yale Myers; that analyzer clock ran 22 min ahead of local time) |
| `Window_start`, `Window_end` | Measurement window on the analyzer clock |
| `DBH_cm` | Diameter (cm) of the stem or branch at the chamber (field-sheet column name) |
| `Distance_cm` | Distance (cm) as recorded on the field sheet; 0 for every stem chamber |
| `Stem_diameter_cm` | Stem diameter at the chamber (cm; Yale Myers) |
| `Air_Temp_C`, `Humidity_RH`, `Soil_Temp_C`, `Stem_Temp_C` | Handheld readings at the time of measurement |
| `VWC_1`–`VWC_3` | Soil volumetric water content, three readings, as recorded (`sat` = saturated) |
| `Position`, `Chamber_type` | Collar position and chamber type (Black Gum Swamp) |
| `Notes` | Field notes |

**Measurement windows.** Harvard Forest and Yale Myers windows were selected by eye on each concentration trace
(goFlux `click.peak2`, displaying 6 min of record before and after the field-recorded closure) and are recorded here. Black Gum Swamp
windows are the field-recorded start and end times on the analyzer clock, with no further selection; two start times
were adjusted in the field record. The Black Gum Swamp field sheet lists analyzer LGR3, but the concentration records
for these closures are the LGR2 files included here.

**Weather matching.** Each closure uses the Fisher station record nearest to its `Field_start`.

### leaf_areas.csv, leaf_area_mapping.csv

`leaf_areas.csv` holds the measured one-sided areas (cm²) of the enclosed shoots by tree species and field sample
number (`Flux_Rep`). `leaf_area_mapping.csv` links each foliage closure to a shoot (`ShootID`) and, where the sample was
identified, to its measured area; `Mapping_basis` gives the basis. Unidentified shoots use the mean measured area of
their species. Closures of the same shoot (for example the ambient and shaded closures of one shoot) share an area.

### black_gum_swamp_chambers.csv

Geometry of the Black Gum Swamp collars and rigid chamber: four collar heights and diameter (cm), cylinder (or
container) volume, cap volume and tubing volume (L), their total, and enclosed area (cm²).

## meteorology/

`hf001-10-15min-m_campaigns.csv` is the subset of Harvard Forest dataset HF001 (Fisher Meteorological Station,
15-minute data, metric units; Boose & VanScoy, Harvard Forest Data Archive) covering 1–7 Oct 2022, 15 Jul–20 Aug 2023
and 25 Aug–3 Sep 2024, with all columns as published. Timestamps are EST. The analysis uses air temperature (`airt`,
°C) and barometric pressure (`bar`, mbar). HF001 is released under CC0 1.0; please cite the dataset if you reuse it.
Source: https://harvardforest1.fas.harvard.edu/data/p00/hf001/hf001-10-15min-m.csv

## photos/

`field/`: canopy lift (Harvard Forest) and arborist climbing (Yale Myers) photographs used in Figure 1.
`chambers/`: stem, branch and foliage chamber photographs used in Figure S1.
