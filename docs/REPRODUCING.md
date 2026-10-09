# Reproducing the ground-to-canopy analysis

The runner is the authoritative processing order. The full workflow recalculates fluxes from saved measurement windows, rebuilds detection flags, and then regenerates the scaling analyses, figures and numerical summaries. It preserves the manually selected windows and does not edit original field records.

## Environment

Use R 4.4.3 and restore `renv.lock` as shown in the root README. See the renv documentation for [restoration](https://rstudio.github.io/renv/reference/restore.html) and [activation](https://rstudio.github.io/renv/reference/activate.html). Activation creates a local `.Rprofile` and `renv/` directory; these are ignored by Git. On Linux, some packages require system libraries for compilation and image rendering. Install the system dependencies reported by the package installer before retrying restoration.

The goFlux revision is `aee8456e62a016b6b496eb669e7445b160e65a9d`. Preflight checks compare directly used package versions and that revision against the lockfile. A different R version produces a warning. Restoring the environment requires network access; normal reproduction with all inputs present does not.

## Inputs and what is regenerated

| Input | Purpose | Included? |
|---|---|---|
| `data processing/input/` | Harvard Forest raw analyzer files and timing keys | Yes |
| `data processing/Field Data Entry - Clean Canopy Lift Total.csv` | Tree, chamber and measurement metadata | Yes |
| `data processing/leaf_areas.csv` | Measured leaf areas (cm²); authoritative area input | Yes |
| `data processing/leaf_area_mapping.csv` | Closure-to-shoot identities and measured-area links | Yes |
| `data processing/YMF Black Oak/` | Yale Myers raw records and field workbook; other ancillary field files | Yes |
| `data processing/goFlux_reprocessing/RData/manID_LGR*.RData` | Selected Harvard Forest windows and their concentration records | Yes; preserve these selections |
| `data processing/goFlux_reprocessing/RData/imp_LGR*_combined.RData` | Imported Harvard Forest series for legacy precision comparisons | Yes |
| `data processing/goFlux_reprocessing/ymf_black_oak/RData/manID_YMF.RData` and `imp_YMF_combined.RData` | Yale Myers window selections and concentration series | Yes |
| `data processing/goFlux_reprocessing/diurnal_blackgum/raw/diurnal/diurnal_final/` | Black gum analyzer files, field windows and chamber volumes | Yes |
| `data processing/goFlux_reprocessing/hf001-10-15min-m.csv` | Fisher station meteorology for 2022–2024 campaigns | Existing cached input included |
| `scaling/soil_jevon2023/fluxes.csv` | Independent soil comparison | No; download separately |
| `IMG_5926_edited.jpg`, `IMG_6437.jpg` | Field photographs in Figure 1 | Yes |

The [soil folder README](../scaling/soil_jevon2023/README.md) identifies the deposited input and required columns. It is needed for full and scaling-only runs, but not flux-only runs.

The meteorological source is the Harvard Forest Fisher Meteorological Station, HF001, [15-minute metric data](https://harvardforest1.fas.harvard.edu/data/p00/hf001/hf001-10-15min-m.csv). Preserve the cached file to reproduce this analysis. Replacing it with a later revision may change matched temperature and pressure. The processing requires the `datetime`, `airt`, and `bar` columns and coverage of October 2022, July–August 2023, and August 2024. If running an auxiliary script directly with no cache, it may download the station file; the main runner requires the input before starting.

Flat import directories, individual-file import caches, optional QC plots and black gum intermediate R objects are regenerated locally and ignored. Combined imports and selected windows remain included because the routine full run deliberately starts after manual window selection.

## Processing order

1. Rebuild Harvard Forest and Yale Myers auxiliary metadata (area, volume, temperature and pressure).
2. Synchronize saved window objects with that metadata and per-day empirical noise thresholds, fit goFlux models, and export compiled fluxes.
3. Rebuild empirical minimum-detectable-flux flags for Harvard Forest and Yale Myers.
4. Import black gum raw files, apply recorded start/end times, fit fluxes and calculate detection flags.
5. Load retained field measurements; test basal extrapolations on held-out upper-stem observations.
6. Integrate measured stem geometry; calculate component means, stand scenarios and uncertainty.
7. Generate the prediction/observation figure, global sensitivity table, main figures and profile SI.
8. Export tree-weighting, wetland-reference and summary statistics; run analysis checks.

`--scaling-only` starts at step 5 using committed compiled tables; it does not rebuild analyzer fluxes. `--flux-only` stops after step 4 and does not claim to refresh downstream scaling results. Preflight (`--check`, optionally combined with either mode) writes no outputs.

## Empirical noise and detection limits

`data processing/goFlux_reprocessing/precision_helpers.R` implements the calculation within this repository. For each selected closure, noise is the normal-consistent MAD of second concentration differences divided by √6. This removes a linear concentration trend without pooling different closure slopes. Differences across missing observations, gaps or logging-interval changes are excluded. The median closure estimate supplies each analyzer × field-day × recorded-interval group. Harvard Forest LGR1/LGR3 and the swamp records have 1 s intervals, LGR2 has 10 s, and Yale Myers has 5 s. These intervals do not establish whether the analyzer internally averaged or subsampled.

The fitting threshold is 1.96 times group noise. MDF uses that threshold, chamber geometry, pressure, temperature and closure duration (timestamp span plus one recorded interval). The pinned goFlux function assumes an extra 1 s in two duration expressions; `goFlux_at_interval()` adapts those expressions in a local copy so the curvature constraint and MDF use the same physical duration at all logging intervals. It leaves the installed package unchanged and stops if the pinned expressions change. Synthetic tests verify agreement with unmodified goFlux at 1 Hz and the expected MDF at 5 s and 10 s.

Each campaign exports `precision_closures_*.csv` and `precision_groups_*.csv` alongside its flux tables. Group diagnostics flag lag-one correlation of second differences above −0.5, a fraction of repeated values at least 0.3, or a centered-first-difference/second-difference noise ratio outside 0.8–1.2. These are review flags, not automatic exclusions. All CH₄ groups pass; the 10 s CO₂ groups warrant caution. Difference-based noise does not bound slow drift, leakage or total flux uncertainty. Signed below-MDF CH₄ observations remain in the main analysis.

The `wass95` output suffix is retained for compatibility; it now denotes this empirical reference MDF. The legacy `MDF_goflux` columns hold the manufacturer comparison, whereas `MDF` is the threshold actually used in the current fits. Whole-record first-difference MAD, per-closure difference SD, datasheet-based thresholds and alternative filters remain comparison diagnostics only. Datasheet scaling by √interval is hypothetical unless interval averaging is independently established. No unreleased package or external local repository is required.

## Selecting windows again

Routine reproduction does not require RStudio or interactive input. If raw data or window selections must change, use the numbered preparation/import scripts (`01`, `02`, and `03`) in `data processing/goFlux_reprocessing/`, then `04_manual_id.R` interactively in RStudio. The Yale Myers equivalents are prefixed `ymf_`. Rerun the full runner afterwards. Preserve previous `manID` files before replacing selections. The `03b`/`03c` scripts are optional metadata repair utilities; the current flux-fitting steps synchronize those fields automatically.

The black gum reference uses field-recorded start/end windows in `diurnal_updated.csv`; it does not use interactive peak selection. Its 105 closures are repeated observations of one tree.

## Units and interpretation

- Positive CH₄ flux means emission; negative means uptake. Per-surface fluxes are nmol CH₄ m⁻² surface s⁻¹. Integrated stand contributions are nmol CH₄ m⁻² ground s⁻¹.
- Soil input is µmol CH₄ m⁻² ground s⁻¹ and is multiplied by 1,000 before comparison.
- Auxiliary chamber area is cm², volume is liters, pressure is kPa and temperature is °C. Harvard/Yale system volume adds 0.028 L analyzer plus 0.029 L tubing. Black gum chamber-volume records already include their tubing, so only 0.028 L analyzer volume is added.
- Fisher timestamps use fixed EST. They are converted to civil Eastern time and relabeled for matching the analyzer's naive local clock. Closure duration is in seconds, including one logging interval.
- Sampled profile integration uses measured diameters, physical height intervals and frustum lateral areas. Boundary values use the nearest observation. It ends at the highest stem chamber.
- Complete stand-area and wetland scenarios use literature-based height proxies and assumed cone taper. Applying sampled compartment means to unmeasured surfaces is an extrapolation.
- Cross-form sensitivity uses one five-tree cohort; the six-tree Harvard Forest baseline and seven-tree descriptive summaries are distinct, labeled outputs. Yale Myers does not enter the Harvard Forest soil comparison.
- External area indices, height proxies, leaf-area imputation and chamber systematics are not resampled. Bootstrap intervals therefore describe flux sampling uncertainty conditional on those assumptions.

## Validation and repeatability

`tests/precision_checks.R` checks noise recovery, invariance to linear trends, gap/missing-data boundaries, day grouping and fitting durations. `tests/analysis_checks.R` also checks consistency between fitting and reported MDFs, the inverse rate–area calculation and component-budget shares, along with clock conversion, chamber durations and synchronized metadata, analytical cylinder/frustum areas, irregular-height integration, cone consistency, extrapolation behavior, cohort integrity, sample counts and stand accounting. When a locally generated black gum cache is absent, it checks committed durations and reports that the window cross-check will run after a flux rebuild.

Fixed seeds and 2,000 draws define the reported bootstrap results. `FLUX_N_BOOT` can be changed for diagnostic runs (minimum 100); such output does not reproduce the committed confidence intervals. Per-run `run_inputs.csv` records input paths and MD5 checksums before execution, `run_metadata.txt` records mode/time/draw count, and `sessionInfo.txt` records the R environment. Input/run files are local provenance, not replacement reference results.

Sparse-profile LOESS fits can issue near-singularity warnings, and the recorded ggplot version may warn about deprecated horizontal error bars. The LOESS curves are descriptive and do not determine flux estimates or bootstrap intervals. Font/rendering differences across systems can change figure pixels without changing tabular results; compare numeric tables separately from image files.

`tests/public_files.R` checks the current public Git index, including staged additions and removals, for manuscript text, Word files, editorial scripts and internal notes. It does not audit past commits or guarantee removal from remote history. Keep private editorial material outside this public index.

## Leaf areas and corrected field metadata

`leaf_areas.csv` is the authoritative table of measured areas in cm²; the XLSX is a matching convenience copy. `leaf_area_mapping.csv` assigns every foliage closure to a shoot and, where identifiable, a measured sample. Missing or unresolved assignments use the mean of measured areas within the same species. The four recovered July black-gum areas are 262.8, 319.8, 334.2 and 228.8 cm² (mean 286.4). Separate closures do not necessarily represent separate shoots: the 13:24 and 13:28:30 closures on July 18 share one shoot in light and dark and therefore share one area. Their conflicting field sample numbers are retained in Notes; the mapping documents use of the species mean. The 12:09 and 12:27 closures are treated as separate shoots. No shared-area assignment is inferred solely from repetition numbers.

Field tissue/chamber labels for the July 18 12:56 stem and 13:15 foliage records were corrected following review of field records and author confirmation. Timing-key heights now match the field-height column. Original UniqueIDs are retained as stable join keys; their embedded heights are not authoritative. The shared area helper validates timing/field metadata, exports the exact leaf-area assignments used, and prevents scaling from using compiled tables with outdated metadata. Saved concentration windows are retained.
