# Whole-Tree Methane Flux

This repository contains data, processing code, and figure-generating scripts for a study of methane (CH₄) fluxes across the full vertical profile of trees — from stem base to canopy.


## Repository Structure

```
whole_tree_flux/
├── data processing/
│   ├── goFlux_reprocessing/          # Harvard Forest flux pipeline (see below)
│   │   ├── ymf_black_oak/            # Yale Myers Forest flux pipeline
│   │   ├── RData/                    # Intermediate R objects
│   │   ├── results/                  # Compiled flux outputs (CSV, XLSX)
│   │   └── plots/                    # Quality-check plots
│   ├── input/                        # Raw LGR data and timing keys
│   ├── functions/                    # Legacy pre-goFlux flux functions
│   ├── YMF Black Oak/                # YMF raw data and gas chromatography
│   ├── Field Data Entry - Clean Canopy Lift Total.csv
│   ├── leaf_areas.csv                # Measured leaf surface areas
│   └── leaf_areas.xlsx
│
├── figures/                          # Figure scripts and rendered outputs
│   ├── figure1_*.R                   # Flux profiles (Harvard Forest, composite)
│   ├── figure2_yale_forest.R         # Yale Myers Forest profile
│   ├── figure3_wu_reanalysis.R       # Wu et al. (2024) synthesis reanalysis
│   ├── figure_sensitivity_breakeven.R
│   ├── figure_wu_below2m.R
│   ├── stats_wu_above2m_slopes.R
│   └── summary_statistics_output.txt # Key manuscript numbers
│
├── scaling/                          # Ground-to-canopy analyses for the GRL manuscript (see below)
├── truncation/                       # Truncation analysis and tree geometry
│   ├── figure_truncation.R           # Standalone truncation figure
│   ├── figure_truncation_combined.R  # Combined 7-panel truncation + sensitivity
│   └── math.R                        # Cone taper geometry calculations
│
├── deprecated/                       # Superseded scripts and figures
│
├── IMG_5926_edited.jpg               # Field photo: canopy lift
├── IMG_6437.jpg                      # Field photo: arborist climbing
└── whole_tree_flux.Rproj
```

## Flux Processing Pipeline

Gas fluxes are calculated using the [goFlux](https://github.com/Qepanna/goFlux) R package. The pipeline is a numbered sequence of scripts in `data processing/goFlux_reprocessing/`:

| Step | Script | Description |
|------|--------|-------------|
| 00 | `00_setup.R` | Load packages, define paths and constants (chamber volume, observation length, best.flux criteria) |
| 01 | `01_prepare_raw_data.R` | Copy raw LGR files into flat staging directories |
| 02 | `02_import.R` | Import raw data via `goFlux::import2RData()` |
| 03 | `03_build_auxfiles.R` | Build auxiliary files (UniqueID, Area, Vtot, Tcham, Pcham) from timing keys and Harvard Forest met data |
| 03b | `03b_patch_manID_leaf_areas.R` | Patch leaf-type measurements with actual measured leaf areas (replaces chamber areas) |
| 03c | `03c_patch_manID_vtot.R` | Sync `Vtot` in existing manID objects with the rebuilt auxfiles (when the analyzer/tubing volume in `00_setup.R` changes) |
| 04 | `04_manual_id.R` | **Interactive** — manually identify gas concentration peaks in RStudio via `click.peak2()` |
| 05 | `05_flux_calculation.R` | Synchronize current metadata/precision, calculate fluxes with `goFlux()` and select best estimate |
| 06 | `06_compile_results.R` | Merge flux results with field metadata; output `canopy_flux_goFlux_compiled.csv` |
| 07 | `07_quality_plots.R` | Generate quality-check plots |
| 08 | `08_ch4_height_plot.R` | CH₄ flux vs. height faceted by tree |
| 09 | `09_mdf_lod_comparison.R` | Rebuild empirical MDF flags; required before current scaling/figures |

A parallel pipeline exists for Yale Myers Forest data in `ymf_black_oak/` (scripts prefixed `ymf_`).

**Note:** Step 04 is interactive and requires RStudio. Steps 03b and 03c patch existing manID objects so that step 04 does not need to be re-run when leaf areas or the system volume change (`ymf_03b_patch_manID_vtot.R` does the same for the Yale Myers pipeline).

**System volume:** `Vtot` = chamber volume + 0.057 L (0.028 L analyzer internal volume + 0.029 L tubing; `vtot_addition` in `00_setup.R`).

## Ground-to-canopy analyses for the GRL manuscript (`scaling/`)

Run every script from the repository root (paths are relative). Inputs are the compiled flux files above, the black gum
results below, and the Harvard Forest soil fluxes of Jevon (2023), which are not redistributed here (see
`scaling/soil_jevon2023/README.md`).

| Script | What it does | Main outputs |
|---|---|---|
| `00_load_field.R` | Loads the Harvard Forest and Yale Myers fluxes (one row per measurement) | object `F` |
| `01_flux_form_test.R` | Tests basal extrapolation using measured stem geometry, a common eligible cohort and hierarchical bootstrap | `out/form_test_*.csv` |
| `03_tree_component.R` | Stem area from measured diameters; tree-weighted rates for the six Harvard Forest trees; woody-surface scenarios vs same-month soil flux; leaves; mixed model | `out/stand_rates_HF.csv`, `out/tree_component.csv`, `fig_SI_extrapolation_fits.png` |
| `04_main_figures.R` | Predicted/measured ratio figure (SI) and global scenario table | `Figure2_v3.png`, `out/global_demo_v3.csv` |
| `08_figures_v3c.R` | Main Figures 1–3 and the SI profile figure | `v3c/Fig1_main_raw.png` / `v3c/Fig1_main_asinh.png` (Fig. 1, raw or arcsinh flux axis), `v3c/Figure2_v3c.png` (Fig. 2), `v3c/Figure3_v3c.png` (Fig. 3), `v3c/FigS1_grid_asinh.png` / `v3c/FigS1_grid_raw.png` (Fig. S1), `out/blackgum_collar_diel_means.csv` |
| `09_tree_weighting.R` | Pooled vs tree-weighted means; leave-one-tree-out | `out/tree_weighting.csv` |
| `11_manuscript_statistics.R` | Current descriptive statistics and filter sensitivity | `out/manuscript_*.csv`, `out/filter_sensitivity_HF.csv` |
| `10_basal_by_site_and_blackgum.R` | Basal vs upper-stem fluxes by tree; black gum share of stem flux above 2 m and the wetland scenario | `out/tree_basal_upper_branch.csv`, `out/blackgum_scaling.csv` |

Global numbers produced by these scripts are scenarios that show sensitivity to assumptions, not estimates. Before multiplying by global woody surface area, woody flux per unit ground area is divided by the woody-area index to obtain flux per unit woody surface area.

## Black Gum Swamp reference tree (`data processing/goFlux_reprocessing/diurnal_blackgum/`)

One *Nyssa sylvatica* in a saturated peat swamp at Harvard Forest, 28–29 August 2024, 105 closures at 0.25–3.6 m.
`bg_run.R` processes the raw analyzer files (`raw/diurnal/diurnal_final/LGR2/`) with goFlux using the field start and end
times, the chamber volumes and the Fisher station met data (HF001), and writes
`results/blackgum_flux_compiled_with_mdf.csv`.

## Key Outputs

- `data processing/goFlux_reprocessing/results/canopy_flux_goFlux_compiled.csv` — Harvard Forest compiled fluxes (136 rows; 134 retained measurements, or 141 combined with Yale Myers)
- `data processing/goFlux_reprocessing/ymf_black_oak/results/ymf_black_oak_flux_compiled.csv` — Yale Myers Forest compiled fluxes
- `scaling/out/manuscript_statistics.csv` and `manuscript_component_summary.csv` — Current manuscript statistics
- `figures/summary_statistics_output.txt` — Earlier-draft statistics (superseded)

## Figures (earlier draft; superseded by `scaling/08_figures_v3c.R`)

The Wu et al. (2024) reanalysis scripts need the Wu et al. (2024) supplementary table (Agricultural and Forest Meteorology 350:109976, Appendix A), which is not redistributed here.

| Figure | Script(s) | Description |
|--------|-----------|-------------|
| Figure 1 | `figures/figure1_composite.R`, `figure1_harvard_forest.R` | CH₄ flux profiles for all Harvard Forest trees + Yale Forest + field photos |
| Figure 2 | `figures/figure2_yale_forest.R` | Yale Myers Forest Black Oak stem CH₄ profile |
| Figure 3 | `figures/figure3_wu_reanalysis.R` | Reanalysis of Wu et al. (2024) synthesis: studies with measurements at ≥ 2 m |
| Truncation | `truncation/figure_truncation_combined.R` | Combined 7-panel: cone geometry, capture fraction, sensitivity analysis, break-even thresholds |
| Sensitivity | `figures/figure_sensitivity_breakeven.R` | Stand-level CH₄ budget sensitivity to surface area indices (Whittaker & Woodwell 1967) |
| Wu < 2 m | `figures/figure_wu_below2m.R` | Wu et al. studies with only below-2 m measurements |
| Wu slopes | `figures/stats_wu_above2m_slopes.R` | Per-study slope comparison (below-2 m vs. full range) |

## Dependencies

### R

Core packages:
- [goFlux](https://github.com/Qepanna/goFlux) — chamber flux calculation
- tidyverse (`dplyr`, `tidyr`, `purrr`, `readr`, `stringr`)
- `ggplot2`, `patchwork` — figures
- `lubridate` — datetime handling
- `openxlsx`, `readxl` — Excel I/O
- `scales`, `jpeg`, `grid` — figure formatting

Install goFlux from GitHub:
```r
remotes::install_github("Qepanna/goFlux@aee8456e62a016b6b496eb669e7445b160e65a9d")
```

### Python (optional)

- `pandas`, `numpy`, `matplotlib`, `statsmodels`, `Pillow`

Only used for an alternative version of Figure 1 (`figures/figure1_composite.py`).

## Reproducing the current analysis

Run commands from the repository root. All active processing paths are relative to that root.
R 4.4.3 and the package versions used for this rebuild are recorded in `renv.lock`, including
goFlux commit `aee8456e62a016b6b496eb669e7445b160e65a9d`. Restore the environment with
`renv::restore(lockfile = "renv.lock")` before reproducing numerical results.
Download the separately distributed soil input as described in `scaling/soil_jevon2023/README.md`.

```sh
# Rebuild fluxes, empirical detection flags, scaling, current figures and manuscript statistics.
# Existing manually selected measurement windows are preserved.
Rscript run_analysis.R
Rscript tests/analysis_checks.R

# Optional: rerun only flux processing or only downstream scaling and figures.
Rscript run_analysis.R --flux-only
Rscript run_analysis.R --scaling-only
```

The flux-calculation scripts automatically synchronize cached measurement objects with rebuilt
auxfiles (area, volume, temperature, pressure and precision). Precision is also passed explicitly
to goFlux. Old imports cannot silently override the settings in `00_setup.R`. Recalculation runs
`09_mdf_lod_comparison.R` before scaling; the field loader checks that MDF-enriched files match
the base flux outputs. The runner skips superseded quality-control figures; run script 09 with
`FLUX_LEGACY_QC=1` to redraw those exploratory plots.

Both Harvard/Yale and black-gum meteorological timestamps are converted from fixed EST to the
civil Eastern clock before matching the analyzer's local clock. Closure durations are in seconds,
including one logging interval. The selected chamber windows themselves are unchanged.

### Geometry, cohorts and uncertainty

- Stem fluxes are integrated over physical height intervals using interpolated measured diameters
  and frustum surface areas. Flux and geometry use the same profile, ending at the highest stem
  chamber. Missing boundary diameters and fluxes use their nearest observed value; no tip is invented.
- Stand area partitioning and the capture diagram share a 23 m representative cone; the black-gum
  scenario uses 15.8 m. Both are literature-based proxies from Sullivan et al. (2017), Table 6.
  Sampled compartment means are applied to the complete scenario areas. The capture diagram counts
  stem area represented by sampling, not all branches/leaves physically below a height.
- All-form comparisons use one fixed five-tree cohort with positive basal height means. Coverage and
  failures are exported in `form_test_coverage.csv`; `form_test_pooled_all_trees.csv` separately reports
  the other forms on all seven trees. Stand cross-form scenarios use the same five Harvard Forest trees
  for the measured benchmark, predictions, basal contribution and area partition. Figure 2's baseline
  remains the full six-tree Harvard Forest scenario.
- Uncertainty uses 2,000 hierarchical draws: resample trees, then closures within tree/component/height.
  Fits and relative-profile basal denominators are recalculated. Missing exponential fits invalidate the
  whole form/draw, never remove individual trees from its denominator. Valid-draw counts are exported.
- `stand_uncertainty_HF.csv` propagates flux variation to component contributions,
  shares, total flux and cancellation thresholds. It is conditional on the assumed external surface-area
  indices and proxy heights, and does not quantify leaf-area imputation or chamber-systematic uncertainty.
- `manuscript_component_summary.csv` and `manuscript_statistics.csv` provide current numerical summaries;
  the basal–upper comparison is paired across trees. `filter_sensitivity_HF.csv` reports zeroing/filtering
  as diagnostic scenarios only. The primary analysis retains signed values below the MDF.
- The 105 black-gum closures are repeated observations of one tree. They are not 105 independent trees.

`FLUX_N_BOOT` can change the number of draws for diagnostic runs (minimum 100); manuscript results use
2,000 and fixed seeds. Each full run records `scaling/out/sessionInfo.txt`.

### Private manuscript build

The private `manuscript.md` is the source of the combined manuscript and SI. The public code repository
does not contain that text, Word drafts, planning notes or third-party source data. Where the private
manuscript is available, build the revised draft with:

```sh
./build_manuscript.sh              # Word only, using current analysis and figures
./build_manuscript.sh --figures    # recompute figure summaries and redraw before building
./build_manuscript.sh --analysis   # full analysis, figures and Word
```

The outputs are `DRAFT_v3d_Ground-to-Canopy_2026-10-07.docx` (combined) plus files with `_main` and `_si`
suffixes (separate main text and supporting information); the previous v3c draft is preserved. These are
Pandoc drafts, using a generic reference document for black headings and figure/caption pagination,
not a validated AGU submission template. Citations remain manually maintained. The split
occurs at the unique `# Supporting Information` heading. After rerunning analysis, review manuscript numbers and conclusions against
the generated statistics before circulating a new draft. Use `./ms-git.sh` for private manuscript files
and plain `git` for code/data; do not stage private files or third-party inputs in the public repository.

## Data Sources

- **Harvard Forest meteorological data** (`hf001-10-15min-m.csv`): Fisher Meteorological Station 15-minute data, used for chamber temperature and pressure
- **Wu et al. (2024)** synthesis dataset (`1-s2.0-S0168192324000911-mmc2.xlsx`): compiled stem CH₄ flux measurements from 50 studies (1,010 observations), published in *Agricultural and Forest Meteorology*
- **Whittaker & Woodwell (1967)**: Surface area indices (stem = 0.45, branch = 1.70, LAI = 4.5 m² m⁻² ground) used in truncation and sensitivity analyses
