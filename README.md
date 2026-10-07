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
| 05 | `05_flux_calculation.R` | Calculate fluxes with `goFlux()` (linear and Hutchinson-Mosier models) and select best estimate |
| 06 | `06_compile_results.R` | Merge flux results with field metadata; output `canopy_flux_goFlux_compiled.csv` |
| 07 | `07_quality_plots.R` | Generate quality-check plots |
| 08 | `08_ch4_height_plot.R` | CH₄ flux vs. height faceted by tree |

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
| `01_flux_form_test.R` | Extrapolates each tree's basal (< 2 m) stem fluxes upward with six forms and scores them against the measured fluxes ≥ 2 m (bootstrap) | `out/form_test_*.csv` |
| `03_tree_component.R` | Stem area from measured diameters; tree-weighted rates for the six Harvard Forest trees; woody-surface scenarios vs same-month soil flux; leaves; mixed model | `out/stand_rates_HF.csv`, `out/tree_component.csv`, `fig_SI_extrapolation_fits.png` |
| `04_main_figures.R` | Predicted/measured ratio figure (SI) and global scenario table | `Figure2_v3.png`, `out/global_demo_v3.csv` |
| `08_figures_v3c.R` | Main Figures 1–3 and the SI profile figure | `v3c/Fig1_main_raw.png` / `v3c/Fig1_main_asinh.png` (Fig. 1, raw or arcsinh flux axis), `v3c/Figure2_v3c.png` (Fig. 2), `v3c/Figure3_v3c.png` (Fig. 3), `v3c/FigS1_grid_asinh.png` / `v3c/FigS1_grid_raw.png` (Fig. S1), `out/blackgum_collar_diel_means.csv` |
| `09_tree_weighting.R` | Pooled vs tree-weighted means; leave-one-tree-out | `out/tree_weighting.csv` |
| `10_basal_by_site_and_blackgum.R` | Basal vs upper-stem fluxes by tree; black gum share of stem flux above 2 m and the wetland scenario | `out/tree_basal_upper_branch.csv`, `out/blackgum_scaling.csv` |

Global numbers produced by these scripts are scenarios that show sensitivity to assumptions, not estimates.

## Black Gum Swamp reference tree (`data processing/goFlux_reprocessing/diurnal_blackgum/`)

One *Nyssa sylvatica* in a saturated peat swamp at Harvard Forest, 28–29 August 2024, 105 closures at 0.25–3.6 m.
`bg_run.R` processes the raw analyzer files (`raw/diurnal/diurnal_final/LGR2/`) with goFlux using the field start and end
times, the chamber volumes and the Fisher station met data (HF001), and writes
`results/blackgum_flux_compiled_with_mdf.csv`.

## Key Outputs

- `data processing/goFlux_reprocessing/results/canopy_flux_goFlux_compiled.csv` — Harvard Forest compiled fluxes (141 measurements)
- `data processing/goFlux_reprocessing/ymf_black_oak/results/ymf_black_oak_flux_compiled.csv` — Yale Myers Forest compiled fluxes
- `figures/summary_statistics_output.txt` — Key manuscript statistics

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
remotes::install_github("Qepanna/goFlux")
```

### Python (optional)

- `pandas`, `numpy`, `matplotlib`, `statsmodels`, `Pillow`

Only used for an alternative version of Figure 1 (`figures/figure1_composite.py`).

## Reproducing the Analysis

1. Open `whole_tree_flux.Rproj` in RStudio
2. Run the pipeline scripts in order: `00_setup.R` through `06_compile_results.R`
   - Step 04 (`04_manual_id.R`) is interactive and requires manual peak identification in RStudio
   - Pre-existing manID RData files are included, so steps 01–04 can be skipped for flux recalculation
3. Run figure scripts from `figures/` and `truncation/`

To recalculate fluxes without re-running the interactive step:
```r
source("data processing/goFlux_reprocessing/00_setup.R")
source("data processing/goFlux_reprocessing/03_build_auxfiles.R")
source("data processing/goFlux_reprocessing/03b_patch_manID_leaf_areas.R")
source("data processing/goFlux_reprocessing/05_flux_calculation.R")
source("data processing/goFlux_reprocessing/06_compile_results.R")
```

## Data Sources

- **Harvard Forest meteorological data** (`hf001-10-15min-m.csv`): Fisher Meteorological Station 15-minute data, used for chamber temperature and pressure
- **Wu et al. (2024)** synthesis dataset (`1-s2.0-S0168192324000911-mmc2.xlsx`): compiled stem CH₄ flux measurements from 50 studies (1,010 observations), published in *Agricultural and Forest Meteorology*
- **Whittaker & Woodwell (1967)**: Surface area indices (stem = 0.45, branch = 1.70, LAI = 4.5 m² m⁻² ground) used in truncation and sensitivity analyses
