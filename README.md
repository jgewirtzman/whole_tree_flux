# Whole-tree methane flux

Data and R code for ground-to-canopy methane measurements, tests of basal-stem extrapolation, and tree-surface scaling scenarios. The field analysis retains 141 measurements on seven temperate trees; a separate wetland reference comprises 105 repeated measurements of one black gum tree.

**Start with `run_analysis.R`.** It runs the active workflow in order, checks inputs and package versions before changing outputs, and runs the analysis checks after scaling. Current figures and numerical results are included so they can be inspected without rerunning the analysis. Manuscript text and editorial materials are maintained separately and are not needed here.

## Reproduce the results

Use **R 4.4.3**, Git, and the package versions in `renv.lock`. Run these commands from the repository root after cloning it:

```sh
# One-time setup: restore and activate a local project library.
Rscript -e 'if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv", repos = "https://cloud.r-project.org"); renv::restore(prompt = FALSE); renv::activate()'
```

Download the separately distributed soil file as described in [the soil-data instructions](scaling/soil_jevon2023/README.md), then:

```sh
Rscript run_analysis.R --check    # Verify packages and inputs; writes nothing
Rscript run_analysis.R            # Rebuild fluxes, analyses, figures and summaries
```

The full run uses the saved, manually selected chamber windows for Harvard Forest and Yale Myers Forest. It does not ask the user to select them again. The black gum reference is rebuilt from raw analyzer files and field-recorded windows. See [reproduction details](docs/REPRODUCING.md) for the inputs, processing order, units, and window-selection procedure.

For a shorter run using the committed flux tables:

```sh
Rscript run_analysis.R --scaling-only
```

To recalculate fluxes and detection flags without the soil comparison or figures:

```sh
Rscript run_analysis.R --flux-only
```

Use `--help` for all options. `--flux-only` and `--scaling-only` cannot be combined. Analysis scripts do not install packages during a run; restore the environment first.

## Workflow and outputs

| Stage | Code | Outputs |
|---|---|---|
| Chamber fluxes and detection flags | `data processing/goFlux_reprocessing/` | Harvard Forest, Yale Myers and black gum flux tables |
| Basal extrapolation tests | `scaling/01_flux_form_test.R` | `scaling/out/form_test_*.csv` |
| Measured stem geometry and stand scenarios | `scaling/03_tree_component.R` | Geometry, component rates, stand totals, uncertainty and SI fit figure |
| Prediction/observation comparison | `scaling/04_main_figures.R` | SI ratio figure and `global_demo_v3.csv` |
| Main figures and profile SI | `scaling/08_figures_v3c.R` | `scaling/v3c/` |
| Tree weighting sensitivity | `scaling/09_tree_weighting.R` | `tree_weighting.csv` |
| Basal/upper and wetland reference summaries | `scaling/10_basal_by_site_and_blackgum.R` | `tree_basal_upper_branch.csv`, `blackgum_scaling.csv` |
| Summary statistics and filter sensitivity | `scaling/11_summary_statistics.R` | `summary_statistics.csv`, `component_summary.csv`, `filter_sensitivity_HF.csv` |

All tabular scaling outputs are in `scaling/out/`. Historical version suffixes in the current figure paths are retained to keep existing references stable; the following table identifies the current set.

| Figure | File |
|---|---|
| 1 — measured profiles and field photographs | `scaling/v3c/Fig1_main_raw.png` |
| 2 — surface-area coverage and scaling | `scaling/v3c/Figure2_v3c.png` |
| 3 — extrapolation assumptions and scenarios | `scaling/v3c/Figure3_v3c.png` |
| S1 — profiles on an arcsinh axis | `scaling/v3c/FigS1_grid_asinh.png` |
| S2 — fitted extrapolations by tree | `scaling/fig_SI_extrapolation_fits.png` |
| S3 — predicted/measured upper-stem ratios | `scaling/Figure2_v3.png` (also PDF) |

Alternative raw/arcsinh profile renderings are also generated. Superseded figure pipelines and manuscript-building tools are not part of this public workflow.

## Data and assumptions

- **Primary observations:** raw analyzer records, chamber timing keys, field metadata and leaf areas are under `data processing/`. Harvard Forest has 136 compiled rows, of which 134 are retained; Yale Myers has 9 compiled rows, of which 7 are retained. The filtering rules are in `scaling/00_load_field.R`.
- **Meteorology:** the Fisher station input (`hf001-10-15min-m.csv`) is the existing cached Harvard Forest HF001 data used for temperature and pressure. Source and replacement instructions are in [reproduction details](docs/REPRODUCING.md).
- **Soil comparison:** [Jevon (2023), Mendeley Data V2](https://doi.org/10.17632/z6wybrtpyk.2), downloaded separately. July–August soil measurements are from different years than the tree campaign.
- **Geometry:** sampled stem fluxes use measured diameters and frustum areas, ending at each highest stem chamber. Stand scaling and capture use an illustrative 23 m cone; the black gum scenario uses 15.8 m. Sullivan et al. (2017), Table 6 supplies height proxies, not a taper law. Shared definitions are in `scaling/analysis_helpers.R`.
- **Area indices:** stems 0.45, branches 1.70 and leaves 4.5 m² per m² ground are midpoints of the ranges in Whittaker and Woodwell (1967). The alternative woody-area scenario totals 3.07. These inputs are assumptions, not measurements of the sampled trees.
- **Uncertainty:** 2,000 hierarchical bootstrap draws resample trees, then closures within tree/height/component, and refit models and basal denominators. Primary all-form comparisons use five eligible trees; other forms are also summarized separately for all seven. Exponential intervals are conditional on complete valid fits. Stand intervals hold external area indices and height proxies fixed.

Global products are sensitivity scenarios, **not global flux estimates**. Ground-area rates are divided by woody-area index before multiplication by global woody surface area. Signed measurements below detection remain in the primary analysis; zeroing and filtering are diagnostic alternatives.

## Check changes

```sh
Rscript tests/analysis_checks.R   # Units, clocks, geometry, cohorts, metadata and totals
Rscript tests/public_files.R     # No private/editorial files in the public Git index
```

Enable the included privacy check before each local commit with `git config core.hooksPath .githooks`. This setting applies to this checkout only.

Each successful run saves `scaling/out/sessionInfo.txt` and local input checksums/run metadata. Numerical changes should be reviewed against the committed tables. See [reproduction details](docs/REPRODUCING.md) for expected warnings and reproducibility limits.
