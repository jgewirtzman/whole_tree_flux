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
| Main figures and profile SI | `scaling/08_figures_v3c.R` | Figures 2, 3 and S2 |
| Tree illustrations | `scaling/12_tree_illustrations.R` | Figures 1 and S3 |
| Chamber photographs | `scaling/14_chamber_photos.R` | Figure S1 |
| Tree weighting sensitivity | `scaling/09_tree_weighting.R` | `tree_weighting.csv` |
| Basal/upper and wetland reference summaries | `scaling/10_basal_by_site_and_blackgum.R` | `tree_basal_upper_branch.csv`, `blackgum_scaling.csv`, `blackgum_component_budget.csv` and SI component-budget figure |
| Individual-flux uncertainty and model-choice sensitivity (optional separate step) | `scaling/13_flux_uncertainty.R` | `scaling/out/flux_uncertainty/` |
| Summary statistics and filter sensitivity | `scaling/11_summary_statistics.R` | `summary_statistics.csv`, `component_summary.csv`, `filter_sensitivity_HF.csv`, `rate_area_tradeoff.csv` |

All tabular scaling outputs are in `scaling/out/`. All figures are written to `figures/` under their manuscript numbers (PNG; PDF also for Figures 1, S1, S3, S5 and S6):

| Figure | Content | Made by |
|---|---|---|
| `Figure1` | Measured profiles, tree silhouettes, field photographs and height-band summary | `scaling/tree_illustrations/scripts/render_side_silhouettes.R` |
| `Figure2` | Surface-area coverage and scaling | `scaling/08_figures_v3c.R` |
| `Figure3` | Extrapolation assumptions and scenarios | `scaling/08_figures_v3c.R` |
| `FigureS1` | Stem, branch and foliage chambers | `scaling/14_chamber_photos.R` |
| `FigureS2` | All profiles on an arcsinh axis | `scaling/08_figures_v3c.R` |
| `FigureS3` | Tree illustrations with stem flux | `scaling/tree_illustrations/scripts/render_final_si.R` |
| `FigureS4` | Fitted extrapolations by tree | `scaling/03_tree_component.R` |
| `FigureS5` | Predicted/measured upper-stem ratios | `scaling/04_main_figures.R` |
| `FigureS6` | Swamp black gum component-budget scenarios | `scaling/10_basal_by_site_and_blackgum.R` |

The tree drawings in Figures 1 and S3 are illustrations built from explicit geometry inputs; see [the tree-illustration notes](scaling/tree_illustrations/README.md).

## Data and assumptions

- **Primary observations:** raw analyzer records, chamber timing keys, field metadata and leaf areas are under `data processing/`. Harvard Forest has 136 compiled rows, of which 134 are retained; Yale Myers has 9 compiled rows, of which 7 are retained. The filtering rules are in `scaling/00_load_field.R`.
- **Detection limits:** a standalone helper estimates noise from second differences within closures, grouped by analyzer, field day and recorded interval. The same empirical threshold is used in flux fitting and MDF flags; diagnostics are exported alongside the flux tables. See [method details](docs/REPRODUCING.md#empirical-noise-and-detection-limits).
- **Meteorology:** the Fisher station input (`hf001-10-15min-m.csv`) is a cached copy of Harvard Forest dataset HF001 (Boose & VanScoy, Harvard Forest Data Archive), used for temperature and pressure. It is redistributed under the archive's CC0 1.0 release; please cite the dataset if you reuse it. Source and replacement instructions are in [reproduction details](docs/REPRODUCING.md).
- **Soil comparison:** [Jevon (2023), Mendeley Data V2](https://doi.org/10.17632/z6wybrtpyk.2), downloaded separately. July–August soil measurements are from different years than the tree campaign.
- **Leaf areas:** `data processing/leaf_areas.csv` holds measured areas; `leaf_area_mapping.csv` links closures to shoots and identified samples. Unresolved areas use within-species means, and paired light/dark closures share an area. The compiled flux tables retain each assignment and its source. See [the reproduction guide](docs/REPRODUCING.md#leaf-areas-and-corrected-field-metadata).
- **Geometry:** sampled stem fluxes use measured diameters and frustum areas, ending at each highest stem chamber. Stand scaling and capture use an illustrative 23 m cone; the black gum scenario uses 15.8 m. Sullivan et al. (2017), Table 6 supplies height proxies, not a taper law. Shared definitions are in `scaling/analysis_helpers.R`.
- **Area indices:** stems 0.55 and branches 1.55 m² per m² ground are midpoints of the ranges reported together for three eastern North American temperate closed-canopy deciduous forests (Whittaker & Woodwell, 1967, p. 937). Leaves retain 4.5, the midpoint of the reported 3–6 range. This is close to the approximately 4.6 inferred for Harvard Forest by subtracting the reported woody contribution (0.9) from its optical canopy area index (about 5.5; [Urbanski et al., 2007](https://doi.org/10.1029/2006JG000293), section 2.2). The alternative woody-area scenario totals 3.07, with the same stem index and the remaining 2.52 assigned to branches. These inputs are centralized in `scaling/analysis_helpers.R`; they are scenario assumptions, not measurements of the sampled trees.
- **Uncertainty:** 2,000 hierarchical bootstrap draws resample trees, then closures within tree/height/component, and refit models and basal denominators. Primary all-form comparisons use five eligible trees; other forms are also summarized separately for all seven. Exponential intervals are conditional on complete valid fits. Stand intervals hold external area indices and height proxies fixed.

Global products are sensitivity scenarios, **not global flux estimates**. Ground-area rates are divided by woody-area index before multiplication by global woody surface area. Signed measurements below detection remain in the primary analysis; zeroing and filtering are diagnostic alternatives.

## Check changes

```sh
Rscript tests/precision_checks.R  # Noise estimator, trend/gap handling and fitting durations
Rscript tests/analysis_checks.R   # Units, clocks, geometry, cohorts, metadata and totals
Rscript tests/public_files.R     # No private/editorial files in the public Git index
```

Enable the included privacy check before each local commit with `git config core.hooksPath .githooks`. This setting applies to this checkout only.

Each successful run saves `scaling/out/sessionInfo.txt` and local input checksums/run metadata. Numerical changes should be reviewed against the committed tables. See [reproduction details](docs/REPRODUCING.md) for expected warnings and reproducibility limits.

## Flux uncertainty sensitivity

After rebuilding fluxes and scaling, run `Rscript scaling/13_flux_uncertainty.R` to reproduce the model-specific interval classifications, correlated-residual sensitivity, and alternative budget and extrapolation point estimates. It uses archived fit standard errors, preserves the primary signed flux estimates, and writes separate tables to `scaling/out/flux_uncertainty/`. The swamp fit archive is rebuilt by `Rscript run_analysis.R --flux-only` and is not redistributed. Restore the updated `renv.lock` for the `sandwich` dependency. Nonlinear intervals are approximate and conditional on the selected fit; zeroing rules are sensitivity scenarios.

Field photographs for Figure 1 are in `data processing/field_photos/`; chamber photographs for Figure S1 are in `data processing/chamber_photos/`.
