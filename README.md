# Whole-tree methane flux

Data and R code for ground-to-canopy methane (CH₄) flux measurements on stems, branches and leaves of temperate
trees, tests of extrapolation from basal-stem measurements, and tree-surface scaling scenarios. The analysis uses
141 closures on seven upland trees (Harvard Forest, Massachusetts, 2023; Yale Myers Forest, Connecticut, 2022) and a
separate wetland reference of 105 repeated closures on one black gum (*Nyssa sylvatica*) in a peat swamp (2024).

Everything is rebuilt from the raw analyzer files and the closure metadata in `data/` by `run_analysis.R`. The
committed flux tables, results and figures let you inspect the outputs without rerunning anything.

## Layout

```
data/        inputs: raw analyzer files, closure metadata, leaf areas, weather, photographs (see data/README.md)
flux/        flux fitting and detection limits         -> flux/out/   (one table per site)
scaling/     extrapolation tests and scaling scenarios -> scaling/out/ (tables)
figures/     Figures 1-3 and S1-S6
tests/       checks run by the workflow
scripts/     input and package checks used by run_analysis.R
docs/        reproduction details and methods notes
```

## Reproduce the results

Use **R 4.4.3** and the package versions in `renv.lock` (goFlux 0.4.0, commit `aee8456`). From the repository root:

```sh
# One-time setup: restore and activate a project library
Rscript -e 'if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv", repos = "https://cloud.r-project.org"); renv::restore(prompt = FALSE); renv::activate()'
```

Download the soil CH₄ file as described in [data/soil_jevon2023/README.md](data/soil_jevon2023/README.md), then:

```sh
Rscript run_analysis.R --check          # verify packages and inputs; writes nothing
Rscript run_analysis.R                  # fluxes, analyses, figures and summaries (about 3 min)
Rscript run_analysis.R --flux-only      # fluxes and detection limits only (no soil file needed)
Rscript run_analysis.R --scaling-only   # analyses and figures from the committed flux tables
Rscript scaling/13_flux_uncertainty.R   # optional: fit-uncertainty sensitivity (after a flux run)
```

No step is interactive. The measurement window of each closure is recorded in the metadata (see
[data/README.md](data/README.md)).

## Workflow

| Step | Code | Output |
|---|---|---|
| Harvard Forest fluxes | `flux/01_harvard_forest.R` | `flux/out/harvard_forest_fluxes_fitted.csv` |
| Yale Myers fluxes | `flux/02_yale_myers.R` | `flux/out/yale_myers_fluxes_fitted.csv` |
| Black Gum Swamp fluxes and detection limits | `flux/03_black_gum_swamp.R` | `flux/out/black_gum_swamp_fluxes.csv` |
| Detection limits (Harvard Forest, Yale Myers) | `flux/04_detection_limits.R` | `flux/out/harvard_forest_fluxes.csv`, `yale_myers_fluxes.csv`, `precision_comparison.csv` |
| Basal extrapolation tests | `scaling/01_flux_form_test.R` | `scaling/out/form_test_*.csv` |
| Stem geometry, component rates and stand scenarios | `scaling/03_tree_component.R` | component and stand tables; Figure S4 |
| Predicted/measured comparison; global scenarios | `scaling/04_main_figures.R` | Figure S5; `global_scenarios.csv` |
| Main figures | `scaling/08_figures.R` | Figures 2, 3 and S2 |
| Tree weighting sensitivity | `scaling/09_tree_weighting.R` | `tree_weighting.csv` |
| Basal/upper summaries; wetland reference | `scaling/10_basal_by_site_and_blackgum.R` | `blackgum_*.csv`; Figure S6 |
| Summary statistics and filter sensitivity | `scaling/11_summary_statistics.R` | `summary_statistics.csv`, `component_summary.csv`, `filter_sensitivity_HF.csv` |
| Tree illustrations | `scaling/12_tree_illustrations.R` | Figures 1 and S3 |
| Chamber photographs | `scaling/14_chamber_photos.R` | Figure S1 |

Each flux step imports the site's analyzer files, cuts each closure's window, attaches chamber area, system volume
and station temperature and pressure, estimates analyzer noise, and fits linear and Hutchinson–Mosier models with
goFlux, keeping the better-supported model (`best.flux`). The per-site flux tables in `flux/out/` hold one row per
closure: its metadata, CO₂ and CH₄ fluxes, fit diagnostics and detection limits. Fluxes are nmol CH₄ (µmol CO₂)
m⁻² s⁻¹ of enclosed surface; positive values are emission.

## Figures

| File | Content |
|---|---|
| `Figure1` | Measured profiles, tree silhouettes, field photographs and height-band summary |
| `Figure2` | Surface-area coverage and scaling |
| `Figure3` | Extrapolation assumptions and scenarios |
| `FigureS1` | Stem, branch and foliage chambers |
| `FigureS2` | All profiles on an arcsinh axis |
| `FigureS3` | Tree illustrations with stem flux (see [scaling/tree_illustrations/](scaling/tree_illustrations/README.md)) |
| `FigureS4` | Fitted extrapolations by tree |
| `FigureS5` | Predicted/measured upper-stem ratios |
| `FigureS6` | Swamp black gum component-budget scenarios |

## Data and assumptions

- **Closures used:** 133 of the 136 Harvard Forest closures (black gum tag 3, with three basal closures only, is not
  one of the six profiled trees) and 8 of the 9 Yale Myers closures (the closure recorded as "Seam" has no height).
  These rules are in `scaling/00_load_field.R`.
- **Detection limits:** analyzer noise is estimated from second differences within each closure, grouped by
  analyzer, field day and logging interval (`flux/precision_helpers.R`). The same empirical threshold (1.96 σ) is used
  in flux fitting and as the minimum detectable flux (MDF). See
  [docs/REPRODUCING.md](docs/REPRODUCING.md#empirical-noise-and-detection-limits).
- **Weather:** air temperature and pressure from Harvard Forest dataset HF001 (Fisher station; Boose & VanScoy,
  Harvard Forest Data Archive), included for the campaign periods under its CC0 1.0 release. Please cite the dataset
  if you reuse it.
- **Soil comparison:** [Jevon (2023), Mendeley Data V2](https://doi.org/10.17632/z6wybrtpyk.2), downloaded separately.
  July–August soil measurements are from different years than the tree campaign.
- **Leaf areas:** measured shoot areas where the sample was identified; otherwise the species mean. Closures of one
  shoot share an area (`flux/leaf_area_helpers.R`).
- **Geometry:** sampled stem fluxes use measured diameters and frustum areas, ending at each highest stem chamber.
  Stand scaling uses an illustrative 23 m cone; the black gum scenario uses 15.8 m (height proxies from Sullivan et al.
  2017, Table 6). Shared definitions are in `scaling/analysis_helpers.R`.
- **Area indices:** stems 0.55 and branches 1.55 m² per m² ground, midpoints of the ranges reported for three eastern
  North American temperate deciduous forests (Whittaker & Woodwell, 1967, p. 937); leaves 4.5, the midpoint of the
  reported 3–6 range. The alternative woody-area scenario totals 3.07 (stem 0.55, branches 2.52). These are scenario
  assumptions, not measurements of the sampled trees.
- **Uncertainty:** 2,000 hierarchical bootstrap draws resample trees, then closures within tree, height and tissue,
  and refit models and basal denominators. Stand intervals hold area indices and height proxies fixed.

Global products are sensitivity scenarios, **not global flux estimates**. Signed measurements below detection remain
in the primary analysis; zeroing and filtering are diagnostic alternatives.

## Checks

```sh
Rscript tests/precision_checks.R   # noise estimator, trend and gap handling, fitting durations
Rscript tests/analysis_checks.R    # units, clocks, windows, geometry, cohorts, metadata and totals
Rscript tests/public_files.R       # only project files are tracked
```

`run_analysis.R` runs the first two automatically. Each run writes `scaling/out/sessionInfo.txt` and local input
checksums. See [docs/REPRODUCING.md](docs/REPRODUCING.md) for expected warnings and reproducibility limits.
