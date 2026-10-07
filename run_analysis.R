#!/usr/bin/env Rscript
# Run from the repository root. Preserves all manually selected chamber windows.
stopifnot(file.exists("whole_tree_flux.Rproj"))
Sys.setenv(FLUX_LEGACY_QC = "0")
steps <- c(
  "data processing/goFlux_reprocessing/03_build_auxfiles.R",
  "data processing/goFlux_reprocessing/ymf_black_oak/ymf_03_build_auxfile.R",
  "data processing/goFlux_reprocessing/05_flux_calculation.R",
  "data processing/goFlux_reprocessing/06_compile_results.R",
  "data processing/goFlux_reprocessing/ymf_black_oak/ymf_05_flux_calculation.R",
  "data processing/goFlux_reprocessing/ymf_black_oak/ymf_06_compile_results.R",
  "data processing/goFlux_reprocessing/09_mdf_lod_comparison.R",
  "data processing/goFlux_reprocessing/diurnal_blackgum/bg_run.R")
args <- commandArgs(trailingOnly = TRUE)
if (!"--flux-only" %in% args && !file.exists("scaling/soil_jevon2023/fluxes.csv"))
  stop("Download the soil input first; see scaling/soil_jevon2023/README.md")
if (!"--flux-only" %in% args) steps <- c(steps,
  "scaling/01_flux_form_test.R", "scaling/03_tree_component.R",
  "scaling/04_main_figures.R", "scaling/08_figures_v3c.R",
  "scaling/09_tree_weighting.R", "scaling/10_basal_by_site_and_blackgum.R",
  "scaling/11_manuscript_statistics.R")
if ("--scaling-only" %in% args) steps <- steps[startsWith(steps, "scaling/")]
for (f in steps) {
  message("\nRUNNING ", f)
  sys.source(f, envir = new.env(parent = globalenv()))
}
dir.create("scaling/out", showWarnings = FALSE)
writeLines(trimws(capture.output(sessionInfo()), which = "right"), "scaling/out/sessionInfo.txt")
message("Analysis complete. Build the Word draft with ./build_manuscript.sh")
