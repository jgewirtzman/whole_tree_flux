#!/usr/bin/env Rscript
# Reproduce the analysis from the repository root; retain selected chamber windows.
args <- commandArgs(trailingOnly = TRUE)
allowed <- c("--help", "--check", "--flux-only", "--scaling-only")
if (any(!args %in% allowed)) stop("Unknown option: ", paste(setdiff(args, allowed), collapse = ", "), ". Use --help.")
if (all(c("--flux-only", "--scaling-only") %in% args)) stop("Choose only one of --flux-only and --scaling-only.")
if ("--help" %in% args) {
  cat("Usage: Rscript run_analysis.R [--flux-only | --scaling-only] [--check]\n",
      "  no options      rebuild fluxes, analyses, figures and summaries; run checks\n",
      "  --flux-only     recalculate fluxes and detection flags from saved windows\n",
      "  --scaling-only  use committed flux tables to rebuild analyses and figures\n",
      "  --check         check dependencies and inputs without writing outputs\n",
      "See README.md and docs/REPRODUCING.md for environment setup and soil data.\n", sep = "")
  quit(status = 0)
}
if (!file.exists("whole_tree_flux.Rproj")) stop("Run from the whole_tree_flux repository root.")
source("scripts/workflow.R")
mode <- if ("--flux-only" %in% args) "flux" else if ("--scaling-only" %in% args) "scaling" else "full"
inputs <- check_workflow(mode)
if ("--check" %in% args) quit(status = 0)
sys.source("tests/precision_checks.R", envir = new.env(parent = globalenv()))
input_md5 <- tools::md5sum(inputs)
Sys.setenv(FLUX_LEGACY_QC = "0")
flux_steps <- c(
  "data processing/goFlux_reprocessing/03_build_auxfiles.R",
  "data processing/goFlux_reprocessing/ymf_black_oak/ymf_03_build_auxfile.R",
  "data processing/goFlux_reprocessing/05_flux_calculation.R",
  "data processing/goFlux_reprocessing/06_compile_results.R",
  "data processing/goFlux_reprocessing/ymf_black_oak/ymf_05_flux_calculation.R",
  "data processing/goFlux_reprocessing/ymf_black_oak/ymf_06_compile_results.R",
  "data processing/goFlux_reprocessing/09_mdf_lod_comparison.R",
  "data processing/goFlux_reprocessing/diurnal_blackgum/bg_run.R")
scaling_steps <- c("scaling/01_flux_form_test.R", "scaling/03_tree_component.R",
  "scaling/04_main_figures.R", "scaling/08_figures_v3c.R",
  "scaling/09_tree_weighting.R", "scaling/10_basal_by_site_and_blackgum.R",
  "scaling/11_summary_statistics.R")
steps <- switch(mode, flux = flux_steps, scaling = scaling_steps, full = c(flux_steps, scaling_steps))
for (i in seq_along(steps)) {
  message(sprintf("\n[%d/%d] %s", i, length(steps), steps[i]))
  sys.source(steps[i], envir = new.env(parent = globalenv()))
}
if (mode != "flux") sys.source("tests/analysis_checks.R", envir = new.env(parent = globalenv()))
dir.create("scaling/out", showWarnings = FALSE)
writeLines(trimws(capture.output(sessionInfo()), which = "right"), "scaling/out/sessionInfo.txt")
write.csv(data.frame(path = names(input_md5), md5 = unname(input_md5)), "scaling/out/run_inputs.csv", row.names = FALSE)
writeLines(c(paste("mode:", mode), paste("completed_utc:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("bootstrap_draws:", Sys.getenv("FLUX_N_BOOT", "2000"))), "scaling/out/run_metadata.txt")
message("Analysis complete. Outputs: scaling/out/ and scaling/v3c/. See README.md for the figure map.")
