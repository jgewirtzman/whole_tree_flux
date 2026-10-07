# Read-only checks used before the runner can change outputs.
check_workflow <- function(mode) {
  stopifnot(mode %in% c("full", "flux", "scaling"))
  pkgs <- c("goFlux", "dplyr", "purrr", "readr", "lubridate", "openxlsx", "tidyr", "stringr", "readxl", "jsonlite")
  if (mode != "flux") pkgs <- c(pkgs, "ggplot2", "patchwork", "scales", "jpeg", "lme4", "Matrix")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "), ". Restore renv.lock as described in README.md.")
  lock <- jsonlite::fromJSON("renv.lock", simplifyVector = FALSE)
  mismatch <- pkgs[vapply(pkgs, function(p) packageVersion(p) != package_version(lock$Packages[[p]]$Version), logical(1))]
  if (length(mismatch)) stop("Package versions differ from renv.lock: ", paste(mismatch, collapse = ", "), ". Restore the recorded environment first.")
  sha <- packageDescription("goFlux")$RemoteSha
  if (is.null(sha) || sha != lock$Packages$goFlux$RemoteSha) stop("goFlux commit differs from renv.lock; restore the recorded environment.")
  if (as.character(getRversion()) != lock$R$Version) warning("Recorded R version is ", lock$R$Version, "; current version is ", getRversion(), ".")
  nboot <- suppressWarnings(as.numeric(Sys.getenv("FLUX_N_BOOT", "2000")))
  if (!is.finite(nboot) || nboot < 100 || nboot != floor(nboot)) stop("FLUX_N_BOOT must be an integer of at least 100; reported results use 2000.")
  rd <- "data processing/goFlux_reprocessing"
  inputs <- character()
  if (mode != "scaling") {
    inputs <- c("data processing/Field Data Entry - Clean Canopy Lift Total.csv", "data processing/leaf_areas.csv",
      "data processing/YMF Black Oak/Black Oak Tree Project Data.xlsx",
      paste0("data processing/input/times_key_tree - Canopy Lift_LGR", 1:3, " (2).csv"),
      file.path(rd, "hf001-10-15min-m.csv"),
      file.path(rd, "RData", paste0("manID_LGR", 1:3, ".RData")),
      file.path(rd, "RData", paste0("imp_LGR", 1:3, "_combined.RData")),
      file.path(rd, "ymf_black_oak/RData", c("manID_YMF.RData", "imp_YMF_combined.RData")),
      file.path(rd, "diurnal_blackgum/raw/diurnal/diurnal_final", c("diurnal_updated.csv", "diurnal_volumes.csv")))
    raw <- list.files(file.path(rd, "diurnal_blackgum/raw/diurnal/diurnal_final/LGR2"),
      pattern = "_f[0-9]+\\.txt(\\.zip)?$", recursive = TRUE, full.names = TRUE)
    raw <- raw[!dir.exists(raw) & file.size(raw) > 0]
    if (!length(raw)) stop("Missing blackgum raw analyzer files.")
    inputs <- c(inputs, raw)
  }
  if (mode != "flux") {
    inputs <- c(inputs, "scaling/soil_jevon2023/fluxes.csv", "IMG_5926_edited.jpg", "IMG_6437.jpg")
    if (mode == "scaling") inputs <- c(inputs,
      file.path(rd, "results", c("canopy_flux_goFlux_compiled.csv", "canopy_flux_goFlux_compiled_with_mdf.csv")),
      file.path(rd, "ymf_black_oak/results", c("ymf_black_oak_flux_compiled.csv", "ymf_black_oak_flux_compiled_with_mdf.csv", "ymf_field_data.csv")),
      file.path(rd, "diurnal_blackgum/results/blackgum_flux_compiled_with_mdf.csv"))
  }
  missing <- inputs[!file.exists(inputs)]
  if (length(missing)) stop("Missing required input(s):\n", paste(missing, collapse = "\n"),
    "\nFor the separately distributed soil data, see scaling/soil_jevon2023/README.md.")
  if (mode != "flux") {
    soil <- read.csv("scaling/soil_jevon2023/fluxes.csv")
    if (!all(c("date", "CH4.flux") %in% names(soil))) stop("Soil input must contain date and CH4.flux columns.")
  }
  message("Preflight passed: ", mode, " mode; recorded package versions and required inputs present.")
  unique(inputs)
}
