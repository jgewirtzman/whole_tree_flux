# Read-only checks used before the runner can change outputs.
check_workflow <- function(mode) {
  stopifnot(mode %in% c("full", "flux", "scaling"))
  pkgs <- c("goFlux", "dplyr", "purrr", "readr", "lubridate", "tidyr", "stringr", "jsonlite")
  if (mode != "flux") pkgs <- c(pkgs, "ggplot2", "patchwork", "scales", "jpeg", "png", "lme4", "Matrix")
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
  meta <- file.path("data/metadata", c("harvard_forest_closures.csv", "yale_myers_closures.csv",
    "black_gum_swamp_closures.csv", "black_gum_swamp_chambers.csv", "leaf_areas.csv", "leaf_area_mapping.csv"))
  inputs <- meta
  if (mode != "scaling") {
    raw <- list.files("data/analyzer", pattern = "_f[0-9]+\\.txt$", recursive = TRUE, full.names = TRUE)
    for (site in c("harvard_forest/LGR1", "harvard_forest/LGR2", "harvard_forest/LGR3", "yale_myers", "black_gum_swamp"))
      if (!any(startsWith(raw, file.path("data/analyzer", site)))) stop("Missing analyzer files for ", site)
    inputs <- c(inputs, "data/meteorology/hf001-10-15min-m_campaigns.csv", raw)
  }
  if (mode != "flux") {
    inputs <- c(inputs, "data/soil_jevon2023/fluxes.csv",
      file.path("data/photos/field", c("IMG_5926_edited.jpg", "IMG_6437.jpg")),
      file.path("data/photos/chambers", c("IMG_5018.png", "IMG_5103.png", "IMG_5115.png", "IMG_5186.png")))
    if (mode == "scaling") inputs <- c(inputs,
      file.path("flux/out", c("harvard_forest_fluxes.csv", "yale_myers_fluxes.csv", "black_gum_swamp_fluxes.csv")))
  }
  missing <- inputs[!file.exists(inputs)]
  if (length(missing)) stop("Missing required input(s):\n", paste(missing, collapse = "\n"),
    "\nFor the separately distributed soil data, see data/soil_jevon2023/README.md.")
  source("flux/leaf_area_helpers.R")
  leaf_area_assignments()
  if (mode != "flux") {
    soil <- read.csv("data/soil_jevon2023/fluxes.csv")
    if (!all(c("date", "CH4.flux") %in% names(soil))) stop("Soil input must contain date and CH4.flux columns.")
  }
  message("Preflight passed: ", mode, " mode; recorded package versions and required inputs present.")
  unique(inputs)
}
