# =============================================================================
# 00_setup.R
# Packages, paths, and constants for goFlux reprocessing workflow
# Source this script at the top of every other script.
# =============================================================================

# --- Package Installation and Loading ----------------------------------------

pkgs <- c("goFlux", "dplyr", "purrr", "readr", "lubridate", "openxlsx", "tidyr", "stringr")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "), ". Restore renv.lock; see README.md.")
for (pkg in pkgs) library(pkg, character.only = TRUE)

# --- Path Definitions --------------------------------------------------------

base_dir <- normalizePath(".", mustWork = TRUE)
stopifnot(file.exists(file.path(base_dir, "whole_tree_flux.Rproj")))
input_dir <- file.path(base_dir, "data processing", "input")
reprocess_dir <- file.path(base_dir, "data processing", "goFlux_reprocessing")

# Original raw data directories (read-only)
lgr1_raw <- file.path(input_dir, "LGR1")
lgr2_raw <- file.path(input_dir, "LGR2")
lgr3_raw <- file.path(input_dir, "LGR3")

# Timing key CSVs
tk_lgr1 <- file.path(input_dir, "times_key_tree - Canopy Lift_LGR1 (2).csv")
tk_lgr2 <- file.path(input_dir, "times_key_tree - Canopy Lift_LGR2 (2).csv")
tk_lgr3 <- file.path(input_dir, "times_key_tree - Canopy Lift_LGR3 (2).csv")

# Field data entry (for final merge)
field_data_path <- file.path(base_dir, "data processing",
                             "Field Data Entry - Clean Canopy Lift Total.csv")

# Working/output directories (created by scripts as needed)
goflux_import_dir <- file.path(reprocess_dir, "import")
rdata_dir         <- file.path(reprocess_dir, "RData")
results_dir       <- file.path(reprocess_dir, "results")
plots_dir         <- file.path(reprocess_dir, "plots")

# --- Constants ---------------------------------------------------------------

# Volume addition: analyzer internal volume (0.028 L) + tubing (0.029 L)
# Lab convention shared across projects (see ch4-data-filtering 00_setup.R):
# 0.028 L analyzer internal volume for both the ABB GLA131-GGA microportable
# and the LI-7810.
# NOTE (2026-09): earlier runs used 0.070 L, which is goFlux's example-auxfile
# value for the larger LGR UGGA, not the GLA131 microportable units used here.
analyzer_vol  <- 0.028  # Liters
tubing_vol    <- 0.029  # Liters
vtot_addition <- analyzer_vol + tubing_vol  # 0.057 L

# goFlux instrument precision for the ABB/LGR GLA131-GGA Microportable UGGA
# (1σ at 1 s; ABB GLA131-GGA datasheet 3KXG167001R1001 Rev. J)
# c(CO2dry_ppm, CH4dry_ppb, H2O_ppm)
# Precision also controls HM curvature constraints and model selection. Pass it
# as import defaults only. Before fitting, apply_precision() replaces the CH4
# and CO2 fields with 1.96 x empirical group sigma; H2O keeps its import default.
ugga_prec <- c(0.35, 0.9, 200)

# Date format in raw LGR files (mm/dd/yyyy)
lgr_date_format <- "mdy"

# Default observation length (seconds) for obs.win
obs_length <- 180

# best.flux selection criteria
flux_criteria <- c("MAE", "AICc", "g.factor", "MDF")

# Harvard Forest Fisher met data URL (15-minute intervals, metric)
met_data_url <- "https://harvardforest1.fas.harvard.edu/data/p00/hf001/hf001-10-15min-m.csv"

# --- Create output directories -----------------------------------------------

for (d in c(goflux_import_dir, rdata_dir, results_dir, plots_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

message("Setup complete. Base directory: ", base_dir)

# Refresh metadata without changing manually selected measurement windows.
# Called by every flux-calculation step, so rebuilding auxfiles cannot leave stale
# areas, volumes, weather or precision in the cached manID objects.
sync_manID <- function(manID, aux) {
  stopifnot(!anyDuplicated(aux$UniqueID))
  i <- match(manID$UniqueID, aux$UniqueID)
  if (anyNA(i)) stop("Unmatched manID UniqueID in auxfile")
  for (nm in c("Area", "Vtot", "Tcham", "Pcham")) {
    stopifnot(all(is.finite(aux[[nm]])))
    manID[[nm]] <- aux[[nm]][i]
  }
  stopifnot(all(manID$Area > 0), all(manID$Vtot > 0), all(manID$Pcham > 0), all(manID$Tcham > -273.15))
  manID$CO2_prec <- ugga_prec[1]
  manID$CH4_prec <- ugga_prec[2]
  manID$H2O_prec <- ugga_prec[3]
  manID
}
