# =============================================================================
# 00_setup.R
# Packages, paths, and constants for goFlux reprocessing workflow
# Source this script at the top of every other script.
# =============================================================================

# --- Package Installation and Loading ----------------------------------------

if (!require("remotes", quietly = TRUE)) install.packages("remotes")

# Install goFlux from GitHub if not already installed
if (!require("goFlux", quietly = TRUE)) {
  remotes::install_github("Qepanna/goFlux")
}
library(goFlux)

# CRAN packages
pkgs <- c("dplyr", "purrr", "readr", "lubridate", "openxlsx", "tidyr", "stringr")
for (pkg in pkgs) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# --- Path Definitions --------------------------------------------------------

base_dir <- "/Users/jongewirtzman/My Drive/Research/whole_tree_flux"
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

# Volume addition: analyzer cell (0.070 L) + tubing (0.029 L)
vtot_addition <- 0.099  # Liters

# goFlux instrument precision for UGGA (GLA132 series)
# c(CO2dry_ppm, CH4dry_ppb, H2O_ppm)
ugga_prec <- c(0.2, 1.4, 50)

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
