# =============================================================================
# ymf_00_setup.R
# Packages, paths, and constants for YMF Black Oak goFlux reprocessing.
# Source this script at the top of every other ymf_* script.
# =============================================================================

# --- Load shared packages and constants from parent setup --------------------
parent_setup <- file.path(
  "/Users/jongewirtzman/My Drive/Research/whole_tree_flux",
  "data processing", "goFlux_reprocessing", "00_setup.R")
source(parent_setup)

# Also need readxl for Excel input
if (!require("readxl", quietly = TRUE)) install.packages("readxl")
library(readxl)

# --- YMF-specific paths ------------------------------------------------------

ymf_dir       <- file.path(base_dir, "data processing", "YMF Black Oak")
ymf_reprocess <- file.path(reprocess_dir, "ymf_black_oak")

# Raw LGR data
ymf_lgr_raw <- file.path(ymf_dir, "2022-10-04")

# Excel file with timing and field metadata
ymf_excel <- file.path(ymf_dir, "Black Oak Tree Project Data.xlsx")
ymf_sheet <- "Tree Stem Fluxes"

# Working/output directories
ymf_import_dir  <- file.path(ymf_reprocess, "import")
ymf_rdata_dir   <- file.path(ymf_reprocess, "RData")
ymf_results_dir <- file.path(ymf_reprocess, "results")
ymf_plots_dir   <- file.path(ymf_reprocess, "plots")

# --- YMF-specific constants ---------------------------------------------------

# Chamber: Large stem (s6)
ymf_chamber_vol  <- 2.30    # Liters
ymf_chamber_area <- 446     # cm² (0.0446 m²)
ymf_vtot         <- ymf_chamber_vol + vtot_addition  # 2.30 + 0.099 = 2.399 L

# Measurement date
ymf_date <- as.Date("2022-10-04")

# Time offset: system time = real time + 22 min
ymf_time_offset_min <- 22

# --- Create output directories -----------------------------------------------

for (d in c(ymf_import_dir, ymf_rdata_dir, ymf_results_dir, ymf_plots_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

message("YMF setup complete.")
message("  Excel: ", ymf_excel)
message("  LGR raw: ", ymf_lgr_raw)
message("  Output: ", ymf_reprocess)
