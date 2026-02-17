# =============================================================================
# ymf_04_manual_id.R
# Interactive manual identification of measurement start/end times using
# goFlux click.peak2().
#
# *** THIS SCRIPT MUST BE RUN IN AN INTERACTIVE R SESSION (RStudio). ***
#
# How it works:
#   1. For each measurement, a pop-up plot shows the gas concentration trace
#      with a "shoulder" of extra data before and after the expected window.
#   2. Blue vertical lines show the original start and end times from the auxfile.
#   3. Click ONCE to mark the true start of the concentration change.
#   4. Click ONCE to mark the true end of the concentration change.
#   5. A validation plot briefly appears showing your selection.
#   6. The corrected times are saved automatically.
#
# Only 9 measurements — all done in a single batch.
# =============================================================================

source(file.path(
  "/Users/jongewirtzman/My Drive/Research/whole_tree_flux",
  "data processing", "goFlux_reprocessing", "ymf_black_oak", "ymf_00_setup.R"))

# --- Create output directory for click.peak plots ----------------------------
dir.create(file.path(ymf_plots_dir, "click_peak"), recursive = TRUE,
           showWarnings = FALSE)

# --- Load imported data and auxfile -------------------------------------------

load(file.path(ymf_rdata_dir, "imp_YMF_combined.RData"))
load(file.path(ymf_rdata_dir, "aux_YMF.RData"))

# =============================================================================
# YMF: 9 measurements, single batch
# =============================================================================

message("\n========================================")
message("YMF Black Oak: Creating observation windows...")
message("========================================")

# obs.length is not passed — obs.win will calculate it from start.time and end.time
ow.YMF <- obs.win(inputfile = imp_YMF, auxfile = aux.YMF, shoulder = 360)
message("YMF: ", length(ow.YMF), " observation windows created")

# --- Manual identification using CO2 trace ------------------------------------
# Identify start/end times once using CO2; the same times apply to both gases.

message("\n--- Manual identification using CO2 (9 measurements) ---")
message("Click start and end for each measurement...")
manID.YMF <- click.peak2(
  ow.list = ow.YMF,
  seq = seq(1, length(ow.YMF)),
  plot.lim = c(0, 2000),
  save.plots = file.path(ymf_plots_dir, "click_peak", "YMF")
)
save(manID.YMF, file = file.path(ymf_rdata_dir, "manID_YMF.RData"))
message("\nYMF manual ID complete: ", length(unique(manID.YMF$UniqueID)),
        " measurements saved")
message("Saved: manID_YMF.RData")
message("Proceed to ymf_05_flux_calculation.R")
