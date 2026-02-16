# =============================================================================
# 04_manual_id.R
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
# Tips:
#   - If you make a mistake, you can re-run a specific batch.
#   - Results are saved after each instrument, so you won't lose work.
#   - Batches are kept to <=18 measurements to avoid fatigue errors.
#   - If a measurement looks bad (no clear signal), click start and end
#     at roughly the same location -- you can flag it later.
# =============================================================================

# Source setup (works both when source()'d and run interactively in RStudio)
setup_path <- file.path(
  "/Users/jongewirtzman/My Drive/Research/whole_tree_flux",
  "data processing", "goFlux_reprocessing", "00_setup.R")
source(setup_path)

# --- Create output directory for click.peak plots ----------------------------
dir.create(file.path(plots_dir, "click_peak"), recursive = TRUE,
           showWarnings = FALSE)

# --- Load imported data and auxfiles ------------------------------------------

load(file.path(rdata_dir, "imp_LGR1_combined.RData"))
load(file.path(rdata_dir, "imp_LGR2_combined.RData"))
load(file.path(rdata_dir, "imp_LGR3_combined.RData"))
load(file.path(rdata_dir, "aux_LGR1.RData"))
load(file.path(rdata_dir, "aux_LGR2.RData"))
load(file.path(rdata_dir, "aux_LGR3.RData"))

# =============================================================================
# LGR1: 26 measurements, 2 batches of 13
# =============================================================================

message("\n========================================")
message("LGR1: Creating observation windows...")
message("========================================")
ow.LGR1 <- obs.win(inputfile = imp.LGR1, auxfile = aux.LGR1,
                     obs.length = obs_length, shoulder = 360)
message("LGR1: ", length(ow.LGR1), " observation windows created")

message("\nLGR1 Batch 1 of 2 (measurements 1-13)")
message("Click start and end for each measurement...")
manID.LGR1_a <- click.peak2(
  ow.list = ow.LGR1,
  seq = seq(1, 13),
  plot.lim = c(0, 2000),
  save.plots = file.path(plots_dir, "click_peak", "LGR1_batch1")
)

message("\nLGR1 Batch 2 of 2 (measurements 14-26)")
manID.LGR1_b <- click.peak2(
  ow.list = ow.LGR1,
  seq = seq(14, length(ow.LGR1)),
  plot.lim = c(0, 2000),
  save.plots = file.path(plots_dir, "click_peak", "LGR1_batch2")
)

manID.LGR1 <- rbind(manID.LGR1_a, manID.LGR1_b)
save(manID.LGR1, file = file.path(rdata_dir, "manID_LGR1.RData"))
message("LGR1 manual ID complete: ", length(unique(manID.LGR1$UniqueID)),
        " measurements saved")

# =============================================================================
# LGR2: 41 measurements, 3 batches (14, 14, 13)
# =============================================================================

message("\n========================================")
message("LGR2: Creating observation windows...")
message("========================================")
ow.LGR2 <- obs.win(inputfile = imp.LGR2, auxfile = aux.LGR2,
                     obs.length = obs_length, shoulder = 360)
message("LGR2: ", length(ow.LGR2), " observation windows created")

message("\nLGR2 Batch 1 of 3 (measurements 1-14)")
manID.LGR2_a <- click.peak2(
  ow.list = ow.LGR2,
  seq = seq(1, 14),
  plot.lim = c(0, 2000),
  save.plots = file.path(plots_dir, "click_peak", "LGR2_batch1")
)

message("\nLGR2 Batch 2 of 3 (measurements 15-28)")
manID.LGR2_b <- click.peak2(
  ow.list = ow.LGR2,
  seq = seq(15, 28),
  plot.lim = c(0, 2000),
  save.plots = file.path(plots_dir, "click_peak", "LGR2_batch2")
)

message("\nLGR2 Batch 3 of 3 (measurements 29-41)")
manID.LGR2_c <- click.peak2(
  ow.list = ow.LGR2,
  seq = seq(29, length(ow.LGR2)),
  plot.lim = c(0, 2000),
  save.plots = file.path(plots_dir, "click_peak", "LGR2_batch3")
)

manID.LGR2 <- rbind(manID.LGR2_a, manID.LGR2_b, manID.LGR2_c)
save(manID.LGR2, file = file.path(rdata_dir, "manID_LGR2.RData"))
message("LGR2 manual ID complete: ", length(unique(manID.LGR2$UniqueID)),
        " measurements saved")

# =============================================================================
# LGR3: 69 measurements, 4 batches (18, 18, 18, 15)
# =============================================================================

message("\n========================================")
message("LGR3: Creating observation windows...")
message("========================================")
ow.LGR3 <- obs.win(inputfile = imp.LGR3, auxfile = aux.LGR3,
                     obs.length = obs_length, shoulder = 360)
message("LGR3: ", length(ow.LGR3), " observation windows created")

message("\nLGR3 Batch 1 of 4 (measurements 1-18)")
manID.LGR3_a <- click.peak2(
  ow.list = ow.LGR3,
  seq = seq(1, 18),
  plot.lim = c(0, 2000),
  save.plots = file.path(plots_dir, "click_peak", "LGR3_batch1")
)

message("\nLGR3 Batch 2 of 4 (measurements 19-36)")
manID.LGR3_b <- click.peak2(
  ow.list = ow.LGR3,
  seq = seq(19, 36),
  plot.lim = c(0, 2000),
  save.plots = file.path(plots_dir, "click_peak", "LGR3_batch2")
)

message("\nLGR3 Batch 3 of 4 (measurements 37-54)")
manID.LGR3_c <- click.peak2(
  ow.list = ow.LGR3,
  seq = seq(37, 54),
  plot.lim = c(0, 2000),
  save.plots = file.path(plots_dir, "click_peak", "LGR3_batch3")
)

message("\nLGR3 Batch 4 of 4 (measurements 55-69)")
manID.LGR3_d <- click.peak2(
  ow.list = ow.LGR3,
  seq = seq(55, length(ow.LGR3)),
  plot.lim = c(0, 2000),
  save.plots = file.path(plots_dir, "click_peak", "LGR3_batch4")
)

manID.LGR3 <- rbind(manID.LGR3_a, manID.LGR3_b, manID.LGR3_c, manID.LGR3_d)
save(manID.LGR3, file = file.path(rdata_dir, "manID_LGR3.RData"))
message("LGR3 manual ID complete: ", length(unique(manID.LGR3$UniqueID)),
        " measurements saved")

# =============================================================================
# Summary
# =============================================================================

message("\n========================================")
message("Manual identification complete!")
message("========================================")
message("LGR1: ", length(unique(manID.LGR1$UniqueID)), " measurements")
message("LGR2: ", length(unique(manID.LGR2$UniqueID)), " measurements")
message("LGR3: ", length(unique(manID.LGR3$UniqueID)), " measurements")
message("Total: ", length(unique(manID.LGR1$UniqueID)) +
                   length(unique(manID.LGR2$UniqueID)) +
                   length(unique(manID.LGR3$UniqueID)), " measurements")
message("\nClick peak plots saved to: ", file.path(plots_dir, "click_peak"))
message("Proceed to 05_flux_calculation.R")
