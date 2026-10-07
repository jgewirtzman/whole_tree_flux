# =============================================================================
# 03c_patch_manID_vtot.R
# Patch the existing manID RData files so that Vtot matches the auxfiles
# rebuilt by 03_build_auxfiles.R.
#
# This avoids re-running the interactive 04_manual_id.R (click.peak2) when
# vtot_addition (analyzer + tubing volume, 00_setup.R) changes. Only the Vtot
# column is replaced; the manually identified start/end times are untouched.
#
# Run this AFTER 03_build_auxfiles.R and BEFORE 05_flux_calculation.R.
# =============================================================================

# Source setup
setup_path <- file.path(
  ".",
  "data processing", "goFlux_reprocessing", "00_setup.R")
source(setup_path)

# --- Patch function ----------------------------------------------------------

patch_vtot <- function(manID, aux, label) {
  idx <- match(manID$UniqueID, aux$UniqueID)
  if (any(is.na(idx))) {
    stop(label, ": UniqueIDs in manID with no auxfile match: ",
         paste(unique(manID$UniqueID[is.na(idx)]), collapse = ", "))
  }
  old_vtot <- manID$Vtot
  manID$Vtot <- aux$Vtot[idx]

  delta <- unique(round(manID$Vtot - old_vtot, 6))
  message(label, ": ", length(unique(manID$UniqueID)), " measurements; Vtot ",
          round(min(old_vtot), 3), "-", round(max(old_vtot), 3), " L → ",
          round(min(manID$Vtot), 3), "-", round(max(manID$Vtot), 3),
          " L (change: ", paste(delta, collapse = ", "), " L)")
  manID
}

# --- Load, patch, and save ---------------------------------------------------

load(file.path(rdata_dir, "aux_LGR1.RData"))
load(file.path(rdata_dir, "aux_LGR2.RData"))
load(file.path(rdata_dir, "aux_LGR3.RData"))

load(file.path(rdata_dir, "manID_LGR1.RData"))
manID.LGR1 <- patch_vtot(manID.LGR1, aux.LGR1, "LGR1")
save(manID.LGR1, file = file.path(rdata_dir, "manID_LGR1.RData"))

load(file.path(rdata_dir, "manID_LGR2.RData"))
manID.LGR2 <- patch_vtot(manID.LGR2, aux.LGR2, "LGR2")
save(manID.LGR2, file = file.path(rdata_dir, "manID_LGR2.RData"))

load(file.path(rdata_dir, "manID_LGR3.RData"))
manID.LGR3 <- patch_vtot(manID.LGR3, aux.LGR3, "LGR3")
save(manID.LGR3, file = file.path(rdata_dir, "manID_LGR3.RData"))

message("\n=== Vtot patching complete ===")
message("Updated manID files saved to: ", rdata_dir)
message("Proceed to 05_flux_calculation.R")
