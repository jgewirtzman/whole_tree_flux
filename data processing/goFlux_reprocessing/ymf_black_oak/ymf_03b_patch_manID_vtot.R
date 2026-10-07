# =============================================================================
# ymf_03b_patch_manID_vtot.R
# Patch the existing YMF manID RData file so that Vtot matches the auxfile
# rebuilt by ymf_03_build_auxfile.R.
#
# This avoids re-running the interactive ymf_04_manual_id.R (click.peak2) when
# vtot_addition (analyzer + tubing volume, 00_setup.R) changes. Only the Vtot
# column is replaced; the manually identified start/end times are untouched.
#
# Run this AFTER ymf_03_build_auxfile.R and BEFORE ymf_05_flux_calculation.R.
# =============================================================================

source(file.path(
  ".",
  "data processing", "goFlux_reprocessing", "ymf_black_oak", "ymf_00_setup.R"))

load(file.path(ymf_rdata_dir, "aux_YMF.RData"))
load(file.path(ymf_rdata_dir, "manID_YMF.RData"))

idx <- match(manID.YMF$UniqueID, aux.YMF$UniqueID)
if (any(is.na(idx))) {
  stop("UniqueIDs in manID.YMF with no auxfile match: ",
       paste(unique(manID.YMF$UniqueID[is.na(idx)]), collapse = ", "))
}

old_vtot <- unique(manID.YMF$Vtot)
manID.YMF$Vtot <- aux.YMF$Vtot[idx]
message("YMF: ", length(unique(manID.YMF$UniqueID)), " measurements; Vtot ",
        paste(old_vtot, collapse = ", "), " L → ",
        paste(unique(manID.YMF$Vtot), collapse = ", "), " L")

save(manID.YMF, file = file.path(ymf_rdata_dir, "manID_YMF.RData"))
message("Saved: manID_YMF.RData")
message("Proceed to ymf_05_flux_calculation.R")
