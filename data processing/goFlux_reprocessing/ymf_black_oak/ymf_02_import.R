# =============================================================================
# ymf_02_import.R
# Import staged LGR files using goFlux import.UGGA(), then combine.
#
# Note: import2RData() throws "no restart 'muffleError' found" with these
# files due to a progress-bar bug in goFlux v0.2.0. We call import.UGGA()
# directly instead, which works correctly.
# =============================================================================

source(file.path(
  "/Users/jongewirtzman/My Drive/Research/whole_tree_flux",
  "data processing", "goFlux_reprocessing", "ymf_black_oak", "ymf_00_setup.R"))

# --- Import with import.UGGA() directly --------------------------------------

staging_dir <- file.path(ymf_import_dir, "YMF")
txt_files <- list.files(staging_dir, pattern = "[.]txt$", full.names = TRUE)

message("=== Importing YMF LGR data ===")
message("Source: ", staging_dir)
message("Files to import: ", length(txt_files))

imp_list <- list()
for (f in txt_files) {
  message("  Importing: ", basename(f))
  imp_list[[basename(f)]] <- import.UGGA(
    inputfile   = f,
    date.format = lgr_date_format,
    timezone    = "UTC"
  )
  message("    -> ", nrow(imp_list[[basename(f)]]), " rows")
}

# --- Combine all imported data ------------------------------------------------

imp_YMF <- bind_rows(imp_list)

message("\nCombined data: ", nrow(imp_YMF), " rows, ",
        ncol(imp_YMF), " columns")
message("Time range: ", min(imp_YMF$POSIX.time), " to ", max(imp_YMF$POSIX.time))

# Save combined
save(imp_YMF, file = file.path(ymf_rdata_dir, "imp_YMF_combined.RData"))
message("Saved: imp_YMF_combined.RData")
