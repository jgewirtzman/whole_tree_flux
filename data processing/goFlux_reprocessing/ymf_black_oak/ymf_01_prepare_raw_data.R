# =============================================================================
# ymf_01_prepare_raw_data.R
# Stage LGR f-files into a flat directory for goFlux import.
# =============================================================================

source(file.path(
  "/Users/jongewirtzman/My Drive/Research/whole_tree_flux",
  "data processing", "goFlux_reprocessing", "ymf_black_oak", "ymf_00_setup.R"))

# --- Stage f-files ------------------------------------------------------------

staging_dir <- file.path(ymf_import_dir, "YMF")
dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)

# Find all f-files (flame data) > 0 bytes
all_files <- list.files(ymf_lgr_raw, pattern = "_f\\d+\\.txt$",
                        full.names = TRUE, recursive = TRUE)
f_files <- all_files[file.size(all_files) > 0]

message("Found ", length(f_files), " non-empty f-files in ", ymf_lgr_raw)

# Copy to flat staging directory
n_copied <- 0
for (f in f_files) {
  dest <- file.path(staging_dir, basename(f))
  file.copy(f, dest, overwrite = TRUE)
  n_copied <- n_copied + 1
}

message("Staged ", n_copied, " f-files to ", staging_dir)

# Verify
staged <- list.files(staging_dir, pattern = "\\.txt$")
message("Files in staging directory:")
for (s in staged) message("  ", s)
