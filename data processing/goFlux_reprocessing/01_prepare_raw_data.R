# =============================================================================
# 01_prepare_raw_data.R
# Stage raw LGR f-files into flat directories for goFlux import.
#
# Problem: goFlux import2RData() for UGGA uses list.files() without
# recursive=TRUE, so it only finds .txt files in the top-level path dir.
# Raw data is organized in date subdirectories and also contains non-data
# files (b, l, p, lims). Some zip extractions created directories ending
# in .txt that contain the actual file inside.
#
# Solution: Recursively find only the concentration data f-files (> 0 bytes)
# and copy them into flat staging directories.
# =============================================================================

# Source setup (works both when source()'d and run interactively in RStudio)
setup_path <- file.path(
  "/Users/jongewirtzman/My Drive/Research/whole_tree_flux",
  "data processing", "goFlux_reprocessing", "00_setup.R")
source(setup_path)

# --- Helper function ---------------------------------------------------------

stage_f_files <- function(src_dir, dest_dir) {
  # Create destination

  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

  # Find ALL .txt files recursively (list.files only returns actual files,

  # not directories, even if dir name ends in .txt)
  all_txt <- list.files(src_dir, pattern = "\\.txt$",
                        recursive = TRUE, full.names = TRUE)

  # Keep only f-files: filename contains _f followed by digits before .txt
  f_files <- all_txt[grepl("_f\\d+\\.txt$", basename(all_txt))]

  # Exclude empty files (0 bytes)
  f_files <- f_files[file.size(f_files) > 0]

  if (length(f_files) == 0) {
    warning("No f-files found in: ", src_dir)
    return(invisible(NULL))
  }

  # Copy each file to the flat staging directory
  copied <- 0
  for (f in f_files) {
    dest_file <- file.path(dest_dir, basename(f))
    if (file.exists(dest_file)) {
      message("  Skipping duplicate: ", basename(f))
    } else {
      file.copy(f, dest_file)
      copied <- copied + 1
    }
  }

  n_staged <- length(list.files(dest_dir, pattern = "\\.txt$"))
  message("Staged ", n_staged, " f-files to ", dest_dir,
          " (", copied, " newly copied)")
}

# --- Stage files for each instrument -----------------------------------------

message("=== Staging LGR1 raw f-files ===")
stage_f_files(lgr1_raw, file.path(goflux_import_dir, "LGR1"))

message("=== Staging LGR2 raw f-files ===")
stage_f_files(lgr2_raw, file.path(goflux_import_dir, "LGR2"))

message("=== Staging LGR3 raw f-files ===")
stage_f_files(lgr3_raw, file.path(goflux_import_dir, "LGR3"))

# --- Verification ------------------------------------------------------------

for (lgr in c("LGR1", "LGR2", "LGR3")) {
  staged <- list.files(file.path(goflux_import_dir, lgr), pattern = "\\.txt$")
  message(lgr, ": ", length(staged), " files staged -> ",
          paste(staged, collapse = ", "))
}

message("\nFile staging complete. Proceed to 02_import.R")
