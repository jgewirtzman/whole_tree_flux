# =============================================================================
# 02_import.R
# Import raw LGR UGGA data using goFlux import2RData().
#
# import2RData() saves RData files into a "RData/" folder relative to the
# working directory (not the path argument). We setwd() to per-instrument
# subdirectories so that each instrument's RData is kept separate.
#
# After import, we load and combine all RData per instrument, then save
# combined files for use in subsequent scripts.
# =============================================================================

# Source setup (works both when source()'d and run interactively in RStudio)
setup_path <- file.path(
  "/Users/jongewirtzman/My Drive/Research/whole_tree_flux",
  "data processing", "goFlux_reprocessing", "00_setup.R")
source(setup_path)

original_wd <- getwd()

# --- Create per-instrument working directories --------------------------------

for (lgr in c("LGR1", "LGR2", "LGR3")) {
  dir.create(file.path(rdata_dir, lgr), recursive = TRUE, showWarnings = FALSE)
}

# --- Import LGR1 -------------------------------------------------------------
message("\n=== Importing LGR1 ===")
setwd(file.path(rdata_dir, "LGR1"))
import2RData(
  path = file.path(goflux_import_dir, "LGR1"),
  instrument = "UGGA",
  date.format = lgr_date_format,
  prec = ugga_prec
)

# --- Import LGR2 -------------------------------------------------------------
message("\n=== Importing LGR2 ===")
setwd(file.path(rdata_dir, "LGR2"))
import2RData(
  path = file.path(goflux_import_dir, "LGR2"),
  instrument = "UGGA",
  date.format = lgr_date_format,
  prec = ugga_prec
)

# --- Import LGR3 -------------------------------------------------------------
message("\n=== Importing LGR3 ===")
setwd(file.path(rdata_dir, "LGR3"))
import2RData(
  path = file.path(goflux_import_dir, "LGR3"),
  instrument = "UGGA",
  date.format = lgr_date_format,
  prec = ugga_prec
)

# Restore working directory
setwd(original_wd)

# --- Load and combine imported data per instrument ----------------------------

message("\n=== Loading and combining imported data ===")

# LGR1
lgr1_files <- list.files(file.path(rdata_dir, "LGR1", "RData"),
                          pattern = "imp\\.RData$", full.names = TRUE)
imp.LGR1 <- lgr1_files %>% map_df(~ get(load(.x)))
message("LGR1: ", nrow(imp.LGR1), " rows from ", length(lgr1_files), " files")
message("  Date range: ", min(imp.LGR1$POSIX.time), " to ", max(imp.LGR1$POSIX.time))

# LGR2
lgr2_files <- list.files(file.path(rdata_dir, "LGR2", "RData"),
                          pattern = "imp\\.RData$", full.names = TRUE)
imp.LGR2 <- lgr2_files %>% map_df(~ get(load(.x)))
message("LGR2: ", nrow(imp.LGR2), " rows from ", length(lgr2_files), " files")
message("  Date range: ", min(imp.LGR2$POSIX.time), " to ", max(imp.LGR2$POSIX.time))

# LGR3
lgr3_files <- list.files(file.path(rdata_dir, "LGR3", "RData"),
                          pattern = "imp\\.RData$", full.names = TRUE)
imp.LGR3 <- lgr3_files %>% map_df(~ get(load(.x)))
message("LGR3: ", nrow(imp.LGR3), " rows from ", length(lgr3_files), " files")
message("  Date range: ", min(imp.LGR3$POSIX.time), " to ", max(imp.LGR3$POSIX.time))

# --- Save combined per-instrument data ----------------------------------------

save(imp.LGR1, file = file.path(rdata_dir, "imp_LGR1_combined.RData"))
save(imp.LGR2, file = file.path(rdata_dir, "imp_LGR2_combined.RData"))
save(imp.LGR3, file = file.path(rdata_dir, "imp_LGR3_combined.RData"))

message("\nImport complete. Combined RData saved to: ", rdata_dir)
message("Proceed to 03_build_auxfiles.R")
