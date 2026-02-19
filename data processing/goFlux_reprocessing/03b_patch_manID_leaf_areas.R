# =============================================================================
# 03b_patch_manID_leaf_areas.R
# Patch the existing manID RData files to replace chamber-based surface areas
# with actual measured leaf areas for leaf-type measurements.
#
# This avoids re-running the interactive 04_manual_id.R (click.peak2).
# It loads each manID file, updates the Area column for leaf rows using the
# same lookup as 03_build_auxfiles.R, and re-saves.
#
# Run this AFTER 03_build_auxfiles.R and BEFORE 05_flux_calculation.R.
# =============================================================================

# Source setup
setup_path <- file.path(
  "/Users/jongewirtzman/My Drive/Research/whole_tree_flux",
  "data processing", "goFlux_reprocessing", "00_setup.R")
source(setup_path)

# --- Load measured leaf areas ------------------------------------------------

leaf_areas_path <- file.path(base_dir, "data processing", "leaf_areas.csv")
leaf_areas <- read.csv(leaf_areas_path, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
# Normalise column names (handle BOM and spaces)
names(leaf_areas) <- trimws(gsub("[^A-Za-z0-9_]", ".", names(leaf_areas)))
names(leaf_areas) <- sub("^\\.+", "", names(leaf_areas))
leaf_areas$Tree <- trimws(leaf_areas$Tree)

# Build explicit UniqueID -> measured area (cm²) lookup for August trees
leaf_area_lookup <- c(
  "321902 11 leaf 1" = leaf_areas$Area[leaf_areas$Tree == "Hemlock" & leaf_areas$Flux.Rep == 1],
  "321902 15 leaf 1" = leaf_areas$Area[leaf_areas$Tree == "Hemlock" & leaf_areas$Flux.Rep == 2],
  "321902 6 leaf 1"  = leaf_areas$Area[leaf_areas$Tree == "Hemlock" & leaf_areas$Flux.Rep == 3],
  "321071 16 leaf 1" = leaf_areas$Area[leaf_areas$Tree == "Maple"   & leaf_areas$Flux.Rep == 1],
  "321071 18 leaf 2" = leaf_areas$Area[leaf_areas$Tree == "Maple"   & leaf_areas$Flux.Rep == 2],
  "321071 22 leaf 1" = leaf_areas$Area[leaf_areas$Tree == "Maple"   & leaf_areas$Flux.Rep == 3],
  "300607 16 leaf 1" = leaf_areas$Area[leaf_areas$Tree == "Oak"     & leaf_areas$Flux.Rep == 1],
  "300607 12 leaf 1" = leaf_areas$Area[leaf_areas$Tree == "Oak"     & leaf_areas$Flux.Rep == 2],
  "300607 8 leaf 1"  = leaf_areas$Area[leaf_areas$Tree == "Oak"     & leaf_areas$Flux.Rep == 3]
)

# Species averages for July trees without measured leaf areas
hem_avg <- mean(leaf_areas$Area[leaf_areas$Tree == "Hemlock"])
maple_avg <- mean(leaf_areas$Area[leaf_areas$Tree == "Maple"])
overall_avg <- mean(leaf_areas$Area)

message("Leaf area values for patching:")
message("  Direct lookups: ", length(leaf_area_lookup), " UniqueIDs")
message("  Hemlock avg: ", round(hem_avg, 1), " cm²")
message("  Maple avg:   ", round(maple_avg, 1), " cm²")
message("  Overall avg:  ", round(overall_avg, 1), " cm²")

# --- Patch function ----------------------------------------------------------

# Map Tree_tag (extracted from UniqueID) to species average
tag_to_area <- function(uid) {
  tag <- as.integer(sub(" .*", "", uid))
  if (tag == 4) return(hem_avg)       # Tree 4 = hemlock (Swamp Rd)
  if (tag == 5) return(maple_avg)     # Tree 5 = red maple (Swamp Rd)
  if (tag == 2) return(overall_avg)   # Tree 2 = black gum (no species data)
  return(NA_real_)
}

patch_manID <- function(manID, label) {
  leaf_rows <- grepl("leaf", manID$UniqueID)
  n_leaf <- sum(leaf_rows)
  if (n_leaf == 0) {
    message(label, ": No leaf measurements found, skipping.")
    return(manID)
  }

  leaf_uids <- unique(manID$UniqueID[leaf_rows])
  message("\n", label, ": Patching ", length(leaf_uids),
          " leaf UniqueIDs (", n_leaf, " rows total)")

  for (uid in leaf_uids) {
    uid_rows <- manID$UniqueID == uid
    old_area <- manID$Area[uid_rows][1]

    if (uid %in% names(leaf_area_lookup)) {
      new_area <- leaf_area_lookup[uid]
    } else {
      new_area <- tag_to_area(uid)
    }

    if (!is.na(new_area)) {
      manID$Area[uid_rows] <- new_area
      message("  ", uid, ": ", round(old_area, 1), " → ", round(new_area, 1), " cm²")
    } else {
      warning("  ", uid, ": No mapping found, keeping ", round(old_area, 1), " cm²")
    }
  }

  manID
}

# --- Load, patch, and save ---------------------------------------------------

# LGR1 (July: Tree 2 bg, Tree 4 hem)
message("\n=== LGR1 ===")
load(file.path(rdata_dir, "manID_LGR1.RData"))
manID.LGR1 <- patch_manID(manID.LGR1, "LGR1")
save(manID.LGR1, file = file.path(rdata_dir, "manID_LGR1.RData"))

# LGR2 (July: Tree 4 hem, Tree 5 rm)
message("\n=== LGR2 ===")
load(file.path(rdata_dir, "manID_LGR2.RData"))
manID.LGR2 <- patch_manID(manID.LGR2, "LGR2")
save(manID.LGR2, file = file.path(rdata_dir, "manID_LGR2.RData"))

# LGR3 (August: all EMS trees with direct lookups)
message("\n=== LGR3 ===")
load(file.path(rdata_dir, "manID_LGR3.RData"))
manID.LGR3 <- patch_manID(manID.LGR3, "LGR3")
save(manID.LGR3, file = file.path(rdata_dir, "manID_LGR3.RData"))

message("\n=== Patching complete ===")
message("Updated manID files saved to: ", rdata_dir)
message("Proceed to 05_flux_calculation.R")
