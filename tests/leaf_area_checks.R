# Leaf-area assignments: measured shoots, within-species means and shared shoots.
source("flux/leaf_area_helpers.R")
a <- leaf_area_assignments()
i <- match(c("2 10.5 leaf 1", "2 10.5 leaf 2", "2 11.7 leaf 1",
             "2 11.7 leaf 2", "2 11.7 leaf_s 1", "2 9.1 leaf 1", "2 7 leaf 1"), a$UniqueID)
stopifnot(!anyNA(i), max(abs(a$Leaf_area_cm2[i] - c(286.4, 334.2, 286.4, 286.4, 286.4, 228.8, 286.4))) < 1e-8,
          a$ShootID[i[1]] != a$ShootID[i[2]], a$ShootID[i[4]] == a$ShootID[i[5]])
# The flux table uses these areas for every foliage closure.
d <- read.csv("flux/out/harvard_forest_fluxes.csv")
stopifnot(max(abs(d$Leaf_area_cm2[match(a$UniqueID, d$UniqueID)] - a$Leaf_area_cm2)) < 1e-8)
m <- read.csv("data/metadata/harvard_forest_closures.csv")
stopifnot(m$Type[m$UniqueID == "2 11.7 stem 1"] == "stem", m$Type[m$UniqueID == "2 11.7 leaf 1"] == "leaf")
cat("PASS: tissue and height metadata, within-species leaf areas and shared-shoot assignments\n")
