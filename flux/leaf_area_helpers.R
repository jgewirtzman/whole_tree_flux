# Leaf areas for the Harvard Forest foliage closures. Each closure is linked to a shoot
# (leaf_area_mapping.csv); a shoot uses its measured area when the sample is identified
# (leaf_areas.csv) and otherwise the mean measured area of its species. Closures of one
# shoot share an area.
leaf_area_assignments <- function(meta = file.path("data", "metadata")) {
  f <- read.csv(file.path(meta, "harvard_forest_closures.csv"), stringsAsFactors = FALSE)
  a <- read.csv(file.path(meta, "leaf_areas.csv"), stringsAsFactors = FALSE)
  map <- read.csv(file.path(meta, "leaf_area_mapping.csv"), stringsAsFactors = FALSE)
  stopifnot(!anyDuplicated(f$UniqueID), !anyDuplicated(map$UniqueID),
            !anyDuplicated(paste(a$Tree, a$Flux_Rep)), all(is.finite(a$Area_cm2)), all(a$Area_cm2 > 0))
  l <- f[grepl("leaf", f$Type, ignore.case = TRUE), ]
  stopifnot(setequal(map$UniqueID, l$UniqueID), all(nzchar(map$ShootID)))
  map <- map[match(l$UniqueID, map$UniqueID), ]
  sp <- c(bg = "Blackgum", hem = "Hemlock", rm = "Maple", ro = "Oak")
  species <- unname(sp[l$Species]); stopifnot(!anyNA(species), all(species == map$Tree))
  av <- tapply(a$Area_cm2, a$Tree, mean)
  area <- unname(av[species]); method <- rep("species_mean", nrow(l))
  for (i in which(!is.na(map$Flux_Rep))) {
    k <- which(a$Tree == map$Tree[i] & a$Flux_Rep == map$Flux_Rep[i])
    stopifnot(length(k) == 1L)
    area[i] <- a$Area_cm2[k]; method[i] <- "measured"
  }
  stopifnot(all(is.finite(area)), all(area > 0))
  for (s in unique(map$ShootID)) {
    j <- which(map$ShootID == s)
    if (length(unique(paste(species[j], map$Flux_Rep[j], area[j], method[j]))) != 1L)
      stop("Inconsistent area assignment within shoot: ", s)
  }
  data.frame(UniqueID = l$UniqueID, ShootID = map$ShootID,
             Leaf_area_cm2 = area, Leaf_area_source = method, Light_condition = map$Light_condition)
}
