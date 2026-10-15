# Harvard Forest canopy-lift campaign (July-August 2023; three analyzers, 136 closures):
# stem, branch and leaf CO2 and CH4 fluxes from the raw analyzer files and closure metadata.
# Output: flux/out/harvard_forest_fluxes_fitted.csv (detection limits are added by 04_detection_limits.R).
source(file.path("flux", "00_setup.R"))
source(file.path("flux", "leaf_area_helpers.R"))

cl <- read_closures("harvard_forest")
stopifnot(!anyDuplicated(cl$UniqueID), all(cl$Analyzer %in% c("LGR1", "LGR2", "LGR3")))
leaf <- leaf_area_assignments()
i <- match(cl$UniqueID, leaf$UniqueID)
stopifnot(identical(!is.na(i), grepl("leaf", cl$Type, ignore.case = TRUE)))
# Chamber area: enclosed bark area for stems and branches, measured or species-mean leaf area for foliage.
cl$Area <- ifelse(is.na(i), cl$Chamber_area_cm2, leaf$Leaf_area_cm2[i])
cl$Vtot <- cl$Chamber_volume_L + (cl$Tubing_volume_L + analyzer_vol)
cl <- cbind(cl, match_met(as_clock(cl$Field_start), read_met()))

fits <- list()
for (lgr in c("LGR1", "LGR2", "LGR3")) {
  imp <- import_analyzer(file.path(ANALYZER, "harvard_forest", lgr))
  fits[[lgr]] <- fit_fluxes(window_traces(imp, cl[cl$Analyzer == lgr, ]), lgr)
  fits[[lgr]]$imp <- imp
  message(lgr, ": ", nrow(fits[[lgr]]$CH4), " closures fitted")
}
saveRDS(fits, file.path(CACHE, "harvard_forest.rds"))

out <- cl %>% select(-Area, -Vtot, -Tcham, -Pcham) %>%
  left_join(bind_rows(lapply(fits, function(f) flux_columns(f$CO2, "CO2"))), by = "UniqueID") %>%
  left_join(bind_rows(lapply(fits, function(f) flux_columns(f$CH4, "CH4"))), by = "UniqueID") %>%
  left_join(leaf, by = "UniqueID")
stopifnot(nrow(out) == nrow(cl), all(is.finite(out$CH4_best.flux)), all(is.finite(out$CO2_best.flux)))
write.csv(out, file.path(FLUX_OUT, "harvard_forest_fluxes_fitted.csv"), row.names = FALSE)
message("Harvard Forest: ", nrow(out), " closures -> ", file.path(FLUX_OUT, "harvard_forest_fluxes_fitted.csv"))
