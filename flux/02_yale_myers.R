# Yale Myers Forest black oak (Quercus velutina), 4 October 2022: nine stem closures at 0.5-10 m.
# The analyzer clock ran 22 min ahead of local time; windows are on the analyzer clock and weather is
# matched at the field (local) start time.
# Output: flux/out/yale_myers_fluxes_fitted.csv (detection limits are added by 04_detection_limits.R).
source(file.path("flux", "00_setup.R"))

cl <- read_closures("yale_myers")
cl$Area <- cl$Chamber_area_cm2
cl$Vtot <- cl$Chamber_volume_L + (cl$Tubing_volume_L + analyzer_vol)
cl <- cbind(cl, match_met(as_clock(cl$Field_start), read_met()))

imp <- import_analyzer(file.path(ANALYZER, "yale_myers"))
fit <- fit_fluxes(window_traces(imp, cl), "YMF")
fit$imp <- imp
saveRDS(list(YMF = fit), file.path(CACHE, "yale_myers.rds"))

out <- cl %>% select(-Area, -Vtot, -Tcham, -Pcham) %>%
  left_join(flux_columns(fit$CO2, "CO2"), by = "UniqueID") %>%
  left_join(flux_columns(fit$CH4, "CH4"), by = "UniqueID")
stopifnot(nrow(out) == nrow(cl), all(is.finite(out$CH4_best.flux)))
write.csv(out, file.path(FLUX_OUT, "yale_myers_fluxes_fitted.csv"), row.names = FALSE)
message("Yale Myers: ", nrow(out), " closures -> ", file.path(FLUX_OUT, "yale_myers_fluxes_fitted.csv"))
