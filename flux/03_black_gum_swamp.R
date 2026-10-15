# Black Gum Swamp reference tree (Nyssa sylvatica in saturated peat, Harvard Forest), 28-29 August 2024:
# 105 repeated stem closures at 0.25, 1, 2, 3 and 3.6 m (three collars per height and one rigid chamber).
# Windows are the field-recorded start and end times converted to the analyzer clock (no interactive
# selection). The MDF uses the same empirical threshold as the fits.
# Output: flux/out/black_gum_swamp_fluxes.csv
source(file.path("flux", "00_setup.R"))

cl <- read_closures("black_gum_swamp")
cl$Area <- cl$Chamber_area_cm2
cl$Vtot <- cl$Chamber_volume_L + (cl$Tubing_volume_L + analyzer_vol)
cl <- cbind(cl, match_met(as_clock(cl$Field_start), read_met()))

imp <- import_analyzer(file.path(ANALYZER, "black_gum_swamp")) %>%
  distinct(POSIX.time, .keep_all = TRUE) %>% arrange(POSIX.time)
fit <- fit_fluxes(window_traces(imp, cl), "BG_LGR2")
saveRDS(list(BG = fit), file.path(CACHE, "black_gum_swamp.rds"))

noise <- read.csv(file.path(FLUX_OUT, "precision_closures_BG_LGR2.csv")) %>%
  filter(gas == "CH4") %>% select(UniqueID, field_day, dt_s, t_sec, sigma_CH4_ppb = sigma_group)
out <- cl %>% select(UniqueID, ChamberID, Position, height_m = Height_m, chamber = Chamber_type, start_field = Field_start, Notes) %>%
  left_join(fit$CH4 %>% transmute(UniqueID, CH4_best.flux = best.flux, CH4_model = model, CH4_quality.check = quality.check,
                                  CH4_flux.term = flux.term, nb.obs, CH4_LM.r2 = LM.r2), by = "UniqueID") %>%
  left_join(fit$CO2 %>% transmute(UniqueID, CO2_best.flux = best.flux, CO2_model = model), by = "UniqueID") %>%
  left_join(noise, by = "UniqueID") %>%
  mutate(CH4_MDF_emp95 = qnorm(.975) * sigma_CH4_ppb / t_sec * CH4_flux.term,
         CH4_below_MDF = abs(CH4_best.flux) < CH4_MDF_emp95)
stopifnot(nrow(out) == nrow(cl), all(is.finite(out$CH4_best.flux)))
write.csv(out, file.path(FLUX_OUT, "black_gum_swamp_fluxes.csv"), row.names = FALSE)
message("Black Gum Swamp: ", nrow(out), " closures, ", sum(out$CH4_below_MDF), " below MDF -> ",
        file.path(FLUX_OUT, "black_gum_swamp_fluxes.csv"))
