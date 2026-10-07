# =============================================================================
# bg_run.R — Black gum (Nyssa sylvatica), Black Gum Swamp next to Swamp Rd, Harvard Forest (fully saturated peat swamp):
# intensive stem CH4/CO2 fluxes, 28-29 Aug 2024, heights 0.25, 1, 2, 3, 3.6 m (three collars per height plus a
# Tupperware chamber). Reprocessed with the same goFlux workflow as the canopy and YMF data, except:
#   - measurement windows come from the field sheet's start/end times on the analyzer clock (diurnal_updated.csv:
#     updated_start/updated_end; 2 of 105 starts were adjusted in the field record) instead of interactive click.peak2
#     (Jon, 2026-10-06: "start with a"; click.peak can be revisited). No additional trimming.
#   - raw files are the LGR2 f-files in the folder (the field sheet lists LGR3; Jon: use the files in the folder).
# Volume convention (lab convention, 2026-09-30): Vtot = chamber volume (diurnal_volumes.csv, which already includes cap
#   and tubing: "Two.tube.Volume" 0.063 L) + 0.028 L analyzer cell. Area = surfarea (m2) x 1e4 cm2.
# Tcham/Pcham: Fisher station (HF001, hf001-10-15min-m.csv; EST -> EDT) at the field (iPad) time.
# MDF: empirical 95 % reference method of 09_mdf_lod_comparison.R: MDF = 1.96 x sigma / t x flux.term, sigma = median closure second-difference
#   precision per analyzer x field day x logging interval (precision_helpers.R).
# Run: LANG=en_US.UTF-8 Rscript "data processing/goFlux_reprocessing/diurnal_blackgum/bg_run.R"  (from whole_tree_flux/)
# =============================================================================
suppressPackageStartupMessages({library(goFlux); library(dplyr); library(readr); library(lubridate); library(tidyr)})
source("data processing/goFlux_reprocessing/precision_helpers.R")
BG <- "data processing/goFlux_reprocessing/diurnal_blackgum"   # run from the repository root
RAW <- file.path(BG, "raw/diurnal/diurnal_final"); STAGE <- file.path(BG, "import"); RES <- file.path(BG, "results")
for (d in c(STAGE, RES)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
analyzer_vol <- 0.028; ugga_prec <- c(0.35, 0.9, 200); flux_criteria <- c("MAE", "AICc", "g.factor", "MDF")

## 1. stage raw f-files (unzip nested .txt.zip; skip empty files and folders)
fz <- list.files(file.path(RAW, "LGR2"), pattern = "_f[0-9]+\\.txt(\\.zip)?$", recursive = TRUE, full.names = TRUE)
fz <- fz[!dir.exists(fz) & file.size(fz) > 0]
for (f in fz) {
  if (grepl("\\.zip$", f)) { tmp <- tempfile(); unzip(f, exdir = tmp); g <- list.files(tmp, pattern = "\\.txt$", recursive = TRUE, full.names = TRUE)
    for (x in g) if (file.size(x) > 0) file.copy(x, file.path(STAGE, basename(x)), overwrite = TRUE)
  } else file.copy(f, file.path(STAGE, basename(f)), overwrite = TRUE)
}
txt <- list.files(STAGE, pattern = "\\.txt$", full.names = TRUE); cat("staged files:", length(txt), "\n")

## 2. import
imp <- bind_rows(lapply(txt, function(f) import.UGGA(inputfile = f, date.format = "mdy", timezone = "UTC", prec = ugga_prec))) %>%
  distinct(POSIX.time, .keep_all = TRUE) %>% arrange(POSIX.time)
cat("imported rows:", nrow(imp), "range:", format(range(imp$POSIX.time)), "\n")

## 3. auxfile
fd <- read.csv(file.path(RAW, "diurnal_updated.csv"), check.names = FALSE, stringsAsFactors = FALSE)
vol <- read.csv(file.path(RAW, "diurnal_volumes.csv"), stringsAsFactors = FALSE) %>% transmute(VolumeID, Vch = Volume, Area_m2 = surfarea)
met <- read.csv("data processing/goFlux_reprocessing/hf001-10-15min-m.csv", stringsAsFactors = FALSE) %>%
  mutate(t = parse_date_time(datetime, c("ymd HM", "ymd_HM", "Ymd HM", "Ymd HMS"), tz = "EST")) %>% filter(!is.na(t)) %>%
  filter(t >= as.POSIXct("2024-08-27", tz = "EST"), t <= as.POSIXct("2024-08-31", tz = "EST"))
met$t <- force_tz(with_tz(met$t, "America/New_York"), "UTC")
aux <- fd %>% left_join(vol, by = "VolumeID") %>%
  mutate(start.time = as.POSIXct(updated_start, tz = "UTC"), end.time = as.POSIXct(updated_end, tz = "UTC"),
         field_time = as.POSIXct(`Start (iPad)`, tz = "UTC"))
stopifnot(!any(is.na(aux$Vch)), !any(is.na(aux$start.time)), !any(is.na(aux$end.time)))
aux$Tcham <- sapply(aux$field_time, function(x) met$airt[which.min(abs(as.numeric(difftime(met$t, x, units = "mins"))))])
aux$Pcham <- sapply(aux$field_time, function(x) met$bar[which.min(abs(as.numeric(difftime(met$t, x, units = "mins"))))]) / 10
auxfile <- aux %>% transmute(UniqueID, start.time, end.time, Vtot = Vch + analyzer_vol, Area = Area_m2 * 1e4, Tcham, Pcham)
cat("aux:", nrow(auxfile), "| Vtot", paste(round(range(auxfile$Vtot), 3), collapse = "-"), "L | Area", paste(round(unique(auxfile$Area), 2), collapse = ", "),
    "cm2 | Tcham", paste(round(range(auxfile$Tcham), 1), collapse = "-"), "| Pcham", paste(round(range(auxfile$Pcham), 1), collapse = "-"), "\n")

## 4. windows from the field times (non-interactive replacement for click.peak2)
ow <- obs.win(inputfile = imp, auxfile = auxfile, shoulder = 0)
manID <- bind_rows(lapply(ow, function(w) {
  w %>% mutate(flag = as.numeric(POSIX.time >= start.time & POSIX.time <= end.time),
               start.time_corr = start.time, end.time_corr = end.time,
               obs.length_corr = as.numeric(difftime(end.time, start.time, units = "secs")),
               Etime = as.numeric(difftime(POSIX.time, start.time, units = "secs")))
}))
nwin <- manID %>% group_by(UniqueID) %>% summarise(n = sum(flag == 1), .groups = "drop")
cat("windows with data:", sum(nwin$n > 10), "of", nrow(auxfile), "| points per window", paste(range(nwin$n), collapse = "-"), "\n")
miss <- setdiff(auxfile$UniqueID, nwin$UniqueID[nwin$n > 10]); if (length(miss)) cat("WINDOWS WITHOUT DATA:", paste(miss, collapse = ", "), "\n")

## 5. Empirical precision and fluxes, using every closure at each analyzer-day interval
manID <- apply_precision(manID, "BG_LGR2", RES)
CH4 <- best.flux(goFlux_at_interval(manID, "CH4dry_ppb"), flux_criteria)
CO2 <- best.flux(goFlux_at_interval(manID, "CO2dry_ppm"), flux_criteria)
noise <- read.csv(file.path(RES, "precision_closures_BG_LGR2.csv")) %>%
  filter(gas == "CH4") %>% select(UniqueID, field_day, dt_s, t_sec, sigma_CH4_ppb = sigma_group)

## 6. Empirical MDF, identical to the threshold used inside the fits
out <- fd %>% select(UniqueID, VolumeID, Position, height_m = `Height (m)`, chamber = `Chamber Type`, start_field = `Start (iPad)`, Notes) %>%
  left_join(CH4 %>% transmute(UniqueID, CH4_best.flux = best.flux, CH4_model = model, CH4_quality.check = quality.check,
                              CH4_flux.term = flux.term, nb.obs, CH4_LM.r2 = LM.r2), by = "UniqueID") %>%
  left_join(CO2 %>% transmute(UniqueID, CO2_best.flux = best.flux, CO2_model = model), by = "UniqueID") %>%
  left_join(noise, by = "UniqueID") %>%
  mutate(CH4_MDF_emp95 = qnorm(.975) * sigma_CH4_ppb / t_sec * CH4_flux.term,
         CH4_below_MDF = abs(CH4_best.flux) < CH4_MDF_emp95)
write_csv(out, file.path(RES, "blackgum_flux_compiled_with_mdf.csv")); save(manID, CH4, CO2, auxfile, file = file.path(RES, "blackgum_goflux.RData"))
cat("sigma CH4 (ppb):", paste(round(unique(out$sigma_CH4_ppb), 3), collapse = ", "), "| fluxes:", sum(!is.na(out$CH4_best.flux)), "| below MDF:", sum(out$CH4_below_MDF, na.rm = TRUE), "\n")
print(out %>% group_by(height_m) %>% summarise(n = n(), mean = mean(CH4_best.flux, na.rm = TRUE), median = median(CH4_best.flux, na.rm = TRUE),
  min = min(CH4_best.flux, na.rm = TRUE), neg = sum(CH4_best.flux < 0, na.rm = TRUE), below_mdf = sum(CH4_below_MDF, na.rm = TRUE)))
