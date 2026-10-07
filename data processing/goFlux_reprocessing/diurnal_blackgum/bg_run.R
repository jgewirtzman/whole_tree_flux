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
# MDF: empirical 95 % reference method of 09_mdf_lod_comparison.R: MDF = 1.96 x sigma / t x flux.term, sigma = campaign
#   MAD of first differences / sqrt(2) (per run of constant logging interval).
# Run: LANG=en_US.UTF-8 Rscript "data processing/goFlux_reprocessing/diurnal_blackgum/bg_run.R"  (from whole_tree_flux/)
# =============================================================================
suppressPackageStartupMessages({library(goFlux); library(dplyr); library(readr); library(lubridate); library(tidyr)})
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
  mutate(t = parse_date_time(datetime, c("ymd HM", "ymd_HM", "Ymd HM", "Ymd HMS"), tz = "EST") + hours(1)) %>% filter(!is.na(t)) %>%
  filter(t >= as.POSIXct("2024-08-27", tz = "EST"), t <= as.POSIXct("2024-08-31", tz = "EST"))
attr(met$t, "tzone") <- "UTC"
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

## 5. fluxes
CH4 <- best.flux(goFlux(manID, "CH4dry_ppb"), flux_criteria); CO2 <- best.flux(goFlux(manID, "CO2dry_ppm"), flux_criteria)

## 6. empirical MDF (09_mdf_lod_comparison.R reference method)
# sigma: the canopy campaigns used the campaign MAD of first differences of the whole record, but here closures with
# concentration rises of up to ~50 ppb s-1 dominate the record, so raw first differences measure the trend, not noise
# (whole-record value reported for comparison). Instead: within each closure, first differences of the residuals from a
# linear fit (trend removed), MAD / sqrt(2). Campaign sigma = median across the low-flux closures (>= 3 m), because
# curvature in high-flux closures leaves structure in linear residuals; the all-closure median is reported too.
dtv <- diff(as.numeric(imp$POSIX.time)); runs <- abs(dtv - median(dtv)) < 0.2
sigma_all <- mad(diff(imp$CH4dry_ppb)[runs], constant = 1.4826) / sqrt(2)
sig_cl <- manID %>% filter(flag == 1) %>% group_by(UniqueID) %>%
  summarise(s = mad(diff(resid(lm(CH4dry_ppb ~ Etime))), constant = 1.4826) / sqrt(2), .groups = "drop")
hts <- setNames(fd$`Height (m)`, fd$UniqueID); sigma_allcl <- median(sig_cl$s)
sigma_ch4 <- median(sig_cl$s[hts[sig_cl$UniqueID] >= 3])
cat("sigma CH4: whole-record first differences", round(sigma_all, 3), "ppb; detrended, all closures", round(sigma_allcl, 3), "ppb; detrended, closures >= 3 m (used)", round(sigma_ch4, 3),
    "ppb (closure range", paste(round(range(sig_cl$s), 2), collapse = "-"), ")\n")
out <- fd %>% select(UniqueID, VolumeID, Position, height_m = `Height (m)`, chamber = `Chamber Type`, start_field = `Start (iPad)`, Notes) %>%
  left_join(CH4 %>% transmute(UniqueID, CH4_best.flux = best.flux, CH4_model = model, CH4_quality.check = quality.check,
                              CH4_flux.term = flux.term, nb.obs, CH4_LM.r2 = LM.r2), by = "UniqueID") %>%
  left_join(CO2 %>% transmute(UniqueID, CO2_best.flux = best.flux, CO2_model = model), by = "UniqueID") %>%
  left_join(manID %>% filter(flag == 1) %>% group_by(UniqueID) %>% summarise(t_sec = as.numeric(diff(range(POSIX.time))), .groups = "drop"), by = "UniqueID") %>%
  mutate(CH4_MDF_emp95 = 1.96 * sigma_ch4 / t_sec * CH4_flux.term, CH4_below_MDF = abs(CH4_best.flux) < CH4_MDF_emp95,
         sigma_CH4_ppb = sigma_ch4)
write_csv(out, file.path(RES, "blackgum_flux_compiled_with_mdf.csv")); save(manID, CH4, CO2, auxfile, file = file.path(RES, "blackgum_goflux.RData"))
cat("sigma CH4 (ppb):", round(sigma_ch4, 3), "| fluxes:", sum(!is.na(out$CH4_best.flux)), "| below MDF:", sum(out$CH4_below_MDF, na.rm = TRUE), "\n")
print(out %>% group_by(height_m) %>% summarise(n = n(), mean = mean(CH4_best.flux, na.rm = TRUE), median = median(CH4_best.flux, na.rm = TRUE),
  min = min(CH4_best.flux, na.rm = TRUE), neg = sum(CH4_best.flux < 0, na.rm = TRUE), below_mdf = sum(CH4_below_MDF, na.rm = TRUE)))
