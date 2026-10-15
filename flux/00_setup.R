# Shared setup for the flux workflow: packages, paths, constants and helpers.
# Every flux script sources this file. Run all scripts from the repository root.
#
# Inputs (data/):
#   analyzer/<site>/            raw concentration files (ABB/LGR GLA131-GGA "f" files), as logged
#   metadata/<site>_closures.csv one row per chamber closure: chamber geometry, field times and the
#                                selected measurement window on the analyzer clock
#   meteorology/                Harvard Forest Fisher station (HF001) air temperature and pressure
# Outputs: flux/out/ (tables). Fitted models and windowed traces are cached in flux/cache/ (not tracked).

pkgs <- c("goFlux", "dplyr", "purrr", "readr", "lubridate", "tidyr", "stringr")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "), ". Restore renv.lock; see README.md.")
suppressPackageStartupMessages(for (pkg in pkgs) library(pkg, character.only = TRUE))
stopifnot(file.exists("whole_tree_flux.Rproj"))

DATA      <- "data"
META      <- file.path(DATA, "metadata")
ANALYZER  <- file.path(DATA, "analyzer")
MET_FILE  <- file.path(DATA, "meteorology", "hf001-10-15min-m_campaigns.csv")
FLUX_OUT  <- file.path("flux", "out")
CACHE     <- file.path("flux", "cache")
for (d in c(FLUX_OUT, CACHE)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# System volume = chamber + tubing (per closure, in the metadata) + analyzer cell.
analyzer_vol <- 0.028  # L, ABB GLA131-GGA microportable

# Datasheet precision of the GLA131-GGA (1 sigma at 1 s): CO2 ppm, CH4 ppb, H2O ppm. These are import
# defaults only; apply_precision() replaces CO2 and CH4 with 1.96 x the empirical group noise before fitting.
ugga_prec <- c(0.35, 0.9, 200)

# goFlux best.flux() selection criteria
flux_criteria <- c("MAE", "AICc", "g.factor", "MDF")

source(file.path("flux", "precision_helpers.R"))

read_closures <- function(site) {
  read.csv(file.path(META, paste0(site, "_closures.csv")), stringsAsFactors = FALSE, check.names = FALSE,
           na.strings = c("", "NA"))
}

# Timestamps in the metadata and analyzer files are local clock times without a time zone. They are
# read as UTC throughout so that no daylight-saving conversion is applied.
as_clock <- function(x) as.POSIXct(x, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# Import every concentration file of one analyzer, in file-name order.
import_analyzer <- function(dir) {
  files <- sort(list.files(dir, pattern = "_f[0-9]+\\.txt$", full.names = TRUE))
  stopifnot(length(files) > 0)
  bind_rows(lapply(files, function(f)
    import.UGGA(inputfile = f, date.format = "mdy", timezone = "UTC", prec = ugga_prec, save = FALSE)))
}

# Fisher station records are in fixed EST. Convert the instant to civil Eastern time and relabel that
# clock as UTC, matching the naive local clock of the field records.
read_met <- function() {
  m <- read.csv(MET_FILE, stringsAsFactors = FALSE)
  t <- parse_date_time(m$datetime, orders = c("ymd HM", "ymd HMS", "ymd_HM", "ymd_HMS"), tz = "EST")
  m$clock <- force_tz(with_tz(t, "America/New_York"), "UTC")
  m <- m[!is.na(m$clock), ]
  stopifnot(all(c("airt", "bar") %in% names(m)))
  m
}

# Air temperature (deg C) and pressure (kPa) of the station record nearest to each field start time.
match_met <- function(times, met) {
  i <- vapply(seq_along(times), function(k)
    which.min(abs(as.numeric(difftime(met$clock, times[k], units = "mins")))), integer(1))
  data.frame(Tcham = met$airt[i], Pcham = met$bar[i] / 10)
}

# Cut each closure's measurement window from the analyzer record. Timestamps are rounded to the
# nearest second (halves up); every observation inside [Window_start, Window_end] is used.
# `cl` must hold UniqueID, Window_start, Window_end, Area (cm2), Vtot (L), Tcham (deg C), Pcham (kPa).
window_traces <- function(imp, cl) {
  imp$POSIX.time <- as.POSIXct(floor(as.numeric(imp$POSIX.time) + 0.5), origin = "1970-01-01", tz = "UTC")
  s <- as_clock(cl$Window_start); e <- as_clock(cl$Window_end)
  stopifnot(!anyNA(s), !anyNA(e), all(e > s), !anyDuplicated(cl$UniqueID))
  bind_rows(lapply(seq_len(nrow(cl)), function(i) {
    w <- imp[imp$POSIX.time >= s[i] & imp$POSIX.time <= e[i], ]
    if (nrow(w) < 10) stop("Fewer than 10 analyzer records in the window of ", cl$UniqueID[i])
    w$UniqueID <- cl$UniqueID[i]
    for (nm in c("Area", "Vtot", "Tcham", "Pcham")) w[[nm]] <- cl[[nm]][i]
    w$flag <- 1
    w$Etime <- as.numeric(difftime(w$POSIX.time, s[i], units = "secs"))
    w
  }))
}

# Fit CO2 and CH4 fluxes for one analyzer record with empirical precision (see precision_helpers.R).
fit_fluxes <- function(traces, label) {
  stopifnot(all(traces$Area > 0), all(traces$Vtot > 0), all(traces$Pcham > 0))
  traces <- apply_precision(traces, label, FLUX_OUT)
  list(traces = traces,
       CO2 = best.flux(goFlux_at_interval(traces, "CO2dry_ppm"), flux_criteria),
       CH4 = best.flux(goFlux_at_interval(traces, "CH4dry_ppb"), flux_criteria))
}

# One row per closure: selected model, fit diagnostics and goFlux MDF, prefixed by gas.
flux_columns <- function(best, gas) {
  best %>% transmute(UniqueID = trimws(UniqueID), best.flux, model, quality.check, LM.flux, LM.r2, LM.p.val, LM.RMSE,
                     HM.flux, HM.r2, HM.RMSE, LM.diagnose, HM.diagnose, nb.obs, flux.term, MDF, g.fact) %>%
    rename_with(~ paste0(gas, "_", .x), -UniqueID)
}
