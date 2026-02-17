# =============================================================================
# ymf_03_build_auxfile.R
# Read timing and field data from Excel, match Fisher met data, and build
# the goFlux auxfile for the YMF black oak measurements.
#
# System times are used for subsetting the LGR data (start.time in auxfile).
# Real times (system − 22 min offset) are used for matching Fisher met data.
# =============================================================================

source(file.path(
  "/Users/jongewirtzman/My Drive/Research/whole_tree_flux",
  "data processing", "goFlux_reprocessing", "ymf_black_oak", "ymf_00_setup.R"))

# =============================================================================
# 1. Read Excel timing data
# =============================================================================

message("=== Reading Excel timing data ===")
tk <- read_xlsx(ymf_excel, sheet = ymf_sheet)
message("Read ", nrow(tk), " rows from '", ymf_sheet, "'")

# --- Parse times --------------------------------------------------------------
# readxl reads Excel fractional-day times as POSIXct with date 1899-12-31.
# Extract HH:MM:SS and combine with the actual measurement date (2022-10-04).

extract_time <- function(posix_col, date) {
  time_str <- format(posix_col, "%H:%M:%S")
  as.POSIXct(paste(date, time_str), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
}

tk$sys_start <- extract_time(tk$`Start time (system)`, ymf_date)
tk$sys_end   <- extract_time(tk$`end time (system)`, ymf_date)
tk$real_start <- extract_time(tk$`Start Time (real)`, ymf_date)
tk$real_end   <- extract_time(tk$`End Time (real)`, ymf_date)

# --- Build UniqueID -----------------------------------------------------------
# Heights: 10, 8, 6, 4, "4 (restarted)", 2, 1.25, 0.5, Seam

height_labels <- c("10m", "8m", "6m", "4m", "4m_restart", "2m", "1.25m", "0.5m", "Seam")
tk$UniqueID <- paste0("YMF_", height_labels)

message("\nMeasurements:")
for (i in seq_len(nrow(tk))) {
  message(sprintf("  %s: system %s - %s | real %s - %s",
                  tk$UniqueID[i],
                  format(tk$sys_start[i], "%H:%M"),
                  format(tk$sys_end[i], "%H:%M"),
                  format(tk$real_start[i], "%H:%M"),
                  format(tk$real_end[i], "%H:%M")))
}

# --- Observation length -------------------------------------------------------
# Use actual duration from system times (some measurements are 10 min, some 5 min)
tk$obs.length <- as.numeric(difftime(tk$sys_end, tk$sys_start, units = "secs"))
message("\nObs lengths (seconds): ", paste(tk$obs.length, collapse = ", "))

# =============================================================================
# 2. Download and process Harvard Forest Fisher met data
# =============================================================================

met_file <- file.path(reprocess_dir, "hf001-10-15min-m.csv")

if (!file.exists(met_file)) {
  message("\nDownloading Harvard Forest Fisher 15-min met data...")
  tryCatch({
    download.file(met_data_url, destfile = met_file, mode = "wb")
    message("Downloaded successfully to: ", met_file)
  }, error = function(e) {
    stop("Could not download met data. Please manually download:\n",
         met_data_url, "\nand save to: ", met_file)
  })
} else {
  message("\nUsing existing met data file: ", met_file)
}

met <- read.csv(met_file, stringsAsFactors = FALSE)

# Parse Fisher datetime (EST)
met$datetime_est <- lubridate::parse_date_time(met$datetime,
                                                orders = c("ymd HM", "ymd HMS",
                                                           "ymd_HM", "ymd_HMS"),
                                                tz = "EST")

# Convert EST to EDT (add 1 hour) then label as UTC to match our "fake UTC"
# timestamps. EDT = EST + 1h. Our raw data is in local EDT but stored as UTC.
met$datetime_utc <- met$datetime_est + lubridate::hours(1)
attr(met$datetime_utc, "tzone") <- "UTC"

# Filter to Oct 2022
met_filtered <- met %>%
  filter(!is.na(datetime_utc)) %>%
  filter(as.Date(datetime_utc) >= as.Date("2022-10-03") &
         as.Date(datetime_utc) <= as.Date("2022-10-05"))

message("Met data filtered to ", nrow(met_filtered), " rows for Oct 3-5, 2022")
stopifnot("airt" %in% names(met_filtered))
stopifnot("bar" %in% names(met_filtered))

# =============================================================================
# 3. Match met data to each measurement using REAL times
# =============================================================================

message("\nMatching met data to measurements using real times...")

tk <- tk %>%
  rowwise() %>%
  mutate(
    met_idx = which.min(abs(difftime(met_filtered$datetime_utc, real_start,
                                      units = "mins"))),
    Tcham = met_filtered$airt[met_idx],
    Pcham = met_filtered$bar[met_idx] / 10  # mbar -> kPa
  ) %>%
  ungroup() %>%
  select(-met_idx)

n_na_t <- sum(is.na(tk$Tcham))
n_na_p <- sum(is.na(tk$Pcham))
if (n_na_t > 0 | n_na_p > 0) {
  warning(n_na_t, " NA Tcham, ", n_na_p, " NA Pcham values")
}

# =============================================================================
# 4. Build final auxfile
# =============================================================================

# start.time and end.time = system times (for subsetting the LGR data).
# obs.win() will calculate obs.length per measurement from start.time and end.time.
aux.YMF <- tk %>%
  transmute(
    UniqueID   = UniqueID,
    start.time = sys_start,
    end.time   = sys_end,
    Vtot       = ymf_vtot,
    Area       = ymf_chamber_area,
    Tcham      = Tcham,
    Pcham      = Pcham
  )

message("\n=== Auxfile summary ===")
message("YMF: ", nrow(aux.YMF), " measurements")
message("  Vtot: ", ymf_vtot, " L (all same)")
message("  Area: ", ymf_chamber_area, " cm2 (all same)")
message("  Tcham range: ",
        round(min(aux.YMF$Tcham, na.rm = TRUE), 1), " - ",
        round(max(aux.YMF$Tcham, na.rm = TRUE), 1), " C")
message("  Pcham range: ",
        round(min(aux.YMF$Pcham, na.rm = TRUE), 1), " - ",
        round(max(aux.YMF$Pcham, na.rm = TRUE), 1), " kPa")
message("  Time range: ",
        format(min(aux.YMF$start.time), "%H:%M"), " - ",
        format(max(aux.YMF$start.time), "%H:%M"))

# =============================================================================
# 5. Save auxfile
# =============================================================================

save(aux.YMF, file = file.path(ymf_rdata_dir, "aux_YMF.RData"))
write.csv(aux.YMF, file.path(ymf_results_dir, "aux_YMF.csv"), row.names = FALSE)

# Also save the full timing key with all metadata for later merging
ymf_field_data <- tk %>%
  transmute(
    UniqueID       = UniqueID,
    Height_m       = `Height (m)`,
    Chamber        = Chamber,
    Stem_Temp_C    = `Stem Temperature`,
    Stem_Diam_mm   = `Stem Diameter`,
    Air_Temp_C     = `Air Temp`,
    Notes          = Notes,
    sys_start      = sys_start,
    sys_end        = sys_end,
    real_start     = real_start,
    real_end       = real_end,
    obs_length_sec = obs.length
  )

save(ymf_field_data, file = file.path(ymf_rdata_dir, "ymf_field_data.RData"))
write.csv(ymf_field_data, file.path(ymf_results_dir, "ymf_field_data.csv"),
          row.names = FALSE)

message("\nAuxfile saved to: ", ymf_rdata_dir)
message("CSV copies saved to: ", ymf_results_dir)
message("Proceed to ymf_04_manual_id.R")
