# =============================================================================
# 03_build_auxfiles.R
# Build goFlux auxiliary files from timing key CSVs and Harvard Forest
# Fisher meteorological station data (hf001, 15-min intervals).
#
# goFlux auxfile required columns:
#   UniqueID, start.time (POSIXct), Area (cm2), Vtot (L), Tcham (C), Pcham (kPa)
#
# Tcham and Pcham come from Fisher station 15-min data:
#   airt = air temperature (C), bar = barometric pressure (millibar)
#   Pcham = bar / 10  (mbar -> kPa)
#
# Timezone: All times are treated as UTC consistently (goFlux default).
# Raw LGR timestamps are local EDT but imported as UTC. Auxfile start.times
# are also built in UTC (treating local time as UTC). Fisher met data is in
# EST; we add 1 hour to convert EST->EDT, matching the raw data's local clock.
# =============================================================================

# Source setup (works both when source()'d and run interactively in RStudio)
setup_path <- file.path(
  "/Users/jongewirtzman/My Drive/Research/whole_tree_flux",
  "data processing", "goFlux_reprocessing", "00_setup.R")
source(setup_path)

# =============================================================================
# 1. Read and parse timing key CSVs
# =============================================================================

tk1 <- read.csv(tk_lgr1, stringsAsFactors = FALSE)
tk2 <- read.csv(tk_lgr2, stringsAsFactors = FALSE)
tk3 <- read.csv(tk_lgr3, stringsAsFactors = FALSE)

# Parse dates (all are mdy format; LGR3 uses 2-digit year which mdy() handles)
tk1$date_parsed <- lubridate::mdy(tk1$Date)
tk2$date_parsed <- lubridate::mdy(tk2$Date)
tk3$date_parsed <- lubridate::mdy(tk3$Date)

# Build start.time as POSIXct in UTC (treating local time as UTC)
# Handle comp_start_time with possibly missing seconds (e.g., "10:43" vs "10:43:00")
build_start_time <- function(df) {
  df %>%
    mutate(
      time_parsed = lubridate::parse_date_time(comp_start_time, orders = c("HMS", "HM")),
      start.time = as.POSIXct(
        paste(date_parsed, format(time_parsed, "%H:%M:%S")),
        format = "%Y-%m-%d %H:%M:%S",
        tz = "UTC"
      )
    ) %>%
    select(-time_parsed)
}

tk1 <- build_start_time(tk1)
tk2 <- build_start_time(tk2)
tk3 <- build_start_time(tk3)

# Check for any failed start.time parsing
for (nm in c("tk1", "tk2", "tk3")) {
  df <- get(nm)
  n_na <- sum(is.na(df$start.time))
  if (n_na > 0) {
    warning(nm, ": ", n_na, " rows with NA start.time:\n",
            paste(df$UniqueID[is.na(df$start.time)], collapse = ", "))
  }
}

# =============================================================================
# 2. Calculate Vtot and Area
# =============================================================================

calc_vtot_area <- function(df) {
  df %>%
    mutate(
      Vtot = Rep_vol_L + vtot_addition,
      Area = surface_area_m2 * 10000  # m2 -> cm2
    )
}

tk1 <- calc_vtot_area(tk1)
tk2 <- calc_vtot_area(tk2)
tk3 <- calc_vtot_area(tk3)

# =============================================================================
# 3. Download and process Harvard Forest Fisher met data
# =============================================================================

met_file <- file.path(reprocess_dir, "hf001-10-15min-m.csv")

if (!file.exists(met_file)) {
  message("Downloading Harvard Forest Fisher 15-min met data...")
  tryCatch({
    download.file(met_data_url, destfile = met_file, mode = "wb")
    message("Downloaded successfully to: ", met_file)
  }, error = function(e) {
    stop("Could not download met data. Please manually download:\n",
         met_data_url, "\nand save to: ", met_file)
  })
} else {
  message("Using existing met data file: ", met_file)
}

met <- read.csv(met_file, stringsAsFactors = FALSE)

# Parse Fisher datetime. Fisher timestamps are in EST (Eastern Standard Time).
# The format is typically "2023-07-18T12:00" or "2023-07-18 12:00"
# Try multiple formats
met$datetime_est <- lubridate::parse_date_time(met$datetime,
                                                orders = c("ymd HM", "ymd HMS",
                                                           "ymd_HM", "ymd_HMS"),
                                                tz = "EST")

# Convert EST to EDT (add 1 hour) then label as UTC to match our "fake UTC"
# timestamps. EDT = EST + 1h. Our raw data is in local EDT but stored as UTC.
met$datetime_utc <- met$datetime_est + lubridate::hours(1)
attr(met$datetime_utc, "tzone") <- "UTC"

# Filter to relevant date ranges to speed up matching
# LGR1/2: July 18-19, 2023; LGR3: August 16-17, 2023
relevant_dates <- as.Date(c("2023-07-17", "2023-07-20",
                             "2023-08-15", "2023-08-18"))
met_filtered <- met %>%
  filter(!is.na(datetime_utc)) %>%
  filter(
    (as.Date(datetime_utc) >= relevant_dates[1] & as.Date(datetime_utc) <= relevant_dates[2]) |
    (as.Date(datetime_utc) >= relevant_dates[3] & as.Date(datetime_utc) <= relevant_dates[4])
  )

message("Met data filtered to ", nrow(met_filtered), " rows for relevant dates")
message("Met data columns: ", paste(names(met_filtered), collapse = ", "))

# Verify expected columns exist
stopifnot("airt" %in% names(met_filtered))
stopifnot("bar" %in% names(met_filtered))

# =============================================================================
# 4. Match met data to each measurement
# =============================================================================

match_met <- function(df, met_data) {
  df %>%
    rowwise() %>%
    mutate(
      met_idx = which.min(abs(difftime(met_data$datetime_utc, start.time,
                                        units = "mins"))),
      Tcham = met_data$airt[met_idx],
      Pcham = met_data$bar[met_idx] / 10  # mbar -> kPa
    ) %>%
    ungroup() %>%
    select(-met_idx)
}

message("Matching met data to LGR1 measurements...")
tk1 <- match_met(tk1, met_filtered)

message("Matching met data to LGR2 measurements...")
tk2 <- match_met(tk2, met_filtered)

message("Matching met data to LGR3 measurements...")
tk3 <- match_met(tk3, met_filtered)

# Handle any NA Tcham/Pcham (missing met data)
for (nm in c("tk1", "tk2", "tk3")) {
  df <- get(nm)
  n_na_t <- sum(is.na(df$Tcham))
  n_na_p <- sum(is.na(df$Pcham))
  if (n_na_t > 0 | n_na_p > 0) {
    warning(nm, ": ", n_na_t, " NA Tcham, ", n_na_p, " NA Pcham values. ",
            "These measurements may have missing met data.")
  }
}

# =============================================================================
# 5. Build final auxfiles
# =============================================================================

build_auxfile <- function(df) {
  df %>%
    mutate(UniqueID = trimws(UniqueID)) %>%
    select(UniqueID, start.time, Area, Vtot, Tcham, Pcham) %>%
    filter(!is.na(start.time))
}

aux.LGR1 <- build_auxfile(tk1)
aux.LGR2 <- build_auxfile(tk2)
aux.LGR3 <- build_auxfile(tk3)

# Print summaries
message("\n=== Auxfile summaries ===")
for (nm in c("aux.LGR1", "aux.LGR2", "aux.LGR3")) {
  df <- get(nm)
  message(nm, ": ", nrow(df), " measurements")
  message("  Area range: ", min(df$Area), " - ", max(df$Area), " cm2")
  message("  Vtot range: ", round(min(df$Vtot), 3), " - ",
          round(max(df$Vtot), 3), " L")
  message("  Tcham range: ", round(min(df$Tcham, na.rm = TRUE), 1), " - ",
          round(max(df$Tcham, na.rm = TRUE), 1), " C")
  message("  Pcham range: ", round(min(df$Pcham, na.rm = TRUE), 1), " - ",
          round(max(df$Pcham, na.rm = TRUE), 1), " kPa")
}

# =============================================================================
# 6. Save auxfiles
# =============================================================================

save(aux.LGR1, file = file.path(rdata_dir, "aux_LGR1.RData"))
save(aux.LGR2, file = file.path(rdata_dir, "aux_LGR2.RData"))
save(aux.LGR3, file = file.path(rdata_dir, "aux_LGR3.RData"))

# Also save as CSV for manual inspection
write.csv(aux.LGR1, file.path(results_dir, "aux_LGR1.csv"), row.names = FALSE)
write.csv(aux.LGR2, file.path(results_dir, "aux_LGR2.csv"), row.names = FALSE)
write.csv(aux.LGR3, file.path(results_dir, "aux_LGR3.csv"), row.names = FALSE)

message("\nAuxfiles saved to: ", rdata_dir)
message("CSV copies saved to: ", results_dir)
message("Proceed to 04_manual_id.R")
