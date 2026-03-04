# =============================================================================
# 09_mdf_lod_comparison.R
# Compute MDF/LOD using three published approaches, flag below-detection fluxes,
# and save updated compiled results.
#
# Three approaches:
#   1. goFlux (Christiansen et al. 2015 / goFlux package):
#      MDF = manufacturer_precision / t × flux.term
#   2. Wassmann et al. 2018:
#      LOD_flux = (3 × SD_empirical_global) / t × flux.term
#   3. Christiansen et al. 2015 (MQL variant):
#      MDF = (SD_empirical_per_meas × 3 × t₉₉%) / t × flux.term
#
# Empirical precision is derived two ways:
#   - Allan deviation within each measurement window (per-measurement)
#   - Rolling window minimum-slope-and-noise scan (global per-instrument)
#
# Instrument: ABB GLA131-GGA (Microportable UGGA)
#   CO₂ precision: 0.35 ppm (1σ, 1s)
#   CH₄ precision: 0.9 ppb (1σ, 1s)
# =============================================================================

# Source setup
setup_path <- file.path(
  "/Users/jongewirtzman/My Drive/Research/whole_tree_flux",
  "data processing", "goFlux_reprocessing", "00_setup.R")
source(setup_path)

# =============================================================================
# Step A: Load data
# =============================================================================

message("\n=== Loading compiled results ===")

hf_compiled <- read.csv(
  file.path(results_dir, "canopy_flux_goFlux_compiled.csv"),
  stringsAsFactors = FALSE)
hf_compiled$UniqueID <- trimws(hf_compiled$UniqueID)

ymf_compiled <- read.csv(
  file.path(reprocess_dir, "ymf_black_oak", "results",
            "ymf_black_oak_flux_compiled.csv"),
  stringsAsFactors = FALSE)
ymf_compiled$UniqueID <- trimws(ymf_compiled$UniqueID)

message("HF compiled: ", nrow(hf_compiled), " rows")
message("YMF compiled: ", nrow(ymf_compiled), " rows")

# Load manID data (per-measurement concentration timeseries)
message("\n=== Loading manID data ===")
load(file.path(rdata_dir, "manID_LGR1.RData"))
load(file.path(rdata_dir, "manID_LGR2.RData"))
load(file.path(rdata_dir, "manID_LGR3.RData"))
load(file.path(reprocess_dir, "ymf_black_oak", "RData", "manID_YMF.RData"))

manID_all_hf <- rbind(manID.LGR1, manID.LGR2, manID.LGR3)
message("HF manID rows: ", nrow(manID_all_hf),
        " (", length(unique(manID_all_hf$UniqueID)), " measurements)")
message("YMF manID rows: ", nrow(manID.YMF),
        " (", length(unique(manID.YMF$UniqueID)), " measurements)")

# Load full continuous imported data (for rolling window)
message("\n=== Loading continuous imported data ===")
load(file.path(rdata_dir, "imp_LGR1_combined.RData"))
load(file.path(rdata_dir, "imp_LGR2_combined.RData"))
load(file.path(rdata_dir, "imp_LGR3_combined.RData"))
load(file.path(reprocess_dir, "ymf_black_oak", "RData", "imp_YMF_combined.RData"))

# Tag each with instrument name
imp.LGR1$instrument <- "LGR1"
imp.LGR2$instrument <- "LGR2"
imp.LGR3$instrument <- "LGR3"

# Determine which instrument each HF measurement came from
# UniqueIDs in manID data tell us the mapping
lgr1_uids <- unique(manID.LGR1$UniqueID)
lgr2_uids <- unique(manID.LGR2$UniqueID)
lgr3_uids <- unique(manID.LGR3$UniqueID)

hf_compiled$instrument <- NA_character_
hf_compiled$instrument[hf_compiled$UniqueID %in% lgr1_uids] <- "LGR1"
hf_compiled$instrument[hf_compiled$UniqueID %in% lgr2_uids] <- "LGR2"
hf_compiled$instrument[hf_compiled$UniqueID %in% lgr3_uids] <- "LGR3"
ymf_compiled$instrument <- "YMF"

# =============================================================================
# Step B, Method 1: Allan deviation within measurement windows
# =============================================================================

message("\n=== Computing Allan deviation per measurement ===")

# Function: compute Allan deviation (τ=1s) from a concentration timeseries.
#
# During a flux measurement, concentrations change monotonically (the trend IS
# the flux signal).  We need to estimate instrument noise *without* that trend
# contaminating the estimate.  First-differencing achieves this:
#
#   δ_i = C(t_{i+1}) - C(t_i)
#
# If the underlying signal is any smooth trend f(t), the first difference is
# approximately f'(t)*Δt (a constant for linear trends, slowly varying for
# nonlinear ones).  Instrument white noise ε_i, by contrast, does not cancel:
# Δε_i = ε_{i+1} - ε_i has variance 2*σ^2.  So SD(δ) ≈ sqrt(2)*σ, and
#
#   σ_Allan = SD(δ) / sqrt(2)
#
# recovers the 1-sigma instrument noise regardless of the trend shape --
# no model fitting required.  This is the standard Allan deviation at τ = 1 s.
allan_sd <- function(x) {
  diffs <- diff(x)
  if (length(diffs) < 2) return(NA_real_)
  sd(diffs) / sqrt(2)
}

# Compute for each UniqueID in HF manID data
# Use flag == 1 rows (within the identified measurement window)
compute_allan_per_uid <- function(manID_df) {
  manID_df %>%
    filter(flag == 1) %>%
    group_by(UniqueID) %>%
    summarise(
      allan_sd_CO2 = allan_sd(CO2dry_ppm),
      allan_sd_CH4 = allan_sd(CH4dry_ppb),
      n_meas_pts   = n(),
      .groups = "drop"
    )
}

allan_hf  <- compute_allan_per_uid(manID_all_hf)
allan_ymf <- compute_allan_per_uid(manID.YMF)

message("HF Allan deviation computed for ", nrow(allan_hf), " measurements")
message("  CO2 Allan SD range: ", round(min(allan_hf$allan_sd_CO2, na.rm = TRUE), 4),
        " - ", round(max(allan_hf$allan_sd_CO2, na.rm = TRUE), 4), " ppm")
message("  CH4 Allan SD range: ", round(min(allan_hf$allan_sd_CH4, na.rm = TRUE), 4),
        " - ", round(max(allan_hf$allan_sd_CH4, na.rm = TRUE), 4), " ppb")

message("YMF Allan deviation computed for ", nrow(allan_ymf), " measurements")
message("  CO2 Allan SD range: ", round(min(allan_ymf$allan_sd_CO2, na.rm = TRUE), 4),
        " - ", round(max(allan_ymf$allan_sd_CO2, na.rm = TRUE), 4), " ppm")
message("  CH4 Allan SD range: ", round(min(allan_ymf$allan_sd_CH4, na.rm = TRUE), 4),
        " - ", round(max(allan_ymf$allan_sd_CH4, na.rm = TRUE), 4), " ppb")

# =============================================================================
# Step B, Method 2: Rolling window minimum-slope-and-noise scan
# =============================================================================

message("\n=== Rolling window scan for stable periods ===")

# Function: scan a continuous timeseries with rolling windows, scoring each
# by |slope| + residual_SD (both normalized), and return the SD from the
# quietest/flattest bottom 5% of windows.
rolling_window_precision <- function(df, col_name, window_sec = 30,
                                     quantile_cutoff = 0.05) {
  x <- df[[col_name]]
  t_sec <- as.numeric(difftime(df$POSIX.time, df$POSIX.time[1], units = "secs"))
  n <- length(x)

  if (n < window_sec + 1) {
    warning("Timeseries too short for rolling window (", n, " points)")
    return(list(sd = NA_real_, n_windows = 0, n_selected = 0))
  }

  # Pre-allocate
  n_windows <- n - window_sec + 1
  slopes    <- numeric(n_windows)
  resid_sds <- numeric(n_windows)

  for (i in seq_len(n_windows)) {
    idx <- i:(i + window_sec - 1)
    t_win <- t_sec[idx] - t_sec[idx[1]]
    x_win <- x[idx]

    # Quick linear regression via lm.fit for speed
    X <- cbind(1, t_win)
    fit <- .lm.fit(X, x_win)
    slopes[i] <- abs(fit$coefficients[2])
    resid_sds[i] <- sd(fit$residuals)
  }

  # Normalize both to [0, 1]
  slope_range <- range(slopes, na.rm = TRUE)
  sd_range    <- range(resid_sds, na.rm = TRUE)

  if (slope_range[2] - slope_range[1] == 0) {
    norm_slopes <- rep(0, n_windows)
  } else {
    norm_slopes <- (slopes - slope_range[1]) / (slope_range[2] - slope_range[1])
  }
  if (sd_range[2] - sd_range[1] == 0) {
    norm_sds <- rep(0, n_windows)
  } else {
    norm_sds <- (resid_sds - sd_range[1]) / (sd_range[2] - sd_range[1])
  }

  # Combined score: low = flat AND quiet
  combined_score <- norm_slopes + norm_sds

  # Select bottom 5%
  threshold <- quantile(combined_score, probs = quantile_cutoff, na.rm = TRUE)
  selected  <- which(combined_score <= threshold)
  n_selected <- length(selected)

  # Use the median residual SD of the selected windows as the precision estimate
  precision_sd <- median(resid_sds[selected])

  list(sd = precision_sd, n_windows = n_windows, n_selected = n_selected)
}

# Run for each HF instrument
message("Processing LGR1...")
rw_lgr1_co2 <- rolling_window_precision(imp.LGR1, "CO2dry_ppm")
rw_lgr1_ch4 <- rolling_window_precision(imp.LGR1, "CH4dry_ppb")

message("Processing LGR2...")
rw_lgr2_co2 <- rolling_window_precision(imp.LGR2, "CO2dry_ppm")
rw_lgr2_ch4 <- rolling_window_precision(imp.LGR2, "CH4dry_ppb")

message("Processing LGR3...")
rw_lgr3_co2 <- rolling_window_precision(imp.LGR3, "CO2dry_ppm")
rw_lgr3_ch4 <- rolling_window_precision(imp.LGR3, "CH4dry_ppb")

# YMF
ymf_imp_name <- ls(pattern = "^imp\\.")
ymf_imp_name <- ymf_imp_name[grepl("YMF", ymf_imp_name, ignore.case = TRUE)]
if (length(ymf_imp_name) == 0) {
  # Try alternate naming
  ymf_imp_name <- ls(pattern = "^imp_YMF|^imp\\.YMF")
}
ymf_imp <- get(ymf_imp_name[1])

message("Processing YMF...")
rw_ymf_co2 <- rolling_window_precision(ymf_imp, "CO2dry_ppm")
rw_ymf_ch4 <- rolling_window_precision(ymf_imp, "CH4dry_ppb")

# Collect rolling window results
rw_results <- data.frame(
  instrument = c("LGR1", "LGR2", "LGR3", "YMF"),
  flat_sd_CO2 = c(rw_lgr1_co2$sd, rw_lgr2_co2$sd, rw_lgr3_co2$sd, rw_ymf_co2$sd),
  flat_sd_CH4 = c(rw_lgr1_ch4$sd, rw_lgr2_ch4$sd, rw_lgr3_ch4$sd, rw_ymf_ch4$sd),
  n_windows   = c(rw_lgr1_co2$n_windows, rw_lgr2_co2$n_windows,
                  rw_lgr3_co2$n_windows, rw_ymf_co2$n_windows),
  n_selected  = c(rw_lgr1_co2$n_selected, rw_lgr2_co2$n_selected,
                  rw_lgr3_co2$n_selected, rw_ymf_co2$n_selected),
  stringsAsFactors = FALSE
)

message("\n--- Rolling window precision estimates (bottom 5% quietest windows) ---")
print(rw_results)

message("\n--- Comparison: Manufacturer vs Empirical precision ---")
message("Manufacturer (GLA131): CO2 = 0.35 ppm, CH4 = 0.9 ppb")
for (i in seq_len(nrow(rw_results))) {
  message(rw_results$instrument[i], ": CO2 = ",
          round(rw_results$flat_sd_CO2[i], 4), " ppm, CH4 = ",
          round(rw_results$flat_sd_CH4[i], 4), " ppb")
}

# =============================================================================
# Step C: Compute three MDF approaches
# =============================================================================

message("\n=== Computing three MDF approaches ===")

# GLA131 manufacturer precision (1σ, 1s)
prec_co2 <- 0.35  # ppm
prec_ch4 <- 0.9   # ppb

# --- Helper: compute MDFs at multiple confidence levels ---
#
# Wassmann multipliers (z-scores for normal distribution):
#   99%: 2.576,  95%: 1.960,  90%: 1.645
# Christiansen multipliers (SD × 3 × t_α):
#   99%: qt(0.995, df),  95%: qt(0.975, df),  90%: qt(0.95, df)

z_mult <- c("99" = qnorm(0.995), "95" = qnorm(0.975), "90" = qnorm(0.95))

compute_mdf_all <- function(compiled, allan_df, rw_df) {

  # Merge Allan deviation per measurement
  compiled <- merge(compiled, allan_df[, c("UniqueID", "allan_sd_CO2",
                                            "allan_sd_CH4", "n_meas_pts")],
                    by = "UniqueID", all.x = TRUE)

  # Merge rolling window precision per instrument
  compiled <- merge(compiled,
                    rw_df[, c("instrument", "flat_sd_CO2", "flat_sd_CH4")],
                    by = "instrument", all.x = TRUE)

  compiled <- compiled %>%
    mutate(
      # Measurement time in seconds (= nb.obs at 1 Hz)
      t_sec = CO2_nb.obs,

      # --- Approach 1: goFlux / manufacturer precision (no confidence level) ---
      CO2_MDF_goflux = prec_co2 / t_sec * CO2_flux.term,
      CH4_MDF_goflux = prec_ch4 / t_sec * CH4_flux.term,

      # --- Approach 2: Wassmann (z × global empirical SD) at 99/95/90% ---
      CO2_MDF_wass99 = (z_mult["99"] * flat_sd_CO2) / t_sec * CO2_flux.term,
      CH4_MDF_wass99 = (z_mult["99"] * flat_sd_CH4) / t_sec * CH4_flux.term,
      CO2_MDF_wass95 = (z_mult["95"] * flat_sd_CO2) / t_sec * CO2_flux.term,
      CH4_MDF_wass95 = (z_mult["95"] * flat_sd_CH4) / t_sec * CH4_flux.term,
      CO2_MDF_wass90 = (z_mult["90"] * flat_sd_CO2) / t_sec * CO2_flux.term,
      CH4_MDF_wass90 = (z_mult["90"] * flat_sd_CH4) / t_sec * CH4_flux.term,

      # --- Approach 3: Christiansen (per-meas SD × 3 × t_α) at 99/95/90% ---
      df_meas = pmax(n_meas_pts - 2, 1),
      t99 = qt(0.995, df = df_meas),
      t95 = qt(0.975, df = df_meas),
      t90 = qt(0.95,  df = df_meas),
      CO2_MDF_chr99 = (allan_sd_CO2 * 3 * t99) / t_sec * CO2_flux.term,
      CH4_MDF_chr99 = (allan_sd_CH4 * 3 * t99) / t_sec * CH4_flux.term,
      CO2_MDF_chr95 = (allan_sd_CO2 * 3 * t95) / t_sec * CO2_flux.term,
      CH4_MDF_chr95 = (allan_sd_CH4 * 3 * t95) / t_sec * CH4_flux.term,
      CO2_MDF_chr90 = (allan_sd_CO2 * 3 * t90) / t_sec * CO2_flux.term,
      CH4_MDF_chr90 = (allan_sd_CH4 * 3 * t90) / t_sec * CH4_flux.term
    )

  # --- Step D: Flag below-detection fluxes ---
  compiled <- compiled %>%
    mutate(
      # CO2 flags — all methods and confidence levels
      CO2_below_MDF_goflux = abs(CO2_best.flux) < CO2_MDF_goflux,
      CO2_below_MDF_wass99 = abs(CO2_best.flux) < CO2_MDF_wass99,
      CO2_below_MDF_wass95 = abs(CO2_best.flux) < CO2_MDF_wass95,
      CO2_below_MDF_wass90 = abs(CO2_best.flux) < CO2_MDF_wass90,
      CO2_below_MDF_chr99  = abs(CO2_best.flux) < CO2_MDF_chr99,
      CO2_below_MDF_chr95  = abs(CO2_best.flux) < CO2_MDF_chr95,
      CO2_below_MDF_chr90  = abs(CO2_best.flux) < CO2_MDF_chr90,
      # Summary flags use 99% (most conservative)
      CO2_below_any_MDF = CO2_below_MDF_goflux | CO2_below_MDF_wass99 |
                           CO2_below_MDF_chr99,
      CO2_below_all_MDF = CO2_below_MDF_goflux & CO2_below_MDF_wass99 &
                           CO2_below_MDF_chr99,

      # CH4 flags — all methods and confidence levels
      CH4_below_MDF_goflux = abs(CH4_best.flux) < CH4_MDF_goflux,
      CH4_below_MDF_wass99 = abs(CH4_best.flux) < CH4_MDF_wass99,
      CH4_below_MDF_wass95 = abs(CH4_best.flux) < CH4_MDF_wass95,
      CH4_below_MDF_wass90 = abs(CH4_best.flux) < CH4_MDF_wass90,
      CH4_below_MDF_chr99  = abs(CH4_best.flux) < CH4_MDF_chr99,
      CH4_below_MDF_chr95  = abs(CH4_best.flux) < CH4_MDF_chr95,
      CH4_below_MDF_chr90  = abs(CH4_best.flux) < CH4_MDF_chr90,
      # Summary flags use 99% (most conservative)
      CH4_below_any_MDF = CH4_below_MDF_goflux | CH4_below_MDF_wass99 |
                           CH4_below_MDF_chr99,
      CH4_below_all_MDF = CH4_below_MDF_goflux & CH4_below_MDF_wass99 &
                           CH4_below_MDF_chr99
    )

  # Clean up intermediate columns
  compiled <- compiled %>% select(-df_meas, -t99, -t95, -t90)

  compiled
}

# Apply to HF
hf_out <- compute_mdf_all(hf_compiled, allan_hf, rw_results)

# Apply to YMF
ymf_out <- compute_mdf_all(ymf_compiled, allan_ymf, rw_results)

# --- Add best-model R² and SNR ---
# SNR uses per-measurement Allan deviation as empirical noise floor:
#   noise_floor_flux = allan_sd / t_sec × flux.term  (1σ in flux units)
#   SNR = |flux| / noise_floor_flux
add_quality_cols <- function(df) {
  df %>%
    mutate(
      CH4_best.r2 = ifelse(CH4_model == "LM", CH4_LM.r2, CH4_HM.r2),
      CO2_best.r2 = ifelse(CO2_model == "LM", CO2_LM.r2, CO2_HM.r2),
      CH4_noise_floor = allan_sd_CH4 / t_sec * CH4_flux.term,
      CH4_SNR     = abs(CH4_best.flux) / CH4_noise_floor
    )
}
hf_out  <- add_quality_cols(hf_out)
ymf_out <- add_quality_cols(ymf_out)

# =============================================================================
# Step E: Summary and save
# =============================================================================

# Helper to print detection summary for a dataset
print_mdf_summary <- function(out, dataset_name) {
  message("\n=== MDF Summary: ", dataset_name, " ===")

  for (gas in c("CO2", "CH4")) {
    n_valid <- sum(!is.na(out[[paste0(gas, "_best.flux")]]))
    message("\n  ", gas, " (n=", n_valid, " measurements):")

    # MDF ranges
    for (method in c("goflux", "wass99", "wass95", "wass90",
                     "chr99", "chr95", "chr90")) {
      col <- paste0(gas, "_MDF_", method)
      if (col %in% names(out)) {
        label <- switch(method,
          goflux = "goFlux (manufacturer)",
          wass99 = "Wassmann 99%", wass95 = "Wassmann 95%", wass90 = "Wassmann 90%",
          chr99  = "Christiansen 99%", chr95 = "Christiansen 95%", chr90 = "Christiansen 90%")
        message(sprintf("    MDF range (%-20s): %s - %s", label,
                round(min(out[[col]], na.rm = TRUE), 6),
                round(max(out[[col]], na.rm = TRUE), 6)))
      }
    }

    # Detection rates
    message("    --- % below detection ---")
    below_g <- paste0(gas, "_below_MDF_goflux")
    message(sprintf("    goFlux (manufacturer):    %d / %d (%.1f%%)",
            sum(out[[below_g]], na.rm = TRUE), n_valid,
            100 * mean(out[[below_g]], na.rm = TRUE)))

    for (conf in c("99", "95", "90")) {
      wass_col <- paste0(gas, "_below_MDF_wass", conf)
      chr_col  <- paste0(gas, "_below_MDF_chr", conf)
      if (wass_col %in% names(out)) {
        message(sprintf("    Wassmann %s%%:              %d / %d (%.1f%%)", conf,
                sum(out[[wass_col]], na.rm = TRUE), n_valid,
                100 * mean(out[[wass_col]], na.rm = TRUE)))
      }
      if (chr_col %in% names(out)) {
        message(sprintf("    Christiansen %s%%:          %d / %d (%.1f%%)", conf,
                sum(out[[chr_col]], na.rm = TRUE), n_valid,
                100 * mean(out[[chr_col]], na.rm = TRUE)))
      }
    }

    below_any <- paste0(gas, "_below_any_MDF")
    below_all <- paste0(gas, "_below_all_MDF")
    message(sprintf("    Below ANY (99%%):          %d / %d (%.1f%%)",
            sum(out[[below_any]], na.rm = TRUE), n_valid,
            100 * mean(out[[below_any]], na.rm = TRUE)))
    message(sprintf("    Below ALL (99%%):          %d / %d (%.1f%%)",
            sum(out[[below_all]], na.rm = TRUE), n_valid,
            100 * mean(out[[below_all]], na.rm = TRUE)))
  }
}

print_mdf_summary(hf_out, "Harvard Forest")
print_mdf_summary(ymf_out, "Yale-Myers Forest")

# --- Save updated compiled results ---

message("\n=== Saving results ===")

write.csv(hf_out,
          file.path(results_dir, "canopy_flux_goFlux_compiled_with_mdf.csv"),
          row.names = FALSE)
write.xlsx(hf_out,
           file.path(results_dir, "canopy_flux_goFlux_compiled_with_mdf.xlsx"))

ymf_out_dir <- file.path(reprocess_dir, "ymf_black_oak", "results")
write.csv(ymf_out,
          file.path(ymf_out_dir, "ymf_black_oak_flux_compiled_with_mdf.csv"),
          row.names = FALSE)
write.xlsx(ymf_out,
           file.path(ymf_out_dir, "ymf_black_oak_flux_compiled_with_mdf.xlsx"))

# --- Save precision summary ---

precision_summary <- data.frame(
  source = c("GLA131 manufacturer (1σ, 1s)",
             rw_results$instrument),
  CO2_precision_ppm = c(prec_co2, rw_results$flat_sd_CO2),
  CH4_precision_ppb = c(prec_ch4, rw_results$flat_sd_CH4),
  method = c("Datasheet",
             rep("Rolling window (bottom 5%)", nrow(rw_results))),
  stringsAsFactors = FALSE
)

# Add median Allan deviation per instrument
allan_median_hf <- manID_all_hf %>%
  filter(flag == 1) %>%
  group_by(UniqueID) %>%
  summarise(
    allan_co2 = allan_sd(CO2dry_ppm),
    allan_ch4 = allan_sd(CH4dry_ppb),
    .groups = "drop"
  )

# Map UniqueID → instrument
allan_median_hf$instrument <- NA_character_
allan_median_hf$instrument[allan_median_hf$UniqueID %in% lgr1_uids] <- "LGR1"
allan_median_hf$instrument[allan_median_hf$UniqueID %in% lgr2_uids] <- "LGR2"
allan_median_hf$instrument[allan_median_hf$UniqueID %in% lgr3_uids] <- "LGR3"

allan_by_inst <- allan_median_hf %>%
  group_by(instrument) %>%
  summarise(
    median_allan_CO2 = median(allan_co2, na.rm = TRUE),
    median_allan_CH4 = median(allan_ch4, na.rm = TRUE),
    .groups = "drop"
  )

# YMF Allan median
allan_ymf_median <- data.frame(
  instrument = "YMF",
  median_allan_CO2 = median(allan_ymf$allan_sd_CO2, na.rm = TRUE),
  median_allan_CH4 = median(allan_ymf$allan_sd_CH4, na.rm = TRUE)
)
allan_by_inst <- rbind(allan_by_inst, allan_ymf_median)

allan_summary <- data.frame(
  source = allan_by_inst$instrument,
  CO2_precision_ppm = allan_by_inst$median_allan_CO2,
  CH4_precision_ppb = allan_by_inst$median_allan_CH4,
  method = "Allan deviation (median)",
  stringsAsFactors = FALSE
)

precision_summary <- rbind(precision_summary, allan_summary)

# Rename for display
precision_summary$source <- sub("^GLA131.*", "GLA131 spec", precision_summary$source)
precision_summary$source <- sub("^YMF$", "LGR3-YMF", precision_summary$source)

write.csv(precision_summary,
          file.path(results_dir, "precision_comparison.csv"),
          row.names = FALSE)

message("\n--- Precision comparison table ---")
print(precision_summary)

message("\nResults saved to:")
message("  HF: ", file.path(results_dir, "canopy_flux_goFlux_compiled_with_mdf.csv"))
message("  YMF: ", file.path(ymf_out_dir, "ymf_black_oak_flux_compiled_with_mdf.csv"))
message("  Precision: ", file.path(results_dir, "precision_comparison.csv"))

# =============================================================================
# Step F: Quality criteria summary table + heatmap
# =============================================================================

message("\n=== Computing quality criteria summary ===")

# Select only columns needed for quality criteria (avoids type mismatches in metadata)
quality_cols <- c("UniqueID", "CH4_best.flux", "CO2_best.flux",
                  "CH4_below_MDF_goflux", "CH4_below_MDF_wass99", "CH4_below_MDF_wass95",
                  "CH4_below_MDF_wass90", "CH4_below_MDF_chr99", "CH4_below_MDF_chr95",
                  "CH4_below_MDF_chr90", "CH4_best.r2", "CO2_best.r2",
                  "CH4_noise_floor", "CH4_SNR")
all_out <- bind_rows(hf_out[, quality_cols], ymf_out[, quality_cols])
n_total <- nrow(all_out)

message("Combined dataset: ", n_total, " measurements (HF + YMF)")

# --- Build summary table: % of CH4 fluxes passing each criterion ---
# (CO2 MDF criteria are also included for completeness)

criteria <- data.frame(
  criterion = c(
    # MDF-based criteria (CH4)
    "CH4 |flux| > MDF (manufacturer)",
    "CH4 |flux| > MDF (Wassmann 99%)",
    "CH4 |flux| > MDF (Wassmann 95%)",
    "CH4 |flux| > MDF (Wassmann 90%)",
    "CH4 |flux| > MDF (Christiansen 99%)",
    "CH4 |flux| > MDF (Christiansen 95%)",
    "CH4 |flux| > MDF (Christiansen 90%)",
    # R-squared criteria
    paste0("CH4 R\u00B2 > 0.7"),
    paste0("CH4 R\u00B2 > 0.9"),
    paste0("CO2 R\u00B2 > 0.7 (same meas.)"),
    paste0("CO2 R\u00B2 > 0.9 (same meas.)"),
    # SNR criteria (empirical, Allan deviation per measurement)
    "CH4 SNR > 2 (empirical)",
    "CH4 SNR > 3 (empirical)"
  ),
  n_pass = c(
    sum(!all_out$CH4_below_MDF_goflux, na.rm = TRUE),
    sum(!all_out$CH4_below_MDF_wass99, na.rm = TRUE),
    sum(!all_out$CH4_below_MDF_wass95, na.rm = TRUE),
    sum(!all_out$CH4_below_MDF_wass90, na.rm = TRUE),
    sum(!all_out$CH4_below_MDF_chr99,  na.rm = TRUE),
    sum(!all_out$CH4_below_MDF_chr95,  na.rm = TRUE),
    sum(!all_out$CH4_below_MDF_chr90,  na.rm = TRUE),
    sum(all_out$CH4_best.r2 > 0.7, na.rm = TRUE),
    sum(all_out$CH4_best.r2 > 0.9, na.rm = TRUE),
    sum(all_out$CO2_best.r2 > 0.7, na.rm = TRUE),
    sum(all_out$CO2_best.r2 > 0.9, na.rm = TRUE),
    sum(all_out$CH4_SNR > 2, na.rm = TRUE),
    sum(all_out$CH4_SNR > 3, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

criteria$n_total   <- n_total
criteria$pct_pass  <- round(100 * criteria$n_pass / n_total, 1)
criteria$label     <- paste0(criteria$n_pass, "/", n_total,
                             " (", criteria$pct_pass, "%)")

# Add a category column for the heatmap
criteria$category <- c(
  rep("MDF threshold", 7),
  rep(paste0("R\u00B2 threshold"), 4),
  rep("SNR threshold", 2)
)

write.csv(criteria,
          file.path(results_dir, "quality_criteria_summary.csv"),
          row.names = FALSE)

message("\n--- Quality criteria summary (CH4 fluxes, all sites combined) ---")
print(criteria[, c("criterion", "n_pass", "n_total", "pct_pass")])


# =============================================================================
# Step G: Visualize empirical precision estimates
# =============================================================================

library(ggplot2)
library(tidyr)

prec_plot_dir <- file.path(plots_dir, "precision")
dir.create(prec_plot_dir, recursive = TRUE, showWarnings = FALSE)

plot_theme <- theme_classic(base_size = 11) +
  theme(
    strip.text       = element_text(size = 10),
    strip.background = element_blank(),
    legend.position  = "bottom",
    axis.line        = element_line(linewidth = 0.3),
    axis.ticks       = element_line(linewidth = 0.3)
  )

# Instrument display labels (YMF used the same LGR3 unit on a different day)
inst_display <- c("LGR1" = "LGR1", "LGR2" = "LGR2",
                  "LGR3" = "LGR3", "YMF" = "LGR3-YMF")

# --- Plot 1: Allan deviation distributions per instrument ---
# Jittered points + boxplot, with manufacturer spec (dashed line) and
# rolling-window estimate (diamond) as reference. Both in legend.

allan_hf_long <- allan_hf %>%
  mutate(instrument = case_when(
    UniqueID %in% lgr1_uids ~ "LGR1",
    UniqueID %in% lgr2_uids ~ "LGR2",
    UniqueID %in% lgr3_uids ~ "LGR3"
  )) %>%
  pivot_longer(cols = c(allan_sd_CO2, allan_sd_CH4),
               names_to = "gas", values_to = "allan_sd") %>%
  mutate(gas = ifelse(gas == "allan_sd_CO2", "CO2", "CH4"),
         dataset = "Harvard Forest")

allan_ymf_long <- allan_ymf %>%
  mutate(instrument = "YMF") %>%
  pivot_longer(cols = c(allan_sd_CO2, allan_sd_CH4),
               names_to = "gas", values_to = "allan_sd") %>%
  mutate(gas = ifelse(gas == "allan_sd_CO2", "CO2", "CH4"),
         dataset = "Yale-Myers Forest")

allan_long <- rbind(allan_hf_long, allan_ymf_long)
allan_long$instrument <- inst_display[allan_long$instrument]

manuf_ref <- data.frame(
  gas = c("CO2", "CH4"),
  manuf_prec = c(prec_co2, prec_ch4)
)

rw_long <- rw_results %>%
  pivot_longer(cols = c(flat_sd_CO2, flat_sd_CH4),
               names_to = "gas", values_to = "rw_sd") %>%
  mutate(gas = ifelse(gas == "flat_sd_CO2", "CO2", "CH4"))
rw_long$instrument <- inst_display[rw_long$instrument]

# Helper: build Allan deviation plot with rolling-window diamonds in legend
make_allan_plot <- function(data, rw_data, manuf_val, y_lab, title_expr) {
  ggplot(data, aes(x = instrument, y = allan_sd)) +
    geom_boxplot(outlier.shape = NA, width = 0.5, fill = "grey90") +
    geom_jitter(aes(shape = "Allan deviation"),
                width = 0.15, alpha = 0.5, size = 1.5, color = "#2166AC") +
    geom_hline(yintercept = manuf_val, linetype = "dashed",
               color = "red", linewidth = 0.6) +
    geom_point(data = rw_data,
               aes(x = instrument, y = rw_sd, shape = "Rolling window"),
               size = 4, color = "#B2182B") +
    scale_shape_manual(
      values = c("Allan deviation" = 16, "Rolling window" = 18),
      name = NULL
    ) +
    labs(y = y_lab, x = "Instrument", title = title_expr) +
    annotate("text", x = 0.55, y = manuf_val * 1.15,
             label = "GLA131 spec", color = "red", hjust = 0, size = 3) +
    plot_theme +
    guides(shape = guide_legend(override.aes = list(
      color = c("#2166AC", "#B2182B"), size = c(2, 4), alpha = c(0.7, 1)
    )))
}

p_allan_co2 <- make_allan_plot(
  allan_long %>% filter(gas == "CO2"),
  rw_long %>% filter(gas == "CO2"),
  prec_co2,
  expression(Allan~deviation~(ppm~CO[2])),
  expression(CO[2]~empirical~precision~per~measurement)
)

p_allan_ch4 <- make_allan_plot(
  allan_long %>% filter(gas == "CH4"),
  rw_long %>% filter(gas == "CH4"),
  prec_ch4,
  expression(Allan~deviation~(ppb~CH[4])),
  expression(CH[4]~empirical~precision~per~measurement)
)

ggsave(file.path(prec_plot_dir, "allan_deviation_CO2.pdf"),
       p_allan_co2, width = 5, height = 4)
ggsave(file.path(prec_plot_dir, "allan_deviation_CO2.png"),
       p_allan_co2, width = 5, height = 4, dpi = 300)
ggsave(file.path(prec_plot_dir, "allan_deviation_CH4.pdf"),
       p_allan_ch4, width = 5, height = 4)
ggsave(file.path(prec_plot_dir, "allan_deviation_CH4.png"),
       p_allan_ch4, width = 5, height = 4, dpi = 300)

message("Saved Allan deviation plots to: ", prec_plot_dir)

# --- Plot 2: Precision comparison bar chart ---

prec_bar_data <- precision_summary %>%
  pivot_longer(cols = c(CO2_precision_ppm, CH4_precision_ppb),
               names_to = "gas", values_to = "precision") %>%
  mutate(
    gas = ifelse(gas == "CO2_precision_ppm", "CO2", "CH4"),
    method = factor(method,
                    levels = c("Datasheet", "Rolling window (bottom 5%)",
                               "Allan deviation (median)"))
  )

for (g in c("CO2", "CH4")) {
  # Use if/else (not ifelse) so expression() works properly
  if (g == "CO2") {
    gas_expr <- expression(CO[2]~precision~(ppm))
  } else {
    gas_expr <- expression(CH[4]~precision~(ppb))
  }

  p_bar <- ggplot(prec_bar_data %>% filter(gas == g),
                  aes(x = source, y = precision, fill = method)) +
    geom_col(position = "dodge", width = 0.7, alpha = 0.85) +
    scale_fill_manual(values = c("Datasheet" = "#D6604D",
                                 "Rolling window (bottom 5%)" = "#4393C3",
                                 "Allan deviation (median)" = "#2166AC")) +
    labs(y = gas_expr, x = NULL, fill = "Method",
         title = paste0(g, " precision: manufacturer vs. empirical")) +
    plot_theme +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(file.path(prec_plot_dir, paste0("precision_comparison_", g, ".pdf")),
         p_bar, width = 6, height = 4)
  ggsave(file.path(prec_plot_dir, paste0("precision_comparison_", g, ".png")),
         p_bar, width = 6, height = 4, dpi = 300)
}

message("Saved precision comparison bar charts to: ", prec_plot_dir)

# --- Plot 3: Rolling window scan — selected quiet windows by period ---
# Faceted by recording period (gaps > 60 min define separate periods).
# Y-axis zoomed to ambient range to show detail.

plot_rolling_window_selection <- function(imp_df, col_name, instrument_name,
                                          gas_name, window_sec = 30,
                                          quantile_cutoff = 0.05,
                                          ylimits = NULL,
                                          gap_threshold_min = 60) {
  x <- imp_df[[col_name]]
  t_sec <- as.numeric(difftime(imp_df$POSIX.time, imp_df$POSIX.time[1],
                                units = "secs"))
  n <- length(x)
  if (n < window_sec + 1) return(NULL)

  n_windows <- n - window_sec + 1
  slopes    <- numeric(n_windows)
  resid_sds <- numeric(n_windows)

  for (i in seq_len(n_windows)) {
    idx <- i:(i + window_sec - 1)
    t_win <- t_sec[idx] - t_sec[idx[1]]
    x_win <- x[idx]
    X <- cbind(1, t_win)
    fit <- .lm.fit(X, x_win)
    slopes[i] <- abs(fit$coefficients[2])
    resid_sds[i] <- sd(fit$residuals)
  }

  slope_range <- range(slopes, na.rm = TRUE)
  sd_range    <- range(resid_sds, na.rm = TRUE)
  norm_slopes <- if (diff(slope_range) == 0) rep(0, n_windows) else
    (slopes - slope_range[1]) / diff(slope_range)
  norm_sds <- if (diff(sd_range) == 0) rep(0, n_windows) else
    (resid_sds - sd_range[1]) / diff(sd_range)

  combined_score <- norm_slopes + norm_sds
  threshold <- quantile(combined_score, probs = quantile_cutoff, na.rm = TRUE)
  selected <- which(combined_score <= threshold)

  selected_flag <- rep(FALSE, n)
  for (s in selected) {
    selected_flag[s:(s + window_sec - 1)] <- TRUE
  }

  # Detect recording periods (gaps > threshold)
  time_gaps <- c(0, diff(t_sec))
  period_id <- cumsum(time_gaps > gap_threshold_min * 60) + 1
  period_labels <- paste("Period", period_id)

  plot_df <- data.frame(
    time_min = t_sec / 60,
    conc = x,
    selected = selected_flag,
    period = factor(period_labels, levels = unique(period_labels))
  )

  # Y-axis label
  if (gas_name == "CO2") {
    y_lab <- expression(CO[2]~concentration~(ppm))
  } else {
    y_lab <- expression(CH[4]~concentration~(ppb))
  }

  p <- ggplot(plot_df, aes(x = time_min, y = conc)) +
    geom_point(aes(color = selected), size = 0.5, alpha = 0.6) +
    scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "#D6604D"),
                       labels = c("FALSE" = "Normal",
                                  "TRUE" = "Selected (bottom 5%)"),
                       name = NULL) +
    facet_wrap(~ period, scales = "free_x") +
    labs(x = "Time (minutes)", y = y_lab,
         title = paste0(instrument_name,
                        ": rolling window quiet-period selection")) +
    plot_theme

  if (!is.null(ylimits)) {
    p <- p + coord_cartesian(ylim = ylimits)
  }

  p
}

# LGR3 (longest HF record)
p_rw_co2 <- plot_rolling_window_selection(
  imp.LGR3, "CO2dry_ppm", "LGR3", "CO2", ylimits = c(400, 500))
p_rw_ch4 <- plot_rolling_window_selection(
  imp.LGR3, "CH4dry_ppb", "LGR3", "CH4", ylimits = c(1000, 2500))

# LGR3-YMF (same instrument, different day)
p_rw_ymf_co2 <- plot_rolling_window_selection(
  ymf_imp, "CO2dry_ppm", "LGR3-YMF", "CO2", ylimits = c(400, 500))
p_rw_ymf_ch4 <- plot_rolling_window_selection(
  ymf_imp, "CH4dry_ppb", "LGR3-YMF", "CH4", ylimits = c(1000, 2500))

# Save all rolling window plots
rw_plots <- list(
  list(p_rw_co2,     "rolling_window_LGR3_CO2"),
  list(p_rw_ch4,     "rolling_window_LGR3_CH4"),
  list(p_rw_ymf_co2, "rolling_window_YMF_CO2"),
  list(p_rw_ymf_ch4, "rolling_window_YMF_CH4")
)

for (pl in rw_plots) {
  if (!is.null(pl[[1]])) {
    n_periods <- length(unique(ggplot_build(pl[[1]])$data[[1]]$PANEL))
    w <- max(6, 4 * n_periods)
    ggsave(file.path(prec_plot_dir, paste0(pl[[2]], ".pdf")),
           pl[[1]], width = w, height = 4)
    ggsave(file.path(prec_plot_dir, paste0(pl[[2]], ".png")),
           pl[[1]], width = w, height = 4, dpi = 300)
  }
}

message("Saved rolling window visualizations to: ", prec_plot_dir)

# --- Plot 4: Quality criteria heatmap ---
# Horizontal bar-style heatmap with % passing each criterion, values in cells.

criteria$criterion <- factor(criteria$criterion,
                              levels = rev(criteria$criterion))

p_heatmap <- ggplot(criteria,
                     aes(x = 1, y = criterion, fill = pct_pass)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = label), size = 3.2, color = "black") +
  scale_fill_gradient2(low = "#D73027", mid = "#FEE08B", high = "#1A9850",
                       midpoint = 50, limits = c(0, 100),
                       name = "% passing") +
  facet_grid(category ~ ., scales = "free_y", space = "free_y",
             switch = "y") +
  labs(x = NULL, y = NULL,
       title = expression(CH[4]~flux~quality~criteria~(all~sites~combined))) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    panel.grid       = element_blank(),
    strip.placement  = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 9),
    legend.position  = "bottom",
    plot.title       = element_text(hjust = 0.5, size = 12)
  )

ggsave(file.path(prec_plot_dir, "quality_criteria_heatmap.pdf"),
       p_heatmap, width = 7, height = 5.5)
ggsave(file.path(prec_plot_dir, "quality_criteria_heatmap.png"),
       p_heatmap, width = 7, height = 5.5, dpi = 300)

message("Saved quality criteria heatmap to: ", prec_plot_dir)

# =============================================================================
# Step H: Component flux distributions under different MDF filters (HF only)
# =============================================================================

message("\n=== Component flux distributions under MDF filtering ===")

# Merge leaf types into a single canopy component
hf_out$component <- case_when(
  grepl("leaf", hf_out$Type, ignore.case = TRUE) ~ "Leaf",
  hf_out$Type == "branch" ~ "Branch",
  TRUE ~ "Stem"
)
hf_out$component <- factor(hf_out$component, levels = c("Stem", "Branch", "Leaf"))

# Define filter scenarios as functions returning logical "fail" vectors
# TRUE = fails criterion (flux set to 0); FALSE = passes (flux kept)
r2_sym <- "\u00B2"
quality_filters <- list()
quality_filters[["No filter"]]                          <- function(df) rep(FALSE, nrow(df))
quality_filters[["Manufacturer MDF"]]                   <- function(df) df$CH4_below_MDF_goflux
quality_filters[[paste0("CH4 R", r2_sym, " > 0.7")]]   <- function(df) df$CH4_best.r2 <= 0.7
quality_filters[["CH4 SNR > 2 (empirical)"]]             <- function(df) df$CH4_SNR <= 2
quality_filters[["Wassmann 90%"]]                       <- function(df) df$CH4_below_MDF_wass90
quality_filters[["Wassmann 95%"]]                       <- function(df) df$CH4_below_MDF_wass95
quality_filters[["CH4 SNR > 3 (empirical)"]]             <- function(df) df$CH4_SNR <= 3
quality_filters[["Wassmann 99%"]]                       <- function(df) df$CH4_below_MDF_wass99
quality_filters[[paste0("CH4 R", r2_sym, " > 0.9")]]   <- function(df) df$CH4_best.r2 <= 0.9
quality_filters[[paste0("CO2 R", r2_sym, " > 0.7")]]   <- function(df) df$CO2_best.r2 <= 0.7
quality_filters[["Christiansen 90%"]]                   <- function(df) df$CH4_below_MDF_chr90
quality_filters[["Christiansen 95%"]]                   <- function(df) df$CH4_below_MDF_chr95
quality_filters[["Christiansen 99%"]]                   <- function(df) df$CH4_below_MDF_chr99

# Sort by number of retained measurements (most to fewest)
n_pass <- sapply(quality_filters, function(fn) sum(!fn(hf_out)))
quality_filters <- quality_filters[order(-n_pass)]
filter_levels <- names(quality_filters)

# Build long-format data: for each filter, set failing fluxes to 0
plot_data <- do.call(rbind, lapply(filter_levels, function(fname) {
  df <- hf_out
  fails <- quality_filters[[fname]](df)
  df$CH4_filtered <- ifelse(fails, 0, df$CH4_best.flux)
  df$zeroed <- fails
  df$filter <- fname
  df[, c("UniqueID", "component", "CH4_best.flux", "CH4_filtered", "zeroed", "filter")]
}))

# Compute pass stats per filter (HF only) for y-axis labels and fill colors
n_hf <- nrow(hf_out)
pass_pct <- sapply(filter_levels, function(fname) {
  fails <- quality_filters[[fname]](hf_out)
  100 * sum(!fails) / n_hf
})
pass_labels <- sapply(filter_levels, function(fname) {
  fails <- quality_filters[[fname]](hf_out)
  n_pass <- sum(!fails)
  pct <- round(100 * n_pass / n_hf, 1)
  paste0(fname, "\n(", n_pass, "/", n_hf, ", ", pct, "%)")
})

# Create labelled factor levels
plot_data$filter <- factor(plot_data$filter, levels = filter_levels,
                           labels = pass_labels)
labelled_levels <- pass_labels

# Compute % retained per filter × component for fill coloring
comp_stats <- plot_data %>%
  group_by(filter, component) %>%
  summarise(mean_flux = mean(CH4_filtered, na.rm = TRUE),
            n = n(),
            n_zeroed = sum(zeroed),
            pct_retained = 100 * (1 - sum(zeroed) / n()),
            .groups = "drop")

# Merge pct_retained back into plot_data for fill mapping
plot_data <- plot_data %>%
  left_join(comp_stats[, c("filter", "component", "pct_retained")],
            by = c("filter", "component"))

# Plot: ggridges — one facet per component, ridgelines by MDF filter
library(ggridges)

# Reverse filter order so "No filter" is at top of y-axis
plot_data$filter <- factor(plot_data$filter, levels = rev(labelled_levels))
comp_stats$filter <- factor(comp_stats$filter, levels = rev(labelled_levels))

# Helper function to build the rainfall ridge plot
# data_in: plot_data (possibly subsetted), stats_in: comp_stats (possibly subsetted)
make_ridge_plot <- function(data_in, stats_in) {
  set.seed(42)
  data_in$y_num <- as.numeric(data_in$filter)
  data_in$y_jitter <- data_in$y_num - runif(nrow(data_in), 0.08, 0.22)

  # asinh scale breaks at nice flux values
  asinh_trans <- scales::trans_new(
    name = "asinh",
    transform = asinh,
    inverse = sinh,
    breaks = scales::extended_breaks()
  )

  p <- ggplot(data_in, aes(x = CH4_best.flux, y = filter)) +
    geom_density_ridges(
      aes(fill = pct_retained, group = filter),
      scale = 0.45, alpha = 0.7, bandwidth = 0.05,
      rel_min_height = 0.005, color = "grey30", linewidth = 0.3
    ) +
    geom_point(aes(x = CH4_filtered, y = y_jitter, color = zeroed),
               shape = 16, size = 1.8, alpha = 0.7) +
    geom_segment(data = stats_in,
                 aes(x = mean_flux, xend = mean_flux,
                     y = as.numeric(filter) - 0.25, yend = as.numeric(filter) + 0.4),
                 color = "#D73027", linewidth = 0.7, linetype = "solid") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    scale_x_continuous(trans = asinh_trans, n.breaks = 5) +
    scale_fill_gradient2(low = "#D73027", mid = "#F7F7F7", high = "#2166AC",
                         midpoint = 50, limits = c(0, 100),
                         name = "% retained") +
    scale_color_manual(values = c("FALSE" = "grey40", "TRUE" = "#D73027"),
                       labels = c("Retained", "Set to 0"),
                       name = NULL) +
    facet_wrap(~ component, nrow = 1, scales = "free_x") +
    labs(x = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
         y = NULL) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    plot_theme +
    theme(strip.text    = element_text(size = 13, face = "bold"),
          strip.background = element_rect(fill = "grey95", color = "grey50"),
          axis.text.x   = element_text(size = 10),
          axis.text.y   = element_text(size = 9, lineheight = 0.85),
          axis.title.x  = element_text(size = 12),
          panel.border  = element_rect(color = "grey50", fill = NA, linewidth = 0.5),
          legend.position = "bottom",
          legend.box = "horizontal")
  p
}

# --- Version 1: All data (with outliers) ---
p_comp_all <- make_ridge_plot(plot_data, comp_stats)

ggsave(file.path(prec_plot_dir, "component_flux_by_mdf_filter.pdf"),
       p_comp_all, width = 13, height = 11)
ggsave(file.path(prec_plot_dir, "component_flux_by_mdf_filter.png"),
       p_comp_all, width = 13, height = 11, dpi = 300)

# --- Version 2: Outliers removed (|flux| > 3× IQR above Q3, per component) ---
outlier_ids <- plot_data %>%
  filter(filter == levels(plot_data$filter)[length(levels(plot_data$filter))]) %>%  # use "No filter" row
  group_by(component) %>%
  mutate(q1 = quantile(CH4_best.flux, 0.25),
         q3 = quantile(CH4_best.flux, 0.75),
         iqr = q3 - q1,
         is_outlier = CH4_best.flux > q3 + 3 * iqr | CH4_best.flux < q1 - 3 * iqr) %>%
  filter(is_outlier) %>%
  pull(UniqueID) %>%
  unique()

message("Outlier measurements removed: ", length(outlier_ids),
        " (UniqueIDs: ", paste(outlier_ids, collapse = ", "), ")")

plot_data_no_outlier <- plot_data %>% filter(!UniqueID %in% outlier_ids)

# Recompute stats without outliers
comp_stats_no_outlier <- plot_data_no_outlier %>%
  group_by(filter, component) %>%
  summarise(mean_flux = mean(CH4_filtered, na.rm = TRUE),
            n = n(),
            n_zeroed = sum(zeroed),
            pct_retained = 100 * (1 - sum(zeroed) / n()),
            .groups = "drop")

# Update pct_retained in the subsetted data
plot_data_no_outlier <- plot_data_no_outlier %>%
  select(-pct_retained) %>%
  left_join(comp_stats_no_outlier[, c("filter", "component", "pct_retained")],
            by = c("filter", "component"))

p_comp_no_outlier <- make_ridge_plot(plot_data_no_outlier, comp_stats_no_outlier)

ggsave(file.path(prec_plot_dir, "component_flux_by_mdf_filter_no_outliers.pdf"),
       p_comp_no_outlier, width = 13, height = 11)
ggsave(file.path(prec_plot_dir, "component_flux_by_mdf_filter_no_outliers.png"),
       p_comp_no_outlier, width = 13, height = 11, dpi = 300)

message("Saved component flux distribution plots to: ", prec_plot_dir)

# Print mean flux table
message("\n--- Mean CH4 flux by component under each MDF filter ---")
mean_wide <- comp_stats %>%
  mutate(label = paste0(round(mean_flux, 3), " (n=", n, ", zeroed=", n_zeroed, ")")) %>%
  select(filter, component, label) %>%
  tidyr::pivot_wider(names_from = component, values_from = label)
print(as.data.frame(mean_wide))

# =============================================================================
# Step I: Sensitivity of per-tissue flux rates to quality filtering
# =============================================================================
# Recreates the per-tissue flux rate bar chart (cf. figure_truncation_combined
# bottom-row panel) under each quality filter criterion, in a single faceted
# plot.  Also quantifies the uncertainty (range) of each component's mean flux
# across filter definitions.

message("\n=== Sensitivity of per-tissue flux rates to quality filtering ===")

# --- Tree geometry (from figure_truncation_combined.R) ---
tree_height_m     <- 26
dbh_m             <- 0.40
dbh_height_m      <- 1.37
base_radius_m     <- dbh_m / 2
height_above_dbh_m <- tree_height_m - dbh_height_m

radius_at_h <- function(h) pmax(0, base_radius_m * (tree_height_m - h) / height_above_dbh_m)
cone_lat_area <- function(h1, h2) {
  r1 <- radius_at_h(h1); r2 <- radius_at_h(h2)
  slant <- sqrt((h2 - h1)^2 + (r1 - r2)^2)
  pi * (r1 + r2) * slant
}

stem_bark_index   <- 0.45
branch_bark_index <- 1.70
leaf_area_index   <- 4.5

total_cone_area <- cone_lat_area(0, tree_height_m)
cone_below_2m   <- cone_lat_area(0, 2)
frac_stem_lt2   <- cone_below_2m / total_cone_area

A_stem_lt2  <- stem_bark_index * frac_stem_lt2
A_stem_ge2  <- stem_bark_index * (1 - frac_stem_lt2)
A_branch    <- branch_bark_index
A_leaf      <- leaf_area_index

area_indices <- c("Stem < 2 m" = A_stem_lt2, "Stem > 2 m" = A_stem_ge2,
                  "Branches" = A_branch, "Leaves" = A_leaf)

comp_labels_4 <- c("Stem < 2 m", "Stem > 2 m", "Branches", "Leaves")

# --- Assign 4-compartment labels (matching truncation figure) ---
hf_out$compartment <- case_when(
  hf_out$Type == "stem"   & hf_out$Height_m <  2 ~ "Stem < 2 m",
  hf_out$Type == "stem"   & hf_out$Height_m >= 2 ~ "Stem > 2 m",
  hf_out$Type == "branch"                        ~ "Branches",
  grepl("leaf", hf_out$Type, ignore.case = TRUE)  ~ "Leaves"
)

message("Compartment counts: ",
        paste(comp_labels_4,
              sapply(comp_labels_4, function(c) sum(hf_out$compartment == c, na.rm = TRUE)),
              sep = "=", collapse = ", "))

# --- Compute per-compartment mean flux under each filter ---
# filter_levels and quality_filters already defined in Step H (ordered by retention)

rate_by_filter <- do.call(rbind, lapply(filter_levels, function(fname) {
  fails <- quality_filters[[fname]](hf_out)
  filtered_flux <- ifelse(fails, 0, hf_out$CH4_best.flux)

  data.frame(
    filter = fname,
    compartment = comp_labels_4,
    flux_rate = sapply(comp_labels_4, function(comp) {
      idx <- hf_out$compartment == comp & !is.na(hf_out$compartment)
      if (sum(idx) == 0) return(NA_real_)
      mean(filtered_flux[idx], na.rm = TRUE)
    }),
    n_meas = sapply(comp_labels_4, function(comp) {
      sum(hf_out$compartment == comp & !is.na(hf_out$compartment))
    }),
    n_pass = sapply(comp_labels_4, function(comp) {
      idx <- hf_out$compartment == comp & !is.na(hf_out$compartment)
      sum(!fails[idx])
    }),
    area_index = as.numeric(area_indices[comp_labels_4]),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}))

rate_by_filter$integrated <- rate_by_filter$flux_rate * rate_by_filter$area_index
rate_by_filter$pct_pass   <- round(100 * rate_by_filter$n_pass / rate_by_filter$n_meas, 1)

# Use the labelled filter levels (with stats) consistent with ridge plot
rate_by_filter$filter <- factor(rate_by_filter$filter,
                                levels = filter_levels,
                                labels = pass_labels)
rate_by_filter$compartment <- factor(rate_by_filter$compartment,
                                     levels = comp_labels_4)

# Colors matching truncation figure
compartment_colors_4 <- c(
  "Stem < 2 m" = "#8B4513",
  "Stem > 2 m" = "#D4A76A",
  "Branches"   = "#4682B4",
  "Leaves"     = "#2E8B57"
)

# --- Side-by-side stacked budget bars (one per filter), two rows ---
# Row 1: outliers retained.  Row 2: outliers removed.
# Mirrors figure_sensitivity_breakeven panel (c): integrated flux = rate × area

# Outlier detection (3× IQR per compartment) — done here before the uncertainty
# section so both plots can use the same outlier IDs.
outlier_ids_4comp <- hf_out %>%
  filter(!is.na(compartment)) %>%
  group_by(compartment) %>%
  mutate(q1 = quantile(CH4_best.flux, 0.25),
         q3 = quantile(CH4_best.flux, 0.75),
         iqr = q3 - q1,
         is_outlier = CH4_best.flux > q3 + 3 * iqr |
                      CH4_best.flux < q1 - 3 * iqr) %>%
  filter(is_outlier) %>%
  pull(UniqueID) %>%
  unique()

message("Outlier measurements (4-comp, 3\u00d7IQR): ", length(outlier_ids_4comp),
        " (", paste(outlier_ids_4comp, collapse = ", "), ")")

hf_no_outlier <- hf_out %>% filter(!UniqueID %in% outlier_ids_4comp)

# Compute rate_by_filter WITHOUT outliers
rate_by_filter_no_out <- do.call(rbind, lapply(filter_levels, function(fname) {
  fails <- quality_filters[[fname]](hf_no_outlier)
  filtered_flux <- ifelse(fails, 0, hf_no_outlier$CH4_best.flux)
  data.frame(
    filter = fname,
    compartment = comp_labels_4,
    flux_rate = sapply(comp_labels_4, function(comp) {
      idx <- hf_no_outlier$compartment == comp & !is.na(hf_no_outlier$compartment)
      if (sum(idx) == 0) return(NA_real_)
      mean(filtered_flux[idx], na.rm = TRUE)
    }),
    area_index = as.numeric(area_indices[comp_labels_4]),
    stringsAsFactors = FALSE, row.names = NULL
  )
}))
rate_by_filter_no_out$integrated <- rate_by_filter_no_out$flux_rate *
                                     rate_by_filter_no_out$area_index
rate_by_filter_no_out$filter <- factor(rate_by_filter_no_out$filter,
                                        levels = filter_levels,
                                        labels = pass_labels)
rate_by_filter_no_out$compartment <- factor(rate_by_filter_no_out$compartment,
                                             levels = comp_labels_4)

# --- Helper: build one budget stacked-bar panel ---
make_budget_panel <- function(rbf, show_x_labels = TRUE, title_label = NULL) {
  bdata <- rbf
  bdata$compartment <- factor(bdata$compartment, levels = rev(comp_labels_4))

  # Short x-axis labels — map from the labelled filter levels already on the data
  n_hf_total <- nrow(hf_out)
  short_labels <- sapply(filter_levels, function(fname) {
    fails <- quality_filters[[fname]](hf_out)
    pct <- round(100 * sum(!fails) / n_hf_total, 0)
    paste0(fname, "\n(", pct, "%)")
  })
  # rbf$filter already uses pass_labels as levels, so match on those
  bdata$filter_short <- factor(bdata$filter, levels = levels(bdata$filter),
                                labels = short_labels)

  totals <- bdata %>%
    group_by(filter_short) %>%
    summarise(total = sum(integrated, na.rm = TRUE), .groups = "drop")

  p <- ggplot(bdata, aes(x = filter_short, y = integrated, fill = compartment)) +
    geom_col(width = 0.7, color = "grey30", linewidth = 0.2) +
    geom_text(data = totals,
              aes(x = filter_short, y = total, label = sprintf("%.3f", total)),
              inherit.aes = FALSE, vjust = -0.3, size = 2.8, color = "grey20") +
    scale_fill_manual(values = compartment_colors_4,
                      breaks = comp_labels_4, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = NULL,
         y = expression("Integrated flux (nmol m"^{-2}*"ground s"^{-1}*")"))

  if (!is.null(title_label)) {
    p <- p + ggtitle(title_label)
  }

  x_theme <- if (show_x_labels) {
    theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1, lineheight = 0.85))
  } else {
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }

  p <- p +
    theme_classic(base_size = 11) +
    theme(
      axis.text.y     = element_text(size = 10),
      axis.title      = element_text(size = 11),
      legend.position = "none",
      plot.title      = element_text(size = 11, face = "bold"),
      plot.margin     = margin(5, 15, 5, 10)
    ) + x_theme

  p
}

p_budget_retained <- make_budget_panel(rate_by_filter,
                                        show_x_labels = FALSE,
                                        title_label = "Outliers retained")
p_budget_removed  <- make_budget_panel(rate_by_filter_no_out,
                                        show_x_labels = TRUE,
                                        title_label = "Outliers removed")

# Extract legend from one panel
p_legend_src <- ggplot(rate_by_filter, aes(x = filter, y = integrated,
                                            fill = factor(compartment,
                                              levels = rev(comp_labels_4)))) +
  geom_col() +
  scale_fill_manual(values = compartment_colors_4, breaks = comp_labels_4, name = NULL) +
  theme(legend.position = "top", legend.text = element_text(size = 9))

library(patchwork)

# Force shared y-axis range
y_max <- max(
  tapply(rate_by_filter$integrated, rate_by_filter$filter, sum),
  tapply(rate_by_filter_no_out$integrated, rate_by_filter_no_out$filter, sum),
  na.rm = TRUE
) * 1.1

p_budget_retained <- p_budget_retained +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, y_max))
p_budget_removed  <- p_budget_removed +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, y_max)) +
  theme(legend.position = "bottom")

p_budget_combined <- p_budget_retained / p_budget_removed

ggsave(file.path(prec_plot_dir, "budget_stacked_by_filter.pdf"),
       p_budget_combined, width = 12, height = 8)
ggsave(file.path(prec_plot_dir, "budget_stacked_by_filter.png"),
       p_budget_combined, width = 12, height = 8, dpi = 300)

message("Saved integrated budget stacked bar chart (both outlier versions)")

# --- Uncertainty quantification across filters ---
# (outlier_ids_4comp, hf_no_outlier, rate_by_filter_no_out already computed above)

# --- Combine both versions for the uncertainty plot ---

rate_by_filter$version <- "Outliers retained"
rate_by_filter_no_out$version <- "Outliers removed"

# Normalize each version to its OWN no-filter baseline
baseline_retained <- rate_by_filter %>%
  filter(filter == levels(rate_by_filter$filter)[1]) %>%
  select(compartment, baseline_flux = flux_rate)

baseline_removed <- rate_by_filter_no_out %>%
  filter(filter == levels(rate_by_filter_no_out$filter)[1]) %>%
  select(compartment, baseline_flux = flux_rate)

rate_by_filter <- rate_by_filter %>%
  left_join(baseline_retained, by = "compartment") %>%
  mutate(pct_of_baseline_flux = 100 * flux_rate / baseline_flux)

rate_by_filter_no_out <- rate_by_filter_no_out %>%
  left_join(baseline_removed, by = "compartment") %>%
  mutate(pct_of_baseline_flux = 100 * flux_rate / baseline_flux)

# Combine — keep both absolute flux_rate and pct columns
keep_cols <- c("compartment", "flux_rate", "pct_of_baseline_flux", "version")
uncert_combined <- bind_rows(
  rate_by_filter[, keep_cols],
  rate_by_filter_no_out[, keep_cols]
)
uncert_combined$version <- factor(uncert_combined$version,
                                   levels = c("Outliers retained", "Outliers removed"))

# Sensitivity summary for both versions
sens_both <- uncert_combined %>%
  group_by(compartment, version) %>%
  summarise(
    min_pct   = min(pct_of_baseline_flux),
    max_pct   = max(pct_of_baseline_flux),
    range_pct = max_pct - min_pct,
    cv_pct    = 100 * sd(pct_of_baseline_flux) / mean(pct_of_baseline_flux),
    .groups   = "drop"
  )

message("\n--- Relative sensitivity (% of own no-filter baseline) ---")
print(as.data.frame(sens_both))

# --- Two-row uncertainty plot: absolute (top) + relative (bottom) ---

# Shared theme elements
uncert_theme <- theme_classic(base_size = 12) +
  theme(
    axis.text    = element_text(size = 11),
    axis.title   = element_text(size = 12),
    plot.title   = element_text(size = 11, face = "bold"),
    legend.position = "none"
  )

# Precompute medians, split by version for separate geom_point calls
medians_all <- uncert_combined %>%
  group_by(compartment, version) %>%
  summarise(median_flux = median(flux_rate, na.rm = TRUE),
            median_pct  = median(pct_of_baseline_flux, na.rm = TRUE),
            .groups = "drop")
medians_retained <- medians_all %>% filter(version == "Outliers retained")
medians_removed  <- medians_all %>% filter(version == "Outliers removed")

# Dodge offset: retained left (-0.125), removed right (+0.125)
dodge_w <- 0.5

# Row 1: Absolute flux rates
p_abs <- ggplot(uncert_combined,
                 aes(x = compartment, y = flux_rate,
                     color = compartment, shape = version)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  stat_summary(fun.min = min, fun.max = max, geom = "linerange",
               linewidth = 1.2, alpha = 0.3,
               position = position_dodge(width = dodge_w)) +
  geom_point(size = 2.5, alpha = 0.7,
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = dodge_w)) +
  # Median outlines: open circle for retained, open triangle for removed
  geom_point(data = medians_retained,
             aes(x = compartment, y = median_flux),
             shape = 1, color = "black", size = 5, stroke = 1,
             inherit.aes = FALSE,
             position = position_nudge(x = -dodge_w / 4)) +
  geom_point(data = medians_removed,
             aes(x = compartment, y = median_flux),
             shape = 2, color = "black", size = 5, stroke = 1,
             inherit.aes = FALSE,
             position = position_nudge(x = dodge_w / 4)) +
  scale_color_manual(values = compartment_colors_4, guide = "none") +
  scale_shape_manual(values = c("Outliers retained" = 16, "Outliers removed" = 17),
                     name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  labs(x = NULL,
       y = expression(CH[4]~flux~rate~(nmol~m^{-2}~s^{-1})),
       title = "Absolute") +
  uncert_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Row 2: Relative (% of no-filter baseline)
p_rel <- ggplot(uncert_combined,
                 aes(x = compartment, y = pct_of_baseline_flux,
                     color = compartment, shape = version)) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  stat_summary(fun.min = min, fun.max = max, geom = "linerange",
               linewidth = 1.2, alpha = 0.3,
               position = position_dodge(width = dodge_w)) +
  geom_point(size = 2.5, alpha = 0.7,
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = dodge_w)) +
  geom_point(data = medians_retained,
             aes(x = compartment, y = median_pct),
             shape = 1, color = "black", size = 5, stroke = 1,
             inherit.aes = FALSE,
             position = position_nudge(x = -dodge_w / 4)) +
  geom_point(data = medians_removed,
             aes(x = compartment, y = median_pct),
             shape = 2, color = "black", size = 5, stroke = 1,
             inherit.aes = FALSE,
             position = position_nudge(x = dodge_w / 4)) +
  scale_color_manual(values = compartment_colors_4, guide = "none") +
  scale_shape_manual(values = c("Outliers retained" = 16, "Outliers removed" = 17),
                     name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)),
                     labels = function(x) paste0(round(x), "%")) +
  labs(x = NULL,
       y = "% of no-filter baseline",
       title = "Relative",
       caption = "Filled = individual filters; outline = median; bar = range") +
  uncert_theme +
  theme(plot.caption = element_text(size = 9, color = "grey40"))

# Set legend on bottom panel only (top has legend.position = "none" from uncert_theme)
p_rel <- p_rel + theme(legend.position = "bottom")

p_uncertainty <- p_abs / p_rel

ggsave(file.path(prec_plot_dir, "flux_rate_uncertainty_by_component.pdf"),
       p_uncertainty, width = 8, height = 8)
ggsave(file.path(prec_plot_dir, "flux_rate_uncertainty_by_component.png"),
       p_uncertainty, width = 8, height = 8, dpi = 300)

message("Saved two-row uncertainty plot (absolute + relative)")

# --- Print integrated budget sensitivity ---
budget_by_filter <- rate_by_filter %>%
  group_by(filter) %>%
  summarise(total_integrated = sum(integrated, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(pct_of_unfiltered = 100 * total_integrated /
           total_integrated[filter == levels(rate_by_filter$filter)[1]])

message("\n--- Total integrated CH4 budget by filter ---")
print(as.data.frame(budget_by_filter))

message("\nStep I complete.")
