# Detection limits for the Harvard Forest and Yale Myers closures, added to the fitted flux tables.
#
# Approaches (all kept as columns; the empirical 95 % MDF is the reference used in the analysis):
#   1. Manufacturer comparison (column suffix `goflux`): MDF = datasheet precision / t x flux.term (k = 1;
#      Christiansen et al. 2015).
#   2. Empirical (REFERENCE; column suffix `wass*` retained for compatibility): MDF = z x sigma / t x flux.term,
#      z = 1.96 for 95 %; sigma = median over closures of MAD(second differences) / sqrt(6) (MAD scaled by 1.4826),
#      per analyzer x field day x logging interval. This is the same threshold used inside the fits.
#   3. x3 t-critical variant (suffix `chr*`): MDF = (per-closure first-difference SD x 3 x t) / t_closure x flux.term.
#      This is a quantification-limit style comparison, not the Christiansen et al. (2015) MDF.
# t is the closure length in seconds (span of the window + one logging interval), not the number of observations.
# Logging intervals: LGR1 1 s, LGR2 10 s, LGR3 1 s, Yale Myers 5 s.
# Whole-record first-difference MAD and a rolling-window scan of the continuous record are reported as
# precision comparisons only (precision_comparison.csv).
# Instrument: ABB GLA131-GGA (microportable UGGA); datasheet precision CO2 0.35 ppm, CH4 0.9 ppb (1 sigma, 1 s).
#
# Inputs: flux/out/*_fluxes_fitted.csv and flux/cache/ from 01_harvard_forest.R and 02_yale_myers.R.
# Outputs: flux/out/harvard_forest_fluxes.csv, flux/out/yale_myers_fluxes.csv, flux/out/precision_comparison.csv
source(file.path("flux", "00_setup.R"))

hf_compiled <- read.csv(file.path(FLUX_OUT, "harvard_forest_fluxes_fitted.csv"), stringsAsFactors = FALSE)
ymf_compiled <- read.csv(file.path(FLUX_OUT, "yale_myers_fluxes_fitted.csv"), stringsAsFactors = FALSE)
hf_fits <- readRDS(file.path(CACHE, "harvard_forest.rds"))
ymf_fit <- readRDS(file.path(CACHE, "yale_myers.rds"))$YMF
for (lgr in c("LGR1", "LGR2", "LGR3")) {
  assign(paste0("manID.", lgr), hf_fits[[lgr]]$traces)
  assign(paste0("imp.", lgr), hf_fits[[lgr]]$imp)
}
manID.YMF <- ymf_fit$traces
imp.YMF <- ymf_fit$imp
manID_all_hf <- rbind(manID.LGR1, manID.LGR2, manID.LGR3)
lgr1_uids <- unique(manID.LGR1$UniqueID)
lgr2_uids <- unique(manID.LGR2$UniqueID)
lgr3_uids <- unique(manID.LGR3$UniqueID)
hf_compiled$instrument <- hf_compiled$Analyzer
ymf_compiled$instrument <- "YMF"
stopifnot(!anyNA(hf_compiled$instrument))

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
# estimates white noise around a linear trend. Curvature or variable slopes
# can inflate it; this SD-based first-difference diagnostic is not the reference
# estimator, and the recorded interval is not always 1 s.
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

# --- Closure length (seconds) and logging interval per measurement ---
# t must be the closure length in SECONDS, not the number of observations:
# LGR2 logs at 10 s and the YMF record at 5 s, so nb.obs underestimates t
# there by 10× / 5×.  dt = median spacing of POSIX.time within the flagged
# window; t = span + dt (equals goFlux's max(Etime) + 1 at 1 Hz).
compute_closure_per_uid <- function(manID_df) {
  manID_df %>%
    filter(flag == 1) %>%
    group_by(UniqueID) %>%
    summarise(
      dt_s  = median(diff(as.numeric(POSIX.time))),
      t_sec = diff(range(as.numeric(POSIX.time))) + dt_s,
      .groups = "drop"
    )
}

closure_hf  <- compute_closure_per_uid(manID_all_hf)
closure_ymf <- compute_closure_per_uid(manID.YMF)

# Logging interval per instrument (verify: LGR1 = 1, LGR2 = 10, LGR3 = 1, YMF = 5)
closure_hf$instrument <- NA_character_
closure_hf$instrument[closure_hf$UniqueID %in% lgr1_uids] <- "LGR1"
closure_hf$instrument[closure_hf$UniqueID %in% lgr2_uids] <- "LGR2"
closure_hf$instrument[closure_hf$UniqueID %in% lgr3_uids] <- "LGR3"
closure_ymf$instrument <- "YMF"

dt_by_inst <- bind_rows(closure_hf, closure_ymf) %>%
  group_by(instrument) %>%
  summarise(dt_s     = median(dt_s),
            dt_min   = min(dt_s),
            dt_max   = max(dt_s),
            t_median = median(t_sec),
            n_meas   = n(),
            .groups = "drop")

message("\n--- Logging interval (dt, s) and median closure length per instrument ---")
print(as.data.frame(dt_by_inst))
if (any(dt_by_inst$dt_min != dt_by_inst$dt_max)) {
  warning("Logging interval is not constant within an instrument; check dt_by_inst")
}

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

# YMF (object name differs between setups)
ymf_imp <- imp.YMF

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

# =============================================================================
# Step B, Method 3: legacy whole-record first-difference comparison
# =============================================================================

message("\n=== Campaign precision: MAD of first differences ===")

# Same first-difference logic as allan_sd(), but using the median absolute
# deviation (× 1.4826 for normal consistency) so that real concentration
# jumps (chamber placement, purging, gaps between campaigns) do not inflate
# the estimate.  Computed ONCE over the whole continuous imported record of
# each analyzer/campaign; retained only for comparison with earlier outputs.
precision_mad <- function(x, constant = 1.4826) {
  x <- x[!is.na(x)]
  d <- diff(x)
  if (length(d) < 2) return(NA_real_)
  mad(d, constant = constant) / sqrt(2)
}

# ymf_imp was resolved in the rolling-window step above
imp_list <- list(LGR1 = imp.LGR1, LGR2 = imp.LGR2, LGR3 = imp.LGR3, YMF = ymf_imp)

# Guard: a record whose logging interval changed mid-way would mix two noise
# scales in one MAD. Split the record at interval changes and warn if more
# than one interval is present (a guard against mixed logging intervals).
precision_mad_runs <- function(x, posix, constant = 1.4826, tol = 0.2) {
  ok <- is.finite(x) & !is.na(posix); x <- x[ok]; posix <- posix[ok]
  d <- diff(as.numeric(posix))
  brk <- c(TRUE, abs(diff(d)) > tol * pmax(d[-1], 1e-9))
  run <- c(1, cumsum(brk[-1]) + 1)
  out <- do.call(rbind, lapply(split(seq_along(d), run), function(i) {
    if (length(i) < 3) return(NULL)
    data.frame(dt_s = round(median(d[i]) * 2) / 2, n = length(i) + 1,   # nearest 0.5 s
               mad_sigma = mad(diff(x)[i], constant = constant) / sqrt(2))
  }))
  out <- do.call(rbind, lapply(split(out, out$dt_s), function(r)
    if (sum(r$n) < 4) NULL else                       # an interval seen < 3 times is a gap, not a run
    data.frame(dt_s = r$dt_s[1], n = sum(r$n), mad_sigma = weighted.mean(r$mad_sigma, r$n))))
  if (is.null(out)) return(data.frame(dt_s = NA_real_, n = 0L, mad_sigma = NA_real_))
  if (nrow(out) > 1) warning("record has ", nrow(out), " distinct logging intervals: ",
                             paste(out$dt_s, collapse = ", "), " s")
  out
}

mad_results <- do.call(rbind, lapply(names(imp_list), function(nm) {
  d <- imp_list[[nm]]
  ch4 <- precision_mad_runs(d$CH4dry_ppb, d$POSIX.time)
  co2 <- precision_mad_runs(d$CO2dry_ppm, d$POSIX.time)
  data.frame(instrument = nm,
             mad_sd_CO2 = co2$mad_sigma[match(ch4$dt_s, co2$dt_s)],
             mad_sd_CH4 = ch4$mad_sigma,
             n_record   = ch4$n,
             dt_run     = ch4$dt_s,
             stringsAsFactors = FALSE)
}))
if (any(duplicated(mad_results$instrument)))
  stop("a canopy record has more than one logging interval; handle per-run sigma explicitly")
mad_results$dt_run <- NULL
rownames(mad_results) <- NULL

# Attach logging interval and datasheet spec-equivalent at that interval
# (white-noise assumption: σ_dt = σ_1s / sqrt(dt)) so 5–10 s loggers are
# compared under an explicitly hypothetical interval-averaging assumption.
# Timestamp spacing alone does not establish averaging; empirical sigma needs no rescaling.
prec_co2 <- 0.35  # ppm, GLA131 datasheet (1σ, 1 s)
prec_ch4 <- 0.9   # ppb, GLA131 datasheet (1σ, 1 s)

mad_results <- mad_results %>%
  left_join(dt_by_inst[, c("instrument", "dt_s")], by = "instrument") %>%
  mutate(spec_equiv_CO2 = prec_co2 / sqrt(dt_s),
         spec_equiv_CH4 = prec_ch4 / sqrt(dt_s))

message("\n--- Campaign MAD precision (whole imported record) ---")
print(mad_results)

message("\n--- Comparison: Manufacturer vs Empirical precision ---")
message("Manufacturer (GLA131): CO2 = 0.35 ppm, CH4 = 0.9 ppb (1σ, 1 s)")
for (i in seq_len(nrow(rw_results))) {
  message(rw_results$instrument[i], " (dt = ", mad_results$dt_s[i], " s): ",
          "MAD CO2 = ", round(mad_results$mad_sd_CO2[i], 4), " ppm, ",
          "MAD CH4 = ", round(mad_results$mad_sd_CH4[i], 4), " ppb | ",
          "rolling CO2 = ", round(rw_results$flat_sd_CO2[i], 4), " ppm, ",
          "rolling CH4 = ", round(rw_results$flat_sd_CH4[i], 4), " ppb | ",
          "spec-equiv CO2 = ", round(mad_results$spec_equiv_CO2[i], 4), " ppm, ",
          "spec-equiv CH4 = ", round(mad_results$spec_equiv_CH4[i], 4), " ppb")
}

# Reference group precision is computed within closures, then grouped by day and interval.
noise_all <- bind_rows(lapply(c("LGR1", "LGR2", "LGR3", "YMF"), function(inst) {
  estimate_precision(get(paste0("manID.", inst)), inst, warn = FALSE)$closures
}))
noise_ref <- noise_all %>%
  select(UniqueID, instrument, field_day, dt_s, gas, sigma_group) %>%
  pivot_wider(names_from = gas, values_from = sigma_group, names_prefix = "empirical_sd_")
stopifnot(!anyDuplicated(noise_ref$UniqueID))

# =============================================================================
# Step C: Compute three MDF approaches
# =============================================================================

message("\n=== Computing three MDF approaches ===")

# GLA131 manufacturer precision (1σ, 1s) — prec_co2 / prec_ch4 defined above

# --- Helper: compute MDFs at multiple confidence levels ---
#
# Empirical multipliers (two-sided z-scores for normal distribution):
#   99%: 2.576,  95%: 1.960 (REFERENCE),  90%: 1.645
# ×3 t-crit variant multipliers (per-closure SD × 3 × t_α):
#   99%: qt(0.995, df),  95%: qt(0.975, df),  90%: qt(0.95, df)

z_mult <- c("99" = qnorm(0.995), "95" = qnorm(0.975), "90" = qnorm(0.95))

compute_mdf_all <- function(compiled, allan_df, closure_df, rw_df, mad_df) {

  # Merge Allan deviation per measurement
  compiled <- merge(compiled, allan_df[, c("UniqueID", "allan_sd_CO2",
                                            "allan_sd_CH4", "n_meas_pts")],
                    by = "UniqueID", all.x = TRUE)

  # Merge closure length (s) and logging interval per measurement
  compiled <- merge(compiled, closure_df[, c("UniqueID", "dt_s", "t_sec")],
                    by = "UniqueID", all.x = TRUE)

  # Merge rolling window precision per instrument (comparison only)
  compiled <- merge(compiled,
                    rw_df[, c("instrument", "flat_sd_CO2", "flat_sd_CH4")],
                    by = "instrument", all.x = TRUE)

  # Merge legacy whole-record MAD precision (comparison only)
  compiled <- merge(compiled,
                    mad_df[, c("instrument", "mad_sd_CO2", "mad_sd_CH4")],
                    by = "instrument", all.x = TRUE)

  compiled <- compiled %>%
    left_join(noise_ref, by = c("UniqueID", "instrument", "dt_s"))
  stopifnot(all(is.finite(compiled$empirical_sd_CH4)), all(is.finite(compiled$empirical_sd_CO2)))
  compiled <- compiled %>%
    mutate(
      # --- Approach 1: manufacturer comparison (k = 1, 1 s spec) ---
      CO2_MDF_goflux = prec_co2 / t_sec * CO2_flux.term,
      CH4_MDF_goflux = prec_ch4 / t_sec * CH4_flux.term,

      # --- Approach 2: Empirical (z × median closure second-difference σ) at 99/95/90% ---
      # 95 % is the reference method; `wass` suffix kept for compatibility
      CO2_MDF_wass99 = (z_mult["99"] * empirical_sd_CO2) / t_sec * CO2_flux.term,
      CH4_MDF_wass99 = (z_mult["99"] * empirical_sd_CH4) / t_sec * CH4_flux.term,
      CO2_MDF_wass95 = (z_mult["95"] * empirical_sd_CO2) / t_sec * CO2_flux.term,
      CH4_MDF_wass95 = (z_mult["95"] * empirical_sd_CH4) / t_sec * CH4_flux.term,
      CO2_MDF_wass90 = (z_mult["90"] * empirical_sd_CO2) / t_sec * CO2_flux.term,
      CH4_MDF_wass90 = (z_mult["90"] * empirical_sd_CH4) / t_sec * CH4_flux.term,

      # --- Approach 3: ×3 t-crit variant (per-closure SD × 3 × t_α) ---
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
hf_out <- compute_mdf_all(hf_compiled, allan_hf, closure_hf, rw_results, mad_results)

# Apply to YMF
ymf_out <- compute_mdf_all(ymf_compiled, allan_ymf, closure_ymf, rw_results, mad_results)

# Sanity check: at 1 Hz t_sec should equal nb.obs; at 10 s / 5 s it is 10× / 5×
message("\n--- t_sec / nb.obs ratio by instrument (expect = dt_s) ---")
print(hf_out %>% group_by(instrument) %>%
        summarise(dt_s = median(dt_s),
                  ratio_median = median(t_sec / CH4_nb.obs, na.rm = TRUE),
                  .groups = "drop") %>% as.data.frame())
message("YMF: dt_s = ", median(ymf_out$dt_s), ", t_sec/nb.obs median = ",
        round(median(ymf_out$t_sec / ymf_out$CH4_nb.obs, na.rm = TRUE), 3))

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
          goflux = "Manufacturer (1σ)",
          wass99 = "Empirical 99%", wass95 = "Empirical 95%", wass90 = "Empirical 90%",
          chr99  = "x3 t-crit 99%", chr95 = "x3 t-crit 95%", chr90 = "x3 t-crit 90%")
        message(sprintf("    MDF range (%-20s): %s - %s", label,
                round(min(out[[col]], na.rm = TRUE), 6),
                round(max(out[[col]], na.rm = TRUE), 6)))
      }
    }

    # Detection rates
    message("    --- % below detection ---")
    below_g <- paste0(gas, "_below_MDF_goflux")
    message(sprintf("    Manufacturer (1σ):    %d / %d (%.1f%%)",
            sum(out[[below_g]], na.rm = TRUE), n_valid,
            100 * mean(out[[below_g]], na.rm = TRUE)))

    for (conf in c("99", "95", "90")) {
      wass_col <- paste0(gas, "_below_MDF_wass", conf)
      chr_col  <- paste0(gas, "_below_MDF_chr", conf)
      if (wass_col %in% names(out)) {
        message(sprintf("    Empirical %s%%:             %d / %d (%.1f%%)", conf,
                sum(out[[wass_col]], na.rm = TRUE), n_valid,
                100 * mean(out[[wass_col]], na.rm = TRUE)))
      }
      if (chr_col %in% names(out)) {
        message(sprintf("    x3 t-crit %s%%:             %d / %d (%.1f%%)", conf,
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

write.csv(hf_out, file.path(FLUX_OUT, "harvard_forest_fluxes.csv"), row.names = FALSE)
write.csv(ymf_out, file.path(FLUX_OUT, "yale_myers_fluxes.csv"), row.names = FALSE)

# --- Save precision summary ---

# Rows: datasheet (1 s), datasheet spec-equivalent at each instrument's
# logging interval (hypothetical averaging), legacy campaign MAD,
# rolling window and median Allan deviation (comparison only).
precision_summary <- data.frame(
  source = c("GLA131 manufacturer (1σ, 1s)",
             mad_results$instrument,
             mad_results$instrument,
             rw_results$instrument),
  CO2_precision_ppm = c(prec_co2, mad_results$spec_equiv_CO2,
                        mad_results$mad_sd_CO2, rw_results$flat_sd_CO2),
  CH4_precision_ppb = c(prec_ch4, mad_results$spec_equiv_CH4,
                        mad_results$mad_sd_CH4, rw_results$flat_sd_CH4),
  dt_s   = c(1, mad_results$dt_s, mad_results$dt_s, mad_results$dt_s),
  method = c("Datasheet",
             rep("Datasheet if interval-averaged (spec/sqrt(dt))", nrow(mad_results)),
             rep("Campaign MAD (first differences)", nrow(mad_results)),
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
  dt_s   = mad_results$dt_s[match(allan_by_inst$instrument, mad_results$instrument)],
  method = "Allan deviation (median)",
  stringsAsFactors = FALSE
)

reference_summary <- noise_ref %>% group_by(instrument, field_day, dt_s) %>%
  summarise(CO2_precision_ppm = first(empirical_sd_CO2), CH4_precision_ppb = first(empirical_sd_CH4), .groups = "drop") %>%
  transmute(source = paste(instrument, field_day), CO2_precision_ppm, CH4_precision_ppb, dt_s,
            method = "Reference: median closure second differences")
precision_summary <- rbind(precision_summary, allan_summary, as.data.frame(reference_summary))

# Rename for display
precision_summary$source <- sub("^GLA131.*", "GLA131 spec", precision_summary$source)
precision_summary$source <- sub("^YMF$", "LGR3-YMF", precision_summary$source)

write.csv(precision_summary,
          file.path(FLUX_OUT, "precision_comparison.csv"),
          row.names = FALSE)

message("\n--- Precision comparison table ---")
print(precision_summary)

message("Detection limits written to ", FLUX_OUT, ": harvard_forest_fluxes.csv, yale_myers_fluxes.csv, precision_comparison.csv")
