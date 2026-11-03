# Standalone empirical noise and MDF implementation. No unreleased packages.
# Per closure: normal-consistent MAD of second differences / sqrt(6).
# Group: median over closures of an analyzer x field day x logging interval.
# This robust second-difference estimator cancels a linear concentration trend.
# It estimates short-timescale noise, not drift, leakage, or total flux error.

closure_noise <- function(x, time, tol = 0.2) {
  tm <- as.numeric(time)
  if (length(x) != length(tm) || length(x) < 3) stop("Invalid closure series")
  dt <- diff(tm)
  nominal <- median(dt[is.finite(dt) & dt > 0])
  if (!is.finite(nominal)) stop("No positive logging interval")
  interval <- if (nominal >= .5) round(nominal * 2) / 2 else signif(nominal, 2)
  # Preserve row adjacency: missing values, duplicates, backwards timestamps,
  # gaps and interval changes never create a synthetic difference.
  valid <- is.finite(dt) & dt > 0 & abs(dt - nominal) <= tol * nominal &
    is.finite(x[-1]) & is.finite(x[-length(x)])
  dx <- diff(x)
  good2 <- valid[-1] & valid[-length(valid)] &
    abs(diff(dt)) <= tol * pmax(dt[-1], 1e-9)
  ix <- which(good2)
  d2 <- diff(dx)[ix]
  d1 <- dx[valid]
  sig <- if (length(d2) >= 6) mad(d2, constant = 1.4826) / sqrt(6) else NA_real_
  adjacent <- which(diff(ix) == 1)
  ac <- if (length(adjacent) >= 5 && sd(d2[adjacent]) > 0 && sd(d2[adjacent+1]) > 0)
    cor(d2[adjacent], d2[adjacent+1]) else NA_real_
  first <- if (length(d1) >= 3) mad(d1, constant = 1.4826) / sqrt(2) else NA_real_
  list(sigma = sig, dt_s = interval, n = sum(is.finite(x)), n_diff2 = length(d2),
       n_diff1 = length(d1), n_zero = sum(d1 == 0), ac1 = ac,
       zero_frac = mean(d1 == 0), sigma_d1c = first,
       d1c_ratio = first/sig, d1c = d1 - median(d1),
       allan_sd = if (length(d1) >= 2) sd(d1)/sqrt(2) else NA_real_,
       t_sec = diff(range(tm[is.finite(tm)])) + nominal)
}

estimate_precision <- function(traces, instrument, warn = TRUE) {
  needed <- c("UniqueID", "flag", "POSIX.time", "CH4dry_ppb", "CO2dry_ppm")
  stopifnot(all(needed %in% names(traces)))
  tr <- traces[!is.na(traces$flag) & traces$flag == 1, ]
  ids <- split(seq_len(nrow(tr)), as.character(tr$UniqueID))
  rows <- list(); centered <- list()
  for (id in names(ids)) for (gas in c("CH4", "CO2")) {
    d <- tr[ids[[id]], ]
    col <- if (gas == "CH4") "CH4dry_ppb" else "CO2dry_ppm"
    s <- closure_noise(d[[col]], d$POSIX.time)
    key <- paste(id, gas, sep = "|")
    # Imported timestamps represent the analyzer's displayed local clock.
    day <- format(min(d$POSIX.time), "%Y-%m-%d", tz = "UTC")
    rows[[key]] <- data.frame(UniqueID = id, instrument = instrument, field_day = day,
      gas = gas, as.data.frame(s[setdiff(names(s), "d1c")]), stringsAsFactors = FALSE)
    centered[[key]] <- s$d1c
  }
  closures <- do.call(rbind, rows); rownames(closures) <- NULL
  key <- with(closures, paste(instrument, field_day, dt_s, gas, sep = "|"))
  groups <- lapply(split(seq_len(nrow(closures)), key), function(i) {
    r <- closures[i, ]; d <- unlist(centered[paste(r$UniqueID,r$gas,sep="|")])
    finite_median <- function(x) if (any(is.finite(x))) median(x[is.finite(x)]) else NA_real_
    sig <- finite_median(r$sigma); ac <- finite_median(r$ac1)
    zf <- sum(r$n_zero)/sum(r$n_diff1)
    first <- mad(d, constant = 1.4826)/sqrt(2)
    data.frame(r[1,c("instrument","field_day","dt_s","gas")], sigma = sig,
      n_closures = sum(is.finite(r$sigma)), n_diff2 = sum(r$n_diff2), ac1 = ac,
      zero_frac = zf, sigma_d1c = first, d1c_ratio = first/sig,
      check_ac1 = is.finite(ac) && ac > -.5, check_repeated = zf >= .3,
      check_ratio = is.finite(first/sig) && (first/sig < .8 || first/sig > 1.2))
  })
  groups <- do.call(rbind, groups); rownames(groups) <- NULL
  if (any(!is.finite(groups$sigma) | groups$sigma <= 0)) stop("Unresolved group precision")
  j <- match(key, with(groups,paste(instrument,field_day,dt_s,gas,sep="|")))
  closures$sigma_group <- groups$sigma[j]
  if (warn && any(groups$check_ac1 | groups$check_repeated | groups$check_ratio))
    warning("Review empirical precision diagnostics for ",instrument,"; see precision_groups CSV (no automatic exclusions).")
  list(closures=closures,groups=groups)
}

apply_precision <- function(traces, instrument, out_dir = NULL) {
  p <- estimate_precision(traces,instrument)
  for (gas in c("CH4","CO2")) {
    r <- p$closures[p$closures$gas==gas,]
    traces[[paste0(gas,"_prec")]] <- r$sigma_group[match(traces$UniqueID,r$UniqueID)]
  }
  stopifnot(all(is.finite(traces$CH4_prec)),all(is.finite(traces$CO2_prec)))
  if (!is.null(out_dir)) {
    write.csv(p$closures,file.path(out_dir,paste0("precision_closures_",instrument,".csv")),row.names=FALSE)
    write.csv(p$groups,file.path(out_dir,paste0("precision_groups_",instrument,".csv")),row.names=FALSE)
  }
  traces
}

# goFlux (fork v0.5.0.9003) with the campaign sigma as precision and conf = 0.95: the MDF is
# 1.96 sigma / t x flux term, with t the window span plus one logging interval (closure.time()),
# and it bounds the Hutchinson-Mosier curvature (Hueppi et al. 2018). best.flux() then flags a fit
# as noisy when its residual scatter exceeds sigma.
fit_gas <- function(dataframe, gastype) goFlux::goFlux(dataframe, gastype, conf = 0.95)
