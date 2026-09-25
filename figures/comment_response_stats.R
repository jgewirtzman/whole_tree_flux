#!/usr/bin/env Rscript
# =============================================================================
# comment_response_stats.R
# Compute all statistics needed for coauthor comment responses:
#   - Bootstrapped mean + 95% CI for all flux categories
#   - Bootstrapped median + 95% CI
#   - SD, SE, n
#   - Wilcoxon test: stem <2m vs >=2m
#   - MDF summary by component
#   - Key point consistency check
# =============================================================================

set.seed(42)
n_boot <- 10000

# ─── Load data ───────────────────────────────────────────────────────────────
# Harvard Forest
hf <- read.csv(file.path("data processing", "goFlux_reprocessing",
                          "results", "canopy_flux_goFlux_compiled_with_mdf.csv"),
               stringsAsFactors = FALSE)
hf <- hf[!(hf$Species == "bg" & hf$Tree_Tag == 3), ]
hf$Component <- ifelse(hf$Type == "leaf (shaded)", "leaf", hf$Type)

# Yale Myers Forest
ymf <- read.csv(file.path("data processing", "goFlux_reprocessing",
                           "ymf_black_oak", "results",
                           "ymf_black_oak_flux_compiled_with_mdf.csv"),
                stringsAsFactors = FALSE)
# Parse numeric portion of height: "4 (restarted)" -> 4, "Seam" -> NA (dropped)
# Matches figure1_composite.R approach (str_extract of leading digits)
ymf$Height_m <- suppressWarnings(as.numeric(
  stringr::str_extract(as.character(ymf$Height_m), "^[\\d.]+")))
ymf$Component <- "stem"  # all YMF measurements are stem
ymf$Tree_Tag <- "YMF_1"

# Combine using shared columns
shared_cols <- intersect(names(hf), names(ymf))
dat <- rbind(hf[, shared_cols], ymf[, shared_cols])
dat <- dat[!is.na(dat$CH4_best.flux) & !is.na(dat$Height_m), ]

cat("Total measurements after filtering:", nrow(dat), "\n\n")

# ─── Split into categories ───────────────────────────────────────────────────
cats <- list(
  "Stem < 2 m"  = dat[dat$Component == "stem" & dat$Height_m < 2, ],
  "Stem >= 2 m" = dat[dat$Component == "stem" & dat$Height_m >= 2, ],
  "Branch"       = dat[dat$Component == "branch", ],
  "Leaf"         = dat[dat$Component == "leaf", ],
  "All < 2 m"   = dat[dat$Height_m < 2, ],
  "All >= 2 m"  = dat[dat$Height_m >= 2, ]
)

# ─── Bootstrap function ──────────────────────────────────────────────────────
boot_ci <- function(x, stat_fn, n = n_boot) {
  boots <- replicate(n, stat_fn(sample(x, replace = TRUE)))
  c(estimate = stat_fn(x),
    ci_lo = quantile(boots, 0.025, names = FALSE),
    ci_hi = quantile(boots, 0.975, names = FALSE))
}

# ═══════════════════════════════════════════════════════════════════════════════
# 1. SUMMARY STATISTICS TABLE
# ═══════════════════════════════════════════════════════════════════════════════
cat("════════════════════════════════════════════════════════════════════\n")
cat("BOOTSTRAPPED SUMMARY STATISTICS (CH₄ flux, nmol m⁻² s⁻¹)\n")
cat("════════════════════════════════════════════════════════════════════\n\n")

for (nm in names(cats)) {
  x <- cats[[nm]]$CH4_best.flux
  n_obs <- length(x)

  bmean  <- boot_ci(x, mean)
  bmed   <- boot_ci(x, median)

  cat(sprintf("--- %s (n = %d) ---\n", nm, n_obs))
  cat(sprintf("  Mean:   %.4f  [95%% CI: %.4f, %.4f]\n", bmean[1], bmean[2], bmean[3]))
  cat(sprintf("  Median: %.4f  [95%% CI: %.4f, %.4f]\n", bmed[1], bmed[2], bmed[3]))
  cat(sprintf("  SD:     %.4f\n", sd(x)))
  cat(sprintf("  SE:     %.4f\n", sd(x) / sqrt(n_obs)))
  cat(sprintf("  Range:  [%.4f, %.4f]\n", min(x), max(x)))
  cat(sprintf("  Negative: %d (%.1f%%)\n\n", sum(x < 0), 100 * sum(x < 0) / n_obs))
}

# ═══════════════════════════════════════════════════════════════════════════════
# 2. WILCOXON TEST: stem <2m vs >=2m
# ═══════════════════════════════════════════════════════════════════════════════
cat("════════════════════════════════════════════════════════════════════\n")
cat("WILCOXON RANK-SUM TEST: Stem <2m vs >=2m\n")
cat("════════════════════════════════════════════════════════════════════\n\n")

stem_lt2 <- cats[["Stem < 2 m"]]$CH4_best.flux
stem_ge2 <- cats[["Stem >= 2 m"]]$CH4_best.flux

wt <- wilcox.test(stem_lt2, stem_ge2, alternative = "greater")
cat(sprintf("  W = %.0f, p = %.2e\n", wt$statistic, wt$p.value))
cat(sprintf("  Stem <2m mean: %.4f (n=%d), Stem >=2m mean: %.4f (n=%d)\n",
            mean(stem_lt2), length(stem_lt2), mean(stem_ge2), length(stem_ge2)))
cat(sprintf("  Ratio: %.1fx\n", mean(stem_lt2) / mean(stem_ge2)))

# Effect size: rank-biserial correlation r = 1 - 2U/(n1*n2)
n1 <- length(stem_lt2); n2 <- length(stem_ge2)
r_effect <- 1 - 2 * wt$statistic / (n1 * n2)
# Actually for greater alternative, W = sum of ranks in first group minus n1*(n1+1)/2
# Use the simpler: r = Z / sqrt(N)
# Let's use the direct formula
U <- wt$statistic  # This is actually the W statistic
r_rb <- 2 * U / (n1 * n2) - 1
cat(sprintf("  Rank-biserial r = %.3f\n\n", r_rb))

# Also run a t-test for comparison
tt <- t.test(stem_lt2, stem_ge2, alternative = "greater")
cat(sprintf("  Welch's t-test: t = %.2f, df = %.1f, p = %.2e\n\n",
            tt$statistic, tt$parameter, tt$p.value))

# ═══════════════════════════════════════════════════════════════════════════════
# 3. MDF SUMMARY BY COMPONENT
# ═══════════════════════════════════════════════════════════════════════════════
cat("════════════════════════════════════════════════════════════════════\n")
cat("MDF SUMMARY (Empirical 95% threshold)\n")
cat("════════════════════════════════════════════════════════════════════\n\n")

for (nm in c("Stem < 2 m", "Stem >= 2 m", "Branch", "Leaf")) {
  d <- cats[[nm]]
  if ("CH4_below_MDF_wass95" %in% names(d)) {
    n_below <- sum(d$CH4_below_MDF_wass95, na.rm = TRUE)
    n_total <- sum(!is.na(d$CH4_below_MDF_wass95))
    pct <- 100 * n_below / n_total

    # Median MDF value
    if ("CH4_MDF_wass95" %in% names(d)) {
      med_mdf <- median(d$CH4_MDF_wass95, na.rm = TRUE)
      cat(sprintf("  %s: %d/%d (%.0f%%) below MDF, median MDF = %.4f nmol m⁻² s⁻¹\n",
                  nm, n_below, n_total, pct, med_mdf))
    } else {
      cat(sprintf("  %s: %d/%d (%.0f%%) below MDF\n", nm, n_below, n_total, pct))
    }
  }
}

# Overall MDF
if ("CH4_MDF_wass95" %in% names(dat)) {
  cat(sprintf("\n  Overall median CH₄ MDF (Empirical 95%%): %.4f nmol m⁻² s⁻¹\n",
              median(dat$CH4_MDF_wass95, na.rm = TRUE)))
  cat(sprintf("  Overall mean CH₄ MDF (Empirical 95%%): %.4f nmol m⁻² s⁻¹\n",
              mean(dat$CH4_MDF_wass95, na.rm = TRUE)))
}

cat("\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 4. KEY POINT CONSISTENCY CHECK (#51)
# ═══════════════════════════════════════════════════════════════════════════════
cat("════════════════════════════════════════════════════════════════════\n")
cat("KEY POINT CHECK: % negative above 2m\n")
cat("════════════════════════════════════════════════════════════════════\n\n")

stem_ge2_all <- dat[dat$Component == "stem" & dat$Height_m >= 2, ]
cat(sprintf("  Field data (stem >=2m): %d negative of %d = %.1f%%\n",
            sum(stem_ge2_all$CH4_best.flux < 0), nrow(stem_ge2_all),
            100 * sum(stem_ge2_all$CH4_best.flux < 0) / nrow(stem_ge2_all)))
cat("  Wu et al. (>=2m): 9 negative of 127 = 7.1%\n")
cat("  Key points say '6-7%' -- this covers both (6.0% field + 7.1% lit)\n\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 5. ADDITIONAL STATS FOR TEXT
# ═══════════════════════════════════════════════════════════════════════════════
cat("════════════════════════════════════════════════════════════════════\n")
cat("ADDITIONAL STATISTICS FOR MANUSCRIPT TEXT\n")
cat("════════════════════════════════════════════════════════════════════\n\n")

# Mean ratio with bootstrapped CI
ratio_boot <- replicate(n_boot, {
  s1 <- sample(stem_lt2, replace = TRUE)
  s2 <- sample(stem_ge2, replace = TRUE)
  mean(s1) / mean(s2)
})
cat(sprintf("  Mean flux ratio (stem <2m / >=2m): %.1f [95%% CI: %.1f, %.1f]\n",
            mean(stem_lt2) / mean(stem_ge2),
            quantile(ratio_boot, 0.025),
            quantile(ratio_boot, 0.975)))

# All components combined
all_lt2 <- dat[dat$Height_m < 2, ]$CH4_best.flux
all_ge2 <- dat[dat$Height_m >= 2, ]$CH4_best.flux
bm_lt2 <- boot_ci(all_lt2, mean)
bm_ge2 <- boot_ci(all_ge2, mean)
cat(sprintf("  All <2m: mean = %.4f [%.4f, %.4f] (n=%d)\n",
            bm_lt2[1], bm_lt2[2], bm_lt2[3], length(all_lt2)))
cat(sprintf("  All >=2m: mean = %.4f [%.4f, %.4f] (n=%d)\n",
            bm_ge2[1], bm_ge2[2], bm_ge2[3], length(all_ge2)))

cat("\nDone.\n")
