#!/usr/bin/env Rscript
# =============================================================================
# figure_wu_below2m.R
# Wu et al. (2024) reanalysis: studies NOT in the original figure (i.e. studies
# with NO measurement at >= 2 m).  Faceted by Upland / Wetland with free
# x & y scales.  Also prints vertical-trend statistics (linear & exponential).
# =============================================================================

library(readxl)
library(ggplot2)
library(dplyr)
library(stringr)
library(broom)

# ── Paths ────────────────────────────────────────────────────────────────────
base_dir <- "/Users/jongewirtzman/My Drive/Research/whole_tree_flux"
wu_file  <- file.path(base_dir, "1-s2.0-S0168192324000911-mmc2.xlsx")

# ── Load & clean ─────────────────────────────────────────────────────────────

upland  <- read_xlsx(wu_file, sheet = "upland-stem")  %>% mutate(ecosystem = "Upland")
wetland <- read_xlsx(wu_file, sheet = "wetland-stem") %>% mutate(ecosystem = "Wetland")

colnames(upland)[1]  <- "Reference";  colnames(upland)[9]  <- "Height_raw"
colnames(upland)[10] <- "CH4_flux"
colnames(wetland)[1] <- "Reference";  colnames(wetland)[9] <- "Height_raw"
colnames(wetland)[10] <- "CH4_flux"

dat <- bind_rows(
  upland  %>% select(Reference, Height_raw, CH4_flux, ecosystem),
  wetland %>% select(Reference, Height_raw, CH4_flux, ecosystem)
)

# Parse height (ranges → midpoint)
parse_height <- function(x) {
  if (is.na(x) || x == "None") return(NA_real_)
  s <- str_trim(as.character(x))
  s <- str_replace(s, "^0\\.\\.6", "0.6")
  m <- str_match(s, "^(-?[0-9.]+)\\s*-\\s*([0-9.]+)$")
  if (!is.na(m[1, 1])) return(mean(c(as.numeric(m[1, 2]), as.numeric(m[1, 3]))))
  suppressWarnings(as.numeric(s))
}
dat <- dat %>%
  mutate(Height = sapply(Height_raw, parse_height)) %>%
  filter(!is.na(Height), !is.na(CH4_flux))

# ── Filter: studies with NO measurement >= 2 m (complement of fig 3) ────────
refs_2m   <- dat %>% filter(Height >= 2) %>% distinct(Reference) %>% pull(Reference)
dat_below <- dat %>% filter(!(Reference %in% refs_2m))

cat("=== Data summary ===\n")
cat("Total studies with all measurements < 2 m:", length(unique(dat_below$Reference)), "\n")
cat("Observations:", nrow(dat_below), "\n")
cat("Upland studies:", length(unique(dat_below$Reference[dat_below$ecosystem == "Upland"])), "\n")
cat("Wetland studies:", length(unique(dat_below$Reference[dat_below$ecosystem == "Wetland"])), "\n\n")

# Clean reference labels
dat_below <- dat_below %>%
  mutate(
    ref_short = str_replace(Reference, "_[A-Za-z-]+$", "") %>%
      str_replace_all("_", " "),
    facet_label = paste0(ref_short, " (", ecosystem, ")")
  )

dat_below$facet_label <- factor(dat_below$facet_label,
                                levels = sort(unique(dat_below$facet_label)))

# ── Slope-fitting criteria ────────────────────────────────────────────────────
# A slope is only fit if ALL of:
#   (1) height range >= 0.5 m
#   (2) >= 3 unique heights
#   (3) no single consecutive gap > 80% of total range (heights must be spread)
can_fit <- function(heights) {
  h <- sort(unique(heights))
  if (length(h) < 3) return(FALSE)
  rng <- max(h) - min(h)
  if (rng < 0.5) return(FALSE)
  gaps <- diff(h)
  if (max(gaps) / rng > 0.8) return(FALSE)
  TRUE
}

# ── Normalize within study (relative flux) ───────────────────────────────────
# Divide each observation by its study's mean absolute flux so that
# vertical trends are comparable across studies with different magnitudes.

dat_below <- dat_below %>%
  group_by(Reference) %>%
  mutate(mean_abs_flux = mean(abs(CH4_flux), na.rm = TRUE),
         CH4_rel = CH4_flux / mean_abs_flux) %>%
  ungroup()

# ── Vertical-trend statistics (RELATIVE) ─────────────────────────────────────

cat("==================================================================\n")
cat("VERTICAL TREND STATISTICS — RELATIVE (flux / study-mean |flux|)\n")
cat("==================================================================\n\n")

for (eco in c("Upland", "Wetland")) {
  sub <- dat_below %>% filter(ecosystem == eco)
  cat("--- ", eco, " ---\n", sep = "")
  cat("  N studies:", length(unique(sub$Reference)),
      " | N obs:", nrow(sub), "\n")

  if (!can_fit(sub$Height)) {
    cat("  Insufficient data (fails fit criteria).\n\n")
    next
  }

  # Linear: CH4_rel ~ Height (pooled, normalized)
  lm_fit <- lm(CH4_rel ~ Height, data = sub)
  lm_sum <- summary(lm_fit)
  cat("\n  POOLED LINEAR (relative_flux ~ Height):\n")
  cat("    Intercept:", round(coef(lm_fit)[1], 4), "\n")
  cat("    Slope    :", round(coef(lm_fit)[2], 4),
      " (relative units per m height)\n")
  cat("    R²       :", round(lm_sum$r.squared, 4), "\n")
  cat("    p(slope) :", format.pval(lm_sum$coefficients[2, 4], digits = 3), "\n")

  # Exponential: log(CH4_rel) ~ Height  (only for positive relative fluxes)
  sub_pos <- sub %>% filter(CH4_rel > 0)
  if (nrow(sub_pos) >= 3) {
    exp_fit <- lm(log(CH4_rel) ~ Height, data = sub_pos)
    exp_sum <- summary(exp_fit)
    cat("\n  POOLED EXPONENTIAL (log(relative_flux) ~ Height, positive only, N =",
        nrow(sub_pos), "):\n")
    cat("    log-Intercept:", round(coef(exp_fit)[1], 4), "\n")
    cat("    Decay rate   :", round(coef(exp_fit)[2], 4), " per m height\n")
    cat("    R²           :", round(exp_sum$r.squared, 4), "\n")
    cat("    p(slope)     :", format.pval(exp_sum$coefficients[2, 4], digits = 3), "\n")
    if (coef(exp_fit)[2] < 0)
      cat("    Half-life    :", round(-log(2) / coef(exp_fit)[2], 2), " m\n")
  } else {
    cat("\n  EXPONENTIAL: Too few positive fluxes (N =", nrow(sub_pos),
        ") for log-linear fit.\n")
  }

  # AIC comparison (on shared positive-flux subset for fair comparison)
  if (nrow(sub_pos) >= 3) {
    lm_pos  <- lm(CH4_rel ~ Height, data = sub_pos)
    aic_lin <- AIC(lm_pos)
    aic_exp <- AIC(exp_fit) + 2 * sum(log(sub_pos$CH4_rel))
    cat("\n  AIC comparison (positive relative fluxes only):\n")
    cat("    Linear AIC     :", round(aic_lin, 1), "\n")
    cat("    Exponential AIC:", round(aic_exp, 1), "\n")
    cat("    Better fit     :", ifelse(aic_exp < aic_lin, "Exponential", "Linear"), "\n")
  }
  cat("\n")
}

# ── Per-study RELATIVE slopes (linear & exponential) ─────────────────────────

# --- Linear per-study ---
cat("==================================================================\n")
cat("PER-STUDY RELATIVE LINEAR SLOPES (relative_flux ~ Height)\n")
cat("==================================================================\n\n")

lin_slopes <- dat_below %>%
  group_by(Reference, ecosystem) %>%
  filter(can_fit(Height)) %>%
  do(tidy(lm(CH4_rel ~ Height, data = .))) %>%
  filter(term == "Height") %>%
  ungroup()

if (nrow(lin_slopes) > 0) {
  for (eco in c("Upland", "Wetland")) {
    eco_slopes <- lin_slopes %>% filter(ecosystem == eco)
    if (nrow(eco_slopes) == 0) {
      cat(eco, ": No studies with enough data.\n\n")
      next
    }
    cat(eco, " (", nrow(eco_slopes), " studies):\n", sep = "")
    for (i in seq_len(nrow(eco_slopes))) {
      cat("  ", eco_slopes$Reference[i], ": slope = ",
          round(eco_slopes$estimate[i], 4),
          ", p = ", format.pval(eco_slopes$p.value[i], digits = 3), "\n", sep = "")
    }
    cat("\n  Mean  :", round(mean(eco_slopes$estimate), 4), "\n")
    cat("  Median:", round(median(eco_slopes$estimate), 4), "\n")
    if (nrow(eco_slopes) >= 3) {
      tt <- t.test(eco_slopes$estimate, mu = 0)
      cat("  t-test: t =", round(tt$statistic, 3), ", df =", tt$parameter,
          ", p =", format.pval(tt$p.value, digits = 3),
          "  95% CI [", round(tt$conf.int[1], 4), ",",
          round(tt$conf.int[2], 4), "]\n")
    }
    cat("\n")
  }
  cat("All (", nrow(lin_slopes), " studies): mean =",
      round(mean(lin_slopes$estimate), 4),
      ", median =", round(median(lin_slopes$estimate), 4))
  if (nrow(lin_slopes) >= 3) {
    tt <- t.test(lin_slopes$estimate, mu = 0)
    cat(", t =", round(tt$statistic, 3),
        ", p =", format.pval(tt$p.value, digits = 3))
  }
  cat("\n\n")
} else {
  cat("No studies with enough observations for per-study regression.\n\n")
}

# --- Exponential per-study (positive fluxes only) ---
cat("==================================================================\n")
cat("PER-STUDY RELATIVE EXPONENTIAL SLOPES (log(relative_flux) ~ Height)\n")
cat("==================================================================\n\n")

exp_slopes <- dat_below %>%
  filter(CH4_rel > 0) %>%
  group_by(Reference, ecosystem) %>%
  filter(can_fit(Height)) %>%
  do(tidy(lm(log(CH4_rel) ~ Height, data = .))) %>%
  filter(term == "Height") %>%
  ungroup()

if (nrow(exp_slopes) > 0) {
  for (eco in c("Upland", "Wetland")) {
    eco_slopes <- exp_slopes %>% filter(ecosystem == eco)
    if (nrow(eco_slopes) == 0) {
      cat(eco, ": No studies with enough data.\n\n")
      next
    }
    cat(eco, " (", nrow(eco_slopes), " studies):\n", sep = "")
    for (i in seq_len(nrow(eco_slopes))) {
      hl_str <- ""
      if (eco_slopes$estimate[i] < 0)
        hl_str <- paste0(", hl = ", round(-log(2) / eco_slopes$estimate[i], 2), " m")
      cat("  ", eco_slopes$Reference[i], ": decay = ",
          round(eco_slopes$estimate[i], 4),
          ", p = ", format.pval(eco_slopes$p.value[i], digits = 3),
          hl_str, "\n", sep = "")
    }
    cat("\n  Mean decay :", round(mean(eco_slopes$estimate), 4), "\n")
    cat("  Median decay:", round(median(eco_slopes$estimate), 4), "\n")
    if (nrow(eco_slopes) >= 3) {
      tt <- t.test(eco_slopes$estimate, mu = 0)
      cat("  t-test: t =", round(tt$statistic, 3), ", df =", tt$parameter,
          ", p =", format.pval(tt$p.value, digits = 3),
          "  95% CI [", round(tt$conf.int[1], 4), ",",
          round(tt$conf.int[2], 4), "]\n")
    }
    cat("\n")
  }
  cat("All (", nrow(exp_slopes), " studies): mean =",
      round(mean(exp_slopes$estimate), 4),
      ", median =", round(median(exp_slopes$estimate), 4))
  if (nrow(exp_slopes) >= 3) {
    tt <- t.test(exp_slopes$estimate, mu = 0)
    cat(", t =", round(tt$statistic, 3),
        ", p =", format.pval(tt$p.value, digits = 3))
  }
  cat("\n\n")
} else {
  cat("No studies with enough positive observations for exponential fit.\n\n")
}

# ── Plot ─────────────────────────────────────────────────────────────────────

eco_colors <- c(Upland = "#d7191c", Wetland = "#2c7bb6")

n_upland  <- length(unique(dat_below$facet_label[dat_below$ecosystem == "Upland"]))
n_wetland <- length(unique(dat_below$facet_label[dat_below$ecosystem == "Wetland"]))
ncol_facets <- 6

plot_theme <- theme_classic(base_size = 11) +
  theme(
    strip.text       = element_text(size = 7.5, face = "bold"),
    strip.background = element_blank(),
    legend.position  = "none",
    axis.line        = element_line(linewidth = 0.3),
    axis.ticks       = element_line(linewidth = 0.3),
    panel.spacing    = unit(0.6, "lines")
  )

make_panel <- function(eco_type) {
  sub <- dat_below %>% filter(ecosystem == eco_type)
  ggplot(sub, aes(x = Height, y = CH4_flux)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black",
               linewidth = 0.4) +
    geom_point(size = 1.8, alpha = 0.7, shape = 16,
               color = eco_colors[eco_type]) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.6,
                alpha = 0.10, color = eco_colors[eco_type],
                fill = eco_colors[eco_type], na.rm = TRUE) +
    facet_wrap(~ facet_label, scales = "free", ncol = ncol_facets) +
    labs(y = expression(Stem~CH[4]~flux~(nmol~m^{-2}~s^{-1})),
         x = "Height (m)",
         title = eco_type) +
    coord_flip() +
    plot_theme +
    theme(plot.title = element_text(size = 13, face = "bold",
                                    color = eco_colors[eco_type]))
}

p_upland  <- make_panel("Upland")
p_wetland <- make_panel("Wetland")

# Combine with patchwork — top and bottom
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  p_combined <- p_upland / p_wetland +
    plot_layout(heights = c(n_upland, n_wetland))
} else {
  cat("Install 'patchwork' for combined plot. Saving panels separately.\n")
  p_combined <- p_upland
}

# ── Save ─────────────────────────────────────────────────────────────────────

out_dir <- file.path(base_dir, "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Compute dimensions: top/bottom, ~square facets
rows_up  <- ceiling(n_upland / ncol_facets)
rows_wet <- ceiling(n_wetland / ncol_facets)
total_h  <- (rows_up + rows_wet) * 2.2 + 1.5
total_w  <- ncol_facets * 2.2

ggsave(file.path(out_dir, "figure_wu_below2m.pdf"), p_combined,
       width = total_w, height = total_h, units = "in", limitsize = FALSE)
ggsave(file.path(out_dir, "figure_wu_below2m.png"), p_combined,
       width = total_w, height = total_h, units = "in", dpi = 300, limitsize = FALSE)

cat("\nSaved figure_wu_below2m.pdf/.png\n")
cat("Upland facets:", n_upland, "| Wetland facets:", n_wetland, "\n")
