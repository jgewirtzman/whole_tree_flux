#!/usr/bin/env Rscript
# =============================================================================
# stats_wu_above2m_slopes.R
# For the >=2m studies (those in the original Wu reanalysis figure):
# (1) Per-study relative slopes (linear & exponential) in below-2m
#     and full-range windows
# (2) Faceted plot: points colored <2m vs >=2m, with two fit lines:
#     - below-2m fit (extended across full range)
#     - full-range fit
# =============================================================================

library(readxl)
library(ggplot2)
library(dplyr)
library(stringr)
library(broom)
library(patchwork)

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

# ── Keep only studies WITH measurements >= 2 m ──────────────────────────────
refs_2m <- dat %>% filter(Height >= 2) %>% distinct(Reference) %>% pull(Reference)
dat_2m  <- dat %>% filter(Reference %in% refs_2m)

# Normalize within study
dat_2m <- dat_2m %>%
  group_by(Reference) %>%
  mutate(mean_abs_flux = mean(abs(CH4_flux), na.rm = TRUE),
         CH4_rel = CH4_flux / mean_abs_flux) %>%
  ungroup()

# Clean labels
dat_2m <- dat_2m %>%
  mutate(
    ref_short = str_replace(Reference, "_[A-Za-z-]+$", "") %>%
      str_replace_all("_", " "),
    facet_label = paste0(ref_short, " (", ecosystem, ")"),
    height_zone = ifelse(Height < 2, "Below 2 m", "At or above 2 m")
  )
dat_2m$facet_label <- factor(dat_2m$facet_label,
                             levels = sort(unique(dat_2m$facet_label)))

cat("=== >=2m studies ===\n")
cat("N studies:", length(refs_2m), "\n")
for (ref in sort(refs_2m)) {
  sub <- dat_2m %>% filter(Reference == ref)
  cat("  ", ref, " (", unique(sub$ecosystem), "): ",
      nrow(sub), " obs, height ",
      round(min(sub$Height), 2), "-", round(max(sub$Height), 2), " m\n", sep = "")
}
cat("\n")

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

# ── Helper: fit linear & exponential slopes per study ────────────────────────
fit_slopes <- function(data, label) {

  # --- Linear per-study ---
  lin_slopes <- data %>%
    group_by(Reference, ecosystem) %>%
    filter(can_fit(Height)) %>%
    do(tidy(lm(CH4_rel ~ Height, data = .))) %>%
    filter(term == "Height") %>%
    ungroup()

  # --- Exponential per-study (positive fluxes only) ---
  exp_slopes <- data %>%
    filter(CH4_rel > 0) %>%
    group_by(Reference, ecosystem) %>%
    filter(can_fit(Height)) %>%
    do(tidy(lm(log(CH4_rel) ~ Height, data = .))) %>%
    filter(term == "Height") %>%
    ungroup()

  cat("==================================================================\n")
  cat(label, "\n")
  cat("==================================================================\n\n")

  # Print linear
  cat("── LINEAR (relative_flux ~ Height) ──\n\n")
  if (nrow(lin_slopes) == 0) {
    cat("No studies with enough data.\n\n")
  } else {
    for (eco in c("Upland", "Wetland")) {
      eco_sl <- lin_slopes %>% filter(ecosystem == eco)
      if (nrow(eco_sl) == 0) {
        cat(eco, ": No studies with enough data.\n\n")
        next
      }
      cat(eco, " (", nrow(eco_sl), " studies):\n", sep = "")
      for (i in seq_len(nrow(eco_sl))) {
        cat("  ", eco_sl$Reference[i], ": slope = ",
            round(eco_sl$estimate[i], 4),
            ", p = ", format.pval(eco_sl$p.value[i], digits = 3), "\n", sep = "")
      }
      cat("\n  Mean  :", round(mean(eco_sl$estimate), 4), "\n")
      cat("  Median:", round(median(eco_sl$estimate), 4), "\n")
      if (nrow(eco_sl) >= 3) {
        tt <- t.test(eco_sl$estimate, mu = 0)
        cat("  t-test: t =", round(tt$statistic, 3), ", df =", tt$parameter,
            ", p =", format.pval(tt$p.value, digits = 3),
            "  95% CI [", round(tt$conf.int[1], 4), ",",
            round(tt$conf.int[2], 4), "]\n")
      } else {
        cat("  (too few studies for t-test)\n")
      }
      cat("\n")
    }
    # Combined
    cat("All (", nrow(lin_slopes), " studies): mean =",
        round(mean(lin_slopes$estimate), 4),
        ", median =", round(median(lin_slopes$estimate), 4))
    if (nrow(lin_slopes) >= 3) {
      tt <- t.test(lin_slopes$estimate, mu = 0)
      cat(", t =", round(tt$statistic, 3),
          ", p =", format.pval(tt$p.value, digits = 3))
    }
    cat("\n\n")
  }

  # Print exponential
  cat("── EXPONENTIAL (log(relative_flux) ~ Height, positive only) ──\n\n")
  if (nrow(exp_slopes) == 0) {
    cat("No studies with enough data.\n\n")
  } else {
    for (eco in c("Upland", "Wetland")) {
      eco_sl <- exp_slopes %>% filter(ecosystem == eco)
      if (nrow(eco_sl) == 0) {
        cat(eco, ": No studies with enough data.\n\n")
        next
      }
      cat(eco, " (", nrow(eco_sl), " studies):\n", sep = "")
      for (i in seq_len(nrow(eco_sl))) {
        hl_str <- ""
        if (eco_sl$estimate[i] < 0)
          hl_str <- paste0(", hl = ", round(-log(2) / eco_sl$estimate[i], 2), " m")
        cat("  ", eco_sl$Reference[i], ": decay = ",
            round(eco_sl$estimate[i], 4),
            ", p = ", format.pval(eco_sl$p.value[i], digits = 3),
            hl_str, "\n", sep = "")
      }
      cat("\n  Mean decay :", round(mean(eco_sl$estimate), 4), "\n")
      cat("  Median decay:", round(median(eco_sl$estimate), 4), "\n")
      if (nrow(eco_sl) >= 3) {
        tt <- t.test(eco_sl$estimate, mu = 0)
        cat("  t-test: t =", round(tt$statistic, 3), ", df =", tt$parameter,
            ", p =", format.pval(tt$p.value, digits = 3),
            "  95% CI [", round(tt$conf.int[1], 4), ",",
            round(tt$conf.int[2], 4), "]\n")
      } else {
        cat("  (too few studies for t-test)\n")
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
  }

  invisible(list(linear = lin_slopes, exponential = exp_slopes))
}

# ── (A) Below-2m portion only ────────────────────────────────────────────────
dat_low <- dat_2m %>% filter(Height < 2)
fit_slopes(dat_low, "PER-STUDY RELATIVE SLOPES: BELOW 2 m ONLY (>=2m studies)")

# ── (B) Full height range ───────────────────────────────────────────────────
fit_slopes(dat_2m, "PER-STUDY RELATIVE SLOPES: FULL HEIGHT RANGE (>=2m studies)")

# ── Pooled normalized regressions (linear & exponential) ────────────────────

cat("\n==================================================================\n")
cat("POOLED NORMALIZED REGRESSIONS (>=2m studies)\n")
cat("==================================================================\n\n")

for (window_name in c("Below 2m", "Full range")) {
  sub <- switch(window_name,
    "Below 2m"   = dat_low,
    "Full range" = dat_2m
  )

  cat("--- ", window_name, " ---\n", sep = "")

  for (eco in c("Upland", "Wetland")) {
    eco_sub <- sub %>% filter(ecosystem == eco)
    if (!can_fit(eco_sub$Height)) {
      cat("  ", eco, ": insufficient data (fails fit criteria)\n")
      next
    }

    cat("  ", eco, " (N = ", nrow(eco_sub), "):\n", sep = "")

    lm_fit <- lm(CH4_rel ~ Height, data = eco_sub)
    lm_sum <- summary(lm_fit)
    cat("    Linear: slope = ", round(coef(lm_fit)[2], 4),
        ", R² = ", round(lm_sum$r.squared, 4),
        ", p = ", format.pval(lm_sum$coefficients[2, 4], digits = 3), "\n")

    eco_pos <- eco_sub %>% filter(CH4_rel > 0)
    if (nrow(eco_pos) >= 3) {
      exp_fit <- lm(log(CH4_rel) ~ Height, data = eco_pos)
      exp_sum <- summary(exp_fit)
      cat("    Exponential (N+ = ", nrow(eco_pos), "): decay = ",
          round(coef(exp_fit)[2], 4),
          ", R² = ", round(exp_sum$r.squared, 4),
          ", p = ", format.pval(exp_sum$coefficients[2, 4], digits = 3), sep = "")
      if (coef(exp_fit)[2] < 0)
        cat(", half-life = ", round(-log(2) / coef(exp_fit)[2], 2), " m", sep = "")
      cat("\n")

      lm_pos <- lm(CH4_rel ~ Height, data = eco_pos)
      aic_lin <- AIC(lm_pos)
      aic_exp <- AIC(exp_fit) + 2 * sum(log(eco_pos$CH4_rel))
      cat("    AIC: linear = ", round(aic_lin, 1),
          ", exponential = ", round(aic_exp, 1),
          ", better = ", ifelse(aic_exp < aic_lin, "Exponential", "Linear"), "\n")
    }
  }
  cat("\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# PLOT: Faceted by study, points colored <2m vs >=2m,
# four fit lines: linear & exponential × below-2m (extended) & full-range
# ══════════════════════════════════════════════════════════════════════════════

# Natural, earthy color palette
zone_colors <- c("Below 2 m"       = "#c2885c",   # warm sand
                 "At or above 2 m" = "#5b7e5f")   # forest green

# Compute per-study regression lines (linear + exponential, only when criteria met)
line_data <- dat_2m %>%
  group_by(Reference, facet_label, ecosystem) %>%
  do({
    d <- .
    h_range <- seq(min(d$Height), max(d$Height), length.out = 200)
    result <- data.frame()

    # --- Full-range fits ---
    if (can_fit(d$Height)) {
      # Linear
      lm_full <- lm(CH4_flux ~ Height, data = d)
      result <- bind_rows(result, data.frame(
        Height = h_range,
        CH4_flux = predict(lm_full, newdata = data.frame(Height = h_range)),
        fit_type = "Full range (linear)"
      ))
      # Exponential (positive fluxes only)
      d_pos <- d %>% filter(CH4_flux > 0)
      if (nrow(d_pos) >= 3 && length(unique(d_pos$Height)) >= 2) {
        exp_full <- lm(log(CH4_flux) ~ Height, data = d_pos)
        result <- bind_rows(result, data.frame(
          Height = h_range,
          CH4_flux = exp(predict(exp_full, newdata = data.frame(Height = h_range))),
          fit_type = "Full range (exponential)"
        ))
      }
    }

    # --- Below-2m fits, extended across full range ---
    d_low <- d %>% filter(Height < 2)
    if (can_fit(d_low$Height)) {
      # Linear
      lm_low <- lm(CH4_flux ~ Height, data = d_low)
      result <- bind_rows(result, data.frame(
        Height = h_range,
        CH4_flux = predict(lm_low, newdata = data.frame(Height = h_range)),
        fit_type = "Below 2 m (linear, extended)"
      ))
      # Exponential (positive fluxes only)
      d_low_pos <- d_low %>% filter(CH4_flux > 0)
      if (nrow(d_low_pos) >= 3 && length(unique(d_low_pos$Height)) >= 2) {
        exp_low <- lm(log(CH4_flux) ~ Height, data = d_low_pos)
        result <- bind_rows(result, data.frame(
          Height = h_range,
          CH4_flux = exp(predict(exp_low, newdata = data.frame(Height = h_range))),
          fit_type = "Below 2 m (exponential, extended)"
        ))
      }
    }

    if (nrow(result) == 0) {
      data.frame(Height = numeric(0), CH4_flux = numeric(0),
                 fit_type = character(0))
    } else {
      result
    }
  }) %>%
  ungroup()

line_data$facet_label <- factor(line_data$facet_label,
                                levels = levels(dat_2m$facet_label))

# Colors and linetypes for fits
fit_colors <- c(
  "Full range (linear)"              = "#2d4a7a",  # deep slate blue
  "Full range (exponential)"         = "#7a2d4a",  # deep berry
  "Below 2 m (linear, extended)"     = "#c2885c",  # warm sand (matches points)
  "Below 2 m (exponential, extended)" = "#8b6534"   # dark ochre
)
fit_ltypes <- c(
  "Full range (linear)"              = "solid",
  "Full range (exponential)"         = "solid",
  "Below 2 m (linear, extended)"     = "dashed",
  "Below 2 m (exponential, extended)" = "dashed"
)

all_colors <- c(zone_colors, fit_colors)
all_breaks <- c("Below 2 m", "At or above 2 m",
                "Full range (linear)", "Full range (exponential)",
                "Below 2 m (linear, extended)", "Below 2 m (exponential, extended)")

plot_theme <- theme_classic(base_size = 11) +
  theme(
    strip.text       = element_text(size = 8, face = "bold"),
    strip.background = element_blank(),
    legend.position  = "bottom",
    legend.box       = "vertical",
    legend.margin    = margin(0, 0, 0, 0),
    legend.text      = element_text(size = 8),
    legend.title     = element_blank(),
    axis.line        = element_line(linewidth = 0.3),
    axis.ticks       = element_line(linewidth = 0.3),
    panel.spacing    = unit(0.8, "lines")
  )

p <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey60",
             linewidth = 0.3) +
  geom_vline(xintercept = 2, linetype = "dotted", color = "grey60",
             linewidth = 0.3) +
  geom_point(data = dat_2m,
             aes(x = Height, y = CH4_flux, color = height_zone),
             size = 2, alpha = 0.65, shape = 16) +
  geom_line(data = line_data,
            aes(x = Height, y = CH4_flux, linetype = fit_type, color = fit_type),
            linewidth = 0.7) +
  scale_color_manual(
    values = all_colors,
    breaks = all_breaks,
    guide = guide_legend(
      nrow = 3,
      override.aes = list(
        shape    = c(16, 16, NA, NA, NA, NA),
        linetype = c("blank", "blank", "solid", "solid", "dashed", "dashed"),
        linewidth = c(NA, NA, 0.7, 0.7, 0.7, 0.7)
      ))
  ) +
  scale_linetype_manual(
    values = fit_ltypes,
    guide = "none"
  ) +
  facet_wrap(~ facet_label, scales = "free", ncol = 4) +
  labs(y = expression(Stem~CH[4]~flux~(nmol~m^{-2}~s^{-1})),
       x = "Height (m)") +
  coord_flip() +
  plot_theme

# ── Save ─────────────────────────────────────────────────────────────────────
out_dir <- file.path(base_dir, "figures")
n_facets <- length(unique(dat_2m$facet_label))
plot_h <- ceiling(n_facets / 4) * 2.8 + 1.8

ggsave(file.path(out_dir, "figure_wu_above2m_slopes.pdf"), p,
       width = 9, height = plot_h, units = "in", limitsize = FALSE)
ggsave(file.path(out_dir, "figure_wu_above2m_slopes.png"), p,
       width = 9, height = plot_h, units = "in", dpi = 300, limitsize = FALSE)

cat("\nSaved figure_wu_above2m_slopes.pdf/.png\n")
