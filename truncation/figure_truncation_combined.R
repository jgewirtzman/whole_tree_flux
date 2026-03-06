#!/usr/bin/env Rscript
# =============================================================================
# figure_truncation_combined.R
#
# Combined figure: sensitivity/budget row on top, truncation/breakeven on bottom.
#
# ROW 1 (a-d):
#   (a) Per-tissue CH4 flux rates with SE
#   (b) Surface area indices
#   (c) Integrated stand-level budget (stacked bar, flux rates)
#   (d) Dot-and-whisker (% of stand budget +/- SE)
#
# ROW 2 (e-g):
#   (e) 3-D cone schematic
#   (f) Capture-fraction curves
#   (g) Break-even thresholds vs observed, with MDF dotted lines
#
# Output: figure_truncation_combined.pdf/.png
# =============================================================================

library(ggplot2)
library(dplyr)
library(patchwork)

# =============================================================================
# 1. TREE GEOMETRY
# =============================================================================

tree_height      <- 26      # m
dbh              <- 0.40    # m (40 cm)
dbh_height       <- 1.37    # m
base_radius      <- dbh / 2
height_above_dbh <- tree_height - dbh_height

# Whittaker & Woodwell surface area indices
stem_bark_index   <- 0.45   # m2 stem bark / m2 ground
branch_bark_index <- 1.70   # m2 branch bark / m2 ground
leaf_area_index   <- 4.5    # m2 leaf / m2 ground

radius_at <- function(h) {
  pmax(0, base_radius * (tree_height - h) / height_above_dbh)
}

cone_lateral_area <- function(h1, h2) {
  r1 <- radius_at(h1)
  r2 <- radius_at(h2)
  slant <- sqrt((h2 - h1)^2 + (r1 - r2)^2)
  pi * (r1 + r2) * slant
}

total_stem_area   <- cone_lateral_area(0, tree_height)
ground_area       <- total_stem_area / stem_bark_index
total_branch_area <- ground_area * branch_bark_index
total_leaf_area   <- ground_area * leaf_area_index

capture_at <- function(h) {
  stem_below <- cone_lateral_area(0, h)
  c(stem       = 100 * stem_below / total_stem_area,
    stem_br    = 100 * stem_below / (total_stem_area + total_branch_area),
    stem_br_lf = 100 * stem_below / (total_stem_area + total_branch_area +
                                      total_leaf_area))
}

# Cone-based split of stem into <2m / >=2m
total_cone_area <- cone_lateral_area(0, tree_height)
cone_below_2m   <- cone_lateral_area(0, 2)
frac_stem_lt2   <- cone_below_2m / total_cone_area

A_stem_lt2  <- stem_bark_index * frac_stem_lt2
A_stem_ge2  <- stem_bark_index * (1 - frac_stem_lt2)
A_branch    <- branch_bark_index
A_leaf      <- leaf_area_index
A_gt2_total <- A_stem_ge2 + A_branch + A_leaf
A_total     <- A_stem_lt2 + A_gt2_total

# =============================================================================
# 2. EMPIRICAL FLUX DATA
# =============================================================================

dat <- read.csv(file.path("data processing", "goFlux_reprocessing",
                           "results", "canopy_flux_goFlux_compiled_with_mdf.csv"),
                stringsAsFactors = FALSE)

dat <- dat %>%
  filter(!(Species == "bg" & Tree_Tag == 3)) %>%
  mutate(Component = ifelse(Type == "leaf (shaded)", "leaf", Type)) %>%
  filter(!is.na(CH4_best.flux), !is.na(Height_m))

# Per-tissue flux summaries
stem_lt2_dat <- dat %>% filter(Component == "stem", Height_m < 2)
stem_ge2_dat <- dat %>% filter(Component == "stem", Height_m >= 2)
branch_dat   <- dat %>% filter(Component == "branch")
leaf_dat     <- dat %>% filter(Component == "leaf")

F_stem_lt2 <- mean(stem_lt2_dat$CH4_best.flux, na.rm = TRUE)
F_stem_ge2 <- mean(stem_ge2_dat$CH4_best.flux, na.rm = TRUE)
F_branch   <- mean(branch_dat$CH4_best.flux, na.rm = TRUE)
F_leaf     <- mean(leaf_dat$CH4_best.flux, na.rm = TRUE)

# SE of the mean = SD / sqrt(n)
SE_stem_lt2 <- sd(stem_lt2_dat$CH4_best.flux, na.rm = TRUE) / sqrt(nrow(stem_lt2_dat))
SE_stem_ge2 <- sd(stem_ge2_dat$CH4_best.flux, na.rm = TRUE) / sqrt(nrow(stem_ge2_dat))
SE_branch   <- sd(branch_dat$CH4_best.flux, na.rm = TRUE)    / sqrt(nrow(branch_dat))
SE_leaf     <- sd(leaf_dat$CH4_best.flux, na.rm = TRUE)      / sqrt(nrow(leaf_dat))

# Stand-level budget
budget <- tibble(
  compartment = c("Stem < 2 m", "Stem > 2 m", "Branches", "Leaves"),
  flux_rate   = c(F_stem_lt2, F_stem_ge2, F_branch, F_leaf),
  flux_SE     = c(SE_stem_lt2, SE_stem_ge2, SE_branch, SE_leaf),
  area_index  = c(A_stem_lt2, A_stem_ge2, A_branch, A_leaf)
) %>%
  mutate(
    integrated    = flux_rate * area_index,
    integrated_SE = flux_SE * area_index,
    pct           = 100 * integrated / sum(integrated)
  )

# Total uncertainty via quadrature
total_integrated_SE <- sqrt(sum(budget$integrated_SE^2))

gt2_integrated   <- sum(budget$integrated[budget$compartment != "Stem < 2 m"])
total_integrated <- sum(budget$integrated)
F_gt2_observed   <- gt2_integrated / A_gt2_total

# Median MDF (Wassmann 90%) for reference line on breakeven panel
median_MDF <- median(dat$CH4_MDF_wass90, na.rm = TRUE)

# Break-even
baseline <- budget$integrated[budget$compartment == "Stem < 2 m"]
scenarios <- tibble(
  scenario    = c("Double", "Cancel", "Flip"),
  description = c("Net = 2x basal emission",
                   "Net = 0 (break-even)",
                   "Net = -1x basal (equal sink)"),
  multiplier  = c(2, 0, -1)
) %>%
  mutate(
    net_target     = multiplier * baseline,
    F_gt2_required = (net_target - baseline) / A_gt2_total
  )

# Print capture-fraction values for caption
pcts_2  <- capture_at(2)
pcts_10 <- capture_at(10)
cat("Caption values — % captured at 2 m:\n")
cat(sprintf("  Stem: %.1f%%  Stem+branches: %.1f%%  Stem+branches+leaves: %.1f%%\n",
            pcts_2[1], pcts_2[2], pcts_2[3]))
cat("Caption values — % captured at 10 m:\n")
cat(sprintf("  Stem: %.1f%%  Stem+branches: %.1f%%  Stem+branches+leaves: %.1f%%\n\n",
            pcts_10[1], pcts_10[2], pcts_10[3]))

# =============================================================================
# 3. SHARED AESTHETICS
# =============================================================================

compartment_colors <- c(
  "Stem < 2 m"  = "#8B4513",
  "Stem > 2 m"  = "#D4A76A",
  "Branches"    = "#4682B4",
  "Leaves"      = "#2E8B57"
)

comp_labels <- c("Stem < 2 m", "Stem > 2 m", "Branches", "Leaves")

rate_data <- tibble(
  compartment = factor(comp_labels, levels = comp_labels),
  flux_rate   = c(F_stem_lt2, F_stem_ge2, F_branch, F_leaf),
  flux_SE     = c(SE_stem_lt2, SE_stem_ge2, SE_branch, SE_leaf),
  area_index  = c(A_stem_lt2, A_stem_ge2, A_branch, A_leaf)
)

# =============================================================================
# 4. ROW 1: SENSITIVITY PANELS
# =============================================================================

# -- Flux rates with SE --------------------------------------------------------

p_rates <- ggplot(rate_data, aes(x = compartment, y = flux_rate,
                                  fill = compartment)) +
  geom_col(width = 0.6, color = "grey30", linewidth = 0.3) +
  geom_errorbar(aes(ymin = flux_rate - flux_SE, ymax = flux_rate + flux_SE),
                width = 0.2, linewidth = 0.4) +
  geom_text(aes(y = flux_rate + flux_SE, label = sprintf("%.2f", flux_rate)),
            vjust = -0.3, size = 3.5, color = "grey20") +
  scale_fill_manual(values = setNames(unname(compartment_colors), comp_labels),
                    guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    x = NULL,
    y = expression(atop("CH"[4]~"flux rate", "(nmol m"^{-2}~"s"^{-1}*")"))
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.line    = element_line(linewidth = 0.3),
    axis.ticks   = element_line(linewidth = 0.3),
    axis.text    = element_text(size = 10),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title   = element_text(size = 10),
    plot.margin  = margin(10, 10, 2, 10)
  )

# -- Area indices --------------------------------------------------------------

p_area <- ggplot(rate_data, aes(x = compartment, y = area_index,
                                 fill = compartment)) +
  geom_col(width = 0.6, color = "grey30", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", area_index)),
            vjust = -0.3, size = 3.5, color = "grey20") +
  scale_fill_manual(values = setNames(unname(compartment_colors), comp_labels),
                    guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    x = NULL,
    y = expression(atop("Area index", "(m"^{2}~"m"^{-2}*"ground)"))
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.line    = element_line(linewidth = 0.3),
    axis.ticks   = element_line(linewidth = 0.3),
    axis.text    = element_text(size = 10),
    axis.text.x  = element_text(size = 9),
    axis.title   = element_text(size = 10),
    plot.margin  = margin(2, 10, 10, 10)
  )

# -- Integrated budget (stacked bar, flux rates) ------------------------------

budget_plot <- budget
budget_plot$compartment <- factor(budget_plot$compartment,
  levels = c("Leaves", "Branches", "Stem > 2 m", "Stem < 2 m"))

p_budget <- ggplot(budget_plot, aes(x = 1, y = integrated,
                                     fill = compartment)) +
  geom_col(width = 0.5, color = "grey30", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%s\n%.2f", compartment, integrated)),
            position = position_stack(vjust = 0.5),
            size = 3.5, color = "white", fontface = "bold") +
  scale_fill_manual(values = compartment_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = NULL,
    y = expression("Integrated flux (nmol m"^{-2}*"ground s"^{-1}*")")
  ) +
  coord_cartesian(xlim = c(0.4, 1.6)) +
  theme_classic(base_size = 11) +
  theme(
    axis.line.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.line.y  = element_line(linewidth = 0.3),
    axis.ticks.y = element_line(linewidth = 0.3),
    axis.text.y  = element_text(size = 10),
    axis.title   = element_text(size = 10),
    plot.margin  = margin(10, 10, 10, 15)
  )

# -- Dot-and-whisker (% of stand budget +/- SE) -------------------------------

dot_data <- budget %>%
  select(compartment, integrated, integrated_SE) %>%
  mutate(
    pct_budget    = 100 * integrated / total_integrated,
    pct_budget_SE = 100 * integrated_SE / total_integrated
  ) %>%
  mutate(compartment = factor(compartment, levels = comp_labels))

p_dotwhisker <- ggplot(dot_data, aes(x = pct_budget, y = compartment,
                                      color = compartment)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60",
             linewidth = 0.3) +
  geom_pointrange(aes(xmin = pct_budget - pct_budget_SE,
                      xmax = pct_budget + pct_budget_SE),
                  size = 0.5, linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.1f \u00b1 %.1f%%", pct_budget, pct_budget_SE)),
            vjust = -1, hjust = 0.5, size = 3.5, color = "grey20",
            show.legend = FALSE) +
  scale_color_manual(
    values = setNames(unname(compartment_colors), comp_labels),
    guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(
    x = "% of stand budget",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.line  = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    axis.text  = element_text(size = 10),
    axis.title = element_text(size = 10),
    plot.margin = margin(10, 10, 10, 15)
  )

# =============================================================================
# 5. ROW 2: TRUNCATION + BREAK-EVEN PANELS
# =============================================================================

# -- Cone schematic ------------------------------------------------------------

create_cone_panel <- function() {
  stretch <- 5
  heights <- seq(0, tree_height, 0.1)
  radii   <- sapply(heights, radius_at)

  n_seg    <- 50
  h_breaks <- seq(0, tree_height, length.out = n_seg + 1)
  seg_pct  <- numeric(n_seg)
  for (i in seq_len(n_seg)) {
    h_mid <- (h_breaks[i] + h_breaks[i + 1]) / 2
    seg_pct[i] <- if (h_mid == 0) 0 else
      100 * cone_lateral_area(0, h_mid) / total_stem_area
  }
  norm <- (seg_pct - min(seg_pct)) / diff(range(seg_pct))
  pal  <- colorRampPalette(c("#4575b4", "#d73027"))(100)
  seg_col <- pal[pmax(1, pmin(100, round(norm * 99) + 1))]

  max_ext <- max(radii) * stretch * 0.3
  zones <- do.call(rbind, lapply(seq_len(n_seg), function(i) {
    h_seq <- seq(h_breaks[i], h_breaks[i + 1], length.out = 20)
    r_seq <- sapply(h_seq, radius_at) * stretch
    if (h_breaks[i] == 0) {
      h_ext <- c(seq(-max_ext, 0, length.out = 10), h_seq)
      r_ext <- c(rep(r_seq[1], 10), r_seq)
    } else {
      h_ext <- h_seq; r_ext <- r_seq
    }
    data.frame(height = c(h_ext, rev(h_ext)),
               radius = c(r_ext, -rev(r_ext)),
               zone   = paste0("s", i),
               color  = seg_col[i])
  }))

  outline <- data.frame(
    height = c(heights, rev(heights)),
    radius = c(radii * stretch, -rev(radii * stretch))
  )

  make_ellipse <- function(h, n = 100) {
    r <- radius_at(h) * stretch
    theta <- seq(0, 2 * pi, length.out = n)
    data.frame(x = r * cos(theta), y = h + r * 0.3 * sin(theta))
  }
  make_base_arc <- function(h = 0, n = 50) {
    r <- radius_at(h) * stretch
    theta <- seq(pi, 2 * pi, length.out = n)
    data.frame(x = r * cos(theta), y = h + r * 0.3 * sin(theta))
  }

  ell_2  <- make_ellipse(2)
  ell_10 <- make_ellipse(10)
  base_c <- make_base_arc(0)
  xlim   <- max(radii) * stretch * 1.1

  p <- ggplot() +
    geom_polygon(data = zones,
                 aes(x = radius, y = height, fill = color, group = zone),
                 alpha = 0.7) +
    geom_path(data = outline, aes(x = radius, y = height),
              color = "black", linewidth = 1.2) +
    geom_path(data = ell_2,  aes(x = x, y = y),
              color = "black", linewidth = 0.6, linetype = "dashed") +
    geom_path(data = ell_10, aes(x = x, y = y),
              color = "black", linewidth = 0.6, linetype = "dashed") +
    geom_path(data = base_c, aes(x = x, y = y),
              color = "black", linewidth = 1.2) +
    geom_hline(yintercept = c(2, 10), color = "grey60",
               linewidth = 0.4, linetype = "dashed") +
    scale_fill_identity() +
    coord_fixed(ratio = 0.5) +
    scale_x_continuous(limits = c(-xlim, xlim), name = NULL) +
    scale_y_continuous(breaks = seq(0, 25, 5), limits = c(-0.5, 27),
                       name = "Height (m)") +
    theme_minimal(base_size = 11) +
    theme(legend.position    = "none",
          axis.text.x        = element_blank(),
          axis.ticks.x       = element_blank(),
          axis.text.y        = element_text(size = 10),
          axis.title.y       = element_text(size = 11),
          panel.grid.major.x = element_blank(),
          panel.grid.minor   = element_blank(),
          plot.margin = margin(5, 0, 5, 5))
  p
}

# -- Capture curves (clean, no annotations) ------------------------------------

create_capture_panel <- function() {
  h_seq <- seq(0.5, tree_height, 0.5)
  curves <- do.call(rbind, lapply(h_seq, function(h) {
    pcts <- capture_at(h)
    data.frame(height = h,
               category = c("% of stem", "% of stem + branches",
                             "% of stem + branches + leaves"),
               pct = as.numeric(pcts))
  }))
  curves$category <- factor(curves$category,
    levels = c("% of stem", "% of stem + branches",
               "% of stem + branches + leaves"))

  cat_colors <- c("% of stem"                       = "#8B4513",
                   "% of stem + branches"            = "#4682B4",
                   "% of stem + branches + leaves"   = "#2E8B57")

  p <- ggplot(curves, aes(x = pct, y = height, color = category)) +
    geom_line(linewidth = 1.2, alpha = 0.85) +
    geom_hline(yintercept = c(2, 10), color = "grey60",
               linewidth = 0.4, linetype = "dashed") +
    geom_vline(xintercept = 100, color = "grey60",
               linewidth = 0.4, linetype = "dotted") +
    scale_color_manual(values = cat_colors, name = NULL) +
    scale_x_continuous(breaks = seq(0, 100, 25), limits = c(0, 110),
                       name = "% of surface area captured") +
    scale_y_continuous(breaks = seq(0, 25, 5), limits = c(-0.5, 27),
                       name = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position   = c(0.65, 0.15),
          legend.background = element_blank(),
          legend.text       = element_text(size = 9),
          legend.key.width  = unit(1.2, "cm"),
          axis.text         = element_text(size = 10),
          axis.title.x      = element_text(size = 11),
          panel.grid.minor  = element_blank(),
          plot.margin = margin(5, 5, 5, 0))
  p
}

# -- Break-even thresholds with MDF dotted lines -------------------------------

scenarios$scenario <- factor(scenarios$scenario,
  levels = c("Double", "Cancel", "Flip"))
scenarios$scenario_x <- as.numeric(scenarios$scenario)

p_breakeven <- ggplot(scenarios, aes(x = scenario_x, y = F_gt2_required)) +
  geom_hline(yintercept = c(-median_MDF, median_MDF),
             linetype = "dotted", color = "grey50", linewidth = 0.5) +
  annotate("text", x = 3.55, y = median_MDF,
           label = "MDF", size = 3.5, color = "grey50",
           fontface = "italic", vjust = -0.4, hjust = 1) +
  geom_col(width = 0.55, fill = "grey70", color = "grey30", linewidth = 0.3) +
  geom_hline(yintercept = 0, linewidth = 0.6, color = "black") +
  geom_hline(yintercept = F_gt2_observed, linewidth = 0.8,
             linetype = "dashed", color = "#C0392B") +
  annotate("text", x = 3.45, y = F_gt2_observed,
           label = sprintf("Observed: %.2f", F_gt2_observed),
           hjust = 1, vjust = -0.5, size = 3.8, color = "#C0392B",
           fontface = "bold") +
  geom_text(aes(label = sprintf("%.2f", F_gt2_required),
                vjust = ifelse(F_gt2_required >= 0, -0.4, 1.4)),
            size = 3.8, color = "grey20") +
  scale_x_continuous(breaks = 1:3, labels = levels(scenarios$scenario)) +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  labs(
    x = NULL,
    y = expression("Required flux >2 m (nmol m"^{-2}~"s"^{-1}*")")
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.line   = element_line(linewidth = 0.3),
    axis.ticks  = element_line(linewidth = 0.3),
    axis.text   = element_text(size = 10),
    axis.title  = element_text(size = 10),
    plot.margin = margin(10, 15, 10, 10)
  )

# =============================================================================
# 6. ASSEMBLE
# =============================================================================

p_cone   <- create_cone_panel()
p_curves <- create_capture_panel()

# Row 1: flux rates/areas | stacked bar | dot-and-whisker
p_rate_area <- p_rates / p_area + plot_layout(heights = c(1, 1))
row1 <- wrap_elements(p_rate_area) + p_budget + p_dotwhisker +
  plot_layout(widths = c(2, 1, 1))

# Row 2: cone | capture curves | breakeven
row2 <- p_cone + p_curves + p_breakeven +
  plot_layout(widths = c(1, 1.5, 1.5))

# Stack
combined <- row1 / row2 +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(", tag_suffix = ")",
    theme = theme(plot.tag = element_text(size = 12, face = "bold"))
  )

# =============================================================================
# 7. SAVE
# =============================================================================

base_dir <- "/Users/jongewirtzman/My Drive/Research/whole_tree_flux"
out_dir  <- file.path(base_dir, "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "figure_truncation_combined.pdf"), combined,
       width = 13, height = 10, units = "in")
ggsave(file.path(out_dir, "figure_truncation_combined.png"), combined,
       width = 13, height = 10, units = "in", dpi = 300)

cat("Saved: figure_truncation_combined.pdf/.png\n")
