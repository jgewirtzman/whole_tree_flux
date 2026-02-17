#!/usr/bin/env Rscript
# =============================================================================
# figure_truncation.R
# Combined figure: 3D cone schematic + capture curves + bar chart.
#
# Demonstrates the fraction of tree surface area sampled when measuring
# only below 2 m or 10 m, for a representative NE canopy tree.
#
# Assumptions:
#   - 26 m tall, 40 cm DBH (typical canopy tree at Harvard / Yale forests)
#   - Stem tapers linearly from radius at ground to 0 at tree top (cone)
#   - Whittaker & Woodwell (1967) surface area ratios:
#       stem bark   = 0.45 m² m⁻² ground
#       branch bark = 1.70 m² m⁻² ground
#   - LAI = 4.5 m² m⁻² ground
#   - Only stem surface area is measured (no branches or leaves captured),
#     so the numerator is always stem SA below h; the denominator varies
#     to show how little of the whole tree is actually sampled.
#
# Output: 2-panel figure (cone | capture curves)
# =============================================================================

library(ggplot2)
library(patchwork)

# ═══════════════════════════════════════════════════════════════════════════════
# 1. TREE GEOMETRY
# ═══════════════════════════════════════════════════════════════════════════════

tree_height      <- 26      # m
dbh              <- 0.40    # m (40 cm)
dbh_height       <- 1.37    # m
base_radius      <- dbh / 2
height_above_dbh <- tree_height - dbh_height

# Whittaker & Woodwell surface area indices
stem_bark_index   <- 0.45   # m² stem bark / m² ground
branch_bark_index <- 1.70   # m² branch bark / m² ground
leaf_area_index   <- 4.5    # m² leaf / m² ground

# ── Helper functions ──────────────────────────────────────────────────────────

radius_at <- function(h) {
  pmax(0, base_radius * (tree_height - h) / height_above_dbh)
}

cone_lateral_area <- function(h1, h2) {
  r1 <- radius_at(h1)
  r2 <- radius_at(h2)
  slant <- sqrt((h2 - h1)^2 + (r1 - r2)^2)
  pi * (r1 + r2) * slant
}

# ── Totals ────────────────────────────────────────────────────────────────────

total_stem_area   <- cone_lateral_area(0, tree_height)
ground_area       <- total_stem_area / stem_bark_index
total_branch_area <- ground_area * branch_bark_index
total_leaf_area   <- ground_area * leaf_area_index

cat("Tree parameters:\n")
cat("  Height:", tree_height, "m | DBH:", dbh * 100, "cm\n")
cat("  Total stem SA:  ", round(total_stem_area, 1), "m²\n")
cat("  Total branch SA:", round(total_branch_area, 1), "m²\n")
cat("  Total leaf SA:  ", round(total_leaf_area, 1), "m²\n\n")

# ── Capture fraction at height h ──────────────────────────────────────────────
# Only stem surface area below h is measured (no branches or leaves captured).
# Three denominators show the fraction of progressively larger totals.

capture_at <- function(h) {
  stem_below <- cone_lateral_area(0, h)

  c(stem       = 100 * stem_below / total_stem_area,
    stem_br    = 100 * stem_below / (total_stem_area + total_branch_area),
    stem_br_lf = 100 * stem_below / (total_stem_area + total_branch_area +
                                      total_leaf_area))
}

# Print summary table
for (h in c(2, 10)) {
  pcts <- capture_at(h)
  cat(sprintf("At %2d m:  %% of stem = %4.1f%%  |  %% of stem+branch = %4.1f%%  |  %% of stem+branch+leaf = %4.1f%%\n",
              h, pcts["stem"], pcts["stem_br"], pcts["stem_br_lf"]))
}
cat("\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 2. PANEL A: 3-D CONE SCHEMATIC (gradient fill + ellipses)
# ═══════════════════════════════════════════════════════════════════════════════

create_cone_panel <- function() {

  stretch <- 5
  heights <- seq(0, tree_height, 0.1)
  radii   <- sapply(heights, radius_at)

  # ── Gradient-filled segments ──
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

  # Build polygon data for each segment
  max_ext <- max(radii) * stretch * 0.3
  zones <- do.call(rbind, lapply(seq_len(n_seg), function(i) {
    h_seq <- seq(h_breaks[i], h_breaks[i + 1], length.out = 20)
    r_seq <- sapply(h_seq, radius_at) * stretch
    # Extend bottom segment below 0 to cover base ellipse
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

  # Outline
  outline <- data.frame(
    height = c(heights, rev(heights)),
    radius = c(radii * stretch, -rev(radii * stretch))
  )

  # Ellipses (3-D perspective)
  make_ellipse <- function(h, n = 100) {
    r <- radius_at(h) * stretch
    theta <- seq(0, 2 * pi, length.out = n)
    data.frame(x = r * cos(theta),
               y = h + r * 0.3 * sin(theta))
  }
  make_base_arc <- function(h = 0, n = 50) {
    r <- radius_at(h) * stretch
    theta <- seq(pi, 2 * pi, length.out = n)
    data.frame(x = r * cos(theta),
               y = h + r * 0.3 * sin(theta))
  }

  ell_2  <- make_ellipse(2)
  ell_10 <- make_ellipse(10)
  base_c <- make_base_arc(0)

  xlim <- max(radii) * stretch * 1.1

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
    # Height annotations
    annotate("text", x = xlim * 0.85, y = 2,  label = "2 m",
             size = 3.5, fontface = "bold") +
    annotate("text", x = xlim * 0.65, y = 10, label = "10 m",
             size = 3.5, fontface = "bold") +
    scale_fill_identity() +
    coord_fixed(ratio = 0.5) +
    scale_x_continuous(limits = c(-xlim, xlim), name = NULL) +
    scale_y_continuous(breaks = seq(0, 25, 5), limits = c(-0.5, 27),
                       name = "Height (m)") +
    theme_minimal(base_size = 11) +
    theme(legend.position  = "none",
          axis.text.x      = element_blank(),
          axis.ticks.x     = element_blank(),
          axis.text.y      = element_text(size = 10),
          axis.title.y     = element_text(size = 11),
          panel.grid.major.x = element_blank(),
          panel.grid.minor   = element_blank(),
          plot.margin = margin(5, 0, 5, 5))

  list(plot = p, seg_pct_range = range(seg_pct))
}

# ═══════════════════════════════════════════════════════════════════════════════
# 3. PANEL B: CAPTURE CURVES (3 lines: stem, stem+branch, stem+branch+leaf)
# ═══════════════════════════════════════════════════════════════════════════════

create_capture_panel <- function(seg_range) {

  h_seq <- seq(0.5, tree_height, 0.5)
  curves <- do.call(rbind, lapply(h_seq, function(h) {
    pcts <- capture_at(h)
    data.frame(height = h,
               category = c("% of stem", "% of stem + branches", "% of stem + branches + leaves"),
               pct = as.numeric(pcts))
  }))
  curves$category <- factor(curves$category,
    levels = c("% of stem", "% of stem + branches", "% of stem + branches + leaves"))

  cat_colors <- c("% of stem"                        = "#8B4513",
                   "% of stem + branches"             = "#4682B4",
                   "% of stem + branches + leaves"    = "#2E8B57")

  p <- ggplot(curves, aes(x = pct, y = height, color = category, linetype = category)) +
    geom_line(linewidth = 1.2, alpha = 0.85) +
    # Dashed horizontal reference lines
    geom_hline(yintercept = c(2, 10), color = "grey60",
               linewidth = 0.4, linetype = "dashed") +
    geom_vline(xintercept = 100, color = "grey60",
               linewidth = 0.4, linetype = "dotted") +
    scale_color_manual(values = cat_colors, name = NULL) +
    scale_linetype_manual(values = c("solid", "longdash", "dotted"), name = NULL) +
    scale_x_continuous(breaks = seq(0, 100, 25), limits = c(0, 110),
                       name = "% of surface area captured") +
    scale_y_continuous(breaks = seq(0, 25, 5), limits = c(-0.5, 27),
                       name = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position  = c(0.65, 0.15),
          legend.background = element_rect(fill = alpha("white", 0.8),
                                           color = NA),
          legend.text       = element_text(size = 9),
          legend.key.width  = unit(1.2, "cm"),
          axis.text         = element_text(size = 10),
          axis.title.x      = element_text(size = 11),
          panel.grid.minor  = element_blank(),
          plot.margin = margin(5, 5, 5, 0))

  p
}

# ═══════════════════════════════════════════════════════════════════════════════
# 4. PANEL C: BAR CHART (% captured at 2 m and 10 m)
# ═══════════════════════════════════════════════════════════════════════════════

create_bar_panel <- function() {

  bar_data <- do.call(rbind, lapply(c(2, 10), function(h) {
    pcts <- capture_at(h)
    data.frame(
      threshold = paste0(h, " m"),
      category  = c("% of stem", "% of stem + branches", "% of stem + branches + leaves"),
      pct       = as.numeric(pcts)
    )
  }))
  bar_data$category <- factor(bar_data$category,
    levels = c("% of stem", "% of stem + branches", "% of stem + branches + leaves"))

  cat_colors <- c("% of stem"                        = "#8B4513",
                   "% of stem + branches"             = "#4682B4",
                   "% of stem + branches + leaves"    = "#2E8B57")

  p <- ggplot(bar_data, aes(x = threshold, y = pct, fill = category)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65,
             color = "black", linewidth = 0.3) +
    geom_text(aes(label = paste0(round(pct, 1), "%")),
              position = position_dodge(width = 0.75),
              vjust = -0.5, size = 3, fontface = "bold") +
    scale_fill_manual(values = cat_colors, name = NULL) +
    scale_y_continuous(limits = c(0, 75), breaks = seq(0, 75, 25),
                       expand = expansion(mult = c(0, 0.05)),
                       name = "% of surface area captured") +
    labs(x = "Maximum measurement height") +
    theme_classic(base_size = 11) +
    theme(legend.position = "bottom",
          legend.text     = element_text(size = 9),
          axis.line       = element_line(linewidth = 0.3),
          axis.ticks      = element_line(linewidth = 0.3),
          axis.title      = element_text(size = 11),
          axis.text       = element_text(size = 10),
          plot.margin     = margin(5, 10, 5, 5))

  p
}

# ═══════════════════════════════════════════════════════════════════════════════
# 5. COMBINE & SAVE
# ═══════════════════════════════════════════════════════════════════════════════

cone_result <- create_cone_panel()
p_cone    <- cone_result$plot
p_curves  <- create_capture_panel(cone_result$seg_pct_range)
p_bars    <- create_bar_panel()

# Layout: top row = cone + curves (matched y-axes), bottom row = bar chart
combined <- (p_cone | p_curves) / p_bars +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(", tag_suffix = ")")

# ── Save ──────────────────────────────────────────────────────────────────────

base_dir <- "/Users/jongewirtzman/My Drive/Research/whole_tree_flux"
out_dir  <- file.path(base_dir, "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "figure_truncation.pdf"), combined,
       width = 8, height = 9, units = "in")
ggsave(file.path(out_dir, "figure_truncation.png"), combined,
       width = 8, height = 9, units = "in", dpi = 300)

cat("Saved: figure_truncation.pdf/.png\n")
