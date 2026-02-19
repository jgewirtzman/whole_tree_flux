#!/usr/bin/env Rscript
# =============================================================================
# figure_sensitivity_breakeven.R
#
# Sensitivity analysis: stand-level CH₄ budget with vertical truncation at 2 m.
#
# All areas are expressed per m² ground area using Whittaker & Woodwell (1967)
# surface area indices for a temperate deciduous forest:
#   stem bark   = 0.45 m² / m² ground
#   branch bark = 1.70 m² / m² ground
#   LAI         = 4.5  m² / m² ground
#
# The cone taper model (26 m tree, 40 cm DBH) is used ONLY to partition the
# 0.45 m² of stem bark into below-2m vs above-2m fractions. Branches and
# leaves are entirely above 2 m.
#
# Part 1 — STAND-LEVEL BUDGET
#   Multiply observed per-tissue CH₄ flux rates by literature surface area
#   indices to get integrated flux per m² ground for each compartment:
#     (a) stem < 2 m
#     (b) stem ≥ 2 m
#     (c) branches
#     (d) leaves
#
# Part 2 — BREAK-EVEN THRESHOLDS
#   What uniform average flux rate across ALL >2 m surfaces would be needed to:
#     DOUBLE  — double the <2 m integrated emission
#     CANCEL  — exactly offset the <2 m emission (net = 0)
#     FLIP    — make the stand a net sink equal in magnitude to <2 m source
#
# Output: figure_sensitivity_breakeven.pdf/.png + tables to console
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ═══════════════════════════════════════════════════════════════════════════════
# 1. SURFACE AREA INDICES (per m² ground)
# ═══════════════════════════════════════════════════════════════════════════════

# Whittaker & Woodwell (1967) stand-level indices
stem_bark_index   <- 0.45   # m² stem bark  / m² ground
branch_bark_index <- 1.70   # m² branch bark / m² ground
leaf_area_index   <- 4.5    # m² leaf / m² ground

# ── Cone taper model (only used to split stem into <2m / ≥2m) ───────────────
# Representative canopy tree: 26 m tall, 40 cm DBH
tree_height      <- 26      # m
dbh              <- 0.40    # m
dbh_height       <- 1.37    # m
base_radius      <- dbh / 2
height_above_dbh <- tree_height - dbh_height

radius_at <- function(h) {
  pmax(0, base_radius * (tree_height - h) / height_above_dbh)
}

cone_lateral_area <- function(h1, h2) {
  r1 <- radius_at(h1)
  r2 <- radius_at(h2)
  slant <- sqrt((h2 - h1)^2 + (r1 - r2)^2)
  pi * (r1 + r2) * slant
}

# Fraction of stem bark below 2 m (from cone geometry)
total_cone_area <- cone_lateral_area(0, tree_height)
cone_below_2m   <- cone_lateral_area(0, 2)
frac_stem_lt2   <- cone_below_2m / total_cone_area  # ~0.148

# ── Stand-level areas (m² surface per m² ground) ────────────────────────────
A_stem_lt2  <- stem_bark_index * frac_stem_lt2           # stem below 2 m
A_stem_ge2  <- stem_bark_index * (1 - frac_stem_lt2)     # stem above 2 m
A_branch    <- branch_bark_index                          # all branches (>2 m)
A_leaf      <- leaf_area_index                            # all leaves (>2 m)
A_gt2_total <- A_stem_ge2 + A_branch + A_leaf             # total >2 m surface
A_total     <- A_stem_lt2 + A_gt2_total                   # total tree surface

cat("═══════════════════════════════════════════════════════════════════\n")
cat("Stand-level surface areas (m² surface per m² ground)\n")
cat("  Whittaker & Woodwell (1967); cone split at 2 m\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat(sprintf("  Stem < 2 m:     %.4f  (%4.1f%% of total)\n",
            A_stem_lt2, 100 * A_stem_lt2 / A_total))
cat(sprintf("  Stem ≥ 2 m:     %.4f  (%4.1f%%)\n",
            A_stem_ge2, 100 * A_stem_ge2 / A_total))
cat(sprintf("  Branch bark:    %.4f  (%4.1f%%)\n",
            A_branch, 100 * A_branch / A_total))
cat(sprintf("  Leaf area:      %.4f  (%4.1f%%)\n",
            A_leaf, 100 * A_leaf / A_total))
cat(sprintf("  Total:          %.4f\n", A_total))
cat(sprintf("  Total > 2 m:    %.4f  (%4.1f%%)\n\n",
            A_gt2_total, 100 * A_gt2_total / A_total))

# ═══════════════════════════════════════════════════════════════════════════════
# 2. EMPIRICAL FLUX DATA (Harvard Forest)
# ═══════════════════════════════════════════════════════════════════════════════

dat <- read.csv(file.path("data processing", "goFlux_reprocessing",
                           "results", "canopy_flux_goFlux_compiled.csv"),
                stringsAsFactors = FALSE)

dat <- dat %>%
  filter(!(Species == "bg" & Tree_Tag == 3)) %>%
  mutate(Component = ifelse(Type == "leaf (shaded)", "leaf", Type)) %>%
  filter(!is.na(CH4_best.flux), !is.na(Height_m))

# Per-tissue flux summaries
stem_lt2_dat  <- dat %>% filter(Component == "stem", Height_m < 2)
stem_ge2_dat  <- dat %>% filter(Component == "stem", Height_m >= 2)
branch_dat    <- dat %>% filter(Component == "branch")
leaf_dat      <- dat %>% filter(Component == "leaf")

F_stem_lt2   <- mean(stem_lt2_dat$CH4_best.flux, na.rm = TRUE)
F_stem_ge2   <- mean(stem_ge2_dat$CH4_best.flux, na.rm = TRUE)
F_branch     <- mean(branch_dat$CH4_best.flux, na.rm = TRUE)
F_leaf       <- mean(leaf_dat$CH4_best.flux, na.rm = TRUE)

cat("═══════════════════════════════════════════════════════════════════\n")
cat("Empirical per-tissue CH₄ flux rates (nmol m⁻² s⁻¹)\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat(sprintf("  Stem <2 m:   mean = %.4f  median = %.4f  (n = %d)\n",
            F_stem_lt2,
            median(stem_lt2_dat$CH4_best.flux, na.rm = TRUE),
            nrow(stem_lt2_dat)))
cat(sprintf("  Stem ≥2 m:   mean = %.4f  median = %.4f  (n = %d)\n",
            F_stem_ge2,
            median(stem_ge2_dat$CH4_best.flux, na.rm = TRUE),
            nrow(stem_ge2_dat)))
cat(sprintf("  Branch:      mean = %.4f  median = %.4f  (n = %d)\n",
            F_branch,
            median(branch_dat$CH4_best.flux, na.rm = TRUE),
            nrow(branch_dat)))
cat(sprintf("  Leaf:        mean = %.4f  median = %.4f  (n = %d)\n\n",
            F_leaf,
            median(leaf_dat$CH4_best.flux, na.rm = TRUE),
            nrow(leaf_dat)))

# ═══════════════════════════════════════════════════════════════════════════════
# 3. STAND-LEVEL CH₄ BUDGET
#    Flux rate × area index = integrated flux per m² ground
# ═══════════════════════════════════════════════════════════════════════════════

budget <- tibble(
  compartment = c("Stem < 2 m", "Stem ≥ 2 m", "Branches", "Leaves"),
  flux_rate   = c(F_stem_lt2, F_stem_ge2, F_branch, F_leaf),
  area_index  = c(A_stem_lt2, A_stem_ge2, A_branch, A_leaf)
) %>%
  mutate(
    integrated = flux_rate * area_index,              # nmol m⁻² ground s⁻¹
    pct        = 100 * integrated / sum(integrated)   # % of stand total
  )

# >2 m aggregate
gt2_integrated <- sum(budget$integrated[budget$compartment != "Stem < 2 m"])
total_integrated <- sum(budget$integrated)
gt2_pct <- 100 * gt2_integrated / total_integrated

cat("═══════════════════════════════════════════════════════════════════\n")
cat("Stand-level CH₄ budget (flux rate × area index)\n")
cat("  Units: nmol CH₄ per m² ground per second\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat(sprintf("  %-15s  %12s  %12s  %16s  %8s\n",
            "Compartment", "Flux rate", "Area index", "Integrated", "% total"))
cat(strrep("─", 70), "\n")
for (i in seq_len(nrow(budget))) {
  b <- budget[i, ]
  cat(sprintf("  %-15s  %11.4f   %11.4f    %14.4f   %6.1f%%\n",
              b$compartment, b$flux_rate, b$area_index, b$integrated, b$pct))
}
cat(strrep("─", 70), "\n")
cat(sprintf("  %-15s  %11s   %11.4f    %14.4f   %6.1f%%\n",
            "TOTAL", "", A_total, total_integrated, 100))
cat(sprintf("  %-15s  %11s   %11.4f    %14.4f   %6.1f%%\n\n",
            "All > 2 m", "", A_gt2_total, gt2_integrated, gt2_pct))

# Area-weighted mean flux across all >2 m surfaces
F_gt2_observed <- gt2_integrated / A_gt2_total
cat(sprintf("  Observed area-weighted mean flux >2 m: %.4f nmol m⁻² s⁻¹\n\n",
            F_gt2_observed))

# ═══════════════════════════════════════════════════════════════════════════════
# 4. BREAK-EVEN CALCULATION
#    What uniform flux across all >2 m surface (per m² of that surface)
#    would double / cancel / flip the <2 m emission?
# ═══════════════════════════════════════════════════════════════════════════════

baseline <- budget$integrated[budget$compartment == "Stem < 2 m"]  # nmol/m²ground/s

scenarios <- tibble(
  scenario    = c("Double", "Cancel", "Flip"),
  description = c("Net = 2× basal emission",
                   "Net = 0 (break-even)",
                   "Net = −1× basal (equal sink)"),
  multiplier  = c(2, 0, -1)
) %>%
  mutate(
    net_target     = multiplier * baseline,
    F_gt2_required = (net_target - baseline) / A_gt2_total
  )

cat("═══════════════════════════════════════════════════════════════════\n")
cat("Break-even analysis\n")
cat(sprintf("  Baseline (stem <2 m): %.4f nmol m⁻²ground s⁻¹\n", baseline))
cat(sprintf("  Total >2 m area index: %.4f m²surface / m²ground\n", A_gt2_total))
cat("═══════════════════════════════════════════════════════════════════\n")
cat(sprintf("  %-10s  %-35s  %18s\n",
            "Scenario", "Description", "Required F_gt2"))
cat(strrep("─", 70), "\n")
for (i in seq_len(nrow(scenarios))) {
  s <- scenarios[i, ]
  cat(sprintf("  %-10s  %-35s  %15.4f nmol m⁻² s⁻¹\n",
              s$scenario, s$description, s$F_gt2_required))
}
cat(sprintf("\n  Observed area-weighted >2 m flux:   %.4f nmol m⁻² s⁻¹\n",
            F_gt2_observed))
cat(sprintf("  Ratio observed / cancel threshold:  %.1f×\n\n",
            abs(F_gt2_observed / scenarios$F_gt2_required[scenarios$scenario == "Cancel"])))

# ═══════════════════════════════════════════════════════════════════════════════
# 5. FIGURES
# ═══════════════════════════════════════════════════════════════════════════════

# Compartment colors — consistent across panels
compartment_colors <- c(
  "Stem < 2 m"  = "#8B4513",   # saddle brown
  "Stem ≥ 2 m"  = "#D4A76A",   # light tan
  "Branches"    = "#4682B4",   # steel blue
  "Leaves"      = "#2E8B57"    # sea green
)

# Scenario colors — distinct from compartment colors
scenario_colors <- c(
  "Double" = "#E8A87C",    # peach / apricot
  "Cancel" = "#A8A8A8",    # neutral grey
  "Flip"   = "#7B6BA1"     # muted purple
)

# ── Panel A: Per-tissue flux rates with area annotations ─────────────────────

rate_data <- tibble(
  compartment = factor(c("Stem < 2 m", "Stem ≥ 2 m", "Branches", "Leaves"),
                        levels = c("Stem < 2 m", "Stem ≥ 2 m", "Branches", "Leaves")),
  flux_rate   = c(F_stem_lt2, F_stem_ge2, F_branch, F_leaf),
  area_index  = c(A_stem_lt2, A_stem_ge2, A_branch, A_leaf),
  area_pct    = c(100 * A_stem_lt2 / A_total,
                  100 * A_stem_ge2 / A_total,
                  100 * A_branch / A_total,
                  100 * A_leaf / A_total)
)

p_rates <- ggplot(rate_data, aes(x = compartment, y = flux_rate,
                                  fill = compartment)) +
  geom_col(width = 0.6, color = "grey30", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.3f\n(%.1f%% of area)", flux_rate, area_pct)),
            vjust = -0.3, size = 3, color = "grey20") +
  scale_fill_manual(values = compartment_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(
    x = NULL,
    y = expression("Mean CH"[4]~"flux rate (nmol m"^{-2}~"s"^{-1}*")")
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.line  = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    axis.text  = element_text(size = 10),
    axis.text.x = element_text(size = 9),
    axis.title = element_text(size = 10),
    plot.margin = margin(10, 10, 10, 10)
  )

# ── Panel B: Integrated budget (stacked bar) ─────────────────────────────────

budget_plot <- budget
budget_plot$compartment <- factor(budget_plot$compartment,
  levels = c("Leaves", "Branches", "Stem ≥ 2 m", "Stem < 2 m"))

p_budget <- ggplot(budget_plot, aes(x = 1, y = integrated, fill = compartment)) +
  geom_col(width = 0.5, color = "grey30", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%s\n%.1f%%", compartment, pct)),
            position = position_stack(vjust = 0.5),
            size = 3, color = "white", fontface = "bold") +
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

# ── Panel C: Break-even thresholds vs observed ───────────────────────────────

scenarios$scenario <- factor(scenarios$scenario,
  levels = c("Double", "Cancel", "Flip"))

p_breakeven <- ggplot(scenarios, aes(x = scenario, y = F_gt2_required,
                                      fill = scenario)) +
  geom_col(width = 0.55, color = "grey30", linewidth = 0.3) +
  geom_hline(yintercept = 0, linewidth = 0.6, color = "black") +
  # Observed empirical rate
  geom_hline(yintercept = F_gt2_observed, linewidth = 0.8,
             linetype = "dashed", color = "#C0392B") +
  annotate("text", x = 3.45, y = F_gt2_observed,
           label = sprintf("Observed: %.4f", F_gt2_observed),
           hjust = 1, vjust = -0.5, size = 3.2, color = "#C0392B",
           fontface = "bold") +
  geom_text(aes(label = sprintf("%.4f", F_gt2_required),
                vjust = ifelse(F_gt2_required >= 0, -0.4, 1.4)),
            size = 3.2, color = "grey20") +
  scale_fill_manual(values = scenario_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  labs(
    x = NULL,
    y = expression("Required uniform flux over >2 m surface (nmol m"^{-2}~"s"^{-1}*")")
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.line  = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    axis.text  = element_text(size = 10),
    axis.title = element_text(size = 10),
    plot.margin = margin(10, 15, 10, 10)
  )

# ── Combine ──────────────────────────────────────────────────────────────────

combined <- (p_rates | p_budget | p_breakeven) +
  plot_layout(widths = c(2, 1, 2)) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(", tag_suffix = ")",
    theme = theme(plot.tag = element_text(size = 12, face = "bold"))
  )

# ═══════════════════════════════════════════════════════════════════════════════
# 6. SAVE
# ═══════════════════════════════════════════════════════════════════════════════

base_dir <- "/Users/jongewirtzman/My Drive/Research/whole_tree_flux"
out_dir  <- file.path(base_dir, "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "figure_sensitivity_breakeven.pdf"), combined,
       width = 12, height = 5.5, units = "in")
ggsave(file.path(out_dir, "figure_sensitivity_breakeven.png"), combined,
       width = 12, height = 5.5, units = "in", dpi = 300)

cat("Saved: figure_sensitivity_breakeven.pdf/.png\n")
