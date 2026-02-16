# =============================================================================
# 08_ch4_height_plot.R
# Plot CH4 flux vs. tree height, colored by component type (stem, branch, leaf),
# facet-wrapped by tree individual.
#
# Three versions:
#   1. Fixed flux axis across facets
#   2. Free flux axis per facet
#   3. Fixed axis with asinh transform (handles near-zero and negative values)
# =============================================================================

# Source setup (works both when source()'d and run interactively in RStudio)
setup_path <- file.path(
  "/Users/jongewirtzman/My Drive/Research/whole_tree_flux",
  "data processing", "goFlux_reprocessing", "00_setup.R")
source(setup_path)

library(ggplot2)
library(scales)

# --- Load compiled data ------------------------------------------------------

dat <- read.csv(file.path(results_dir, "canopy_flux_goFlux_compiled.csv"),
                stringsAsFactors = FALSE)

# --- Clean up component type -------------------------------------------------
# Combine "leaf (shaded)" into "leaf"
dat$Component <- ifelse(dat$Type == "leaf (shaded)", "leaf", dat$Type)

# --- Species full names and facet ordering -----------------------------------

species_lookup <- c(
  "bg"  = "Nyssa sylvatica",
  "rm"  = "Acer rubrum",
  "ro"  = "Quercus rubra",
  "hem" = "Tsuga canadensis"
)

dat$Species_full <- species_lookup[dat$Species]

# Remove Nyssa sylvatica tag 3
dat <- dat[!(dat$Species == "bg" & dat$Tree_Tag == 3), ]

# Facet label: Site - Species (no tag)
dat$Tree_label <- paste0(dat$Site, " - ", dat$Species_full)

# Order by Site then Species_full so facets group within site
dat$Tree_label <- factor(dat$Tree_label,
                         levels = unique(dat$Tree_label[order(dat$Site, dat$Species_full)]))

# --- Define component aesthetics ---------------------------------------------

component_colors <- c(
  "stem"   = "#D2691E",
  "branch" = "#4682B4",
  "leaf"   = "#2E8B57"
)

# --- Common theme -------------------------------------------------------------

plot_theme <- theme_classic(base_size = 13) +
  theme(
    strip.text       = element_text(face = "italic", size = 11),
    strip.background = element_blank(),
    legend.position  = "bottom",
    legend.text      = element_text(size = 11),
    axis.line        = element_line(linewidth = 0.4),
    axis.ticks       = element_line(linewidth = 0.3),
    panel.spacing    = unit(1, "lines")
  )

# --- Base plot layers ---------------------------------------------------------

base_layers <- list(
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4),
  geom_smooth(aes(fill = Component), method = "loess", se = TRUE,
              alpha = 0.15, linewidth = 0.8),
  geom_point(size = 3, alpha = 0.85),
  scale_color_manual(values = component_colors, name = "Component"),
  scale_fill_manual(values = component_colors, name = "Component"),
  labs(
    y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    x = "Height (m)"
  ),
  plot_theme
)

# =============================================================================
# Plot 1: Fixed flux axis
# =============================================================================

p_fixed <- ggplot(dat, aes(x = Height_m, y = CH4_best.flux, color = Component)) +
  base_layers +
  facet_wrap(~ Tree_label) +
  coord_flip()

print(p_fixed)

ggsave(file.path(plots_dir, "CH4_flux_by_height_fixed.pdf"),
       plot = p_fixed, width = 12, height = 8)
ggsave(file.path(plots_dir, "CH4_flux_by_height_fixed.png"),
       plot = p_fixed, width = 12, height = 8, dpi = 300)

# =============================================================================
# Plot 2: Free flux axis
# =============================================================================

p_free <- ggplot(dat, aes(x = Height_m, y = CH4_best.flux, color = Component)) +
  base_layers +
  facet_wrap(~ Tree_label, scales = "free_x") +
  coord_flip()

print(p_free)

ggsave(file.path(plots_dir, "CH4_flux_by_height_free.pdf"),
       plot = p_free, width = 12, height = 8)
ggsave(file.path(plots_dir, "CH4_flux_by_height_free.png"),
       plot = p_free, width = 12, height = 8, dpi = 300)

# =============================================================================
# Plot 3: Fixed axis with asinh transform
# =============================================================================

p_asinh <- ggplot(dat, aes(x = Height_m, y = CH4_best.flux, color = Component)) +
  base_layers +
  facet_wrap(~ Tree_label) +
  scale_y_continuous(trans = "asinh") +
  coord_flip()

print(p_asinh)

ggsave(file.path(plots_dir, "CH4_flux_by_height_asinh.pdf"),
       plot = p_asinh, width = 12, height = 8)
ggsave(file.path(plots_dir, "CH4_flux_by_height_asinh.png"),
       plot = p_asinh, width = 12, height = 8, dpi = 300)

# =============================================================================
# Summary
# =============================================================================

message("Plots saved to: ", plots_dir)
message("  CH4_flux_by_height_fixed.pdf/png  (fixed axis)")
message("  CH4_flux_by_height_free.pdf/png   (free flux axis)")
message("  CH4_flux_by_height_asinh.pdf/png  (asinh transform)")
