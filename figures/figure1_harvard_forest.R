#!/usr/bin/env Rscript
# =============================================================================
# figure1_harvard_forest.R
# CH4 flux vs height for 6 individual Harvard Forest trees.
# One panel per tree; points colored by component (stem / branch / leaf).
# =============================================================================

library(ggplot2)
library(dplyr)

# ── Paths ────────────────────────────────────────────────────────────────────
base_dir <- "/Users/jongewirtzman/My Drive/Research/whole_tree_flux"
dat <- read.csv(file.path(base_dir, "data processing", "goFlux_reprocessing",
                           "results", "canopy_flux_goFlux_compiled.csv"),
                stringsAsFactors = FALSE)

# ── Clean ────────────────────────────────────────────────────────────────────

# Remove Nyssa sylvatica tag 3
dat <- dat %>% filter(!(Species == "bg" & Tree_Tag == 3))

# Combine "leaf (shaded)" into "leaf"
dat <- dat %>% mutate(Component = ifelse(Type == "leaf (shaded)", "leaf", Type))

# Keep valid fluxes / heights
dat <- dat %>% filter(!is.na(CH4_best.flux), !is.na(Height_m))

# Species lookup
species_lookup <- c(bg = "Nyssa sylvatica", rm = "Acer rubrum",
                    ro = "Quercus rubra",   hem = "Tsuga canadensis")
dat$Species_full <- species_lookup[dat$Species]

# Facet label: italic species (site)
dat <- dat %>%
  mutate(Tree_label = paste0(Species_full, " (", Site, ")"))

# Order: group by site, then alphabetically by species
dat$Tree_label <- factor(
  dat$Tree_label,
  levels = dat %>%
    distinct(Tree_label, Site, Species_full) %>%
    arrange(Site, Species_full) %>%
    pull(Tree_label)
)

# ── Aesthetics ───────────────────────────────────────────────────────────────

component_colors <- c(stem = "#8B4513", branch = "#4682B4", leaf = "#2E8B57")

plot_theme <- theme_classic(base_size = 12) +
  theme(
    strip.text       = element_text(face = "italic", size = 10),
    strip.background = element_blank(),
    legend.position  = "bottom",
    legend.title     = element_blank(),
    legend.text      = element_text(size = 10),
    axis.line        = element_line(linewidth = 0.3),
    axis.ticks       = element_line(linewidth = 0.3),
    panel.spacing    = unit(0.8, "lines")
  )

# ── Plot ─────────────────────────────────────────────────────────────────────

p <- ggplot(dat, aes(x = Height_m, y = CH4_best.flux, color = Component)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black",
             linewidth = 0.4) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(data = dat %>% filter(Component == "stem"),
              method = "loess", se = TRUE, fill = "#8B4513",
              alpha = 0.12, linewidth = 0.7, na.rm = TRUE,
              show.legend = FALSE) +
  scale_color_manual(values = component_colors,
                     labels = c(stem = "Stem", branch = "Branch", leaf = "Leaf")) +
  facet_wrap(~ Tree_label, scales = "free", ncol = 3) +
  labs(y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
       x = "Height (m)") +
  coord_flip() +
  plot_theme

# ── Save ─────────────────────────────────────────────────────────────────────

out_dir <- file.path(base_dir, "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "figure1_harvard_forest.pdf"), p,
       width = 7.5, height = 6, units = "in")
ggsave(file.path(out_dir, "figure1_harvard_forest.png"), p,
       width = 7.5, height = 6, units = "in", dpi = 300)

cat("Saved figure1_harvard_forest.pdf/.png\n")
cat("Trees:", length(unique(dat$Tree_label)), "\n")
cat("Observations:", nrow(dat), "\n")
