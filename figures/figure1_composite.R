#!/usr/bin/env Rscript
# =============================================================================
# figure1_composite.R
# Composite Figure 1 combining:
#   (a) Harvard Forest 6-panel flux profiles (top, full width)
#   (b) Yale Myers Forest single-panel profile (bottom-left)
#   (c) IMG_5926_edited.jpg – DinoLift photo (bottom-right top)
#   (d) IMG_6437.jpg – arborist climbing photo (bottom-right bottom)
#
# Run from the project root (whole_tree_flux/).
# =============================================================================

library(ggplot2)
library(dplyr)
library(stringr)
library(patchwork)
library(jpeg)
library(grid)

# =============================================================================
# Helper: center-crop a JPEG array to a target aspect ratio (w/h)
# =============================================================================
crop_to_aspect <- function(img, target_ratio) {
  h <- dim(img)[1]; w <- dim(img)[2]
  current_ratio <- w / h
  if (current_ratio > target_ratio) {
    # too wide — crop sides
    new_w <- round(h * target_ratio)
    offset <- round((w - new_w) / 2)
    img <- img[, (offset + 1):(offset + new_w), , drop = FALSE]
  } else {
    # too tall — crop top/bottom
    new_h <- round(w / target_ratio)
    offset <- round((h - new_h) / 2)
    img <- img[(offset + 1):(offset + new_h), , , drop = FALSE]
  }
  img
}

# =============================================================================
# Panel (a): Harvard Forest
# =============================================================================

dat <- read.csv(file.path("data processing", "goFlux_reprocessing",
                           "results", "canopy_flux_goFlux_compiled.csv"),
                stringsAsFactors = FALSE)

dat <- dat %>%
  filter(!(Species == "bg" & Tree_Tag == 3)) %>%
  mutate(Component = ifelse(Type == "leaf (shaded)", "leaf", Type)) %>%
  filter(!is.na(CH4_best.flux), !is.na(Height_m))

species_lookup <- c(bg = "Nyssa sylvatica", rm = "Acer rubrum",
                    ro = "Quercus rubra",   hem = "Tsuga canadensis")
dat$Species_full <- species_lookup[dat$Species]
dat <- dat %>% mutate(Tree_label = paste0(Species_full, " (", Site, ")"))

dat$Tree_label <- factor(
  dat$Tree_label,
  levels = dat %>%
    distinct(Tree_label, Site, Species_full) %>%
    arrange(Site, Species_full) %>%
    pull(Tree_label)
)

component_colors <- c(stem = "#8B4513", branch = "#4682B4", leaf = "#2E8B57")

theme_panel <- theme_classic(base_size = 10) +
  theme(
    strip.text       = element_text(face = "italic", size = 8),
    strip.background = element_blank(),
    legend.position  = "none",
    axis.line        = element_line(linewidth = 0.3),
    axis.ticks       = element_line(linewidth = 0.3),
    panel.spacing    = unit(0.6, "lines"),
    plot.tag         = element_text(size = 13, face = "bold"),
    plot.margin      = margin(2, 4, 2, 2)
  )

p_hf <- ggplot(dat, aes(x = Height_m, y = CH4_best.flux, color = Component)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40",
             linewidth = 0.3) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(data = dat %>% filter(Component == "stem"),
              method = "loess", se = TRUE, fill = "#8B4513",
              alpha = 0.10, linewidth = 0.6, na.rm = TRUE,
              show.legend = FALSE) +
  scale_color_manual(values = component_colors,
                     labels = c(stem = "Stem", branch = "Branch",
                                leaf = "Leaf")) +
  facet_wrap(~ Tree_label, scales = "free", ncol = 3) +
  labs(y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
       x = "Height (m)", tag = "a") +
  coord_flip() +
  theme_panel

# Extract legend as a standalone grob
legend_plot <- ggplot(dat, aes(x = Height_m, y = CH4_best.flux,
                               color = Component)) +
  geom_point() +
  scale_color_manual(values = component_colors,
                     labels = c(stem = "Stem", branch = "Branch",
                                leaf = "Leaf")) +
  theme_void() +
  theme(legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.text      = element_text(size = 9),
        legend.key.size  = unit(0.4, "cm"),
        legend.spacing.x = unit(0.2, "cm"))

# Extract legend without cowplot dependency
tmp <- ggplot_gtable(ggplot_build(legend_plot))
legend_idx <- which(grepl("guide-box", sapply(tmp$grobs, function(x) x$name)))
p_legend <- tmp$grobs[[legend_idx[1]]]

# =============================================================================
# Panel (b): Yale Myers Forest
# =============================================================================

ymf <- read.csv(file.path("data processing", "goFlux_reprocessing",
                           "ymf_black_oak", "results",
                           "ymf_black_oak_flux_compiled.csv"),
                stringsAsFactors = FALSE)

ymf <- ymf %>%
  mutate(Height_num = suppressWarnings(
    as.numeric(str_extract(as.character(Height_m), "^[\\d.]+")))
  ) %>%
  filter(!is.na(Height_num), !is.na(CH4_best.flux))

theme_ymf <- theme_classic(base_size = 10) +
  theme(
    plot.title  = element_text(face = "italic", size = 9, hjust = 0.5),
    axis.line   = element_line(linewidth = 0.3),
    axis.ticks  = element_line(linewidth = 0.3),
    plot.tag    = element_text(size = 13, face = "bold"),
    plot.margin = margin(2, 4, 2, 2)
  )

p_ymf <- ggplot(ymf, aes(x = Height_num, y = CH4_best.flux)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40",
             linewidth = 0.3) +
  geom_point(color = "#8B4513", size = 3, alpha = 0.8) +
  geom_smooth(method = "loess", se = TRUE, color = "#8B4513",
              fill = "#8B4513", alpha = 0.10, linewidth = 0.6,
              na.rm = TRUE) +
  labs(y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
       x = "Height (m)",
       title = "Quercus velutina (Yale Myers Forest)",
       tag = "b") +
  coord_flip() +
  theme_ymf

# =============================================================================
# Panels (c) and (d): Field photos — center-cropped to 4:3
# =============================================================================

img_c_raw <- readJPEG("IMG_5926_edited.jpg")
img_d_raw <- readJPEG("IMG_6437.jpg")

# Crop both to 4:3 landscape (w/h = 1.33)
img_c <- crop_to_aspect(img_c_raw, 4 / 3)
img_d <- crop_to_aspect(img_d_raw, 4 / 3)

make_photo <- function(img, tag_label) {
  g <- rasterGrob(img, interpolate = TRUE)
  ggplot() +
    annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_void() +
    theme(plot.tag    = element_text(size = 13, face = "bold"),
          plot.margin = margin(1, 1, 1, 1)) +
    labs(tag = tag_label)
}

p_photo_c <- make_photo(img_c, "c")
p_photo_d <- make_photo(img_d, "d")

# =============================================================================
# Composite layout using patchwork
#
# Design (8 rows × 4 cols):
#   Rows 1–4: Harvard Forest panels (a), full width
#   Row 5:    legend, full width
#   Rows 6–7: Yale (b) left half | photo (c) right half
#   Row 8:    Yale (b) left half | photo (d) right half  [not needed — see below]
#
# Simpler: use nested wrap_plots for the bottom row.
# =============================================================================

# Bottom-right: two photos stacked
right_col <- p_photo_c / p_photo_d + plot_layout(heights = c(1, 1))

# Bottom row: Yale | photos
bottom_row <- (p_ymf | right_col) + plot_layout(widths = c(1, 1))

# Full composite: Harvard / legend / bottom
composite <- p_hf / wrap_elements(p_legend) / bottom_row +
  plot_layout(heights = c(4, 0.3, 3))

# =============================================================================
# Save
# =============================================================================

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

ggsave("figures/figure1_composite.pdf", composite,
       width = 7.5, height = 10, units = "in")
ggsave("figures/figure1_composite.png", composite,
       width = 7.5, height = 10, units = "in", dpi = 300)

cat("Saved figures/figure1_composite.pdf/.png\n")
