#!/usr/bin/env Rscript
# =============================================================================
# figure2_yale_forest.R
# CH4 flux vs height for 1 Quercus velutina at Yale Myers Forest.
# Single panel; stem measurements only.
# =============================================================================

library(ggplot2)
library(dplyr)
library(stringr)

# ── Paths ────────────────────────────────────────────────────────────────────
base_dir <- "/Users/jongewirtzman/My Drive/Research/whole_tree_flux"
ymf <- read.csv(file.path(base_dir, "data processing", "goFlux_reprocessing",
                           "ymf_black_oak", "results",
                           "ymf_black_oak_flux_compiled.csv"),
                stringsAsFactors = FALSE)

# ── Clean ────────────────────────────────────────────────────────────────────

# Parse height: handle entries like "4 (restarted)", drop "Seam"
ymf <- ymf %>%
  mutate(Height_num = suppressWarnings(
    as.numeric(str_extract(as.character(Height_m), "^[\\d.]+")))
  ) %>%
  filter(!is.na(Height_num), !is.na(CH4_best.flux))

# ── Aesthetics ───────────────────────────────────────────────────────────────

plot_theme <- theme_classic(base_size = 12) +
  theme(
    plot.title  = element_text(face = "italic", size = 12, hjust = 0.5),
    axis.line   = element_line(linewidth = 0.3),
    axis.ticks  = element_line(linewidth = 0.3)
  )

# ── Plot ─────────────────────────────────────────────────────────────────────

p <- ggplot(ymf, aes(x = Height_num, y = CH4_best.flux)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black",
             linewidth = 0.4) +
  geom_point(color = "#8B4513", size = 3, alpha = 0.8) +
  geom_smooth(method = "loess", se = TRUE, color = "#8B4513",
              fill = "#8B4513", alpha = 0.12, linewidth = 0.7,
              na.rm = TRUE) +
  labs(y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
       x = "Height (m)",
       title = "Quercus velutina (Yale Myers Forest)") +
  coord_flip() +
  plot_theme

# ── Save ─────────────────────────────────────────────────────────────────────

out_dir <- file.path(base_dir, "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "figure2_yale_forest.pdf"), p,
       width = 4, height = 4, units = "in")
ggsave(file.path(out_dir, "figure2_yale_forest.png"), p,
       width = 4, height = 4, units = "in", dpi = 300)

cat("Saved figure2_yale_forest.pdf/.png\n")
cat("Observations:", nrow(ymf), "\n")
