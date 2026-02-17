#!/usr/bin/env Rscript
# =============================================================================
# figure3_wu_reanalysis.R
# Reanalysis of Wu et al. (2024) compiled stem CH4 flux data.
# Faceted by study; points colored by ecosystem (Upland / Wetland).
# Includes only studies with at least one measurement at >= 2 m.
# =============================================================================

library(readxl)
library(ggplot2)
library(dplyr)
library(stringr)

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

# Filter to studies with any measurement >= 2 m
refs_2m <- dat %>% filter(Height >= 2) %>% distinct(Reference) %>% pull(Reference)
dat_filt <- dat %>% filter(Reference %in% refs_2m)

# Clean reference labels: "Author Year (Ecosystem)"
dat_filt <- dat_filt %>%
  mutate(
    ref_short = str_replace(Reference, "_[A-Za-z-]+$", "") %>%
      str_replace_all("_", " "),
    facet_label = paste0(ref_short, " (", ecosystem, ")")
  )

# Order facets alphabetically
dat_filt$facet_label <- factor(dat_filt$facet_label,
                               levels = sort(unique(dat_filt$facet_label)))

# ── Aesthetics ───────────────────────────────────────────────────────────────

eco_colors <- c(Upland = "#d7191c", Wetland = "#2c7bb6")

n_facets <- length(unique(dat_filt$facet_label))
plot_h   <- ceiling(n_facets / 4) * 2.5 + 0.8

plot_theme <- theme_classic(base_size = 11) +
  theme(
    strip.text       = element_text(size = 8.5, face = "bold"),
    strip.background = element_blank(),
    legend.position  = "bottom",
    legend.title     = element_blank(),
    legend.text      = element_text(size = 10),
    axis.line        = element_line(linewidth = 0.3),
    axis.ticks       = element_line(linewidth = 0.3),
    panel.spacing    = unit(0.6, "lines")
  )

# ── Plot ─────────────────────────────────────────────────────────────────────

p <- ggplot(dat_filt, aes(x = Height, y = CH4_flux,
                          color = ecosystem, fill = ecosystem)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black",
             linewidth = 0.4) +
  geom_point(size = 1.8, alpha = 0.6, shape = 16) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 0.6,
              alpha = 0.10, na.rm = TRUE) +
  scale_color_manual(values = eco_colors) +
  scale_fill_manual(values = eco_colors) +
  facet_wrap(~ facet_label, scales = "free", ncol = 4) +
  labs(y = expression(Stem~CH[4]~flux~(nmol~m^{-2}~s^{-1})),
       x = "Height (m)") +
  coord_flip() +
  plot_theme

# ── Save ─────────────────────────────────────────────────────────────────────

out_dir <- file.path(base_dir, "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "figure3_wu_reanalysis.pdf"), p,
       width = 7.5, height = plot_h, units = "in", limitsize = FALSE)
ggsave(file.path(out_dir, "figure3_wu_reanalysis.png"), p,
       width = 7.5, height = plot_h, units = "in", dpi = 300, limitsize = FALSE)

cat("Saved figure3_wu_reanalysis.pdf/.png\n")
cat("Studies:", n_facets, "| Observations:", nrow(dat_filt), "\n")
