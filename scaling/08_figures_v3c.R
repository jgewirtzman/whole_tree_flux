# =============================================================================
# 08_figures_v3c.R — GRL main figures, third pass (Jon, 2026-10-06): annotations moved to captions; one palette with one
# meaning per colour family, used in every figure:
#   compartments  stem < 2 m  #8B4513 (dark brown) · stem >= 2 m #D4A76A (tan) · branches #4682B4 (blue) · leaves #2E8B57 (green)
#   assumptions   extrapolation forms in an orange-red family (constant = amber, decay = vermilion, linear = dark red),
#                 linetype distinguishes variants; zero = grey dotted
#   measured      black;  detection limit (MDF) = light grey band;  soil = muted purple #7B6A9E
#   branch rules  blues (branch family); global woody area sets = browns (woody family)
# Figure 1: four versions (A free log x / B fixed y + linear free x / C fixed y + fixed linear x / D = A + photos).
# Figure 2: row 1 = flux and area by compartment, share of total; row 2 = model tree, area below height, double/cancel/flip.
# Figure 3: (a) forms on an illustrative basal profile, (b) measured relative fluxes (mean, bootstrap CI) with the forms,
#           (c) crossed assumptions vs soil, (d) global demonstration.
# Run: LANG=en_US.UTF-8 Rscript scaling/08_figures_v3c.R   (from whole_tree_flux/)
# =============================================================================
source("scaling/00_load_field.R")
suppressPackageStartupMessages({library(tidyr); library(ggplot2); library(readr); library(patchwork); library(jpeg); library(grid)})
OUT <- file.path(ROOT, "scaling/out"); FD <- file.path(ROOT, "scaling/v3c"); dir.create(FD, showWarnings = FALSE)
FONT <- "Helvetica"
th <- theme_classic(base_size = 8.5, base_family = FONT) + theme(legend.position = "bottom", strip.background = element_blank(),
  strip.text = element_text(size = 7, hjust = 0, lineheight = 0.9), plot.tag = element_text(size = 10, face = "bold"), legend.text = element_text(size = 7))
COMP <- c("Stem < 2 m", "Stem ≥ 2 m", "Branches", "Leaves"); CC <- setNames(c("#8B4513", "#D4A76A", "#4682B4", "#2E8B57"), COMP)
SOIL <- "#7B6A9E"; MDFFILL <- "grey88"
FORM_COL <- c("Constant (basal mean)" = "#E69F00", "Constant (top chamber)" = "#E69F00", "Exponential decay" = "#D55E00",
              "Linear, floored at zero" = "#8B1A1A", "Linear (unbounded)" = "#8B1A1A", "Zero above 2 m" = "grey45")
FORM_LT <- c("Constant (basal mean)" = "solid", "Constant (top chamber)" = "22", "Exponential decay" = "solid",
             "Linear, floored at zero" = "22", "Linear (unbounded)" = "solid", "Zero above 2 m" = "13")
tr <- scales::pseudo_log_trans(sigma = 0.01); lab_sl <- function(x) sub("-", "−", format(x, drop0trailing = TRUE, trim = TRUE, scientific = FALSE))
SP <- c(hem = "T. canadensis", rm = "A. rubrum", ro = "Q. rubra", bg = "N. sylvatica", qv = "Q. velutina")
micro <- read_csv(file.path(OUT, "microsite.csv"), show_col_types = FALSE) %>% select(tree, microsite)

## ============================ Figure 1 versions
bg <- read.csv(file.path(ROOT, "data processing/goFlux_reprocessing/diurnal_blackgum/results/blackgum_flux_compiled_with_mdf.csv")) %>%
  transmute(site = "Black Gum Swamp", tree = "Black Gum Swamp bg (2024)", species = "bg", component = "stem", height_m, flux = CH4_best.flux, below_mdf = CH4_below_MDF)
LAB1 <- tribble(~tree, ~lab, ~ord,
  "EMS hem 321902", 'italic("T. canadensis")~"1"', 1, "Swamp Rd hem 4", 'italic("T. canadensis")~"2"', 2, "EMS ro 300607", 'italic("Q. rubra")', 3,
  "EMS rm 321071", 'italic("A. rubrum")~"1"', 4, "Swamp Rd rm 5", 'italic("A. rubrum")~"2"', 5, "Swamp Rd bg 2", 'italic("N. sylvatica")', 6,
  "YMF black oak", 'atop(italic("Q. velutina"), "(Yale Myers)")', 7, "Black Gum Swamp bg (2024)", 'atop(italic("N. sylvatica"), "(swamp)")', 8)
P1 <- bind_rows(F, bg) %>% left_join(LAB1, by = "tree") %>%
  mutate(comp = factor(case_when(component == "stem" & height_m < 2 ~ COMP[1], component == "stem" ~ COMP[2], component == "branch" ~ COMP[3], TRUE ~ COMP[4]), COMP))
stopifnot(!any(is.na(P1$lab)))
lev1 <- LAB1$lab; P1$lab <- factor(P1$lab, lev1)
XCAP <- c("EMS rm 321071" = 3.2)   # A. rubrum 1: one basal value (5.4) would compress the whole panel
P1 <- P1 %>% group_by(tree, height_m) %>% mutate(k = rank(flux, ties.method = "first"), n = n(),
  step = ifelse(site == "Black Gum Swamp", 0, pmin(0.25, 0.8 / pmax(n - 1, 1))), yp = height_m + (k - (n + 1) / 2) * step) %>% ungroup() %>%
  mutate(cap = XCAP[tree], clip = !is.na(cap) & flux > cap, xp = ifelse(clip, cap, flux),
         big = site != "Black Gum Swamp") %>% select(-k, -n, -step, -cap)
med1 <- P1 %>% filter(component == "stem") %>% group_by(lab, height_m) %>% summarise(m = median(flux), .groups = "drop")
base1 <- function(d, m) ggplot(d, aes(flux, height_m)) + annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 2, fill = "#F3F3F3") +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.3) +
  # LOESS of stem flux on height, with 95% band, over each tree's measured range, drawn under the points
  geom_smooth(data = ~ filter(.x, component == "stem"), aes(flux, height_m), orientation = "y", method = "loess", formula = y ~ x, se = TRUE,
              colour = "grey25", fill = "grey60", alpha = 0.25, linewidth = 0.5, na.rm = TRUE) +
  # points: x capped where a single basal value would compress a panel (drawn at the cap as ">"); points sharing a tree and
  # chamber height are spread evenly in height (ordered by flux, at most ±0.4 m) rather than jittered at random
  geom_point(data = ~ arrange(filter(.x, !below_mdf, !clip), comp), aes(xp, yp, fill = comp, shape = "above MDF", size = big), colour = "white", stroke = 0.3) +
  geom_point(data = ~ arrange(filter(.x, below_mdf, !clip), comp), aes(xp, yp, colour = comp, shape = "below MDF", size = big), fill = "white", stroke = 0.6) +
  geom_text(data = ~ filter(.x, clip), aes(xp, yp, colour = comp), label = ">", size = 3.2, fontface = "bold", show.legend = FALSE) +
  geom_errorbarh(data = ~ filter(.x, site == "Black Gum Swamp") %>% group_by(lab, height_m) %>% summarise(lo = min(flux), hi = max(flux), .groups = "drop"),
                 aes(xmin = lo, xmax = hi, y = height_m), inherit.aes = FALSE, height = 0, linewidth = 0.4, colour = "grey20") +
  geom_point(data = ~ filter(.x, site == "Black Gum Swamp") %>% group_by(lab, height_m) %>% summarise(m = mean(flux), .groups = "drop"),
             aes(m, height_m), inherit.aes = FALSE, shape = 23, fill = "white", colour = "black", size = 2, stroke = 0.5) +
  scale_size_manual(values = c(`TRUE` = 2.4, `FALSE` = 1.1), guide = "none") +
  scale_shape_manual(values = c("above MDF" = 21, "below MDF" = 21), name = NULL) +
  scale_fill_manual(values = CC, name = NULL, drop = FALSE) + scale_colour_manual(values = CC, guide = "none", drop = FALSE) +
  labs(x = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})), y = "Height (m)") + th +
  guides(fill = guide_legend(nrow = 1, order = 1, override.aes = list(shape = 21, colour = "white", size = 3)),
         shape = guide_legend(nrow = 1, order = 2, override.aes = list(fill = c("grey30", "white"), colour = c("white", "grey30"), size = 2.6, stroke = c(0.35, 0.7)))) +
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.9), axis.line = element_blank(), panel.spacing.x = unit(1.5, "mm"))
vA <- base1(P1, med1) + scale_x_continuous(trans = tr, breaks = c(-0.1, 0, 0.1, 1, 10, 100), labels = lab_sl) + facet_wrap(~lab, nrow = 1, scales = "free", labeller = label_parsed) +
  theme(axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust = 1))
sq <- scales::trans_new("sq15", function(y) ifelse(y <= 15, y, 15 + (y - 15) / 3), function(z) ifelse(z <= 15, z, 15 + (z - 15) * 3))
padx <- P1 %>% group_by(lab) %>% summarise(flux = -0.12 * max(abs(flux)), height_m = 0)
vB <- base1(P1, med1) + geom_blank(data = padx) + geom_hline(yintercept = 15, linetype = "13", colour = "grey55", linewidth = 0.3) +
  scale_y_continuous(trans = sq, limits = c(0, 22.5), breaks = c(0, 5, 10, 15, 20)) + facet_wrap(~lab, nrow = 1, scales = "free_x", labeller = label_parsed) +
  theme(axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust = 1))
P1c <- P1 %>% filter(site != "Black Gum Swamp"); med1c <- med1 %>% filter(!grepl("swamp", as.character(lab))) %>% mutate(lab = droplevels(lab)); P1c$lab <- droplevels(P1c$lab)
vC <- base1(P1c, med1c) + scale_y_continuous(limits = c(0, 22.5)) + scale_x_continuous(limits = c(-0.6, 2.6), breaks = c(0, 1, 2)) + facet_wrap(~lab, nrow = 1, labeller = label_parsed)
for (v in list(list("Fig1_A_freeLogX", vA, 8), list("Fig1_B_fixedY_linearFreeX", vB, 8), list("Fig1_C_fixedY_fixedLinearX_noSwamp", vC, 7)))
  ggsave(file.path(FD, paste0(v[[1]], ".png")), v[[2]], width = 190, height = 140, units = "mm", dpi = 300, bg = "white")
img <- function(f, cx = 0.5, cy = 0.5, asp = 2 / 3) { x <- readJPEG(file.path(ROOT, f)); h <- dim(x)[1]; w <- dim(x)[2]
  ww <- min(w, round(h / asp)); hh <- round(ww * asp)
  c0 <- max(1, min(w - ww + 1, round(cx * w - ww / 2))); r0 <- max(1, min(h - hh + 1, round(cy * h - hh / 2)))
  x <- x[r0:(r0 + hh - 1), c0:(c0 + ww - 1), ]
  ggplot() + annotation_raster(x, 0, 1, 0, 1) + coord_fixed(asp, expand = FALSE, xlim = c(0, 1), ylim = c(0, 1)) + theme_void() + theme(plot.tag = element_text(size = 10, face = "bold")) }
NH <- P1 %>% filter(site != "Black Gum Swamp", component == "stem") %>% mutate(hb = cut(height_m, seq(0, 24, 2), right = FALSE, labels = paste0(seq(0, 22, 2), "–", seq(2, 24, 2))))
NS <- NH %>% group_by(hb) %>% summarise(m = mean(flux), se = sd(flux) / sqrt(n()), n = n(), .groups = "drop"); write_csv(NS, file.path(OUT, "fig1_heightbin_means.csv"))
pN_log <- ggplot(NH, aes(flux, hb)) + geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.3) +
  geom_point(aes(colour = comp, shape = below_mdf), position = position_jitter(height = 0.18, width = 0, seed = 3), size = 0.9, stroke = 0.4, alpha = 0.8) +
  geom_errorbarh(data = NS, aes(xmin = m - se, xmax = m + se, y = hb), inherit.aes = FALSE, height = 0, linewidth = 0.5) +
  geom_point(data = NS, aes(m, hb), inherit.aes = FALSE, shape = 23, fill = "white", colour = "black", size = 1.8, stroke = 0.5) +
  scale_colour_manual(values = CC, guide = "none") + scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), guide = "none") +
  scale_x_continuous(trans = tr, breaks = c(-0.1, 0, 0.1, 1, 5), labels = lab_sl) +
  labs(x = expression(Stem~CH[4]~flux~(nmol~m^{-2}~s^{-1})), y = "Height (m)") + th + theme(axis.text.y = element_text(size = 6.5))
XM <- 3; NHl <- NH %>% mutate(clip = flux > XM, fx = pmin(flux, XM))
pN_lin <- ggplot(NHl, aes(fx, hb)) + geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.3) +
  geom_point(data = ~ filter(.x, !clip), aes(colour = comp, shape = below_mdf), position = position_jitter(height = 0.18, width = 0, seed = 3), size = 1, stroke = 0.4, alpha = 0.8) +
  geom_point(data = ~ filter(.x, clip), aes(colour = comp), shape = 62, size = 2.2) +
  geom_errorbarh(data = NS, aes(xmin = m - se, xmax = m + se, y = hb), inherit.aes = FALSE, height = 0, linewidth = 0.5) +
  geom_point(data = NS, aes(m, hb), inherit.aes = FALSE, shape = 23, fill = "white", colour = "black", size = 1.9, stroke = 0.5) +
  scale_colour_manual(values = CC, guide = "none") + scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), guide = "none") +
  scale_x_continuous(limits = c(-0.2, XM), breaks = seq(0, 3, 0.5), expand = expansion(mult = c(0.02, 0.03))) +
  labs(x = expression(Stem~CH[4]~flux~(nmol~m^{-2}~s^{-1})), y = "Height (m)") + th + theme(axis.text.y = element_text(size = 6.5))
mk_lin <- function(XM) { d <- NH %>% mutate(clip = flux > XM, fx = pmin(flux, XM))
  ggplot(d, aes(fx, hb)) + geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.3) +
    geom_point(data = ~ filter(.x, !clip), aes(colour = comp, shape = below_mdf), position = position_jitter(height = 0.18, width = 0, seed = 3), size = 1, stroke = 0.4, alpha = 0.8) +
    geom_point(data = ~ filter(.x, clip), aes(colour = comp), shape = 62, size = 2.2) +
    geom_errorbarh(data = NS, aes(xmin = m - se, xmax = m + se, y = hb), inherit.aes = FALSE, height = 0, linewidth = 0.5) +
    geom_point(data = NS, aes(m, hb), inherit.aes = FALSE, shape = 23, fill = "white", colour = "black", size = 1.9, stroke = 0.5) +
    scale_colour_manual(values = CC, guide = "none") + scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), guide = "none") +
    scale_x_continuous(limits = c(-0.2, XM), expand = expansion(mult = c(0.02, 0.03))) +
    labs(title = paste("axis to", XM), x = expression(Stem~CH[4]~flux~(nmol~m^{-2}~s^{-1})), y = "Height (m)") + th +
    theme(axis.text.y = element_text(size = 6.5), plot.title = element_text(size = 8)) }
cmp <- wrap_plots(lapply(seq(2, 5.5, 0.5), mk_lin), nrow = 2)
ggsave(file.path(FD, "Fig1d_truncation_options.png"), cmp, width = 240, height = 130, units = "mm", dpi = 250, bg = "white")
for (pv in list(list(pN_lin, "Fig1_D_photos_plus_B.png"), list(pN_log, "Fig1_D_photos_plus_B_logd.png"))) {
bot <- (img("IMG_5926_edited.jpg", 0.5, 0.45) | img("IMG_6437.jpg", 0.5, 0.5) | pv[[1]]) + plot_layout(widths = c(1, 1, 0.9))
vD <- vB / bot + plot_layout(heights = c(1, 0.42)) + plot_annotation(tag_levels = "a")
ggsave(file.path(FD, pv[[2]]), vD, width = 190, height = 200, units = "mm", dpi = 300, bg = "white") }

## ============================ Figure 2
hf <- read.csv(file.path(DATA, "data processing/goFlux_reprocessing/results/canopy_flux_goFlux_compiled_with_mdf.csv")) %>%
  filter(!(Species == "bg" & Tree_Tag == 3), !is.na(CH4_best.flux), !is.na(Height_m)) %>% mutate(component = ifelse(Type == "leaf (shaded)", "leaf", Type))
H <- 23   # canopy height: lidar heights of the dominant species in the HF ForestGEO plot, 20.0–22.9 m (Sullivan et al., 2017)
# stand rates: six HF trees, tree-weighted (03_tree_component.R); stem area split from our own trees' diameters
R <- read_csv(file.path(OUT, "stand_rates_HF.csv"), show_col_types = FALSE) %>% mutate(comp = factor(COMP, COMP), int = mean * area, pct = 100 * int / sum(int))
write_csv(R, file.path(OUT, "fig2_stand_shares.csv")); print(R)
mdf <- median(hf$CH4_MDF_wass95, na.rm = TRUE)
FL <- "atop(Flux~per~unit~surface, (nmol~m^{-2}~s^{-1}))"; AL <- "atop(Surface~area, (m^2~m^{-2}~ground))"
R2 <- R %>% transmute(comp, !!FL := mean, !!AL := area) %>% pivot_longer(-comp) %>%
  left_join(R %>% transmute(comp, name = FL, se), by = c("comp", "name")) %>% mutate(name = factor(name, c(FL, AL)))
p2a <- ggplot(R2, aes(comp, value, fill = comp)) + geom_col(width = 0.65) + geom_errorbar(aes(ymin = value - se, ymax = value + se), width = 0.2, linewidth = 0.3, na.rm = TRUE) +
  geom_text(aes(y = value + ifelse(is.na(se), 0, se), label = ifelse(value < 0.1, sprintf("%.3f", value), sprintf("%.2f", value))), vjust = -0.5, size = 2.5) +
  facet_wrap(~name, ncol = 1, scales = "free", strip.position = "left", labeller = label_parsed) + scale_fill_manual(values = CC, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) + labs(x = NULL, y = NULL) + th +
  theme(strip.placement = "outside", strip.text = element_text(size = 7.5, hjust = 0.5), axis.text.x = element_text(size = 7))
p2b <- ggplot(R, aes(1, pct, fill = forcats::fct_rev(comp))) + geom_col(width = 0.6, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.0f%%", pct)), position = position_stack(vjust = 0.5), size = 2.6, colour = "white", fontface = "bold") +
  scale_fill_manual(values = CC, guide = "none") + scale_y_continuous(breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(x = NULL, y = expression(Share~of~tree~surface~CH[4]~flux~("%"))) + th + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())
cone <- tibble(h = seq(0, H, 0.1)) %>% mutate(r = 0.2 * (H - h) / H)
p2c <- ggplot(cone) + geom_ribbon(aes(y = h, xmin = -r, xmax = r), orientation = "y", fill = CC[2]) +
  geom_ribbon(data = cone %>% filter(h <= 2), aes(y = h, xmin = -r, xmax = r), orientation = "y", fill = CC[1]) +
  geom_hline(yintercept = c(2, 10), linetype = "22", linewidth = 0.3, colour = "grey40") + coord_cartesian(xlim = c(-0.35, 0.35)) +
  labs(x = NULL, y = "Height (m)") + th + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())
cap <- tibble(h = seq(0, H, 0.1)) %>% mutate(s = 1 - ((H - h) / H)^2, `Stem` = 100 * s, `Stem + branches` = 100 * s * 0.45 / 2.15, `All surfaces` = 100 * s * 0.45 / 6.65) %>%
  select(-s) %>% pivot_longer(-h) %>% mutate(name = factor(name, c("Stem", "Stem + branches", "All surfaces")))
p2d <- ggplot(cap, aes(value, h, colour = name)) + geom_hline(yintercept = c(2, 10), linetype = "22", linewidth = 0.3, colour = "grey40") + geom_path(linewidth = 0.8) +
  scale_colour_manual(values = c(Stem = CC[[2]], `Stem + branches` = CC[[3]], `All surfaces` = CC[[4]]), name = NULL) + guides(colour = guide_legend(ncol = 1)) +
  labs(x = "Surface area below height (%)", y = NULL) + th +
  theme(legend.position = c(0.3, 0.98), legend.justification = c(0, 1), legend.background = element_blank())
A_up <- sum(R$area[-1]); base <- R$int[1]; obs <- sum(R$int[-1]) / A_up
TH <- tibble(what = factor(c("Double", "Cancel", "Flip"), c("Double", "Cancel", "Flip")), f = c(base, -base, -2 * base) / A_up)
p2e <- ggplot(TH, aes(as.numeric(what), f)) +
  annotate("rect", xmin = 0.4, xmax = 3.6, ymin = -mdf, ymax = mdf, fill = MDFFILL) +
  geom_hline(yintercept = 0, linewidth = 0.3) + geom_col(width = 0.55, fill = "grey35") +
  geom_hline(yintercept = obs, linetype = "dashed", colour = "black", linewidth = 0.6) +
  annotate("text", x = 3.55, y = obs, label = "measured above 2 m", vjust = -0.5, hjust = 1, size = 2.6) +
  annotate("text", x = 3.55, y = mdf, label = "detection limit (median MDF)", vjust = -0.5, hjust = 1, size = 2.6, colour = "grey30") +
  scale_x_continuous(breaks = 1:3, labels = levels(TH$what), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-0.04, 0.1, 0.02), limits = c(-0.04, 0.1), expand = c(0, 0)) +
  labs(x = NULL, y = expression(atop(Uniform~flux~on~surfaces~above~2~m, (nmol~m^{-2}~s^{-1})))) + th
f2 <- ((p2a | p2b) + plot_layout(widths = c(1.6, 0.6))) / ((p2c | p2d | p2e) + plot_layout(widths = c(0.45, 1, 1))) +
  plot_layout(heights = c(1, 1)) + plot_annotation(tag_levels = "a")
ggsave(file.path(FD, "Figure2_v3c.png"), f2, width = 180, height = 180, units = "mm", dpi = 300, bg = "white")

## ============================ Figure 3
hb <- c(0.5, 1.25); fbas <- c(1.3, 0.7); lin <- coef(lm(fbas ~ hb)); ex <- coef(lm(log(fbas) ~ hb)); top <- 0.7; x <- seq(0.3, 22, 0.05)
FORMS <- bind_rows(tibble(form = "Constant (basal mean)", h = x[x >= 2], y = 1), tibble(form = "Constant (top chamber)", h = x[x >= 2], y = top),
  tibble(form = "Exponential decay", h = x, y = exp(ex[1] + ex[2] * x)), tibble(form = "Linear (unbounded)", h = x, y = lin[1] + lin[2] * x),
  tibble(form = "Linear, floored at zero", h = x[x >= 2], y = pmax(lin[1] + lin[2] * x[x >= 2], 0)), tibble(form = "Zero above 2 m", h = x[x >= 2], y = 0)) %>%
  mutate(form = factor(form, names(FORM_COL)))
formlayers <- list(geom_path(data = FORMS, aes(y, h, colour = form, linetype = form), linewidth = 0.75),
  scale_colour_manual(values = FORM_COL, name = NULL), scale_linetype_manual(values = FORM_LT, name = NULL))
p3a <- ggplot() + annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 2, fill = "#F3F3F3") + geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.3) +
  formlayers + geom_point(aes(fbas, hb), shape = 21, fill = CC[1], colour = "black", size = 2.4) + coord_cartesian(xlim = c(-1.2, 3), ylim = c(0, 22)) +
  labs(x = "Flux relative to basal mean", y = "Height (m)") + th + guides(colour = guide_legend(ncol = 3), linetype = guide_legend(ncol = 3))
nb <- F %>% filter(component == "stem") %>% group_by(tree) %>% mutate(bm = mean(tapply(flux[height_m < 2], height_m[height_m < 2], mean))) %>% ungroup() %>%
  mutate(rel = flux / bm) %>% filter(height_m >= 2) %>% mutate(hbin = cut(height_m, c(2, 4, 7, 10, 14, 23), right = FALSE))
set.seed(1); binm <- nb %>% group_by(hbin) %>% summarise(h = mean(height_m), n = n(), m = mean(rel),
  lo = quantile(replicate(2000, mean(sample(rel, replace = TRUE))), 0.025), hi = quantile(replicate(2000, mean(sample(rel, replace = TRUE))), 0.975), .groups = "drop")
write_csv(binm, file.path(OUT, "fig3_relative_means.csv")); print(binm)
p3b <- ggplot() + annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 2, fill = "#F3F3F3") + geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.3) +
  geom_path(data = FORMS, aes(y, h, colour = form, linetype = form), linewidth = 0.5, alpha = 0.45) + scale_colour_manual(values = FORM_COL, guide = "none") +
  scale_linetype_manual(values = FORM_LT, guide = "none") +
  geom_point(data = nb, aes(rel, height_m), colour = CC[2], size = 0.9, alpha = 0.7) +
  geom_errorbarh(data = binm, aes(xmin = lo, xmax = hi, y = h), height = 0, linewidth = 0.6) + geom_point(data = binm, aes(m, h), size = 2.2) +
  coord_cartesian(xlim = c(-1.2, 3), ylim = c(0, 22)) + labs(x = "Flux relative to basal mean", y = NULL) + th
tc <- read_csv(file.path(OUT, "tree_component.csv"), show_col_types = FALSE)
LAB <- c(measured = "Measured above 2 m", exp_decay = "Exponential decay", linear_zero = "Linear, floored at zero", const_top = "Constant (top chamber)",
         const_mean = "Constant (basal mean)", linear = "Linear (unbounded)", zero = "Zero above 2 m")
soil <- read.csv(file.path(ROOT, "scaling/soil_jevon2023/fluxes.csv")) %>% mutate(m = format(as.Date(date, "%m/%d/%y"), "%m"), f = CH4.flux * 1000) %>% filter(m %in% c("07", "08"))
S_mean <- mean(soil$f); S_q <- quantile(soil$f, c(0.25, 0.75))
tc2 <- tc %>% filter(form %in% names(LAB)) %>% mutate(form = factor(LAB[form], rev(LAB)))
BRC <- c(measured_mean = "#4682B4", measured_median = "#9CC3E4", equal_upper_stem = "#1F3F66")
p3c <- ggplot(tc2, aes(woody_component, form)) +
  geom_rect(aes(xmin = S_q[1], xmax = S_q[2], ymin = -Inf, ymax = Inf), data = data.frame(z = 1), inherit.aes = FALSE, fill = SOIL, alpha = 0.18) +
  geom_vline(xintercept = S_mean, colour = SOIL, linewidth = 0.7) + geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.3) +
  geom_point(aes(colour = branch_rule, shape = area_set), position = position_jitter(height = 0.15, seed = 1), size = 1.6, stroke = 0.6) +
  scale_colour_manual(values = BRC, labels = c(measured_mean = "measured mean", measured_median = "measured median", equal_upper_stem = "same as upper stem"), name = "Branch flux") +
  scale_shape_manual(values = c(16, 2), labels = c("Whittaker & Woodwell (1967)", "Gauci et al. (2024)"), name = "Woody area") +
  guides(colour = guide_legend(ncol = 1, title.position = "top", order = 1), shape = guide_legend(ncol = 1, title.position = "top", order = 2)) +
  labs(x = expression(Woody~surface~CH[4]~(nmol~m^{-2}~ground~s^{-1})), y = NULL) + th
WA <- c("90" = 90.4, "132" = 131.9, "207" = 206.6); WAC <- c("#D9B98A", "#A0703C", "#5C3A1A")
gd <- tc %>% filter(form %in% names(LAB)) %>% mutate(wf = woody_component / (stem + branch)) %>% crossing(tibble(wa = names(WA), A = WA)) %>%
  mutate(Tg = wf * A * 0.506, form = factor(LAB[form], rev(LAB)), wa = factor(wa, names(WA)))
p3d <- ggplot(gd, aes(Tg, form)) + geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.3) +
  geom_point(aes(colour = wa), position = position_jitter(height = 0.15, seed = 2), size = 1.6) +
  scale_colour_manual(values = setNames(WAC, names(WA)), name = expression(Global~woody~area~(10^6~km^2))) + guides(colour = guide_legend(nrow = 1, title.position = "top")) +
  labs(x = expression(CH[4]~(Tg~yr^{-1})), y = NULL) + th + theme(axis.text.y = element_blank())
thL <- theme(legend.key.height = unit(3, "mm"), legend.key.width = unit(5, "mm"), legend.spacing.y = unit(0, "mm"), legend.margin = margin(0, 0, 0, 0),
              legend.box.margin = margin(-6, 0, 0, 0), legend.text = element_text(size = 6.5), legend.title = element_text(size = 7))
top3 <- wrap_plots(p3a + labs(tag = "a"), p3b + labs(tag = "b"), nrow = 1) + plot_layout(guides = "collect") & theme(legend.position = "bottom") & thL
# (A) tightened: forms legend in two rows, shorter figure
f3 <- (wrap_elements(full = top3) / (wrap_plots((p3c + labs(tag = "c")) & thL, (p3d + labs(tag = "d")) & thL, nrow = 1, widths = c(1.35, 0.9)))) + plot_layout(heights = c(1, 0.85))
ggsave(file.path(FD, "Figure3_v3c.png"), f3, width = 180, height = 165, units = "mm", dpi = 300, bg = "white")
# (B) equal-proportion 2 x 2
f3e <- (wrap_elements(full = top3) / (wrap_plots((p3c + labs(tag = "c")) & thL, (p3d + labs(tag = "d")) & thL, nrow = 1, widths = c(1, 1)))) + plot_layout(heights = c(1, 1))
ggsave(file.path(FD, "Figure3_v3c_equal.png"), f3e, width = 180, height = 180, units = "mm", dpi = 300, bg = "white")

## ============================ SI figure: Fig 1 profiles on an arcsinh flux axis (shared across trees), larger panels, true heights
ASC <- 0.01   # asinh(x / 0.01): linear within about ±0.01 nmol m-2 s-1, logarithmic beyond
tr_as <- scales::trans_new("asinh01", function(x) asinh(x / ASC), function(y) sinh(y) * ASC)
PS <- P1 %>% arrange(below_mdf, comp)
pS1 <- ggplot(PS, aes(flux, height_m)) + annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 2, fill = "#F3F3F3") +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.3) +
  geom_point(data = ~ filter(.x, !below_mdf), aes(fill = comp, shape = "above MDF"), colour = "white", size = 2.6, stroke = 0.3) +
  geom_point(data = ~ filter(.x, below_mdf), aes(colour = comp, shape = "below MDF"), fill = "white", size = 2.3, stroke = 0.7) +
  scale_shape_manual(values = c("above MDF" = 21, "below MDF" = 21), name = NULL) +
  scale_fill_manual(values = CC, name = NULL, drop = FALSE) + scale_colour_manual(values = CC, guide = "none", drop = FALSE) +
  scale_x_continuous(trans = tr_as, breaks = c(-0.1, 0, 0.1, 1, 10, 100), minor_breaks = c(-0.01, 0.01), labels = lab_sl) +
  scale_y_continuous(limits = c(0, 22.5), breaks = seq(0, 20, 5)) +
  facet_wrap(~lab, ncol = 4, labeller = label_parsed) +
  labs(x = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})*","~arcsinh~scale), y = "Height (m)") + th +
  guides(fill = guide_legend(nrow = 1, order = 1, override.aes = list(shape = 21, colour = "white", size = 3)),
         shape = guide_legend(nrow = 1, order = 2, override.aes = list(fill = c("grey30", "white"), colour = c("white", "grey30"), size = 2.6, stroke = c(0.3, 0.7)))) +
  theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.7), axis.line = element_blank(),
        panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3), panel.grid.minor.x = element_line(colour = "grey95", linewidth = 0.25), axis.text.x = element_text(size = 6.5))
ggsave(file.path(FD, "FigureS_profiles_asinh.png"), pS1, width = 190, height = 190, units = "mm", dpi = 300, bg = "white")
