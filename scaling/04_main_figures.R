# =============================================================================
# 04_main_figures.R — GRL main-text Figures 2 and 3 (v3) from the scaling outputs.
#   Figure 2: extrapolating basal (< 2 m) stem fluxes upward vs the same trees' measured fluxes >= 2 m: ratio of
#             predicted to measured mean flux over the stem above 2 m (stem-area weighted), bootstrap 95% range, and
#             the number of trees each form predicts to be net sinks above 2 m (measured: none).
#   Figure S4 (was Figure 3 in the first v3 pass; main Figure 3 is now 05_figure3_truncation_v3.R): (a) share of tree surface area below a given sampling height (cone 26 m + Whittaker & Woodwell 1967);
#             (b) tree woody-surface component (stems + branches) per m2 ground vs same-month soil uptake (Jevon 2023);
#             (c) DEMONSTRATION ONLY: the same per-area rates times published global woody surface area.
# Inputs: scaling/out/form_test_bootstrap.csv, form_test_scores.csv, tree_component.csv
# Run: LANG=en_US.UTF-8 Rscript scaling/04_main_figures.R   (from whole_tree_flux/)
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr); library(ggplot2); library(patchwork)})
ROOT <- "."; OUT <- file.path(ROOT, "scaling/out")   # run all scripts from the repository root
th <- theme_classic(base_size = 8.5, base_family = "Helvetica") + theme(plot.title = element_text(size = 8.5, face = "bold"), legend.position = "bottom")
LAB <- c(measured = "Measured above 2 m", exp_decay = "Exponential decay", linear_zero = "Linear, floored at zero", const_top = "Constant (top chamber)",
         const_mean = "Constant (basal mean)", linear = "Linear (unbounded)", zero = "Zero above 2 m")
## ---- Figure 2
bs <- read_csv(file.path(OUT, "form_test_bootstrap.csv"), show_col_types = FALSE) %>% filter(form %in% names(LAB))
pt <- read_csv(file.path(OUT, "form_test_pooled.csv"), show_col_types = FALSE) %>% select(form, ratio = ratio_pred_obs_cone, trees_pred_net_uptake, trees)
d2 <- bs %>% left_join(pt, by = "form") %>% mutate(form = factor(LAB[form], rev(LAB[c("const_mean", "const_top", "exp_decay", "linear_zero", "linear", "zero")])))
f2 <- ggplot(d2, aes(ratio, form)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = "#EAF2F8") +
  geom_vline(xintercept = 1, colour = "black", linewidth = 0.4) + geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = pmax(ratio_lo, -2.9), xmax = ratio_hi), height = 0.25, linewidth = 0.4) +
  geom_segment(data = ~ filter(.x, ratio_lo < -2.9), aes(x = -2.5, xend = -2.95, yend = form), linewidth = 0.4,
               arrow = arrow(length = unit(1.6, "mm"), type = "closed")) +
  geom_text(data = ~ filter(.x, ratio_lo < -2.9), aes(x = -2.9, label = sprintf("to %.0f", ratio_lo)), vjust = -0.7, hjust = 0, size = 2.4) +
  geom_point(size = 2.2) +
  geom_text(aes(x = 3.2, label = sprintf("%d of %d", trees_pred_net_uptake, trees)), hjust = 0, size = 2.7) +
  annotate("text", x = 3.2, y = 6.6, label = "trees predicted\nnet sink above 2 m", hjust = 0, vjust = 0, size = 2.5, lineheight = 0.9) +
  annotate("text", x = 1.05, y = 6.55, label = "= measured", hjust = 0, vjust = 0, size = 2.5) +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2, 3)) + coord_cartesian(xlim = c(-3, 4.3), clip = "off") +
  labs(x = "Predicted / measured mean stem flux above 2 m (stem-area weighted)", y = "Basal flux extrapolated as:") + th +
  theme(plot.margin = margin(14, 30, 4, 4))
ggsave(file.path(ROOT, "scaling/Figure2_v3.png"), f2, width = 120, height = 75, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(ROOT, "scaling/Figure2_v3.pdf"), f2, width = 120, height = 75, units = "mm", device = cairo_pdf)
cat("Figure 2 values:\n"); print(d2 %>% select(form, ratio, ratio_lo, ratio_hi, trees_pred_net_uptake, trees))

## ---- Figure 3a: capture fractions (cone 26 m; Whittaker & Woodwell 1967 stem 0.45, branch 1.70, LAI 4.5; branches and
##      leaves treated as above any sampled stem height)
H <- 26; A <- c(stem = 0.45, branch = 1.70, leaf = 4.5)
cap <- tibble(h = seq(0, H, 0.1)) %>% mutate(stem_frac = 1 - ((H - h) / H)^2,
  `Stem only` = 100 * stem_frac, `Stem + branches` = 100 * stem_frac * A["stem"] / (A["stem"] + A["branch"]),
  `Stem + branches + leaves` = 100 * stem_frac * A["stem"] / sum(A)) %>% select(-stem_frac) %>% pivot_longer(-h, names_to = "denominator", values_to = "pct") %>%
  mutate(denominator = factor(denominator, c("Stem only", "Stem + branches", "Stem + branches + leaves")))
at <- cap %>% filter(abs(h - 2) < 1e-9 | abs(h - 10) < 1e-9); cat("\ncapture fractions:\n"); print(at)
f3a <- ggplot(cap, aes(pct, h, colour = denominator)) + geom_hline(yintercept = 2, linetype = "22", colour = "grey50", linewidth = 0.3) + geom_path(linewidth = 0.7) +
  scale_colour_manual(values = c("#8B6B4A", "#E69F00", "#009E73"), name = NULL) + guides(colour = guide_legend(ncol = 1)) +
  labs(x = "Tree surface area below this height (%)", y = "Sampling height (m)", title = "a  Surface sampled below a given height") + th
## ---- Figure 3b: tree component vs soil (from 03_tree_component.R)
tc <- read_csv(file.path(OUT, "tree_component.csv"), show_col_types = FALSE)
soil <- read.csv(file.path(ROOT, "scaling/soil_jevon2023/fluxes.csv")) %>% mutate(m = format(as.Date(date, "%m/%d/%y"), "%m"), f = CH4.flux * 1000) %>% filter(m %in% c("07", "08"))
S_mean <- mean(soil$f); S_q <- quantile(soil$f, c(0.25, 0.75))
lev <- c("measured", "exp_decay", "linear_zero", "const_top", "const_mean", "linear", "zero")
tc2 <- tc %>% filter(form %in% lev) %>% mutate(form = factor(LAB[form], rev(LAB[lev])))
f3b <- ggplot(tc2, aes(woody_component, form)) +
  annotate("rect", xmin = S_q[1], xmax = S_q[2], ymin = -Inf, ymax = Inf, fill = "#D6E6F2") + geom_vline(xintercept = S_mean, colour = "#0072B2", linewidth = 0.6) +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.3) +
  geom_point(aes(colour = branch_rule, shape = area_set), position = position_jitter(height = 0.15, seed = 1), size = 1.5, alpha = 0.85) +
  annotate("text", x = S_mean, y = 7.45, label = "soil, Jul–Aug", colour = "#0072B2", size = 2.5, vjust = 0) +
  scale_colour_manual(values = c(measured_mean = "#E69F00", measured_median = "#009E73", equal_upper_stem = "#CC79A7"),
                      labels = c(measured_mean = "branch: measured mean", measured_median = "branch: measured median", equal_upper_stem = "branch = upper stem"), name = NULL) +
  scale_shape_manual(values = c(16, 2), labels = c("W&W 1967 areas", "WAI 3.07"), name = NULL) +
  guides(colour = guide_legend(ncol = 1), shape = guide_legend(ncol = 1)) + coord_cartesian(clip = "off") +
  labs(x = expression(CH[4]~(nmol~m^{-2}~ground~s^{-1})), y = "Stem above 2 m:", title = "b  Tree woody surfaces vs soil") + th
## ---- Figure 3c: DEMONSTRATION (woody per-area flux x global woody area; Gauci 2024 Extended Data Table 4)
WA <- c("MODIS VCF (90 M km²)" = 90.4, "Consensus (132 M km²)" = 131.9, "Consensus + shrubs (207 M km²)" = 206.6); K <- 0.506
gd <- tc %>% filter(form %in% lev) %>% mutate(woody_flux = woody_component / (stem + branch)) %>% crossing(tibble(wa = names(WA), A = WA)) %>%
  mutate(Tg = woody_flux * A * K, form = factor(LAB[form], rev(LAB[lev])), wa = factor(wa, names(WA)))
write_csv(gd, file.path(OUT, "global_demo_v3.csv")); cat("\nDEMONSTRATION range (Tg yr-1):", paste(round(range(gd$Tg), 1), collapse = " to "),
  "| measured above 2 m:", paste(round(range(gd$Tg[gd$form == LAB["measured"]]), 1), collapse = " to "), "\n")
f3c <- ggplot(gd, aes(Tg, form)) + geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.3) +
  geom_point(aes(colour = wa), position = position_jitter(height = 0.15, seed = 2), size = 1.5, alpha = 0.85) +
  scale_colour_manual(values = c("#56B4E9", "#0072B2", "#000000"), name = NULL) + guides(colour = guide_legend(ncol = 1)) +
  labs(x = expression(Woody~surface~CH[4]~(Tg~yr^{-1})), y = NULL, title = "c  Demonstration, not an estimate",
       subtitle = "Same per-area rates × global woody area") + th +
  theme(axis.text.y = element_blank(), plot.subtitle = element_text(size = 7, colour = "grey30"))
f3 <- (f3a | f3b | f3c) + plot_layout(widths = c(0.75, 1, 0.85)) & theme(legend.text = element_text(size = 6.5), legend.key.size = unit(3, "mm"))
ggsave(file.path(ROOT, "scaling/FigureS4_v3.png"), f3, width = 190, height = 105, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(ROOT, "scaling/FigureS4_v3.pdf"), f3, width = 190, height = 105, units = "mm", device = cairo_pdf)
