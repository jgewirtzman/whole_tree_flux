# SI prediction/observation figure and global sensitivity table. Run from repository root.
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr); library(ggplot2); library(patchwork)})
source("scaling/analysis_helpers.R")
ROOT <- "."; OUT <- file.path(ROOT, "scaling/out")   # run all scripts from the repository root
th <- theme_classic(base_size = 8.5, base_family = "Helvetica") + theme(plot.title = element_text(size = 8.5, face = "bold"), legend.position = "bottom")
LAB <- c(measured = "Measured above 2 m", exp_decay = "Exponential decay", linear_zero = "Linear, floored at zero", const_top = "Constant (top chamber)",
         const_mean = "Constant (basal mean)", linear = "Linear (unbounded)", zero = "Zero above 2 m")
## ---- Figure 2
bs <- read_csv(file.path(OUT, "form_test_bootstrap.csv"), show_col_types = FALSE) %>% filter(form %in% names(LAB)) %>% select(-trees)
pt <- read_csv(file.path(OUT, "form_test_pooled.csv"), show_col_types = FALSE) %>% select(form, ratio = ratio_pred_obs_area, trees_pred_net_uptake, trees)
d2 <- bs %>% left_join(pt, by = "form") %>% mutate(form = factor(LAB[form], rev(LAB[c("const_mean", "const_top", "exp_decay", "linear_zero", "linear", "zero")])))
label_x <- max(3.2, max(d2$ratio_hi, na.rm=TRUE) + 0.3)
f2 <- ggplot(d2, aes(ratio, form)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = "#EAF2F8") +
  geom_vline(xintercept = 1, colour = "black", linewidth = 0.4) + geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = pmax(ratio_lo, -4.9), xmax = ratio_hi), height = 0.25, linewidth = 0.4) +
  geom_segment(data = ~ filter(.x, ratio_lo < -4.9), aes(x = -4.5, xend = -4.95, yend = form), linewidth = 0.4,
               arrow = arrow(length = unit(1.6, "mm"), type = "closed")) +
  geom_text(data = ~ filter(.x, ratio_lo < -4.9), aes(x = -4.9, label = sprintf("to %.0f", ratio_lo)), vjust = -0.7, hjust = 0, size = 2.4) +
  geom_point(size = 2.2) +
  geom_text(aes(x = label_x, label = sprintf("%d of %d", trees_pred_net_uptake, trees)), hjust = 0, size = 2.7) +
  annotate("text", x = label_x, y = 6.6, label = "trees predicted\nnet sink above 2 m", hjust = 0, vjust = 0, size = 2.5, lineheight = 0.9) +
  annotate("text", x = 1.05, y = 6.55, label = "= measured", hjust = 0, vjust = 0, size = 2.5) +
  scale_x_continuous(breaks = seq(-4, ceiling(label_x), 2)) + coord_cartesian(xlim = c(-5, label_x + 1.1), clip = "off") +
  labs(x = "Predicted / measured upper-stem flux (area weighted)", y = "Basal flux extrapolated as:") + th +
  theme(plot.margin = margin(14, 30, 4, 4))
ggsave(file.path(ROOT, "scaling/Figure2_v3.png"), f2, width = 160, height = 80, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(ROOT, "scaling/Figure2_v3.pdf"), f2, width = 160, height = 80, units = "mm", device = cairo_pdf)
cat("Figure 2 values:\n"); print(d2 %>% select(form, ratio, ratio_lo, ratio_hi, trees_pred_net_uptake, trees))

## Global sensitivity table: convert ground-area rates to woody-surface rates.
tc <- read_csv(file.path(OUT, "tree_component.csv"), show_col_types = FALSE)
lev <- c("measured", "exp_decay", "linear_zero", "const_top", "const_mean", "linear", "zero")
WA <- c("MODIS VCF (90 M km²)" = 90.4, "Consensus (132 M km²)" = 131.9, "Consensus + shrubs (207 M km²)" = 206.6); K <- 0.506
gd <- tc %>% filter(form %in% lev) %>% mutate(woody_flux = woody_component / (stem + branch)) %>% crossing(tibble(wa = names(WA), A = WA)) %>%
  mutate(Tg = woody_flux * A * K, form = factor(LAB[form], rev(LAB[lev])), wa = factor(wa, names(WA)))
write_csv(gd, file.path(OUT, "global_demo_v3.csv")); cat("\nDEMONSTRATION range (Tg yr-1):", paste(round(range(gd$Tg), 1), collapse = " to "),
  "| measured above 2 m:", paste(round(range(gd$Tg[gd$form == LAB["measured"]]), 1), collapse = " to "), "\n")
