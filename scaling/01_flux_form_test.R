# =============================================================================
# 01_flux_form_test.R — GRL: test extrapolation forms against measured fluxes above 2 m.
# For each tree, stem fluxes measured below 2 m (0.5, 1.0-1.7 m) are used the way a basal-only study would use them;
# each candidate form is extrapolated upward and compared with the same tree's measured stem fluxes at >= 2 m.
# Forms (all fitted to the per-height mean flux below 2 m of that tree):
#   const_mean  : constant = mean flux below 2 m (no decline)
#   const_top   : constant = flux at the highest height below 2 m (top-chamber value carried upward)
#   exp_decay   : f = a exp(-b h) through the below-2 m means (only if all means > 0); b < 0 (increase) is bounded at 0,
#                 i.e. no growth with height, so an increasing lower profile is carried upward as constant at the top value
#   exp_pooled  : the tree's top-of-band value decaying at the pooled below-2 m log-linear rate (shape shared across trees)
#   linear      : straight line through the below-2 m means, unbounded (may cross zero into uptake)
#   linear_zero : the same line floored at zero (decline to zero, never uptake)
#   zero        : nothing above 2 m (truncation)
# A rising lower profile (slope > 0) is carried upward as constant at the top value for every decline form.
# Scores per tree and pooled: bias and MAE at each measured >= 2 m point; mean over the >= 2 m profile (per-height means,
# unweighted and cone-weighted with the tallest measured height as the top); predicted vs observed negative share.
# Output: scaling/out/form_test_points.csv, form_test_scores.csv, form_test_pooled.csv; scaling/fig_form_test.png
# Run: LANG=en_US.UTF-8 Rscript scaling/01_flux_form_test.R   (from whole_tree_flux/)
# =============================================================================
source("scaling/00_load_field.R")
suppressPackageStartupMessages({library(tidyr); library(ggplot2); library(readr)})
OUT <- file.path(ROOT, "scaling/out"); dir.create(OUT, showWarnings = FALSE)
S <- F %>% filter(component == "stem")
low <- S %>% filter(height_m < 2) %>% group_by(tree, height_m) %>% summarise(f = mean(flux), .groups = "drop")
high <- S %>% filter(height_m >= 2)
# pooled relative decline below 2 m: slope of log(mean flux) on height across trees with all-positive means (tree intercepts)
pos <- low %>% group_by(tree) %>% filter(all(f > 0), n() >= 2) %>% ungroup()
B_POOL <- if (nrow(pos)) coef(lm(log(f) ~ height_m + tree, data = pos))[["height_m"]] else NA
cat("pooled log-linear slope below 2 m (per m):", round(B_POOL, 3), "\n")
forms <- function(lo) {   # lo: per-height means below 2 m for one tree. No form may increase with height (a rising
  h <- lo$height_m; f <- lo$f; top <- f[which.max(h)]   # lower profile is carried upward as constant at the top value)
  lin <- coef(lm(f ~ h)); rising <- lin[2] > 0
  expf <- if (all(f > 0) && length(unique(h)) >= 2) coef(lm(log(f) ~ h)) else c(NA, NA)
  list(const_mean = function(x) rep(mean(f), length(x)),
       const_top  = function(x) rep(top, length(x)),
       exp_decay  = function(x) if (is.na(expf[1])) rep(NA_real_, length(x)) else if (expf[2] > 0) rep(top, length(x)) else exp(expf[1] + expf[2] * x),
       exp_pooled = function(x) if (is.na(B_POOL) || B_POOL > 0) rep(top, length(x)) else top * exp(B_POOL * (x - max(h))),
       linear     = function(x) if (rising) rep(top, length(x)) else lin[1] + lin[2] * x,
       linear_zero = function(x) if (rising) rep(top, length(x)) else pmax(lin[1] + lin[2] * x, 0),
       zero       = function(x) rep(0, length(x)))
}
P <- bind_rows(lapply(split(high, high$tree), function(hi) {
  fs <- forms(low %>% filter(tree == hi$tree[1]))
  bind_rows(lapply(names(fs), function(k) hi %>% mutate(form = k, pred = fs[[k]](height_m))))
}))
write_csv(P, file.path(OUT, "form_test_points.csv"))
# profile means: per-height means of observed and predicted, then unweighted and cone-weighted (weight = Htop - h) means
prof <- P %>% group_by(tree, form, height_m) %>% summarise(obs = mean(flux), pred = mean(pred), .groups = "drop") %>%
  group_by(tree) %>% mutate(Htop = max(height_m), w = pmax(Htop - height_m, 0) + 0.5) %>% ungroup()  # +0.5 m keeps the top point
sc <- P %>% group_by(tree, form) %>% summarise(n = n(), bias = mean(pred - flux), mae = mean(abs(pred - flux)),
  obs_neg_share = mean(flux < 0), pred_neg_share = mean(pred < 0), .groups = "drop") %>%
  left_join(prof %>% group_by(tree, form) %>% summarise(obs_profile = mean(obs), pred_profile = mean(pred),
    obs_cone = weighted.mean(obs, w), pred_cone = weighted.mean(pred, w), .groups = "drop"), by = c("tree", "form"))
write_csv(sc, file.path(OUT, "form_test_scores.csv"))
pool <- sc %>% group_by(form) %>% summarise(trees = sum(!is.na(pred_profile)), points = sum(n[!is.na(bias)]),
  median_bias = median(bias, na.rm = TRUE), median_mae = median(mae, na.rm = TRUE),
  ratio_pred_obs_profile = sum(pred_profile, na.rm = TRUE) / sum(obs_profile[!is.na(pred_profile)]),
  ratio_pred_obs_cone = sum(pred_cone, na.rm = TRUE) / sum(obs_cone[!is.na(pred_cone)]),
  trees_pred_net_uptake = sum(pred_cone < 0, na.rm = TRUE), trees_obs_net_uptake = sum(obs_cone[!is.na(pred_cone)] < 0), .groups = "drop") %>%
  arrange(abs(log(abs(ratio_pred_obs_cone) + 1e-9)))
write_csv(pool, file.path(OUT, "form_test_pooled.csv"))
cat("below-2 m means per tree:\n"); print(low %>% pivot_wider(names_from = height_m, values_from = f), width = 200)
cat("\npooled scores (ratio = predicted / observed mean flux over the measured >= 2 m profile):\n"); print(pool, width = 200)
cat("\nobserved >= 2 m: points", nrow(high), "negative", sum(high$flux < 0), "\n")
# figure: per tree, observed (points) and forms (lines) over height
grid <- bind_rows(lapply(unique(high$tree), function(t) { fs <- forms(low %>% filter(tree == t)); x <- seq(0.3, max(S$height_m[S$tree == t]), by = 0.1)
  bind_rows(lapply(names(fs), function(k) tibble(tree = t, form = k, height_m = x, pred = fs[[k]](x)))) }))
p <- ggplot() + geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) + geom_hline(yintercept = 2, linetype = "22", colour = "grey60", linewidth = 0.3) +
  geom_path(data = grid %>% filter(!is.na(pred)), aes(pred, height_m, colour = form), linewidth = 0.5) +
  geom_point(data = S, aes(flux, height_m, shape = height_m >= 2), size = 1.3) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), labels = c("used to fit (< 2 m)", "held out (>= 2 m)"), name = NULL) +
  scale_colour_manual(values = c(const_mean = "#E69F00", const_top = "#D55E00", exp_decay = "#009E73", exp_pooled = "#CC79A7", linear = "#0072B2", linear_zero = "#56B4E9", zero = "grey40"), name = "Extrapolation form") +
  facet_wrap(~tree, scales = "free_x", ncol = 4) + coord_cartesian(xlim = NULL) +
  labs(x = expression(Stem~CH[4]~flux~(nmol~m^{-2}~s^{-1})), y = "Height (m)") + theme_classic(base_size = 9) + theme(legend.position = "bottom")
ggsave(file.path(ROOT, "scaling/fig_form_test.png"), p, width = 220, height = 150, units = "mm", dpi = 250, bg = "white")

# ---- bootstrap: resample measurements within tree x height (below and above 2 m), recompute the pooled cone ratio
set.seed(42); NB <- 500
boot <- bind_rows(lapply(seq_len(NB), function(b) {
  Sb <- S %>% group_by(tree, height_m) %>% slice_sample(prop = 1, replace = TRUE) %>% ungroup()
  lo <- Sb %>% filter(height_m < 2) %>% group_by(tree, height_m) %>% summarise(f = mean(flux), .groups = "drop")
  hi <- Sb %>% filter(height_m >= 2)
  bind_rows(lapply(split(hi, hi$tree), function(x) { fs <- forms(lo %>% filter(tree == x$tree[1]))
    bind_rows(lapply(names(fs), function(k) x %>% mutate(form = k, pred = fs[[k]](height_m)))) })) %>%
    group_by(tree, form, height_m) %>% summarise(obs = mean(flux), pred = mean(pred), .groups = "drop") %>%
    group_by(tree) %>% mutate(w = pmax(max(height_m) - height_m, 0) + 0.5) %>% group_by(tree, form) %>%
    summarise(o = weighted.mean(obs, w), p = weighted.mean(pred, w), .groups = "drop") %>% filter(!is.na(p)) %>%
    group_by(form) %>% summarise(ratio = sum(p) / sum(o), uptake_trees = sum(p < 0), .groups = "drop") %>% mutate(b = b)
}))
bs <- boot %>% group_by(form) %>% summarise(ratio_median = median(ratio), ratio_lo = quantile(ratio, 0.025), ratio_hi = quantile(ratio, 0.975),
  uptake_trees_median = median(uptake_trees), .groups = "drop")
write_csv(bs, file.path(OUT, "form_test_bootstrap.csv")); cat("\nbootstrap (", NB, "resamples), cone-weighted ratio:\n"); print(bs)
