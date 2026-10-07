# 10_basal_by_site_and_blackgum.R
# (1) Are Swamp Rd trees higher at the base than the EMS/YMF trees? Tree-level basal (< 2 m) means and their relation to
#     each tree's mean above 2 m and branch flux.
# (2) Black gum (saturated peat swamp, reference only, not in the stand scenarios): share of stem flux above 2 m and
#     absolute flux above 2 m vs the upland trees, under extrapolation forms above the top chamber (3.6 m).
#     Stem geometry: a cone (stem area density proportional to H - h), H = 15.8 m, the mean lidar canopy height of N. sylvatica
#     in the HF ForestGEO plot (Sullivan et al. 2017; a proxy, not measured in Black Gum Swamp).
source("scaling/00_load_field.R")
suppressPackageStartupMessages(library(readr))
TM <- F %>% filter(component == "stem") %>% mutate(z = ifelse(height_m < 2, "basal", "upper")) %>% group_by(tree, z, height_m) %>%
  summarise(f = mean(flux), .groups = "drop") %>% group_by(tree, z) %>% summarise(f = mean(f), .groups = "drop") %>% tidyr::pivot_wider(names_from = z, values_from = f) %>%
  left_join(F %>% filter(component == "branch") %>% group_by(tree) %>% summarise(branch = mean(flux), nbranch = n()), by = "tree") %>%
  mutate(site = sub(" .*", "", tree)) %>% arrange(site, -basal)
print(TM)
cat("\nbasal mean by site:\n"); print(TM %>% group_by(site) %>% summarise(n = n(), basal_mean = mean(basal), basal_min = min(basal), basal_max = max(basal), upper_mean = mean(upper)))
cat("\nSwamp Rd vs other trees, basal (Wilcoxon):", format.pval(wilcox.test(basal ~ site == "Swamp", data = TM)$p.value, 2), "\n")
ct <- cor.test(TM$basal, TM$upper, method = "spearman"); cat("basal vs upper mean across trees: Spearman rho", round(ct$estimate, 2), "p", round(ct$p.value, 2), "n", nrow(TM), "\n")
write_csv(TM, file.path(ROOT, "scaling/out/tree_basal_upper_branch.csv"))

## ---- black gum
bg <- read.csv(file.path(ROOT, "data processing/goFlux_reprocessing/diurnal_blackgum/results/blackgum_flux_compiled_with_mdf.csv"))
fit <- nls(CH4_best.flux ~ a * exp(-b * height_m), data = bg, start = list(a = 200, b = 1)); a <- coef(fit)[1]; b <- coef(fit)[2]
top <- bg %>% filter(height_m == max(height_m)) %>% summarise(f = mean(CH4_best.flux)) %>% pull(f)
prof <- bg %>% group_by(height_m) %>% summarise(f = mean(CH4_best.flux), n = n()); print(prof)
up_upland <- TM %>% filter(site != "Swamp") %>% pull(upper); up_all <- TM$upper
scen <- function(H, form) { h <- seq(0, H, 0.01); w <- (H - h)
  f <- switch(form, exp = a * exp(-b * h), const_top = ifelse(h <= max(bg$height_m), a * exp(-b * h), top), zero_above_top = ifelse(h <= max(bg$height_m), a * exp(-b * h), 0))
  c(share_above2 = sum((w * f)[h >= 2]) / sum(w * f), mean_flux_above2 = sum((w * f)[h >= 2]) / sum(w[h >= 2]), mean_flux_below2 = sum((w * f)[h < 2]) / sum(w[h < 2])) }
BGS <- tidyr::crossing(H = 15.8, form = c("exp", "const_top", "zero_above_top")) %>% rowwise() %>% mutate(as_tibble(as.list(scen(H, form)))) %>% ungroup() %>%
  mutate(ratio_vs_upland_upper = mean_flux_above2 / mean(up_upland), ratio_vs_all_upper = mean_flux_above2 / mean(up_all))
write_csv(BGS, file.path(ROOT, "scaling/out/blackgum_scaling.csv"))
cat(sprintf("\nblack gum fit a = %.1f, b = %.2f; top chamber (%.1f m) mean %.1f; upland trees' mean flux above 2 m %.3f (tree means: %s)\n",
            a, b, max(bg$height_m), top, mean(up_upland), paste(round(up_upland, 3), collapse = ", ")))
print(BGS %>% mutate(across(where(is.numeric), ~ signif(.x, 3))), n = 30)

## ---- global illustration: black gum-like stem flux above 2 m over the forested wetland classes of GLWD v2
# Forested classes 8,10,12,14,16,18,20,22,24,26,28 from Lehner et al. 2025 ESSD Table 3 (10^3 km2); forested = >= 10 % tree cover.
FW <- c(428.8, 378.6, 805.2, 701.2, 72.5, 138.5, 37.3, 1410.4, 431.7, 803.5, 150.8); A_fw <- sum(FW) / 1000   # M km2
f_below <- 0.188   # mean share of stem area below 2 m on our measured trees (03_tree_component.R)
SAI_up <- 0.45 * (1 - f_below)   # Whittaker & Woodwell 1967 stem area 0.45 m2 m-2 (temperate, closed canopy): illustrative
G <- BGS %>% mutate(Tg = mean_flux_above2 * SAI_up * A_fw * 0.506)
cat(sprintf("\nforested wetland area (GLWD v2) %.2f M km2; stem area above 2 m %.3f m2 m-2 ground (illustrative, full stocking)\n", A_fw, SAI_up))
cat(sprintf("black gum-like stem above 2 m: %.1f-%.1f Tg/yr (all H, forms); exp/zero-above-top only %.1f-%.1f; same area at the other trees' mean above 2 m (%.2f): %.2f Tg/yr\n",
  min(G$Tg), max(G$Tg), min(G$Tg[G$form != "const_top"]), max(G$Tg[G$form != "const_top"]), mean(up_all), mean(up_all) * SAI_up * A_fw * 0.506))
write_csv(G, file.path(ROOT, "scaling/out/blackgum_scaling.csv"))

## ---- add branches to the wetland scenario: no wetland branch data, so the HF branch rates (tree-weighted mean and
## observation median, 03_tree_component.R) are used as conservative stand-ins; branch area 1.70 m2 m-2 (W&W 1967)
RT <- read_csv(file.path(ROOT, "scaling/out/stand_rates_HF.csv"), show_col_types = FALSE)
br_mean <- RT$mean[RT$comp == "branch"]; br_med <- median(F$flux[F$component == "branch" & !grepl("YMF", F$tree)])
for (bf in c(br_med, br_mean)) cat(sprintf("branches at %.3f nmol m-2 s-1 x 1.70 x %.2f M km2: %.2f Tg/yr\n", bf, A_fw, bf * 1.70 * A_fw * 0.506))
cat(sprintf("stem above 2 m + branches: %.1f-%.1f Tg/yr\n", min(G$Tg) + br_med * 1.70 * A_fw * 0.506, max(G$Tg) + br_mean * 1.70 * A_fw * 0.506))
