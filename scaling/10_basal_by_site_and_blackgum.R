# 10_basal_by_site_and_blackgum.R
# (1) Are Swamp Rd trees higher at the base than the EMS/YMF trees? Tree-level basal (< 2 m) means and their relation to
#     each tree's mean above 2 m and branch flux.
# (2) Black gum (saturated peat swamp, reference only, not in the stand scenarios): share of stem flux above 2 m and
#     absolute flux above 2 m vs the upland trees, under extrapolation forms above the top chamber (3.6 m).
#     Stem geometry: a cone (stem area density proportional to H - h), H = 15.8 m, the mean maximum lidar canopy height of N. sylvatica
#     in the HF ForestGEO plot (Sullivan et al. 2017; a proxy, not measured in Black Gum Swamp).
source("scaling/00_load_field.R")
source("scaling/analysis_helpers.R")
suppressPackageStartupMessages({library(readr); library(ggplot2)})
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
BGS <- tidyr::crossing(H = BLACKGUM_HEIGHT, form = c("exp", "const_top", "zero_above_top")) %>% rowwise() %>% mutate(as_tibble(as.list(scen(H, form)))) %>% ungroup() %>%
  mutate(ratio_vs_upland_upper = mean_flux_above2 / mean(up_upland), ratio_vs_all_upper = mean_flux_above2 / mean(up_all))
write_csv(BGS, file.path(ROOT, "scaling/out/blackgum_scaling.csv"))
cat(sprintf("\nblack gum fit a = %.1f, b = %.2f; top chamber (%.1f m) mean %.1f; upland trees' mean flux above 2 m %.3f (tree means: %s)\n",
            a, b, max(bg$height_m), top, mean(up_upland), paste(round(up_upland, 3), collapse = ", ")))
print(BGS %>% mutate(across(where(is.numeric), ~ signif(.x, 3))), n = 30)

## ---- global illustration: black gum-like stem flux above 2 m over the forested wetland classes of GLWD v2
# Forested classes 8,10,12,14,16,18,20,22,24,26,28 from Lehner et al. 2025 ESSD Table 3 (10^3 km2); forested = >= 10 % tree cover.
FW <- c(428.8, 378.6, 805.2, 701.2, 72.5, 138.5, 37.3, 1410.4, 431.7, 803.5, 150.8); A_fw <- sum(FW) / 1000   # M km2
f_below <- cone_share_below(2,BLACKGUM_HEIGHT)
SAI_up <- 0.45 * (1 - f_below)   # Whittaker & Woodwell 1967 stem area 0.45 m2 m-2 (temperate, closed canopy): illustrative
G <- BGS %>% mutate(Tg = mean_flux_above2 * SAI_up * A_fw * 0.506)
cat(sprintf("\nforested wetland area (GLWD v2) %.2f M km2; stem area above 2 m %.3f m2 m-2 ground (illustrative, full stocking)\n", A_fw, SAI_up))
cat(sprintf("black gum-like stem above 2 m: %.1f-%.1f Tg/yr (all H, forms); exp/zero-above-top only %.1f-%.1f; same area at the other trees' mean above 2 m (%.2f): %.2f Tg/yr\n",
  min(G$Tg), max(G$Tg), min(G$Tg[G$form != "const_top"]), max(G$Tg[G$form != "const_top"]), mean(up_all), mean(up_all) * SAI_up * A_fw * 0.506))
write_csv(G, file.path(ROOT, "scaling/out/blackgum_scaling.csv"))

## ---- add branches to the wetland scenario: branches were not measured on this tree, so the HF branch rates (tree-weighted mean and
## observation median, 03_tree_component.R) are used as conservative stand-ins; branch area 1.70 m2 m-2 (W&W 1967)
RT <- read_csv(file.path(ROOT, "scaling/out/stand_rates_HF.csv"), show_col_types = FALSE)
br_mean <- RT$mean[RT$comp == "branch"]; br_med <- median(F$flux[F$component == "branch" & !grepl("YMF", F$tree)])
for (bf in c(br_med, br_mean)) cat(sprintf("branches at %.3f nmol m-2 s-1 x 1.70 x %.2f M km2: %.2f Tg/yr\n", bf, A_fw, bf * 1.70 * A_fw * 0.506))
cat(sprintf("stem above 2 m + branches: %.1f-%.1f Tg/yr\n", min(G$Tg) + br_med * 1.70 * A_fw * 0.506, max(G$Tg) + br_mean * 1.70 * A_fw * 0.506))

## ---- SI component budget: same stem scenarios, with explicitly assigned branches.
# No branch or leaf flux was measured on this reference tree. This is a woody
# budget illustration; leaves are omitted, not assumed to have zero exchange.
BGCOMP <- tidyr::crossing(G, branch_rule = c("Upland median", "Upland tree mean")) %>%
  mutate(branch_rate = ifelse(branch_rule == "Upland median", br_med, br_mean),
    stem_lo = mean_flux_below2 * .45 * f_below,
    stem_up = mean_flux_above2 * .45 * (1-f_below), branch = branch_rate * 1.70) %>%
  select(form, branch_rule, stem_lo, stem_up, branch) %>%
  tidyr::pivot_longer(c(stem_lo,stem_up,branch), names_to="component", values_to="integrated") %>%
  group_by(form,branch_rule) %>% mutate(total=sum(integrated),share_pct=100*integrated/total) %>% ungroup()
write_csv(BGCOMP,file.path(ROOT,"scaling/out/blackgum_component_budget.csv"))
BGCOMP <- BGCOMP %>% mutate(form=factor(form,c("zero_above_top","exp","const_top"),
  c("Zero above\n3.6 m","Exponential\ndecay","Constant above\n3.6 m")),
  component=factor(component,c("branch","stem_up","stem_lo"),c("Branches (assigned)","Stem ≥ 2 m","Stem < 2 m")))
cols<-c("Stem < 2 m"="#8B4513","Stem ≥ 2 m"="#D4A76A","Branches (assigned)"="#4682B4")
p <- ggplot(BGCOMP,aes(form,share_pct,fill=component))+
  geom_col(width=.62,colour="white",linewidth=.3)+
  geom_text(aes(label=ifelse(component=="Branches (assigned)","",sprintf("%.1f%%",share_pct))),position=position_stack(vjust=.5),
    size=3,colour="white",fontface="bold")+
  geom_text(data=filter(BGCOMP,component=="Branches (assigned)"),
    aes(y=103,label=sprintf("%.1f%%",share_pct)),colour="#31648C",size=3,fontface="bold")+
  facet_wrap(~branch_rule,nrow=1)+scale_fill_manual(values=cols,breaks=names(cols),name=NULL)+
  scale_y_continuous(limits=c(0,106),breaks=seq(0,100,25),expand=c(0,0))+
  labs(x="Assumed flux above the highest stem chamber",y="Share of woody-surface methane flux (%)")+
  theme_classic(base_size=10,base_family="Helvetica")+
  theme(legend.position="bottom",strip.background=element_blank(),strip.text=element_text(face="bold"),
        axis.text.x=element_text(size=9),legend.text=element_text(size=9))
ggsave(file.path(ROOT,"scaling/v3c/FigS4_blackgum_budget.png"),p,width=7.2,height=4.7,dpi=320,bg="white")
ggsave(file.path(ROOT,"scaling/v3c/FigS4_blackgum_budget.pdf"),p,width=7.2,height=4.7,device=cairo_pdf)
