# =============================================================================
# 03_tree_component.R — GRL: the TREE COMPONENT of the forest CH4 exchange (not an ecosystem budget), per m2 ground,
# compared with same-month soil uptake at Harvard Forest; leaves treated separately; microsite table; black gum
# (saturated peat swamp) as a wetland reference profile; mixed-model extrapolation from basal data.
# Inputs: 00_load_field.R (field fluxes), out/form_test_scores.csv (01), black gum compiled fluxes, Jevon 2023 soil data.
# Sampled stem flux means use measured diameters and physical height intervals.
# Stand area is partitioned with the same Sullivan-height cone as the capture figure.
# Applying sampled compartment rates to that area is the explicit stand scenario.
# Woody area per m2 ground: Whittaker & Woodwell 1967 (stem 0.55, branch 1.55; LAI 4.5) or total woody area 3.07
#   (Gauci 2024 temperate WAI; branch = 3.07 - stem). These are scenario inputs.
# One measured branch rate is applied to branch area; its vertical distribution is not estimated.
# Soil: Jevon et al. 2023 (Prospect Hill, upland), July-August instantaneous fluxes 2016-17 (fluxes.csv, umol m-2 s-1).
# Outputs: scaling/out/tree_component.csv, microsite.csv, mixed_model_pred.csv; figures fig_tree_component.png,
#   fig_SI_extrapolation_fits.png
# Run: LANG=en_US.UTF-8 Rscript scaling/03_tree_component.R   (from whole_tree_flux/)
# =============================================================================
source("scaling/00_load_field.R")
source("scaling/analysis_helpers.R")
suppressPackageStartupMessages({library(tidyr); library(ggplot2); library(readr); library(patchwork); library(lme4)})
OUT <- file.path(ROOT, "scaling/out")
num <- function(x) suppressWarnings(as.numeric(x))

## ---- microsite table (field sheet soil moisture and temperature)
hf <- read.csv(file.path(DATA, "data processing/goFlux_reprocessing/results/canopy_flux_goFlux_compiled_with_mdf.csv"), stringsAsFactors = FALSE) %>%
  filter(!(Species == "bg" & Tree_Tag == 3))
micro <- hf %>% mutate(tree = paste(Site, Species, Tree_Tag)) %>% group_by(tree, Site, Species) %>%
  summarise(date = paste(unique(Date), collapse = ", "), soil_T = paste(unique(Soil_Temp_C), collapse = ", "),
            VWC = paste(unique(paste(VWC_1, VWC_2, VWC_3, sep = "/")), collapse = "; "), .groups = "drop") %>%
  mutate(microsite = case_when(grepl("sat", VWC) ~ "saturated soil (wetland-like)",
                               Site == "Swamp Rd" ~ "wetland margin (Swamp Rd)",
                               TRUE ~ "upland"))
micro <- bind_rows(micro, tibble(tree = "YMF black oak", Site = "Yale Myers", Species = "qv", date = "10-04-2022", soil_T = NA, VWC = NA, microsite = "upland"),
                   tibble(tree = "Black Gum Swamp bg (2024)", Site = "Black Gum Swamp", Species = "bg", date = "08-28/29-2024", soil_T = NA, VWC = "saturated peat (site)", microsite = "wetland (saturated peat swamp)"))
write_csv(micro, file.path(OUT, "microsite.csv")); print(micro %>% select(tree, VWC, microsite), width = 200)

## ---- measured geometry and integrated upper-stem profiles
G <- load_geometry(F)
HFT <- unique(F$tree[F$site != "Yale Myers"])
FH <- F %>% filter(tree %in% HFT)
sc <- read_csv(file.path(OUT, "form_test_scores.csv"), show_col_types=FALSE) %>% filter(tree %in% HFT)
coverage <- read_csv(file.path(OUT, "form_test_coverage.csv"), show_col_types=FALSE)
common_HF <- intersect(HFT, coverage$tree[coverage$exp_eligible])
SA <- bind_rows(lapply(split(filter(G,tree %in% HFT), filter(G,tree %in% HFT)$tree),function(g) {
  tibble(tree=g$tree[1],below2=sum(stem_grid(g,0,2)$area),above2=sum(stem_grid(g,2,g$top[1])$area),top=g$top[1])
})) %>% mutate(share_below2=below2/(below2+above2),
  scenario_height=CANOPY_HEIGHT,scenario_share_below2=cone_share_below(2,CANOPY_HEIGHT))
write_csv(SA,file.path(OUT,"stem_geometry.csv"))
component_rows <- function(d) {
  bind_rows(lapply(split(d,d$tree),function(x) {
    original <- if("source_tree" %in% names(x)) x$source_tree[1] else x$tree[1]
    g <- filter(G,tree==original)
    rate <- function(component,lower=NA,upper=NA) {
      y <- x[x$component==component,,drop=FALSE]
      if(is.finite(lower)) y<-y[y$height_m>=lower & y$height_m<upper,,drop=FALSE]
      if(!nrow(y)) return(NA_real_)
      h<-height_means(y)
      if(component!="stem") return(mean(h$flux))
      gr<-stem_grid(g,lower,min(upper,g$top[1]),h$height_m)
      sum(profile_weights(h$height_m,gr)*h$flux)
    }
    tibble(tree=x$tree[1],source_tree=original,stem_lo=rate("stem",0,2),
      stem_up=rate("stem",2,Inf),branch=rate("branch"),leaf=rate("leaf"),
      share_below2=cone_share_below(2,CANOPY_HEIGHT))
  }))
}
CR <- component_rows(FH)
write_csv(CR,file.path(OUT,"tree_component_rates.csv"))
comp_names<-c("stem_lo","stem_up","branch","leaf")
summarize_stand <- function(cr) {
  f<-vapply(comp_names,function(k)mean(cr[[k]],na.rm=TRUE),numeric(1))
  below<-mean(cr$share_below2)
  area<-unname(c(STAND_AREA["stem"]*below,STAND_AREA["stem"]*(1-below),STAND_AREA["branch"],STAND_AREA["leaf"]))
  val<-f*area
  c(setNames(f,paste0("rate_",comp_names)),setNames(area,paste0("area_",comp_names)),
    setNames(val,paste0("int_",comp_names)),setNames(100*val/sum(val),paste0("pct_",comp_names)),
    total=sum(val),woody=sum(val[1:3]),cancel=unname(-val[1]/sum(area[-1])),flip=unname(-2*val[1]/sum(area[-1])),
    upper_mean=sum(val[-1])/sum(area[-1]))
}
base<-summarize_stand(CR)
set.seed(43)
stand_boot<-t(replicate(N_BOOT,summarize_stand(component_rows(resample_trees(FH,HFT)))))
# Some draws select no branch-measured tree. Retain NA for those terms rather than
# treating unavailable branches as zero; report completeness for every interval.
stand_ci<-bind_rows(lapply(names(base),function(k)tibble(metric=k,estimate=base[k],lo=ci95(stand_boot[,k])[1],
  hi=ci95(stand_boot[,k])[2],valid_draws=sum(is.finite(stand_boot[,k])),total_draws=N_BOOT)))
write_csv(stand_ci,file.path(OUT,"stand_uncertainty_HF.csv"))
RATES<-tibble(comp=comp_names,mean=unname(base[paste0("rate_",comp_names)]),
  se=vapply(comp_names,function(k)sd(CR[[k]],na.rm=TRUE)/sqrt(sum(is.finite(CR[[k]]))),numeric(1)),
  ntree=vapply(comp_names,function(k)sum(is.finite(CR[[k]])),integer(1)),
  area=unname(base[paste0("area_",comp_names)]),
  lo=stand_ci$lo[match(paste0("rate_",comp_names),stand_ci$metric)],
  hi=stand_ci$hi[match(paste0("rate_",comp_names),stand_ci$metric)])
write_csv(RATES,file.path(OUT,"stand_rates_HF.csv")); print(RATES)
# All cross-form scenarios use the same five eligible Harvard Forest trees,
# including the measured benchmark, basal contribution and area partition.
scc<-filter(sc,tree %in% common_HF)
upper<-bind_rows(scc %>% group_by(form) %>% summarise(f_up=mean(pred_area),.groups="drop"),
  scc %>% filter(form=="const_mean") %>% summarise(form="measured",f_up=mean(obs_area)))
stopifnot(all(is.finite(upper$f_up)))
CRc<-filter(CR,source_tree %in% common_HF)
stem_lo<-mean(CRc$stem_lo); f_below<-mean(CRc$share_below2)
br<-FH %>% filter(component=="branch",tree %in% common_HF)
lf<-FH %>% filter(component=="leaf")
tw<-tree_mean
BR<-c(measured_mean=mean(CRc$branch,na.rm=TRUE),measured_median=median(br$flux))
AREA<-tribble(~area_set,~stem,~branch,"Whittaker & Woodwell 1967",STAND_AREA[["stem"]],STAND_AREA[["branch"]],"Woody area index 3.07 (Gauci 2024 temperate)",STAND_AREA[["stem"]],ALTERNATIVE_WOODY_AREA-STAND_AREA[["stem"]])
T<-crossing(upper,AREA,branch_rule=c("measured_mean","measured_median","equal_upper_stem")) %>%
  mutate(f_branch=ifelse(branch_rule=="equal_upper_stem",f_up,BR[branch_rule]),
    stem_lo_g=stem_lo*stem*f_below,stem_up_g=f_up*stem*(1-f_below),branch_g=f_branch*branch,
    woody_component=stem_lo_g+stem_up_g+branch_g,ntree=length(common_HF),cohort=paste(sort(common_HF),collapse="; "))
## ---- soil, same months (July-August), Jevon et al. 2023
soil <- read.csv(file.path(ROOT, "scaling/soil_jevon2023/fluxes.csv")) %>% mutate(m = format(as.Date(date, "%m/%d/%y"), "%m"), f = CH4.flux * 1000) %>% filter(m %in% c("07", "08"))
S_mean <- mean(soil$f); S_q <- quantile(soil$f, c(0.25, 0.75))
T <- T %>% mutate(soil_JulAug = S_mean, offset_pct_of_soil_sink = 100 * woody_component / abs(S_mean))
write_csv(T, file.path(OUT, "tree_component.csv"))
cat(sprintf("\nsoil CH4, Jul-Aug (Jevon 2023, n = %d): mean %.2f (IQR %.2f to %.2f) nmol m-2 s-1\n", nrow(soil), S_mean, S_q[1], S_q[2]))
print(T %>% group_by(form) %>% summarise(woody_min = min(woody_component), woody_max = max(woody_component),
  offset_min = min(offset_pct_of_soil_sink), offset_max = max(offset_pct_of_soil_sink), .groups = "drop"))
## ---- leaves, separately
LV <- tibble(leaf_rule = c("measured mean", "measured median", "zero", "uptake at -1x mean leaf rate"), f = c(mean(lf$flux), median(lf$flux), 0, -mean(lf$flux))) %>%
  mutate(leaf_g = f * STAND_AREA[["leaf"]])
cat("\nleaf term (LAI 4.5), nmol m-2 ground s-1; leaf fluxes below MDF (see summary_statistics.csv):\n"); print(LV)

## ---- mixed model on basal (< 2 m) stem data, extrapolated (asinh scale; tree random intercept and slope)
S <- F %>% filter(component == "stem") %>% mutate(y = asinh(flux / 0.01))
m <- lmer(y ~ height_m + (height_m | tree), data = S %>% filter(height_m < 2), REML = TRUE, control = lmerControl())
fe <- fixef(m); cat("\nmixed model (basal data): fixed slope on asinh scale =", round(fe[2], 3), "per m\n")
grid <- S %>% distinct(tree) %>% crossing(height_m = seq(0.3, 22, by = 0.1)) %>% left_join(S %>% group_by(tree) %>% summarise(hmax = max(height_m)), by = "tree") %>%
  filter(height_m <= hmax)
grid$pred <- 0.01 * sinh(predict(m, newdata = grid, re.form = NULL, allow.new.levels = TRUE))
write_csv(grid, file.path(OUT, "mixed_model_pred.csv"))

## ---- black gum reference (saturated peat swamp)
bg <- read.csv(file.path(ROOT, "data processing/goFlux_reprocessing/diurnal_blackgum/results/blackgum_flux_compiled_with_mdf.csv")) %>%
  transmute(tree = "Black Gum Swamp bg (2024)", height_m, flux = CH4_best.flux, below_mdf = CH4_below_MDF)
bgfit <- nls(flux ~ a * exp(-b * height_m), data = bg, start = list(a = 200, b = 1))
cat("black gum: exponential fit a =", round(coef(bgfit)[1], 1), "b =", round(coef(bgfit)[2], 3), "per m (e-folding", round(1 / coef(bgfit)[2], 2), "m); all", nrow(bg), "fluxes > 0:", all(bg$flux > 0), "\n")

## ---- figures
th <- theme_classic(base_size = 9, base_family = "Helvetica") + theme(strip.background = element_blank(), strip.text = element_text(face = "bold", hjust = 0), legend.position = "bottom")
lab_form <- c(const_mean = "Constant (basal mean)", const_top = "Constant (top chamber)", exp_decay = "Exponential decay", linear = "Linear (unbounded)",
              linear_zero = "Linear, floored at zero", zero = "Zero above 2 m", exp_pooled = "Pooled decline")
fpts <- read_csv(file.path(OUT, "form_test_points.csv"), show_col_types = FALSE)
lines <- bind_rows(lapply(unique(fpts$tree), function(t) {
  lo <- F %>% filter(component == "stem", tree == t, height_m < 2) %>% group_by(height_m) %>% summarise(f = mean(flux), .groups = "drop")
  h <- lo$height_m; f <- lo$f; top <- f[which.max(h)]; lin <- coef(lm(f ~ h)); rising <- lin[2] > 0
  ex <- if (all(f > 0)) coef(lm(log(f) ~ h)) else c(NA, NA); x <- seq(2, max(F$height_m[F$tree == t & F$component == "stem"]), by = 0.1)
  bind_rows(tibble(form = "const_mean", y = mean(f)), tibble(form = "exp_decay", y = if (is.na(ex[1])) NA else if (ex[2] > 0) top else NA),
            tibble(form = "zero", y = 0)) %>% select(-y) %>% {NULL}
  bind_rows(tibble(tree = t, height_m = x, form = "const_mean", pred = mean(f)),
            tibble(tree = t, height_m = x, form = "exp_decay", pred = if (is.na(ex[1])) NA else if (ex[2] > 0) top else exp(ex[1] + ex[2] * x)),
            tibble(tree = t, height_m = x, form = "linear", pred = if (rising) top else lin[1] + lin[2] * x),
            tibble(tree = t, height_m = x, form = "zero", pred = 0)) }))
mic <- micro %>% select(tree, microsite)
Sx <- F %>% filter(component == "stem") %>% left_join(mic, by = "tree") %>% mutate(panel = paste0(tree, "\n", microsite),
  pt = factor(ifelse(below_mdf, "below detection", ifelse(height_m < 2, "basal (< 2 m), used to fit", "held out (>= 2 m)")),
              c("basal (< 2 m), used to fit", "held out (>= 2 m)", "below detection")))
bgx <- bg %>% mutate(panel = "Black Gum Swamp bg (2024)\nwetland reference (saturated peat)",
  pt = factor(ifelse(below_mdf, "below detection", "wetland reference"), c("basal (< 2 m), used to fit", "held out (>= 2 m)", "below detection", "wetland reference")))
bgline <- tibble(height_m = seq(0.2, 3.6, 0.05)) %>% mutate(pred = predict(bgfit, newdata = data.frame(height_m = height_m)), panel = bgx$panel[1])
Sx$pt <- factor(Sx$pt, levels(bgx$pt))
PLEV <- c(sort(unique(Sx$panel)), bgx$panel[1])
lines <- lines %>% left_join(mic, by = "tree") %>% mutate(panel = paste0(tree, "\n", microsite))
mm <- grid %>% left_join(mic, by = "tree") %>% mutate(panel = paste0(tree, "\n", microsite))
tr <- scales::pseudo_log_trans(sigma = 0.02)
fix <- function(d) { d$panel <- factor(d$panel, PLEV); d }
Sx <- fix(Sx); bgx <- fix(bgx); bgline <- fix(bgline); lines <- fix(lines); mm <- fix(mm)
pS <- ggplot() + annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 2, fill = "#F2F2F2") +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_path(data = mm, aes(pred, height_m, linetype = "Mixed model (all trees' basal data)"), colour = "black", linewidth = 0.6) +
  geom_path(data = bgline, aes(pred, height_m, linetype = "Exponential fit (all heights)"), colour = "#7F3C8D", linewidth = 0.6) +
  geom_path(data = lines %>% filter(!is.na(pred)), aes(pred, height_m, colour = form), linewidth = 0.55, alpha = 0.9) +
  geom_point(data = bind_rows(Sx, bgx), aes(flux, height_m, shape = pt), size = 1.6, stroke = 0.5) +
  scale_shape_manual(values = c("basal (< 2 m), used to fit" = 21, "held out (>= 2 m)" = 16, "below detection" = 4, "wetland reference" = 17), name = NULL, drop = FALSE) +
  scale_linetype_manual(values = c("Mixed model (all trees' basal data)" = "solid", "Exponential fit (all heights)" = "22"), name = NULL) +
  scale_colour_manual(values = c(const_mean = "#E69F00", exp_decay = "#009E73", linear = "#0072B2", zero = "grey45"), labels = lab_form, name = "Per-tree extrapolation:") +
  scale_x_continuous(trans = tr, breaks = c(-1, 0, 0.1, 1, 10, 100), labels = function(x) sub("-", "−", format(x, drop0trailing = TRUE, trim = TRUE))) +
  facet_wrap(~panel, ncol = 4, scales = "free") +
  guides(colour = guide_legend(nrow = 2, order = 1), shape = guide_legend(nrow = 2, order = 2), linetype = guide_legend(nrow = 2, order = 3)) +
  labs(x = expression(Stem~CH[4]~flux~(nmol~m^{-2}~s^{-1})*","~signed~log~scale), y = "Height (m)",
       caption = "Grey band: below 2 m (the usual sampling range). Coloured lines: each tree's basal (< 2 m) data extrapolated upward (forms may not increase with height).\nBlack line: mixed model of all trees' basal data (tree-level intercept and slope) extrapolated upward. Last panel: a saturated peat-swamp black gum measured to 3.6 m (reference only).") +
  th + theme(legend.box = "vertical", plot.caption = element_text(size = 7, hjust = 0))
ggsave(file.path(ROOT, "scaling/fig_SI_extrapolation_fits.png"), pS, width = 230, height = 190, units = "mm", dpi = 300, bg = "white")
