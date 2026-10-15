# Regression checks for physical units, area integration, cohort integrity and outputs.
# Run from the repository root after run_analysis.R.
suppressPackageStartupMessages({library(lubridate); library(dplyr)})
source("scaling/analysis_helpers.R")
source("flux/precision_helpers.R")
source("scaling/00_load_field.R")
# Summer station clock 12:00 EST is 13:00 EDT, relabelled 13:00 for the naive analyzer clock.
t<-force_tz(with_tz(as.POSIXct("2023-07-18 12:00:00",tz="EST"),"America/New_York"),"UTC")
stopifnot(format(t,"%H:%M",tz="UTC")=="13:00")
# Analytic cylinder and linear flux, with intentionally unequal sampling intervals.
g<-data.frame(h=c(0,2,10),d=c(.4,.4,.4))
q<-stem_grid(g,2,10,c(2,3,10)); q2<-stem_grid(g,2,10,c(2,3,5,10))
stopifnot(abs(sum(q$area)-pi*.4*8)<1e-10,
 abs(sum(profile_weights(c(2,3,10),q)*c(2,3,10))-6)<1e-10,
 abs(sum(profile_weights(c(2,3,5,10),q2)*c(2,3,5,10))-6)<1e-10)
# Frustum lateral area has a closed form, including the slant correction.
g<-data.frame(h=c(0,10),d=c(.4,.1));q<-stem_grid(g,0,10)
stopifnot(abs(sum(q$area)-pi*(.2+.05)*sqrt(10^2+.15^2))<1e-10)
# Whole-stand cone partition matches the capture diagram, independently of base diameter.
g<-data.frame(h=c(0,CANOPY_HEIGHT),d=c(.4,0))
stopifnot(abs(sum(stem_grid(g,0,2)$area)/sum(stem_grid(g,0,CANOPY_HEIGHT)$area)-cone_share_below(2,CANOPY_HEIGHT))<1e-10)
# A rising profile is held constant; a declining zero-floor form never predicts uptake.
fs<-extrapolation_forms(data.frame(height_m=c(.5,1.25),flux=c(.1,.2)))
stopifnot(all(fs$linear(c(2,10))==.2),all(fs$exp_decay(c(2,10))==.2))
fs<-extrapolation_forms(data.frame(height_m=c(.5,1.25),flux=c(1,.5)))
stopifnot(fs$linear(10)<0,fs$linear_zero(10)==0)
# No bootstrap draw silently changes the cohort denominator.
p<-read.csv("scaling/out/form_test_pooled.csv");b<-read.csv("scaling/out/form_test_bootstrap.csv")
cov<-read.csv("scaling/out/form_test_coverage.csv")
stopifnot(all(p$trees==sum(cov$exp_eligible)),all(b$trees==sum(cov$exp_eligible)),all(b$valid_draws<=b$total_draws))
tc<-read.csv("scaling/out/tree_component.csv")
stopifnot(length(unique(tc$cohort))==1,all(tc$ntree==5),all(is.finite(tc$woody_component)))
# Verify the metadata and thresholds actually used by goFlux (requires the flux cache from a flux run).
if (all(file.exists(file.path("flux/cache", c("harvard_forest.rds", "yale_myers.rds", "black_gum_swamp.rds"))))) {
 cached <- c(readRDS("flux/cache/harvard_forest.rds"), readRDS("flux/cache/yale_myers.rds"))
 for (inst in c("LGR1","LGR2","LGR3","YMF")) {
  m <- cached[[inst]]$traces
  o <- read.csv(file.path("flux/out", if (inst == "YMF") "yale_myers_fluxes.csv" else "harvard_forest_fluxes.csv"))
  stopifnot(all(m$H2O_prec == 200), all(m$flag == 1))
  p <- estimate_precision(m, inst, warn = FALSE)
  for (gas in c("CH4","CO2")) {
   pc <- p$closures[p$closures$gas == gas, ]
   stopifnot(max(abs(m[[paste0(gas,"_prec")]] - qnorm(.975)*pc$sigma_group[match(m$UniqueID, pc$UniqueID)])) < 1e-10)
  }
  # Windows in the traces are exactly the metadata windows.
  w <- m %>% group_by(UniqueID) %>% summarise(s = min(POSIX.time), e = max(POSIX.time), .groups = "drop")
  j <- match(w$UniqueID, o$UniqueID)
  stopifnot(all(w$s >= as.POSIXct(o$Window_start[j], tz = "UTC")), all(w$e <= as.POSIXct(o$Window_end[j], tz = "UTC")))
  cl <- m %>% group_by(UniqueID) %>% summarise(t = diff(range(as.numeric(POSIX.time))) + median(diff(as.numeric(POSIX.time))), .groups = "drop")
  stopifnot(max(abs(cl$t - o$t_sec[match(cl$UniqueID, o$UniqueID)])) < 1e-8)
  stopifnot(max(abs(o$CH4_MDF - o$CH4_MDF_wass95)) < 1e-8, max(abs(o$CO2_MDF - o$CO2_MDF_wass95)) < 1e-8)
 }
 bg <- readRDS("flux/cache/black_gum_swamp.rds")$BG
 o <- read.csv("flux/out/black_gum_swamp_fluxes.csv")
 stopifnot(max(abs(bg$CH4$MDF - o$CH4_MDF_emp95[match(bg$CH4$UniqueID, o$UniqueID)])) < 1e-8)
 cl <- bg$traces %>% group_by(UniqueID) %>% summarise(t = diff(range(as.numeric(POSIX.time))) + median(diff(as.numeric(POSIX.time))), .groups = "drop")
 stopifnot(max(abs(cl$t - o$t_sec[match(cl$UniqueID, o$UniqueID)])) < 1e-8)
} else message("Flux cache absent: fit-level checks run after a flux rebuild; committed tables checked.")
o <- read.csv("flux/out/black_gum_swamp_fluxes.csv")
stopifnot(all(is.finite(o$t_sec)), all(o$t_sec > 30))
# Basic sampling and stand accounting invariants.
stopifnot(nrow(F)==141,n_distinct(F$tree)==7,sum(F$component=="branch")==11)
r<-read.csv("scaling/out/stand_rates_HF.csv");u<-read.csv("scaling/out/stand_uncertainty_HF.csv")
stopifnot(abs(sum(r$mean*r$area)-u$estimate[u$metric=="total"])<1e-10,all(u$lo<=u$hi))
stopifnot(abs(r$area[r$comp=="stem_lo"]/STAND_AREA[["stem"]]-cone_share_below(2,CANOPY_HEIGHT))<1e-10)
cat("PASS: timestamp, units, geometry, extrapolation, fixed cohorts, metadata, durations and stand accounting\n")

trade<-read.csv("scaling/out/rate_area_tradeoff.csv")
stopifnot(diff(range(trade$rate_matching_basal*trade$area))<1e-10)
budget<-read.csv("scaling/out/blackgum_component_budget.csv")
stopifnot(all(abs(tapply(budget$share_pct,interaction(budget$form,budget$branch_rule),sum)-100)<1e-10))

# Area assumptions propagate consistently to stand, sensitivity, capture and swamp outputs.
stopifnot(abs(sum(r$area[r$comp %in% c("stem_lo","stem_up")])-STAND_AREA[["stem"]])<1e-12,
  abs(r$area[r$comp=="branch"]-STAND_AREA[["branch"]])<1e-12,
  all(abs(tc$stem-STAND_AREA[["stem"]])<1e-12))
alt<-grepl("Gauci",tc$area_set)
stopifnot(all(abs(tc$stem[alt]+tc$branch[alt]-ALTERNATIVE_WOODY_AREA)<1e-12),
  all(abs(tc$branch[!alt]-STAND_AREA[["branch"]])<1e-12))
cap<-read.csv("scaling/out/surface_capture.csv")
for (nm in c("Stem","Stem + branches","All surfaces")) {
  z<-cap[cap$name==nm,]
  denom<-switch(nm,"Stem"=STAND_AREA[["stem"]],"Stem + branches"=sum(STAND_AREA[c("stem","branch")]),"All surfaces"=sum(STAND_AREA))
  stopifnot(max(abs(z$value-100*cone_share_below(z$h,CANOPY_HEIGHT)*STAND_AREA[["stem"]]/denom))<1e-10)
}
bgs<-read.csv("scaling/out/blackgum_scaling.csv")
for (fm in unique(budget$form)) {
  z<-budget[budget$form==fm,]; v<-bgs[bgs$form==fm,]
  expected<-v$mean_flux_above2*STAND_AREA[["stem"]]*(1-cone_share_below(2,BLACKGUM_HEIGHT))
  stopifnot(all(abs(z$integrated[z$component=="stem_up"]-expected)<1e-10))
}
cat("PASS: shared area assumptions across stand, sensitivity, capture and swamp scenarios\n")

source("tests/leaf_area_checks.R")
