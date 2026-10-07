# Regression checks for physical units, area integration, cohort integrity and outputs.
# Run from the repository root after run_analysis.R.
suppressPackageStartupMessages({library(lubridate); library(dplyr)})
source("scaling/analysis_helpers.R")
source("data processing/goFlux_reprocessing/precision_helpers.R")
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
# Verify metadata actually used by goFlux, rather than merely the setup constants.
rd<-"data processing/goFlux_reprocessing"
for(inst in c("LGR1","LGR2","LGR3","YMF")) {
 d<-if(inst=="YMF")file.path(rd,"ymf_black_oak") else rd
 e<-new.env();load(file.path(d,"RData",paste0("manID_",inst,".RData")),e)
 m<-e[[paste0("manID.",inst)]];a<-read.csv(file.path(d,"results",paste0("aux_",inst,".csv")))
 stopifnot(all(m$H2O_prec==200))
 p<-estimate_precision(m,inst,warn=FALSE)
 for(gas in c("CH4","CO2")) {
  pc<-p$closures[p$closures$gas==gas,]
  stopifnot(max(abs(m[[paste0(gas,"_prec")]]-qnorm(.975)*pc$sigma_group[match(m$UniqueID,pc$UniqueID)]))<1e-10)
 }
 for(k in c("Area","Vtot","Tcham","Pcham"))stopifnot(max(abs(m[[k]]-a[[k]][match(m$UniqueID,a$UniqueID)]))<1e-10)
 mf<-if(inst=="YMF")"ymf_black_oak_flux_compiled_with_mdf.csv" else "canopy_flux_goFlux_compiled_with_mdf.csv"
 o<-read.csv(file.path(d,"results",mf))
 cl<-m %>% filter(flag==1) %>% group_by(UniqueID) %>% summarise(t=diff(range(as.numeric(POSIX.time)))+median(diff(as.numeric(POSIX.time))),.groups="drop")
 stopifnot(max(abs(cl$t-o$t_sec[match(cl$UniqueID,o$UniqueID)]))<1e-8)
 stopifnot(max(abs(o$CH4_MDF-o$CH4_MDF_wass95))<1e-8,
   max(abs(o$CO2_MDF-o$CO2_MDF_wass95))<1e-8)
}
o<-read.csv(file.path(rd,"diurnal_blackgum/results/blackgum_flux_compiled_with_mdf.csv"))
stopifnot(all(is.finite(o$t_sec)),all(o$t_sec>30))
bg_cache<-file.path(rd,"diurnal_blackgum/results/blackgum_goflux.RData")
if (file.exists(bg_cache)) {
bg<-new.env();load(file.path(rd,"diurnal_blackgum/results/blackgum_goflux.RData"),bg)
o<-read.csv(file.path(rd,"diurnal_blackgum/results/blackgum_flux_compiled_with_mdf.csv"))
stopifnot(max(abs(bg$CH4$MDF-o$CH4_MDF_emp95[match(bg$CH4$UniqueID,o$UniqueID)]))<1e-8)
cl<-bg$manID %>% filter(flag==1) %>% group_by(UniqueID) %>% summarise(t=diff(range(as.numeric(POSIX.time)))+median(diff(as.numeric(POSIX.time))),.groups="drop")
stopifnot(max(abs(cl$t-o$t_sec[match(cl$UniqueID,o$UniqueID)]))<1e-8,all(o$t_sec>30))
} else {
 message("Blackgum cache absent: duration/window cross-check runs after flux rebuild; committed durations checked.")
}
# Basic sampling and stand accounting invariants.
stopifnot(nrow(F)==141,n_distinct(F$tree)==7,sum(F$component=="branch")==11)
r<-read.csv("scaling/out/stand_rates_HF.csv");u<-read.csv("scaling/out/stand_uncertainty_HF.csv")
stopifnot(abs(sum(r$mean*r$area)-u$estimate[u$metric=="total"])<1e-10,all(u$lo<=u$hi))
stopifnot(abs(r$area[r$comp=="stem_lo"]/.45-cone_share_below(2,CANOPY_HEIGHT))<1e-10)
cat("PASS: timestamp, units, geometry, extrapolation, fixed cohorts, metadata, durations and stand accounting\n")

trade<-read.csv("scaling/out/rate_area_tradeoff.csv")
stopifnot(diff(range(trade$rate_matching_basal*trade$area))<1e-10)
budget<-read.csv("scaling/out/blackgum_component_budget.csv")
stopifnot(all(abs(tapply(budget$share_pct,interaction(budget$form,budget$branch_rule),sum)-100)<1e-10))
