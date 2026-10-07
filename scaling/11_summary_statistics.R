# Current manuscript numbers, with tree-level inference instead of closure-level tests.
source("scaling/00_load_field.R")
source("scaling/analysis_helpers.R")
suppressPackageStartupMessages(library(readr))
OUT<-"scaling/out"
comp<-function(d) d %>% mutate(comp=case_when(component=="stem" & height_m<2~"stem_lo",component=="stem"~"stem_up",TRUE~component))
summary_rates<-function(d) comp(d) %>% group_by(comp,tree,height_m) %>% summarise(f=mean(flux),.groups="drop") %>%
  group_by(comp,tree) %>% summarise(f=mean(f),.groups="drop") %>% group_by(comp) %>% summarise(mean=mean(f),ntree=n(),.groups="drop")
set.seed(45)
bo<-bind_rows(lapply(seq_len(N_BOOT),function(b)summary_rates(resample_trees(F)) %>% mutate(b=b)))
ci<-bo %>% group_by(comp) %>% summarise(lo=ci95(mean)[1],hi=ci95(mean)[2],valid_draws=n(),.groups="drop")
cs<-comp(F) %>% group_by(comp) %>% summarise(n=n(),median=median(flux),pooled_mean=mean(flux),negative=sum(flux<0),below=sum(below_mdf),.groups="drop") %>%
  left_join(summary_rates(F),by="comp") %>% left_join(ci,by="comp")
write_csv(cs,file.path(OUT,"component_summary.csv"))
TM<-comp(F) %>% filter(component=="stem") %>% group_by(tree,comp,height_m) %>% summarise(f=mean(flux),.groups="drop") %>%
  group_by(tree,comp) %>% summarise(f=mean(f),.groups="drop") %>% pivot_wider(names_from=comp,values_from=f)
wt<-wilcox.test(TM$stem_lo,TM$stem_up,paired=TRUE,alternative="two.sided",exact=TRUE)
brat<-bo %>% select(b,comp,mean) %>% filter(comp %in% c("stem_lo","stem_up")) %>% pivot_wider(names_from=comp,values_from=mean) %>% mutate(ratio=stem_lo/stem_up)
hf<-read.csv("data processing/goFlux_reprocessing/results/canopy_flux_goFlux_compiled_with_mdf.csv")
ym<-read.csv("data processing/goFlux_reprocessing/ymf_black_oak/results/ymf_black_oak_flux_compiled_with_mdf.csv")
dat<-bind_rows(hf %>% select(UniqueID,starts_with("CH4_"),starts_with("CO2_")),ym %>% select(UniqueID,starts_with("CH4_"),starts_with("CO2_"))) %>% filter(UniqueID %in% F$uid)
metrics<-c(n=nrow(F),ntree=n_distinct(F$tree),mdf_median=median(dat$CH4_MDF_wass95),
  basal_upper_ratio=mean(TM$stem_lo)/mean(TM$stem_up),ratio_lo=unname(ci95(brat$ratio)[1]),ratio_hi=unname(ci95(brat$ratio)[2]),
  paired_W=unname(wt$statistic),paired_p=wt$p.value,
  upper5_n=sum(F$component=="stem" & F$height_m>=5),upper5_neg=sum(F$component=="stem" & F$height_m>=5 & F$flux<0),
  upper10_n=sum(F$component=="stem" & F$height_m>=10),upper10_neg=sum(F$component=="stem" & F$height_m>=10 & F$flux<0))
# Transparent filter sensitivity on the retained HF cohort. Zeroing is a diagnostic
# scenario, never the default treatment of observations below the MDF.
h<-hf[hf$UniqueID %in% F$uid,]; x<-filter(F,site!="Yale Myers"); h<-h[match(x$uid,h$UniqueID),]
filters<-list(none=rep(FALSE,nrow(h)),manufacturer=h$CH4_below_MDF_goflux,
  r2_07=h$CH4_best.r2<=.7,r2_09=h$CH4_best.r2<=.9,co2_r2=h$CO2_best.r2<=.7,
  snr2=h$CH4_SNR<=2,snr3=h$CH4_SNR<=3,emp90=h$CH4_below_MDF_wass90,emp95=h$CH4_below_MDF_wass95,
  emp99=h$CH4_below_MDF_wass99,chr90=h$CH4_below_MDF_chr90,chr95=h$CH4_below_MDF_chr95,chr99=h$CH4_below_MDF_chr99)
G<-load_geometry(F); areas<-read.csv(file.path(OUT,"stand_rates_HF.csv")); names_area<-setNames(areas$area,areas$comp)
filter_rows<-bind_rows(lapply(names(filters),function(k) {
  d<-comp(x); fail<-filters[[k]]; stopifnot(!anyNA(fail)); d$flux[fail]<-0
  rates<-d %>% group_by(comp,tree) %>% group_modify(function(.x,.y) {
    hm<-height_means(.x)
    if(.y$comp %in% c("stem_lo","stem_up")) {
      g<-filter(G,tree==.y$tree); lo<-if(.y$comp=="stem_lo") 0 else 2; up<-if(lo==0) 2 else g$top[1]
      v<-sum(profile_weights(hm$height_m,stem_grid(g,lo,up,hm$height_m))*hm$flux)
    } else v<-mean(hm$flux)
    tibble(f=v)
  }) %>% ungroup() %>% group_by(comp) %>% summarise(rate=mean(f),.groups="drop") %>% mutate(filter=k,integrated=rate*names_area[comp])
  rates
}))
write_csv(filter_rows,file.path(OUT,"filter_sensitivity_HF.csv"))
filter_totals<-filter_rows %>% group_by(filter) %>% summarise(total=sum(integrated),.groups="drop")
metrics<-c(metrics,filter_total_min=min(filter_totals$total),filter_total_max=max(filter_totals$total))
write_csv(tibble(metric=names(metrics),value=as.numeric(metrics)),file.path(OUT,"summary_statistics.csv"))
print(cs); print(metrics)
