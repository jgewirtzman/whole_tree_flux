# Basal-only extrapolation evaluated on held-out upper stems using measured geometry.
source("scaling/00_load_field.R")
source("scaling/analysis_helpers.R")
suppressPackageStartupMessages({library(readr); library(ggplot2)})
OUT <- "scaling/out"; dir.create(OUT, showWarnings=FALSE)
S <- F %>% filter(component=="stem")
G <- load_geometry(F)
low <- S %>% filter(height_m<2) %>% group_by(tree,height_m) %>% summarise(flux=mean(flux),.groups="drop")
B_POOL <- pooled_slope(low)
# Establish a fixed common cohort for all forms before any comparisons/bootstraps.
coverage <- low %>% group_by(tree) %>% summarise(exp_eligible=all(flux>0),.groups="drop") %>%
  mutate(reason=ifelse(exp_eligible,"eligible","nonpositive basal height mean; log-linear exponential unavailable"))
common <- coverage$tree[coverage$exp_eligible]
write_csv(coverage,file.path(OUT,"form_test_coverage.csv"))
grids <- lapply(split(G,G$tree),function(g) {
  h<-S$height_m[S$tree==g$tree[1] & S$height_m>=2]
  stem_grid(g,2,g$top[1],h)
})
score <- function(d, original, bp) {
  lo<-height_means(d[d$height_m<2,]); hi<-d[d$height_m>=2,]; hm<-height_means(hi)
  fs<-extrapolation_forms(lo,bp); grid<-grids[[original]]
  obs_area<-sum(hm$flux*profile_weights(hm$height_m,grid))
  bind_rows(lapply(names(fs),function(k) {
    pred<-fs[[k]](hi$height_m)
    tibble(tree=d$tree[1],source_tree=original,form=k,n=nrow(hi),bias=mean(pred-hi$flux),mae=mean(abs(pred-hi$flux)),
      obs_neg_share=mean(hi$flux<0),pred_neg_share=mean(pred<0),obs_profile=mean(hm$flux),pred_profile=mean(fs[[k]](hm$height_m)),
      obs_area=obs_area,pred_area=sum(grid$area*fs[[k]](grid$h))/sum(grid$area))
  }))
}
sc <- bind_rows(lapply(split(S,S$tree),function(d)score(d,d$tree[1],B_POOL)))
write_csv(sc,file.path(OUT,"form_test_scores.csv"))
P<-bind_rows(lapply(split(S,S$tree),function(d) {
  fs<-extrapolation_forms(height_means(d[d$height_m<2,]),B_POOL)
  bind_rows(lapply(names(fs),function(k)d[d$height_m>=2,] %>% mutate(form=k,pred=fs[[k]](height_m))))
}))
write_csv(P,file.path(OUT,"form_test_points.csv"))
pool <- function(d) d %>% group_by(form) %>% summarise(trees=n(),points=sum(n),median_bias=median(bias),median_mae=median(mae),
  ratio_pred_obs_profile=sum(pred_profile)/sum(obs_profile),ratio_pred_obs_area=sum(pred_area)/sum(obs_area),
  trees_pred_net_uptake=sum(pred_area<0),trees_obs_net_uptake=sum(obs_area<0),.groups="drop")
write_csv(pool(filter(sc,tree %in% common)),file.path(OUT,"form_test_pooled.csv"))
# All-tree sensitivity is explicitly separate: exponential is unavailable for this cohort.
write_csv(pool(filter(sc,form!="exp_decay")),file.path(OUT,"form_test_pooled_all_trees.csv"))
set.seed(42)
boot <- bind_rows(lapply(seq_len(N_BOOT),function(b) {
  sb<-resample_trees(S,common)
  lb<-sb %>% filter(height_m<2) %>% group_by(tree,height_m) %>% summarise(flux=mean(flux),.groups="drop")
  bp<-pooled_slope(lb) # pooled slope is re-fitted, not held fixed across bootstrap draws
  ss<-bind_rows(lapply(split(sb,sb$tree),function(d)score(d,d$source_tree[1],bp)))
  ss %>% group_by(form) %>% summarise(ratio=if(any(!is.finite(pred_area))) NA_real_ else sum(pred_area)/sum(obs_area),
    uptake_trees=if(anyNA(pred_area)) NA_real_ else sum(pred_area<0),.groups="drop") %>% mutate(b=b)
}))
bs<-boot %>% group_by(form) %>% summarise(ratio_median=median(ratio,na.rm=TRUE),ratio_lo=ci95(ratio)[1],ratio_hi=ci95(ratio)[2],
  uptake_trees_median=median(uptake_trees,na.rm=TRUE),valid_draws=sum(is.finite(ratio)),total_draws=n(),
  invalid_fraction=mean(!is.finite(ratio)),trees=length(common),.groups="drop")
# Missing exponential fits invalidate the entire form/draw; never shrink its denominator.
write_csv(bs,file.path(OUT,"form_test_bootstrap.csv"))
print(coverage); print(pool(filter(sc,tree %in% common))); print(bs)
