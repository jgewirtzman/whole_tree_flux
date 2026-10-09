# Flux uncertainty and model-choice sensitivity; the primary signed estimates remain unchanged.
# Run from the repository root: Rscript scaling/13_flux_uncertainty.R
# Optional: Rscript scaling/13_flux_uncertainty.R PROJECT_ROOT OUTPUT_DIR
args<-commandArgs(TRUE)
if(length(args)==0) args<-c(".","scaling/out/flux_uncertainty")
stopifnot(length(args)==2)
root<-normalizePath(args[1]);out<-normalizePath(args[2],mustWork=FALSE)
dir.create(out,recursive=TRUE,showWarnings=FALSE);setwd(root)
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(sandwich)})
source('scaling/00_load_field.R');source('scaling/analysis_helpers.R')
stopifnot(nrow(F)==141,!anyDuplicated(F$uid))
base<-'data processing/goFlux_reprocessing'
fits<-list();traces<-list();input_paths<-c('scaling/00_load_field.R','scaling/analysis_helpers.R')
for(inst in c('LGR1','LGR2','LGR3','YMF')) {
 sub<-if(inst=='YMF')'ymf_black_oak/RData' else 'RData'
 pf<-file.path(base,sub,paste0('flux_results_',inst,'.RData'))
 pt<-file.path(base,sub,paste0('manID_',inst,'.RData'))
 e<-new.env();load(pf,e);r<-get(paste0('CH4_best.',inst),e);r$instrument<-inst;fits[[inst]]<-r
 e<-new.env();load(pt,e);traces[[inst]]<-get(paste0('manID.',inst),e)
 input_paths<-c(input_paths,pf,pt)
}
fit<-bind_rows(fits);stopifnot(!anyDuplicated(fit$UniqueID))
a<-F %>% left_join(fit,by=c('uid'='UniqueID'))
stopifnot(all(is.finite(a$best.flux)),max(abs(a$flux-a$best.flux))<1e-10,all(a$model %in% c('LM','HM')))
a<-a %>% mutate(comp=case_when(component=='stem' & height_m<2~'stem_lo',component=='stem'~'stem_up',TRUE~component),
 selected_se=ifelse(model=='LM',LM.SE,HM.SE),selected_df=nb.obs-ifelse(model=='LM',2,3),
 selected_half=qt(.975,selected_df)*selected_se,selected_lo=flux-selected_half,selected_hi=flux+selected_half,
 selected_p=2*pt(-abs(flux/selected_se),selected_df),selected_sig=selected_lo>0|selected_hi<0,
 above_mdf=!below_mdf,HM_near_upper=abs(HM.k-k.max*k.mult)<=pmax(1e-10,abs(k.max*k.mult)*.01),
 HM_near_zero=abs(HM.k)<1e-6)
stopifnot(all(is.finite(a$selected_se)),all(a$selected_se>0))
# Conditional HM Wald/delta intervals with n-3 t reference are exploratory,
# not calibrated post-selection intervals. The linear fits have exact usual t
# interpretation only under their model assumptions. Geometry factors fixed.
lin<-lapply(seq_len(nrow(a)),function(i){
 r<-a[i,];d<-traces[[r$instrument]];d<-d[d$UniqueID==r$uid & !is.na(d$flag) & d$flag==1,]
 d<-d[is.finite(d$Etime)&is.finite(d$CH4dry_ppb),];n<-nrow(d)
 stopifnot(n==r$nb.obs)
 # Preserve original archive row order; report timestamp anomalies separately.
 stopifnot(all(diff(d$Etime)>=0))
 f<-lm(CH4dry_ppb~Etime,data=d);ct<-coef(summary(f));G<-r$flux.term
 flux<-unname(coef(f)[2])*G;se<-ct[2,2]*G;p<-ct[2,4]
 stopifnot(abs(flux-r$LM.flux)<1e-9,abs(p-r$LM.p.val)<1e-8)
 # goFlux formats G to 6 decimal places in its delta-method SE expression.
 stopifnot(abs(se-r$LM.SE)<1e-5)
 # Automatic Newey-West bandwidth, no prewhitening, n/(n-k) adjustment.
 # Also double bandwidth as a diagnostic, not a calibrated error envelope.
 bw<-min(n-2L,max(0L,floor(sandwich::bwNeweyWest(f,prewhite=FALSE))))
 bw2<-min(n-2L,max(1L,2L*bw))
 hs<-sqrt(sandwich::NeweyWest(f,lag=bw,prewhite=FALSE,adjust=TRUE)[2,2])*G
 hs2<-sqrt(sandwich::NeweyWest(f,lag=bw2,prewhite=FALSE,adjust=TRUE)[2,2])*G
 dt<-diff(d$Etime);res<-residuals(f)
 data.frame(uid=r$uid,linear_flux=flux,linear_se=se,linear_p=p,
 linear_lo=flux-qt(.975,n-2)*se,linear_hi=flux+qt(.975,n-2)*se,
 hac_se=hs,hac_lag=bw,hac_lag_seconds=bw*median(dt[dt>0]),
 hac_lo=flux-qt(.975,n-2)*hs,hac_hi=flux+qt(.975,n-2)*hs,
 hac2_se=hs2,hac2_lag=bw2,hac2_lo=flux-qt(.975,n-2)*hs2,hac2_hi=flux+qt(.975,n-2)*hs2,
 residual_ac1=cor(head(res,-1),tail(res,-1)),duplicate_times=sum(dt==0),
 irregular_intervals=sum(dt>0 & abs(dt-median(dt[dt>0]))>.2*median(dt[dt>0])))
})
a<-a %>% left_join(bind_rows(lin),by='uid') %>% mutate(linear_sig=linear_p<.05,hac_sig=hac_lo>0|hac_hi<0,hac2_sig=hac2_lo>0|hac2_hi<0, selected_BH=p.adjust(selected_p,'BH'), linear_BH=p.adjust(linear_p,'BH'))
stopifnot(all(a$linear_sig==(a$linear_lo>0|a$linear_hi<0)),all(a$selected_sig==(a$selected_p<.05)))
wr<-function(x,n)write.csv(x,file.path(out,paste0(n,'.csv')),row.names=FALSE)
wr(a,'closure_audit')
wr(a %>% group_by(comp) %>% summarise(n=n(),negative=sum(flux<0),LM=sum(model=='LM'),HM=sum(model=='HM'),
 n_above_mdf=sum(above_mdf),n_selected_sig=sum(selected_sig),n_linear_sig=sum(linear_sig),n_hac_sig=sum(hac_sig),n_hac2_sig=sum(hac2_sig),
 negative_above_mdf=sum(flux<0&above_mdf),negative_selected_sig=sum(flux<0&selected_sig),
 negative_linear_sig=sum(linear_flux<0&linear_sig),negative_hac_sig=sum(linear_flux<0&hac_sig),
 median_hac_se_ratio=median(hac_se/linear_se),median_residual_ac1=median(residual_ac1),.groups='drop'),'component_counts')
wr(a %>% count(comp,model,above_mdf,selected_sig),'classification_crosswalk')
wr(a %>% filter(flux<0) %>% select(uid,tree,comp,height_m,model,flux,MDF,above_mdf,selected_se,selected_lo,selected_hi,selected_p,selected_sig,
 linear_flux,linear_lo,linear_hi,linear_p,hac_lo,hac_hi,hac_sig,hac2_sig,residual_ac1,HM_near_upper,HM_near_zero),'negative_fluxes')
wr(a %>% group_by(model) %>% summarise(n=n(),significant=sum(selected_sig),near_HM_upper=sum(HM_near_upper,na.rm=TRUE),near_HM_zero=sum(HM_near_zero,na.rm=TRUE),.groups='drop'),'model_diagnostics')
# Preserve sampled trees, heights, area weights and all observations when setting
# flagged estimates to zero. These are stress tests, not recommended estimators.
scenarios<-list(selected_all=a$flux,selected_mdf_zero=ifelse(a$above_mdf,a$flux,0),
 selected_interval_zero=ifelse(a$selected_sig,a$flux,0),linear_all=a$linear_flux,
 linear_iid_zero=ifelse(a$linear_sig,a$linear_flux,0),linear_hac_zero=ifelse(a$hac_sig,a$linear_flux,0),
 linear_hac_double_zero=ifelse(a$hac2_sig,a$linear_flux,0))
G<-load_geometry(F);areas<-read.csv('scaling/out/stand_rates_HF.csv');ar<-setNames(areas$area,areas$comp)
component_rates<-function(d) d %>% group_by(comp,tree) %>% group_modify(function(.x,.y){
 h<-height_means(.x)
 if(.y$comp %in% c('stem_lo','stem_up')) {
  g<-filter(G,tree==.y$tree);lo<-if(.y$comp=='stem_lo')0 else 2;up<-if(lo==0)2 else g$top[1]
  v<-sum(profile_weights(h$height_m,stem_grid(g,lo,up,h$height_m))*h$flux)
 }else v<-mean(h$flux)
 tibble(rate=v)
}) %>% ungroup()
stand<-list();paired<-list();forms<-list();prof<-list()
fixed_cohort<-read.csv('scaling/out/form_test_coverage.csv');fixed_cohort<-fixed_cohort$tree[fixed_cohort$exp_eligible]
for(k in names(scenarios)){
 d<-a;d$flux<-scenarios[[k]]
 cr<-component_rates(filter(d,site!='Yale Myers')) %>% group_by(comp) %>% summarise(rate=mean(rate),.groups='drop') %>% mutate(scenario=k,area=ar[comp],integrated=rate*area)
 stand[[k]]<-cr
 # Same paired mean convention as current manuscript (height then tree).
 tm<-d %>% filter(component=='stem') %>% group_by(tree,comp,height_m) %>% summarise(f=mean(flux),.groups='drop') %>%
 group_by(tree,comp) %>% summarise(f=mean(f),.groups='drop') %>% pivot_wider(names_from=comp,values_from=f)
 w<-suppressWarnings(wilcox.test(tm$stem_lo,tm$stem_up,paired=TRUE,exact=TRUE))
 paired[[k]]<-data.frame(scenario=k,basal_mean=mean(tm$stem_lo),upper_mean=mean(tm$stem_up),ratio=mean(tm$stem_lo)/mean(tm$stem_up),paired_wilcoxon_p=w$p.value)
 prof[[k]]<-component_rates(filter(d,component=='stem')) %>% mutate(scenario=k)
 # Point-estimate basal extrapolation audit at FIXED original five-tree cohort.
 ds<-filter(d,component=='stem');low<-ds %>% filter(height_m<2) %>% group_by(tree,height_m) %>% summarise(flux=mean(flux),.groups='drop')
 bp<-pooled_slope(low)
 forms[[k]]<-bind_rows(lapply(split(filter(ds,tree %in% fixed_cohort),filter(ds,tree %in% fixed_cohort)$tree),function(x){
  lo<-height_means(filter(x,height_m<2));hi<-height_means(filter(x,height_m>=2));g<-filter(G,tree==x$tree[1])
  grid<-stem_grid(g,2,g$top[1],hi$height_m);fs<-extrapolation_forms(lo,bp)
  obs<-sum(profile_weights(hi$height_m,grid)*hi$flux)
  bind_rows(lapply(names(fs),function(f)tibble(scenario=k,tree=x$tree[1],form=f,observed=obs,predicted=sum(grid$area*fs[[f]](grid$h))/sum(grid$area))))
 }))
}
st<-bind_rows(stand);wr(st,'stand_components');wr(st %>% group_by(scenario) %>% summarise(total=sum(integrated),woody=sum(integrated[comp!='leaf']),branch_share=100*integrated[comp=='branch']/sum(integrated),leaf_share=100*integrated[comp=='leaf']/sum(integrated),.groups='drop'),'stand_totals')
ref<-read.csv('scaling/out/stand_rates_HF.csv');st0<-filter(st,scenario=='selected_all');stopifnot(max(abs(st0$rate-ref$mean[match(st0$comp,ref$comp)]))<1e-10)
wr(bind_rows(paired),'basal_upper_comparison');wr(bind_rows(prof),'stem_profiles');fo<-bind_rows(forms);wr(fo,'extrapolation_tree_scores')
wr(fo %>% group_by(scenario,form) %>% summarise(n_tree=n(),n_valid=sum(is.finite(predicted)),ratio=if(all(is.finite(predicted)))sum(predicted)/sum(observed) else NA_real_,predicted_negative=sum(predicted<0,na.rm=TRUE),observed_negative=sum(observed<0),.groups='drop'),'extrapolation_summary')
# Wetland reference is separate from the 141-observation upland cohort.
bgp<-file.path(base,'diurnal_blackgum/results/blackgum_goflux.RData')
if(!file.exists(bgp)) stop('Rebuild archived swamp fits first: Rscript run_analysis.R --flux-only')
e<-new.env();load(bgp,e)
bgc<-read.csv(file.path(base,'diurnal_blackgum/results/blackgum_flux_compiled_with_mdf.csv'))
bg<-get('CH4',e);j<-match(bg$UniqueID,bgc$UniqueID)
stopifnot(!anyNA(j),max(abs(bg$best.flux-bgc$CH4_best.flux[j]))<1e-10)
bg<-bg %>% mutate(height_m=bgc$height_m[j],selected_se=ifelse(model=='LM',LM.SE,HM.SE),
 selected_df=nb.obs-ifelse(model=='LM',2,3),selected_lo=best.flux-qt(.975,selected_df)*selected_se,
 selected_hi=best.flux+qt(.975,selected_df)*selected_se,selected_sig=selected_lo>0|selected_hi<0)
wr(bg,'wetland_reference_audit')
wr(bg %>% group_by(height_m) %>% summarise(n=n(),positive=sum(best.flux>0),selected_significant=sum(selected_sig),linear_significant=sum(LM.p.val<.05),selected_mean=mean(best.flux),linear_mean=mean(LM.flux),.groups='drop'),'wetland_reference_summary')
input_paths<-c(input_paths,bgp,file.path(base,'diurnal_blackgum/results/blackgum_flux_compiled_with_mdf.csv'))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
input_paths<-c(input_paths,file.path(base,'results/canopy_flux_goFlux_compiled_with_mdf.csv'),file.path(base,'ymf_black_oak/results/ymf_black_oak_flux_compiled_with_mdf.csv'),'scaling/out/stand_rates_HF.csv','scaling/out/form_test_coverage.csv')
wr(data.frame(path=input_paths,md5=unname(tools::md5sum(input_paths))),'input_checksums')
print(read.csv(file.path(out,'component_counts.csv')));print(read.csv(file.path(out,'stand_totals.csv')));print(read.csv(file.path(out,'negative_fluxes.csv')))
