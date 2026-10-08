# Flux-painted woody models. Interpolation + endpoint extension are visual rules,
# not a new fitted whole-tree budget. No observations or manuscript files changed.
suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(tidyr);library(readr);library(patchwork)})
args <- commandArgs(trailingOnly=TRUE)
PROJECT_ROOT <- normalizePath(if(length(args)) args[1] else '.', mustWork=TRUE)
BUNDLE <- normalizePath(if(length(args)>1) args[2] else file.path(PROJECT_ROOT,'scaling/tree_illustrations'),mustWork=TRUE)
DEST <- file.path(BUNDLE,'presentations')
CANDIDATES <- file.path(BUNDLE,'si_candidates')
dir.create(DEST,recursive=TRUE,showWarnings=FALSE)
dir.create(CANDIDATES,recursive=TRUE,showWarnings=FALSE)
oldwd <- setwd(PROJECT_ROOT)
source('scaling/00_load_field.R')
ROOT <- PROJECT_ROOT
meta <- read_csv(file.path(BUNDLE,'inputs/individual_crown_inputs.csv'),show_col_types=FALSE)
audit <- read_csv(file.path(BUNDLE,'inputs/surface_area_audit.csv'),show_col_types=FALSE)
stems <- read_csv(file.path(BUNDLE,'inputs/modeled_stem_frusta.csv'),show_col_types=FALSE)
branches <- read_csv(file.path(BUNDLE,'inputs/modeled_branch_frusta.csv'),show_col_types=FALSE)
models <- lapply(meta$tree,function(id) list(stem=filter(stems,tree==id),branch=filter(branches,tree==id)))
D <- bind_rows(F,read_csv(file.path(ROOT,'scaling/out/blackgum_collar_diel_means.csv'),show_col_types=FALSE) %>% select(-any_of('n')))
sampling_marks <- FALSE
gray_color <- '#DCDCDC'
M<-D %>% group_by(tree,component,height_m) %>% summarise(flux=mean(flux),n=n(),.groups='drop')
stand<-read_csv(file.path(ROOT,'scaling/out/stand_rates_HF.csv'),show_col_types=FALSE)
BRANCH_MEAN<-stand$mean[stand$comp=='branch'];stopifnot(length(BRANCH_MEAN)==1)
mode<-Sys.getenv('MISSING_BRANCH_MODE','assigned')
stopifnot(mode %in% c('assigned','gray'))
branch_style<-Sys.getenv('BRANCH_STYLE','flux')
stopifnot(branch_style %in% c('flux','light_gray'))
interp<-function(d,h) {if(nrow(d)==0)return(rep(if(mode=='assigned')BRANCH_MEAN else NA_real_,length(h)));if(nrow(d)==1)return(rep(d$flux,length(h)));approx(d$height_m,d$flux,xout=h,rule=2,ties=mean)$y}
# Splitting tapered segments retains their surface area, and gives smoother color
# changes along long rising branches; no geometry/distribution is re-estimated.
split_seg<-function(d) {
 counts<-pmax(1,ceiling(abs(d$y1-d$y0)/.14))
 dd<-d %>% mutate(original=seq_len(n()),pieces=counts) %>% uncount(pieces,.remove=FALSE,.id='piece')
 a<-(dd$piece-1)/dd$pieces;b<-dd$piece/dd$pieces
 old<-dd
 for(v in c('x','y','z','r')) {lo<-paste0(v,'0');hi<-paste0(v,'1');dd[[lo]]<-old[[lo]]+(old[[hi]]-old[[lo]])*a;dd[[hi]]<-old[[lo]]+(old[[hi]]-old[[lo]])*b}
 dd$L<-sqrt((dd$x1-dd$x0)^2+(dd$y1-dd$y0)^2+(dd$z1-dd$z0)^2)
 dd
}
polygons<-function(seg) {
 dx<-seg$x1-seg$x0;dy<-seg$y1-seg$y0;ll<-pmax(sqrt(dx^2+dy^2),1e-9);nx<--dy/ll;ny<-dx/ll
 tibble(id=rep(seq_len(nrow(seg)),each=4),x=as.vector(rbind(seg$x0+nx*seg$r0,seg$x1+nx*seg$r1,seg$x1-nx*seg$r1,seg$x0-nx*seg$r0)),y=as.vector(rbind(seg$y0+ny*seg$r0,seg$y1+ny*seg$r1,seg$y1-ny*seg$r1,seg$y0-ny*seg$r0)),rate=rep(seg$rate,each=4))
}
# Equal absolute fluxes have equal HCL lightness and chroma on both sides of zero.
# The displayed legend is cropped to the retained data range, not independently
# stretched on its negative and positive sides. White is exactly zero.
lim<-c(-.05,160) # Full-data bounds used to validate the initial assignments.
tr<-scales::trans_new('asinh_flux',function(x)asinh(x/.02),function(x)sinh(x)*.02)
flux_colors<-function(rate,reference) {
 intensity<-(abs(asinh(rate/.02))/asinh(reference/.02))^.7
 grDevices::hcl(h=ifelse(rate<0,250,10),c=60*intensity,l=100-65*intensity,fixup=FALSE)
}
make_scale<-function(limits,breaks) {
 reference<-max(abs(limits))
 color_rates<-sort(unique(c(sinh(seq(asinh(limits[1]/.02),asinh(limits[2]/.02),length.out=2049))*.02,0)))
 col<-flux_colors(color_rates,reference)
 vals<-scales::rescale(asinh(color_rates/.02),from=asinh(limits/.02))
 stopifnot(!anyNA(col),flux_colors(0,reference)=='#FFFFFF')
 scale_fill_gradientn(colors=col,values=vals,transform=tr,limits=limits,
  breaks=breaks,labels=function(x)gsub('-', '−', vapply(x,function(z)format(z,trim=TRUE,scientific=FALSE),character(1))),
  na.value='#BCC4BC',name=expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
  guide=guide_colorbar(direction='horizontal',barwidth=unit(105,'mm'),barheight=unit(4,'mm'),title.position='top'))
}
painted<-list();assignment<-list()
for(i in 1:8) {
 m<-meta[i,];st<-M %>% filter(tree==m$tree,component=='stem');br<-M %>% filter(tree==m$tree,component=='branch')
 a<-split_seg(models[[i]]$stem);b<-split_seg(models[[i]]$branch)
 a$rate<-interp(st,(a$y0+a$y1)/2);b$rate<-interp(br,(b$y0+b$y1)/2)
 area<-function(x)sum(pi*(x$r0+x$r1)*sqrt(x$L^2+(x$r0-x$r1)^2))
 stopifnot(abs(area(b)/area(a)-audit$achieved_ratio[i])<1e-8,all(is.finite(a$rate)),all(a$rate>=lim[1]&a$rate<=lim[2]))
 if(mode=='assigned')stopifnot(all(is.finite(b$rate)),all(b$rate>=lim[1]&b$rate<=lim[2]))
 # Interpolation passes through every height mean and holds both endpoints.
 stopifnot(max(abs(interp(st,st$height_m)-st$flux))<1e-12,
           identical(interp(st,c(-1,100)),c(st$flux[which.min(st$height_m)],st$flux[which.max(st$height_m)])))
 painted[[i]]<-list(stem=a,branch=b,assigned=nrow(br)==0)
 assignment[[i]]<-tibble(tree=m$tree,branch_source=if(nrow(br))'own measured height means' else if(mode=='assigned')'assigned measured upland tree-weighted branch mean' else 'unmeasured: gray',branch_measured_heights=nrow(br),assigned_branch_rate=if(nrow(br)||mode=='gray')NA_real_ else BRANCH_MEAN,stem_low=min(st$height_m),stem_high=max(st$height_m),stem_bottom_rate=interp(st,0),stem_top_rate=interp(st,m$H),branch_stem_area_ratio=area(b)/area(a))
}
write_csv(bind_rows(assignment),file.path(DEST,'painted_flux_assignment.csv'))
write_csv(M,file.path(DEST,'painted_input_height_means.csv'))
if(Sys.getenv('EXPORT_PAINTED_SEGMENTS','0')=='1') write_csv(bind_rows(lapply(1:8,function(i)bind_rows(painted[[i]]$stem %>% mutate(component='stem'),painted[[i]]$branch %>% mutate(component='branch')) %>% mutate(tree=meta$tree[i]))),file.path(DEST,'painted_segment_rates.csv'))
# Draw each stem as ONE filled path. Adjacent opaque quadrilaterals can expose
# white antialias seams in PNG/PDF viewers even when their edges coincide.
# A native continuous gradient removes all internal polygon boundaries.
stem_gradient<-function(i,reference) {
 m<-meta[i,];stem<-models[[i]]$stem
 h<-c(stem$y0,tail(stem$y1,1));r<-c(stem$r0,tail(stem$r1,1))
 stopifnot(all(stem$x0==0),all(stem$x1==0),all(diff(h)>0))
 w<-max(r)
 x<-c(-r,rev(r));y<-c(h,rev(h))
 dat<-M %>% filter(tree==m$tree,component=='stem')
 # Dense color stops approximate the same flux interpolation and nonlinear
 # palette continuously; include observation heights and zero crossings exactly.
 crossing<-which(head(dat$flux,-1)*tail(dat$flux,-1)<0)
 zero_h<-vapply(crossing,function(j)dat$height_m[j]-dat$flux[j]*(dat$height_m[j+1]-dat$height_m[j])/(dat$flux[j+1]-dat$flux[j]),numeric(1))
 stops_h<-sort(unique(c(seq(0,m$H,length.out=2049),dat$height_m,zero_h)))
 stops_h<-stops_h[stops_h>=0&stops_h<=m$H]
 color<-flux_colors(interp(dat,stops_h),reference)
 stopifnot(!anyNA(color))
 fill<-grid::linearGradient(colours=color,stops=stops_h/m$H,x1=.5,x2=.5,y1=0,y2=1)
 grob<-grid::polygonGrob(x=(x+w)/(2*w),y=y/m$H,
   gp=grid::gpar(fill=fill,col=NA))
 annotation_custom(grob,xmin=-w,xmax=w,ymin=0,ymax=m$H)
}
plot_tree<-function(i,ncol,color_scale,reference) {
 m<-meta[i,];z<-painted[[i]];b<-polygons(z$branch);a<-polygons(z$stem)
 branch_layer<-if(branch_style=='light_gray')geom_polygon(data=b,aes(x,y,group=id),fill=gray_color,color=NA) else geom_polygon(data=b,aes(x,y,group=id,fill=rate),color=NA)
 p<-ggplot()+branch_layer+geom_blank(data=a,aes(x,y,fill=rate),show.legend=TRUE)+stem_gradient(i,reference)+color_scale+coord_fixed(ratio=1,xlim=c(-6.0,6.0),ylim=c(0,24),expand=FALSE)+
 scale_y_continuous(breaks=c(0,5,10,15,20))+scale_x_continuous(breaks=NULL)+labs(title=paste0(m$label,if(z$assigned && branch_style=='flux')' *' else ''),subtitle=if(m$site=='Harvard Forest')NULL else m$site,x=NULL,y=NULL)+
 theme_void(base_family='Helvetica',base_size=10)+theme(plot.title=element_text(face='italic',size=11,hjust=.5),plot.subtitle=element_text(size=8,color='#60685E',hjust=.5),plot.margin=margin(6,2,4,2),legend.position='bottom')
 if(sampling_marks) {
  observed <- M %>% filter(tree==m$tree,component=='stem')
  observed$r <- approx(c(models[[i]]$stem$y0,tail(models[[i]]$stem$y1,1)),c(models[[i]]$stem$r0,tail(models[[i]]$stem$r1,1)),xout=observed$height_m,rule=2)$y
  p <- p+geom_segment(data=observed,aes(x=r+.06,xend=r+.32,y=height_m,yend=height_m),inherit.aes=FALSE,color='#737373',linewidth=.22)
 }
 p
}
save<-function(p,name,w,h){ggsave(file.path(DEST,paste0(name,'.png')),p,width=w,height=h,dpi=210,bg='white',limitsize=FALSE);ggsave(file.path(DEST,paste0(name,'.pdf')),p,width=w,height=h,device=cairo_pdf,bg='white',limitsize=FALSE)}
tree_layouts <- list()
scale_audit<-list();color_checks<-list()
for(branch_style in c('flux','light_gray')) for(layout in c('upland_row','all_two_rows')) {
 ids<-if(layout=='upland_row')1:7 else 1:8;nc<-if(layout=='upland_row')7 else 4
 rates<-unlist(lapply(ids,function(i)c(painted[[i]]$stem$rate,painted[[i]]$branch$rate)))
 actual<-range(rates,na.rm=TRUE)
 # Round outward to one significant digit, using only trees in this layout.
 outward<-function(x){step<-10^floor(log10(abs(x)));sign(x)*ceiling(abs(x)/step)*step}
 limits<-if(layout=='upland_row')vapply(actual,outward,numeric(1)) else c(-.05,160)
 breaks<-if(layout=='upland_row')sort(unique(c(limits[1],0,.1,.5,1,limits[2]))) else c(-.05,0,.1,1,10,160)
 breaks<-breaks[breaks>=limits[1]&breaks<=limits[2]]
 stopifnot(all(rates[is.finite(rates)]>=limits[1]&rates[is.finite(rates)]<=limits[2]))
 reference<-max(abs(limits))
 checks<-tibble(layout=layout,magnitude=sort(unique(c(.01,.05,.1,1,reference))))
 checks<-checks %>% filter(magnitude<=reference) %>% mutate(negative=flux_colors(-magnitude,reference),positive=flux_colors(magnitude,reference))
 checks$negative_lightness<-farver::decode_colour(checks$negative,to='lab')[,1]
 checks$positive_lightness<-farver::decode_colour(checks$positive,to='lab')[,1]
 stopifnot(max(abs(checks$negative_lightness-checks$positive_lightness))<.5)
 color_checks[[layout]]<-checks
 scale_audit[[layout]]<-tibble(layout=layout,actual_min=actual[1],actual_max=actual[2],legend_min=limits[1],legend_max=limits[2],magnitude_reference=reference)
 p<-wrap_plots(lapply(ids,plot_tree,ncol=nc,color_scale=make_scale(limits,breaks),reference=reference),ncol=nc,guides='collect') & theme(legend.position='bottom')
 p<-p+plot_annotation(theme=theme(plot.margin=margin(8,8,8,8)))
 tree_layouts[[paste(layout,branch_style,sep='_')]] <- p
 save(p,paste0('T_',layout,if(branch_style=='light_gray')'_gray_branches' else ''),if(nc==7)16 else 12,if(nc==7)6.2 else 12.2)
}
cat('PASS: all painted geometries preserve their audited skeleton area; stem interpolation passes through means; color limits retain all rates.\n')

cat('PASS: blue and red have matched lightness/chroma at equal magnitudes; white is zero.\n')

write_csv(bind_rows(scale_audit),file.path(DEST,'color_scale_ranges.csv'))
write_csv(bind_rows(color_checks),file.path(DEST,'color_balance_check.csv'))
print(bind_rows(scale_audit))

# A separate publication candidate: slightly darker gray branches and precise
# sampling-height notches; no assigned branch colors. Presentation exports above
# retain the original uncluttered versions.
branch_style <- 'light_gray'; gray_color <- '#CBCBCB'; sampling_marks <- TRUE
candidate_scale <- make_scale(c(-.04,3),c(-.04,0,.1,.5,1,3))
candidate_scale$name <- expression(Stem~CH[4]~flux~(nmol~m^{-2}~s^{-1}))
upland_candidate <- wrap_plots(lapply(1:7,plot_tree,ncol=7,color_scale=candidate_scale,reference=3),ncol=7,guides='collect') & theme(legend.position='bottom',plot.title=element_text(size=7.5),plot.subtitle=element_text(size=6),legend.text=element_text(size=6),legend.title=element_text(size=8))
upland_candidate <- upland_candidate+plot_annotation(theme=theme(plot.margin=margin(5,5,5,5)))
ggsave(file.path(CANDIDATES,'upland_stem_flux_sampling.png'),upland_candidate,width=190,height=81,units='mm',dpi=300,bg='white')
ggsave(file.path(CANDIDATES,'upland_stem_flux_sampling.pdf'),upland_candidate,width=190,height=81,units='mm',device=cairo_pdf,bg='white')

# Two-row SI alternatives: enlarged trees and sampling marks. The seven-tree
# version reserves the eighth cell for its legend; the eight-tree comparison
# keeps the wider wetland-inclusive range and a bottom legend.
for(include_wetland in c(FALSE,TRUE)) {
 ids <- if(include_wetland)1:8 else 1:7
 limits <- if(include_wetland)c(-.05,160)else c(-.04,3)
 breaks <- if(include_wetland)c(-.05,0,.1,1,10,160)else c(-.04,0,.1,.5,1,3)
 cs <- make_scale(limits,breaks)
 cs$name <- expression(atop(Stem~CH[4]~flux,(nmol~m^{-2}~s^{-1})))
 if(!include_wetland)cs$guide <- guide_colorbar(direction='vertical',barwidth=unit(4,'mm'),barheight=unit(42,'mm'),title.position='top',title.hjust=0)
 plots <- lapply(ids,plot_tree,ncol=4,color_scale=cs,reference=max(abs(limits)))
 if(!include_wetland)plots <- c(plots,list(guide_area()))
 two_rows <- wrap_plots(plots,ncol=4,guides='collect') & theme(legend.position=if(include_wetland)'bottom'else'right',plot.title=element_text(size=9),plot.subtitle=element_text(size=7),legend.text=element_text(size=7),legend.title=element_text(size=8))
 two_rows <- two_rows+plot_annotation(theme=theme(plot.margin=margin(5,5,5,5)))
 name <- if(include_wetland)'all_stem_flux_sampling_two_rows'else'upland_stem_flux_sampling_two_rows'
 for(ext in c('png','pdf')) {
  args_out <- list(filename=file.path(CANDIDATES,paste0(name,'.',ext)),plot=two_rows,width=190,height=if(include_wetland)205 else 190,units='mm',bg='white')
  if(ext=='png')args_out$dpi<-300 else args_out$device<-cairo_pdf
  do.call(ggsave,args_out)
 }
}

# Current Figure S5 and color comparisons are rendered by render_color_options.R.
