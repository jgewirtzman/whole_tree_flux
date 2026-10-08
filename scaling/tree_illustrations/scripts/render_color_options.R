# Color/layout comparisons for Figure S5. Run from the repository root.
args <- commandArgs(trailingOnly=TRUE)
project <- normalizePath(if(length(args))args[1]else'.',mustWork=TRUE)
bundle <- normalizePath(if(length(args)>1)args[2]else file.path(project,'scaling/tree_illustrations'),mustWork=TRUE)
source(file.path(bundle,'scripts/render_trees.R'))
# All comparisons preserve geometry and flux interpolation. Only color mapping,
# labels and surrounding space change. Remove sampling-height marks throughout.
sampling_marks <- FALSE;branch_style <- 'light_gray';gray_color <- '#CBCBCB'
COLOR_CURVE <- 'original';COLOR_HUE <- 10
# A thin neutral edge keeps nearly white, near-zero stems legible without
# changing their flux colors or dimensions.
body_text <- paste(deparse(body(stem_gradient)),collapse='\n')
body_text <- sub('col = NA','col = "#D9D9D9", lwd = .25',body_text,fixed=TRUE)
body(stem_gradient) <- parse(text=body_text)[[1]]
color_position <- function(rate) {
 switch(COLOR_CURVE,linear=rate,gentle=sign(rate)*abs(rate)^.65,
        shared=asinh(rate/.5),original=sign(rate)*abs(asinh(rate/.02))^.7)
}
flux_colors <- function(rate,reference) {
 intensity <- abs(color_position(rate))/abs(color_position(reference))
 stopifnot(all(intensity<=1+1e-10))
 grDevices::hcl(h=ifelse(rate<0,250,COLOR_HUE),c=60*intensity,l=100-65*intensity,fixup=FALSE)
}
option_scale <- function(limits,breaks,label,order=1,shared=FALSE) {
 tr_option <- switch(COLOR_CURVE,
  linear=scales::identity_trans(),
  gentle=scales::trans_new('signed_power_065',function(x)sign(x)*abs(x)^.65,function(x)sign(x)*abs(x)^(1/.65)),
  shared=scales::trans_new('asinh_05',function(x)asinh(x/.5),function(x)sinh(x)*.5),
  original=scales::trans_new('asinh_002',function(x)asinh(x/.02),function(x)sinh(x)*.02))
 rates <- sort(unique(c(seq(limits[1],limits[2],length.out=4097),0)))
 col <- flux_colors(rates,max(abs(limits)));stopifnot(!anyNA(col))
 scale_fill_gradientn(colors=col,values=scales::rescale(tr_option$transform(rates),from=tr_option$transform(limits)),transform=tr_option,limits=limits,breaks=breaks,
 labels=function(x)gsub('-','−',format(x,trim=TRUE,scientific=FALSE)),name=label,
 guide=guide_colorbar(direction='horizontal',barwidth=unit(if(shared)135 else 62,'mm'),barheight=unit(3,'mm'),title.position='top',order=order))
}
# 5.4 m half-width contains all model polygons; 23 m contains all height proxies.
# Fixed coordinates use identical scales for all eight trees in every option.
for(i in 1:8)stopifnot(max(abs(polygons(painted[[i]]$branch)$x))<5.4,meta$H[i]<23)
render_option <- function(name,upland_curve='linear',wetland_hue=10,shared=FALSE) {
 plots <- lapply(1:8,function(i) {
  wet <- i==8
  COLOR_CURVE <<- if(shared)'shared'else if(wet && upland_curve!='original')'linear'else upland_curve
  COLOR_HUE <<- if(wet)wetland_hue else 10
  limits <- if(shared)c(-.04,160)else if(wet)c(0,160)else c(-.04,3)
  breaks <- if(shared)c(0,1,10,160)else if(wet)c(0,40,80,120,160)else if(upland_curve=='original')c(-.04,0,.1,.5,1,3)else c(0,.5,1,2,3)
  if(upland_curve=='original'&&wet)breaks<-c(0,1,10,160)
  label <- if(shared)expression(Stem~CH[4]~flux~(nmol~m^{-2}~s^{-1}))else if(wet)expression(atop(Wetland~reference,CH[4]~flux~(nmol~m^{-2}~s^{-1})))else expression(atop(Upland~trees,CH[4]~flux~(nmol~m^{-2}~s^{-1})))
  cs <- option_scale(limits,breaks,label,order=if(wet&&!shared)2 else 1,shared=shared)
  p <- plot_tree(i,ncol=4,color_scale=cs,reference=max(abs(limits)))
  p <- suppressMessages(p+coord_fixed(ratio=1,xlim=c(-5.4,5.4),ylim=c(0,23),expand=FALSE))+
   labs(subtitle=if(wet)'Wetland reference'else NULL)+
   theme(plot.margin=margin(1,1,1,1),plot.title=element_text(size=9,hjust=.5),plot.subtitle=element_text(size=7,hjust=.5))
  p
 })
 p <- wrap_plots(plots,ncol=4,guides='collect') & theme(legend.position='bottom',legend.box='horizontal',legend.text=element_text(size=7),legend.title=element_text(size=8),legend.spacing.x=unit(5,'mm'),legend.margin=margin(2,1,1,1))
 p <- p+plot_annotation(theme=theme(plot.margin=margin(2,2,2,2)))
 for(ext in c('png','pdf')) {
  out <- list(filename=file.path(CANDIDATES,paste0(name,'.',ext)),plot=p,width=190,height=192,units='mm',bg='white')
  if(ext=='png')out$dpi<-300 else out$device<-cairo_pdf
  do.call(ggsave,out)
 }
}
# Retain the current color choice in the document pending a color decision,
# while applying the requested removal of ticks/site subtitle and tighter layout.
render_option('all_stem_flux_sampling_separate_scales','original',10)
render_option('S5_A_linear_red','linear',10)
render_option('S5_B_linear_purple','linear',295)
render_option('S5_C_gentle_purple','gentle',295)
render_option('S5_D_shared_scale','shared',10,TRUE)
cat('PASS: four color options; no sampling ticks or Yale Myers subtitle; shared spatial limits contain every tree; all palette values valid.\n')
