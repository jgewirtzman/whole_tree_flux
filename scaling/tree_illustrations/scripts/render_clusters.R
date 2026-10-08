# Neutral structural illustrations and simple/detailed genus clusters for slides.
suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(readr);library(patchwork)})
args <- commandArgs(trailingOnly=TRUE)
PROJECT_ROOT <- normalizePath(if(length(args))args[1]else'.',mustWork=TRUE)
BUNDLE <- normalizePath(if(length(args)>1)args[2]else file.path(PROJECT_ROOT,'scaling/tree_illustrations'),mustWork=TRUE)
DEST <- file.path(BUNDLE,'presentations');dir.create(DEST,recursive=TRUE,showWarnings=FALSE)
meta <- read_csv(file.path(BUNDLE,'inputs/individual_crown_inputs.csv'),show_col_types=FALSE)
st <- read_csv(file.path(BUNDLE,'inputs/modeled_stem_frusta.csv'),show_col_types=FALSE)
br <- read_csv(file.path(BUNDLE,'inputs/modeled_branch_frusta.csv'),show_col_types=FALSE)
models <- lapply(meta$tree,function(id)list(stem=filter(st,tree==id),branch=filter(br,tree==id)))
for(e in parse(file.path(BUNDLE,'scripts/render_trees.R')))if(is.call(e)&&identical(e[[1]],as.name('<-'))&&identical(e[[2]],as.name('polygons')))eval(e)
for(e in parse(file.path(BUNDLE,'scripts/render_figure1_options.R')))if(is.call(e)&&identical(e[[1]],as.name('<-'))&&identical(e[[2]],as.name('shape_data')))eval(e)
poly <- function(z,cx=0) {z$rate<-0;p<-polygons(z);p$x<-p$x+cx;p}
neutral <- function(i) {
 m<-meta[i,];z<-models[[i]]
 ggplot()+geom_polygon(data=poly(z$branch),aes(x,y,group=id),fill='#45574A',color=NA)+geom_polygon(data=poly(z$stem),aes(x,y,group=id),fill='#45574A',color=NA)+coord_fixed(xlim=c(-6,6),ylim=c(0,24),expand=FALSE)+labs(title=m$label,subtitle=if(m$site=='Harvard Forest')NULL else m$site)+theme_void(base_family='Helvetica')+theme(plot.title=element_text(face='italic',size=11,hjust=.5),plot.subtitle=element_text(size=8,hjust=.5),plot.margin=margin(8,4,8,4))
}
save <- function(p,name,w,h) {
 ggsave(file.path(DEST,paste0(name,'.png')),p,width=w,height=h,dpi=210,bg='white')
 ggsave(file.path(DEST,paste0(name,'.pdf')),p,width=w,height=h,device=cairo_pdf,bg='white')
}
save(wrap_plots(lapply(1:8,neutral),ncol=4),'T_structure',12,12)
for(style in c('simple','complex')) {
 p<-ggplot();ids<-c(1,3,4,6);cx<-c(0,4,8,12)
 for(j in seq_along(ids)) {
  i<-ids[j];z<-models[[i]]
  if(style=='simple') {q<-shape_data(i,cx[j],1);crown<-q$crown;stem<-q$stem} else {crown<-poly(z$branch,cx[j]);stem<-poly(z$stem,cx[j])}
  p<-p+geom_polygon(data=crown,aes(x,y,group=id),fill=c('#BBC9B6','#CAD3C5','#AEBFA8','#D2DACD')[j],color=NA)+geom_polygon(data=stem,aes(x,y,group=id),fill='#859A7D',color=NA)
 }
 p<-p+coord_fixed(xlim=c(-5,17),ylim=c(0,24),expand=FALSE)+theme_void()+theme(plot.margin=margin(10,10,10,10))
 save(p,paste0('genus_cluster_',style),6,6.5)
}
