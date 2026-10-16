# Data-only Figure 1 with a faint genus cluster beside panel a.
suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(readr);library(patchwork);library(jpeg);library(grid);library(gtable)})
args <- commandArgs(trailingOnly=TRUE)
PROJECT_ROOT <- normalizePath(if(length(args))args[1]else'.',mustWork=TRUE)
BUNDLE <- normalizePath(if(length(args)>1)args[2]else file.path(PROJECT_ROOT,'scaling/tree_illustrations'),mustWork=TRUE)
setwd(PROJECT_ROOT);source('scaling/00_load_field.R');ROOT <- PROJECT_ROOT
grDevices::pdf(file=NULL)
DEST_F1 <- file.path(BUNDLE,'figure1_options')
meta <- read_csv(file.path(BUNDLE,'inputs/individual_crown_inputs.csv'),show_col_types=FALSE)
st <- read_csv(file.path(BUNDLE,'inputs/modeled_stem_frusta.csv'),show_col_types=FALSE)
br <- read_csv(file.path(BUNDLE,'inputs/modeled_branch_frusta.csv'),show_col_types=FALSE)
models <- lapply(meta$tree,function(id)list(stem=filter(st,tree==id),branch=filter(br,tree==id)))
load_assignments <- function(path,names_keep) {
 for(e in parse(path)) if(is.call(e)&&identical(e[[1]],as.name('<-'))&&is.symbol(e[[2]])&&as.character(e[[2]]) %in% names_keep) eval(e,envir=.GlobalEnv)
}
load_assignments(file.path(ROOT,'scaling/08_figures.R'),c('FONT','th','COMP','CC','bg_all','bg_raw','bg_col','LAB1','prep1','P1m','ASC','tr_as','sq','prof','img','pd'))
LAB1$lab[8] <- 'italic("N. sylvatica")~""';P1m <- prep1(bg_col)
NS <- read_csv(file.path(ROOT,'scaling/out/fig1_heightbin_means.csv'),show_col_types=FALSE)
NS$hb <- factor(NS$hb,levels=paste0(seq(0,22,2),'–',seq(2,24,2)))
NH <- P1m %>% filter(site!='Black Gum Swamp',component=='stem') %>% mutate(hb=cut(height_m,seq(0,24,2),right=FALSE,labels=levels(NS$hb)))
# Simple filled crowns retain the latest modeled height and lower crown extent.
# Horizontal positions are decorative within free flux axes, never spatial data.
shape_data <- function(i,cx,hscale) {
 m <- meta[i,]; br <- models[[i]]$branch; st <- models[[i]]$stem
 low <- min(m$modeled_crown_base_m,br$y0,br$y1)
 H <- m$H; w <- m$crown_width_m/2*hscale
 if(m$species=='hem') {
  # A single stepped conifer outline, with no overlapping translucent polygons.
  t <- sort(unique(c(0,seq(.06,.94,length.out=9),seq(.06,.94,length.out=9)+.025,1)))
  t <- t[t<=1]; breadth <- (1-t)^.85*(.86+.14*(seq_along(t)%%2))
  breadth[1] <- .35; breadth[length(breadth)] <- 0
 } else {
  t <- seq(0,1,length.out=220)
  breadth <- sin(pi*t)^.38*(1+.025*sin(11*pi*t))
  if(m$species=='bg') breadth <- breadth*(1-.18*t)
  breadth <- breadth/max(breadth)
 }
 crown <- tibble(x=c(cx-w*breadth,rev(cx+w*breadth)),y=c(low+(H-low)*t,rev(low+(H-low)*t)),id=i)
 h <- c(st$y0,tail(st$y1,1)); r <- c(st$r0,tail(st$r1,1))*hscale
 # Slight common enlargement keeps the illustrative stems visible at panel size.
 stem <- tibble(x=c(cx-r*1.6,rev(cx+r*1.6)),y=c(h,rev(h)),id=i)
 list(crown=crown,stem=stem,lower=low)
}
profile <- function() {
 p <- prof(P1m,'row','raw')
 suppressMessages(p+scale_y_continuous(trans=sq,limits=c(0,24),breaks=c(0,5,10,15,20)))
}
# Brackets span exact facet-panel columns, excluding the y-axis and outer margin.
habitat_header <- function(p) {
 gt <- ggplotGrob(p)
 panels <- gt$layout[grepl('^panel',gt$layout$name),]; panels <- panels[order(panels$l),]
 stopifnot(nrow(panels)==8)
 gt <- gtable_add_rows(gt,unit(7,'mm'),pos=0)
 for(ids in list(1:7,8)) {
  lab <- if(length(ids)==7)'Upland' else 'Wetland'
  g <- grobTree(segmentsGrob(x0=unit(c(0,0,1),'npc'),x1=unit(c(1,0,1),'npc'),y0=unit(c(.18,.18,.18),'npc'),y1=unit(c(.18,.04,.04),'npc'),gp=gpar(col='#737373',lwd=.6)),textGrob(lab,x=.5,y=.64,gp=gpar(fontfamily=FONT,fontsize=8,col='#404040')))
  gt <- gtable_add_grob(gt,g,t=1,l=min(panels$l[ids]),r=max(panels$r[ids]),clip='off',name=paste0('habitat-',lab))
 }
 gt
}
save_figure <- function(p,name,width=190,height=205) {
 ggsave(file.path(DEST_F1,paste0(name,'.png')),p,width=width,height=height,units='mm',dpi=300,bg='white',limitsize=FALSE)
 ggsave(file.path(DEST_F1,paste0(name,'.pdf')),p,width=width,height=height,units='mm',device=cairo_pdf,bg='white',limitsize=FALSE)
}
bottom <- (img('data/photos/field/IMG_5926_edited.jpg',.5,.45)+labs(tag='b') | img('data/photos/field/IMG_6437.jpg',.5,.5)+labs(tag='c') | pd('raw')+labs(tag='d'))+plot_layout(widths=c(1,1,.9))
# Use the exact transformed panel range including its expansion, not an
# independently laid-out ggplot. This aligns y=0 and every other height exactly.
p <- profile()+theme(strip.text=element_text(size=6.2))
built <- ggplot_build(p);yrange <- built$layout$panel_params[[1]]$y.range
stopifnot(all(vapply(built$layout$panel_params,function(z)isTRUE(all.equal(z$y.range,yrange)),logical(1))))
y_npc <- function(y)(sq$transform(y)-yrange[1])/diff(yrange)
ids <- c(1,3,4,6);centers <- c(.235,.43,.63,.80)
children <- list();check <- list()
for(j in seq_along(ids)) {
 i <- ids[j];z <- shape_data(i,centers[j],.052)
 stopifnot(all(z$crown$x>=0&z$crown$x<=1),min(z$stem$y)==0)
 children <- c(children,list(polygonGrob(z$crown$x,y_npc(z$crown$y),default.units='npc',gp=gpar(fill=c('#EEF1EA','#E7ECE3','#F0F2ED','#E4EAE0')[j],col=NA)),polygonGrob(z$stem$x,y_npc(z$stem$y),default.units='npc',gp=gpar(fill='#D6DFD0',col=NA))))
 check[[j]] <- tibble(tree=meta$tree[i],base_m=0,base_panel_fraction=y_npc(0),height_proxy_m=meta$H[i],top_panel_fraction=y_npc(meta$H[i]),height_axis='same compression above 15 m and same expansion as data')
}
cluster <- do.call(grobTree,children)
# Figure 1: silhouettes to the left of the profiles (layout E). Set
# EXPORT_ILLUSTRATION_ALTERNATIVES=1 to also write the right-side and reordered layouts.
ALTERNATIVES <- Sys.getenv('EXPORT_ILLUSTRATION_ALTERNATIVES','0')=='1'
if(ALTERNATIVES) dir.create(DEST_F1,showWarnings=FALSE)
for(side in if(ALTERNATIVES) c('left','right') else 'left') {
 gt <- habitat_header(p)
 # Insert outside the complete axes so tree silhouettes cannot cover the data.
 if(side=='left') {
  gt <- gtable_add_cols(gt,unit(22,'mm'),pos=0);col <- 1
 } else {
  gt <- gtable_add_cols(gt,unit(22,'mm'),pos=-1);col <- length(gt$widths)
 }
 panels <- gt$layout[grepl('^panel',gt$layout$name),]
 stopifnot(length(unique(panels$t))==1,length(unique(panels$b))==1)
 gt <- gtable_add_grob(gt,cluster,t=panels$t[1],b=panels$b[1],l=col,r=col,clip='on',name=paste0('aligned-silhouettes-',side))
 top <- wrap_elements(full=gt)+labs(tag='a')
 full <- (top/bottom)+plot_layout(heights=c(1,.42))+plot_annotation(theme=theme(plot.margin=margin(6,6,6,6)))
 if(side=='left') {
  for(ext in c('png','pdf')) {
   out <- list(filename=file.path(PROJECT_ROOT,'figures',paste0('Figure1.',ext)),plot=full,width=190,height=205,units='mm',bg='white',limitsize=FALSE)
   if(ext=='png') out$dpi <- 300 else out$device <- cairo_pdf
   do.call(ggsave,out)
  }
 }
 if(!ALTERNATIVES) next
 name <- paste0('F1_',if(side=='left')'E'else'F','_faint_side_',side)
 save_figure(full,name)
 save_figure(top,paste0(name,'_panel_a'),height=147)
 # Reorder the second-row panels: summary first, followed by both photographs.
 swapped_bottom <- (pd('raw')+labs(tag='b') | img('data/photos/field/IMG_5926_edited.jpg',.5,.45)+labs(tag='c') | img('data/photos/field/IMG_6437.jpg',.5,.5)+labs(tag='d'))+plot_layout(widths=c(.9,1,1))
 swapped <- (top/swapped_bottom)+plot_layout(heights=c(1,.42))+plot_annotation(theme=theme(plot.margin=margin(6,6,6,6)))
 swapped_name <- paste0('F1_',if(side=='left')'G'else'H','_summary_first_',side)
 save_figure(swapped,swapped_name)
 save_figure(top,paste0(swapped_name,'_panel_a'),height=147)
}
write_csv(bind_rows(check),file.path(BUNDLE,'checks','figure1_silhouette_alignment.csv'))
cat('PASS: both side clusters occupy the exact data-panel row; tree bases map to the same y=0; all heights use the same transformed range. Data-only profiles are reused unchanged.\n')
