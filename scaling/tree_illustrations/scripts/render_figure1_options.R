# Optional Figure 1 layouts, kept separate from the production figure pipeline.
# Run from the repository root; an optional second argument selects another bundle.
args <- commandArgs(trailingOnly=TRUE)
PROJECT_ROOT <- normalizePath(if(length(args)) args[1] else '.',mustWork=TRUE)
BUNDLE <- normalizePath(if(length(args)>1) args[2] else file.path(PROJECT_ROOT,'scaling/tree_illustrations'),mustWork=TRUE)
source(file.path(BUNDLE,'scripts/render_trees.R'))
suppressPackageStartupMessages({library(jpeg);library(grid);library(gtable)})
# Grid needs an active device for layout measurements; do not create Rplots.pdf.
grDevices::pdf(file=NULL)
DEST_F1 <- file.path(BUNDLE,'figure1_options')
dir.create(DEST_F1,recursive=TRUE,showWarnings=FALSE)
# Import only definitions and input preparation. Never run the production
# script's figure writes, bootstraps or other analyses.
load_assignments <- function(path,names_keep) {
 for(e in parse(path)) if(is.call(e)&&identical(e[[1]],as.name('<-'))&&is.symbol(e[[2]])&&as.character(e[[2]]) %in% names_keep) eval(e,envir=.GlobalEnv)
}
load_assignments(file.path(PROJECT_ROOT,'scaling/08_figures_v3c.R'),c('FONT','th','COMP','CC','bg_all','bg_raw','bg_col','LAB1','prep1','P1m','P1s','ASC','tr_as','sq','prof','img','pd'))
# Habitat is now above the panels; remove the redundant '(swamp)' subtitle.
LAB1$lab[8] <- 'italic("N. sylvatica")~""'
P1m <- prep1(bg_col)
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
profile <- function(background=FALSE) {
 p <- prof(P1m,'row','raw')
 p <- suppressMessages(p+scale_y_continuous(trans=sq,limits=c(0,24),breaks=c(0,5,10,15,20)))
 if(background) {
  cr <- list();tk <- list();alignment <- list()
  for(i in 1:8) {
   di <- filter(P1m,tree==meta$tree[i])
   left <- min(-.12*max(abs(di$flux)),min(di$flux));right <- max(di$flux)
   z <- shape_data(i,(left+right)/2,.82*(right-left)/max(meta$crown_width_m))
   z$crown$lab <- factor(LAB1$lab[i],levels=LAB1$lab); z$stem$lab <- z$crown$lab[1]
   cr[[i]] <- z$crown; tk[[i]] <- z$stem
   obs <- filter(D,tree==meta$tree[i],component %in% c('leaf','branch'))
   stopifnot(min(z$crown$y)>=0,abs(max(z$crown$y)-meta$H[i])<1e-10,all(z$crown$x>=left),all(z$crown$x<=right),!nrow(obs)||min(obs$height_m)>=z$lower-1e-8)
   alignment[[i]] <- tibble(tree=meta$tree[i],top_m=meta$H[i],crown_bottom_m=z$lower,lowest_observed_branch_or_leaf_m=if(nrow(obs))min(obs$height_m)else NA_real_,horizontal_scale='schematic within flux panel')
  }
  layers <- list(geom_polygon(data=bind_rows(cr),aes(x,y,group=id),inherit.aes=FALSE,fill='#F0F2EE',color=NA),geom_polygon(data=bind_rows(tk),aes(x,y,group=id),inherit.aes=FALSE,fill='#E0E5DD',color=NA))
  p$layers <- append(p$layers,layers,after=1)
  write_csv(bind_rows(alignment),file.path(DEST_F1,'silhouette_alignment.csv'))
 }
 p
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
bottom <- (img('IMG_5926_edited.jpg',.5,.45)+labs(tag='b') | img('IMG_6437.jpg',.5,.5)+labs(tag='c') | pd('raw')+labs(tag='d'))+plot_layout(widths=c(1,1,.9))
for(bg in c(FALSE,TRUE)) {
 name <- if(bg)'F1_B_faint_aligned' else 'F1_A_data_only'
 pg <- profile(bg)
 top <- wrap_elements(full=habitat_header(pg))+labs(tag='a')
 full <- (top/bottom)+plot_layout(heights=c(1,.42))+plot_annotation(theme=theme(plot.margin=margin(6,6,6,6)))
 save_figure(full,name)
 save_figure(top,paste0(name,'_panel_a'),height=147)
 # A matched third-row option for both backgrounds makes the tradeoff visible.
 third <- wrap_elements(full=patchworkGrob(upland_candidate))+labs(tag='e')
 combined <- (top/bottom/third)+plot_layout(heights=c(1,.42,.60))+plot_annotation(theme=theme(plot.margin=margin(6,6,6,6)))
 extra <- if(bg)'F1_D_faint_plus_tree_row' else 'F1_C_data_plus_tree_row'
 save_figure(combined,extra,height=270)
 save_figure(top,paste0(extra,'_panel_a'),height=147)
}
# Test the actual built panels: silhouettes cannot change flux axes or points.
a <- ggplot_build(profile(FALSE));b <- ggplot_build(profile(TRUE))
stopifnot(isTRUE(all.equal(lapply(a$layout$panel_scales_x,function(x)x$range$range),lapply(b$layout$panel_scales_x,function(x)x$range$range))))
# Two extra background layers follow the basal band; all original layers match.
for(j in seq_along(a$data)) stopifnot(isTRUE(all.equal(a$data[[j]],b$data[[if(j==1)j else j+2]],check.attributes=FALSE)))
cat('PASS: habitat spans 7 upland + 1 wetland panels; simple crowns cover observed branch/leaf heights; raw data, curves and flux scales unchanged.\n')
