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
load_assignments(file.path(ROOT,'scaling/08_figures_v3c.R'),c('FONT','th','COMP','CC','bg_all','bg_raw','bg_col','LAB1','prep1','P1m','ASC','tr_as','sq','prof','img','pd'))
LAB1$lab[8] <- 'italic("N. sylvatica")~""';P1m <- prep1(bg_col)
NS <- read_csv(file.path(ROOT,'scaling/out/fig1_heightbin_means.csv'),show_col_types=FALSE)
NS$hb <- factor(NS$hb,levels=paste0(seq(0,22,2),'–',seq(2,24,2)))
NH <- P1m %>% filter(site!='Black Gum Swamp',component=='stem') %>% mutate(hb=cut(height_m,seq(0,24,2),right=FALSE,labels=levels(NS$hb)))
load_assignments(file.path(BUNDLE,'scripts/render_figure1_options.R'),c('shape_data','profile','habitat_header','save_figure','bottom'))
# Use the exact transformed panel range including its expansion, not an
# independently laid-out ggplot. This aligns y=0 and every other height exactly.
p <- profile(FALSE)+theme(strip.text=element_text(size=6.2))
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
for(side in c('left','right')) {
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
 name <- paste0('F1_',if(side=='left')'E'else'F','_faint_side_',side)
 save_figure(full,name)
 save_figure(top,paste0(name,'_panel_a'),height=147)
 # Reorder the second-row panels: summary first, followed by both photographs.
 swapped_bottom <- (pd('raw')+labs(tag='b') | img('IMG_5926_edited.jpg',.5,.45)+labs(tag='c') | img('IMG_6437.jpg',.5,.5)+labs(tag='d'))+plot_layout(widths=c(.9,1,1))
 swapped <- (top/swapped_bottom)+plot_layout(heights=c(1,.42))+plot_annotation(theme=theme(plot.margin=margin(6,6,6,6)))
 swapped_name <- paste0('F1_',if(side=='left')'G'else'H','_summary_first_',side)
 save_figure(swapped,swapped_name)
 save_figure(top,paste0(swapped_name,'_panel_a'),height=147)
}
write_csv(bind_rows(check),file.path(DEST_F1,'side_silhouette_alignment.csv'))
cat('PASS: both side clusters occupy the exact data-panel row; tree bases map to the same y=0; all heights use the same transformed range. Data-only profiles are reused unchanged.\n')
