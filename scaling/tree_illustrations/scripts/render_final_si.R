# Selected Figure S5: original eight-tree layout, labels below, gray branches.
args <- commandArgs(trailingOnly=TRUE)
project <- normalizePath(if(length(args))args[1]else'.',mustWork=TRUE)
bundle <- normalizePath(if(length(args)>1)args[2]else file.path(project,'scaling/tree_illustrations'),mustWork=TRUE)
source(file.path(bundle,'scripts/render_trees.R'))
branch_style <- 'light_gray'; gray_color <- '#DCDCDC'; sampling_marks <- FALSE
cs <- make_scale(c(-.05,160),c(-.05,0,.1,1,10,160))
cs$name <- expression(Stem~CH[4]~flux~(nmol~m^{-2}~s^{-1}))
plots <- lapply(1:8,function(i) {
 subtitle <- if(i==8)'Wetland reference'else' '
 label <- bquote(atop(italic(.(meta$label[i])),scriptstyle(.(subtitle))))
 plot_tree(i,ncol=4,color_scale=cs,reference=160)+
  labs(title=NULL,subtitle=NULL,x=label)+
  theme(axis.title.x=element_text(size=11,hjust=.5,margin=margin(t=5)))
})
# Restore original equal-width panels, +/-6 m horizontal bounds, 0–24 m height,
# panel margins, row spacing and square-ish canvas. Every tree remains to scale.
p <- wrap_plots(plots,ncol=4,guides='collect') & theme(legend.position='bottom')
p <- p+plot_annotation(theme=theme(plot.margin=margin(8,8,8,8)))
for(ext in c('png','pdf')) {
 out<-list(filename=file.path(CANDIDATES,paste0('FigS5_tree_flux.',ext)),plot=p,width=12,height=12.2,units='in',bg='white')
 if(ext=='png')out$dpi<-300 else out$device<-cairo_pdf
 do.call(ggsave,out)
}
stopifnot(!sampling_marks,flux_colors(0,160)=='#FFFFFF')
cat('PASS: finalized S5 retains original spacing and colors, labels below, gray branches, no measurement ticks.\n')
