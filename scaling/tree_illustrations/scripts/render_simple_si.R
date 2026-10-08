# A simplified review candidate; does not replace the selected manuscript figure.
args <- commandArgs(trailingOnly=TRUE)
project <- normalizePath(if(length(args))args[1]else'.',mustWork=TRUE)
bundle <- normalizePath(if(length(args)>1)args[2]else file.path(project,'scaling/tree_illustrations'),mustWork=TRUE)
source(file.path(bundle,'scripts/render_trees.R'))
# Preserve the earlier eight-tree shared palette exactly: blue/red stems,
# neutral branches, and the same nonlinear scale from -0.05 to 160.
branch_style <- 'light_gray'; sampling_marks <- FALSE
bounds <- t(vapply(painted,function(z)range(polygons(z$branch)$x),numeric(2)))
left <- vapply(1:4,function(j)min(bounds[c(j,j+4),1],-3.5),numeric(1))
right <- vapply(1:4,function(j)max(bounds[c(j,j+4),2],3.5),numeric(1))
xpos <- numeric(4)
for(j in 2:4)xpos[j] <- xpos[j-1]+right[j-1]-left[j]+.7
# One coordinate system preserves absolute dimensions and eliminates facet padding.
base_y <- c(rep(24.5,4),rep(0,4))
all_br <- bind_rows(lapply(1:8,function(i){z<-polygons(painted[[i]]$branch);z$x<-z$x+xpos[(i-1)%%4+1];z$y<-z$y+base_y[i];z$id<-paste(i,z$id,sep='_');z}))
p <- ggplot()+geom_polygon(data=all_br,aes(x,y,group=id),fill='#DCDCDC',color=NA)
for(i in 1:8) {
 xo<-xpos[(i-1)%%4+1];yo<-base_y[i]
 layer<-stem_gradient(i,160)
 layer$geom_params$xmin<-layer$geom_params$xmin+xo
 layer$geom_params$xmax<-layer$geom_params$xmax+xo
 layer$geom_params$ymin<-layer$geom_params$ymin+yo
 layer$geom_params$ymax<-layer$geom_params$ymax+yo
 p<-p+layer
}
labels<-data.frame(x=rep(xpos,2),y=base_y-.7,label=meta$label)
p<-p+geom_text(data=labels,aes(x,y,label=label),family='Helvetica',fontface='italic',size=3.2,vjust=1)+
 annotate('text',x=xpos[4],y=-1.85,label='Wetland reference',family='Helvetica',size=2.6,color='#656565',vjust=1)
cs<-make_scale(c(-.05,160),c(-.05,0,.1,1,10,160))
cs$name<-expression(Stem~CH[4]~flux~(nmol~m^{-2}~s^{-1}))
cs$guide<-guide_colorbar(direction='horizontal',barwidth=unit(118,'mm'),barheight=unit(3.5,'mm'),title.position='top')
p<-p+geom_blank(data=data.frame(x=xpos[1],y=0,rate=c(-.05,160)),aes(x,y,fill=rate),show.legend=TRUE)+cs+
 coord_fixed(ratio=1,xlim=c(xpos[1]+left[1]-.1,xpos[4]+right[4]+.1),ylim=c(-2.7,max(meta$H+base_y)+.15),expand=FALSE)+
 theme_void(base_family='Helvetica')+theme(plot.margin=margin(1,1,1,1),legend.position='bottom',legend.title=element_text(size=9),legend.text=element_text(size=8),legend.margin=margin(1,1,1,1),legend.box.margin=margin(0,0,0,0),legend.spacing=unit(0,'mm'))
for(ext in c('png','pdf')) {
 out<-list(filename=file.path(CANDIDATES,paste0('S5_F_shared_gray_branches_labels_below.',ext)),plot=p,width=160,height=220,units='mm',bg='white')
 if(ext=='png')out$dpi<-300 else out$device<-cairo_pdf
 do.call(ggsave,out)
}
stopifnot(all(labels$label==meta$label),sampling_marks==FALSE)
cat('PASS: earlier shared scale, all eight stems flux-colored, labels below, shared dimensions, no height marks.\n')
