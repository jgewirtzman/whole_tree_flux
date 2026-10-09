# Four unaltered field photographs, arranged to preserve their aspect ratios.
# Run from the repository root: Rscript scaling/14_chamber_photos.R
library(grid)
files <- c('IMG_5103.png','IMG_5018.png','IMG_5186.png','IMG_5115.png')
photos <- lapply(files,function(n) png::readPNG(file.path('data processing/chamber_photos',n)))
draw <- function() {
  grid.newpage()
  # The portrait view fills the left; three landscape details stack at right.
  boxes <- list(c(.31,.5,.576,.96),c(.803,.834,.346,.324),
                c(.803,.5,.346,.324),c(.803,.166,.346,.324))
  for(i in seq_along(photos)) {
    b<-boxes[[i]]
    grid.raster(photos[[i]],x=b[1],y=b[2],width=b[3],height=b[4],interpolate=TRUE)
    grid.rect(x=b[1]-b[3]/2+.019,y=b[2]+b[4]/2-.025,width=.031,height=.043,
              gp=gpar(fill='white',col=NA))
    grid.text(letters[i],x=b[1]-b[3]/2+.019,y=b[2]+b[4]/2-.025,
              gp=gpar(fontsize=15,fontface='bold',col='black'))
  }
}
png('scaling/v3c/FigS1_chambers.png',width=8,height=6.4,units='in',res=300)
draw();invisible(dev.off())
pdf('scaling/v3c/FigS1_chambers.pdf',width=8,height=6.4)
draw();invisible(dev.off())
