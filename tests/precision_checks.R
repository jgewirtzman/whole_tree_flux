# Scientific regression tests: trend invariance, gap boundaries, grouping,
# and consistency of the noise threshold used in fitting and reporting.
source("data processing/goFlux_reprocessing/precision_helpers.R")
set.seed(1207)
t <- 0:9999; e <- rnorm(length(t),sd=2)
a <- closure_noise(e,t); b <- closure_noise(e+20*t,t)
stopifnot(abs(a$sigma/2-1)<.05, abs(a$sigma-b$sigma)<1e-9, abs(a$ac1+2/3)<.03)
# A gap and an arbitrary jump cannot contribute a difference across the break.
g <- closure_noise(c(e[1:5000],1e6+e[5001:10000]),c(t[1:5000],t[5001:10000]+100))
stopifnot(g$n_diff2==length(t)-4,abs(g$sigma/a$sigma-1)<.01)
# Missing observations break adjacency too.
x<-e;x[5000]<-NA;g<-closure_noise(x,t)
stopifnot(g$n_diff2==length(t)-5)
# Large changes between slopes and field days must not inflate pooled precision.
t0<-as.POSIXct("2023-07-18",tz="UTC")
d<-do.call(rbind,lapply(1:8,function(i) {
 z<-0:999; day<-if(i<=4)0 else 1; noise<-if(day==0)1 else 4
 data.frame(UniqueID=paste0("id",i),flag=1,POSIX.time=t0+day*86400+(i%%4)*11000+z*10,
 CH4dry_ppb=2000+100*i*z+rnorm(1000,sd=noise),CO2dry_ppm=400+rnorm(1000,sd=noise))
}))
p<-estimate_precision(d,"test",warn=FALSE)
stopifnot(nrow(p$groups)==4,all(p$groups$dt_s==10),
 all(abs(p$groups$sigma[p$groups$field_day=="2023-07-18"]-1)<.1),
 all(abs(p$groups$sigma[p$groups$field_day=="2023-07-19"]-4)<.3))
# Real fitting smoke check: the adapter matches unmodified goFlux at 1 Hz and
# uses physical seconds for slow logging.
library(goFlux)
for(dt in c(1,5,10)) {
 z<-0:99
 w<-data.frame(UniqueID="test",flag=1,Etime=z*dt,CH4dry_ppb=2000+.05*z*dt+rnorm(100),
 H2O_ppm=10000,CH4_prec=1.96,Vtot=.2,Area=100,Pcham=101.3,Tcham=20)
 r<-suppressWarnings(goFlux_at_interval(w,"CH4dry_ppb"))
 stopifnot(abs(r$MDF-1.96/(100*dt)*r$flux.term)<1e-10)
 if(dt==1) {
  ref<-suppressWarnings(goFlux(w,"CH4dry_ppb"))
  stopifnot(isTRUE(all.equal(ref,r)))
 }
}
cat("PASS: empirical precision, trend/gap/missing-data handling, day grouping and goFlux duration adapter\n")
