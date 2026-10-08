# Shared geometry, extrapolation and resampling rules for the current analysis.
# Run from the repository root. No files are written by this module.
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
N_BOOT <- as.integer(Sys.getenv("FLUX_N_BOOT", "2000"))
stopifnot(is.finite(N_BOOT), N_BOOT >= 100)

# Illustrative eastern North American temperate closed-canopy deciduous forests.
# Whittaker & Woodwell (1967), p. 937: stem 0.5–0.6 and branch 1.5–1.6
# m2 per m2 ground; choose the midpoints of the ranges reported together.
# Leaf index retains the midpoint of their abstract's 3–6 range.
# These are fixed scenario inputs, not measured areas of the sampled trees.
STAND_AREA <- c(stem = 0.55, branch = 1.55, leaf = 4.5)
ALTERNATIVE_WOODY_AREA <- 3.07 # Gauci et al. (2024); partition is our assumption

# Representative whole-stem heights, not measured heights of our sampled trees.
# Sullivan et al. (2017), Table 6: rounded upper-canopy proxy and Nyssa proxy.
CANOPY_HEIGHT <- 23
# Species means of individual-crown maximum lidar heights, uncorrected column.
# Use for individual-tree scenarios when a total height is needed; measured
# chamber positions remain measurements, not estimates of the tree's total height.
SULLIVAN_HEIGHT <- c(hem=20.86, rm=20.03, ro=22.88, bg=15.81, qv=22.64)
BLACKGUM_HEIGHT <- round(unname(SULLIVAN_HEIGHT['bg']),1)
cone_share_below <- function(h, H) {
  stopifnot(all(is.finite(H)), all(H > 0))
  1 - (1 - pmin(pmax(h, 0), H)/H)^2
}

load_geometry <- function(F) {
  hf <- read.csv("data processing/goFlux_reprocessing/results/canopy_flux_goFlux_compiled_with_mdf.csv")
  ym <- read.csv("data processing/goFlux_reprocessing/ymf_black_oak/results/ymf_field_data.csv")
  g <- bind_rows(
    hf %>% filter(Type == "stem") %>% transmute(tree = paste(Site, Species, Tree_Tag), h = Height_m,
      d = suppressWarnings(as.numeric(DBH_cm)) / 100),
    ym %>% transmute(tree = "YMF black oak", h = suppressWarnings(as.numeric(stringr::str_extract(as.character(Height_m), "^[0-9.]+"))),
      d = Stem_Diam_mm / 100)) # original column name says mm; field values are cm
  g <- g %>% filter(tree %in% F$tree, is.finite(h), is.finite(d), d > 0) %>%
    group_by(tree, h) %>% summarise(d = mean(d), .groups = "drop")
  tops <- F %>% filter(component == "stem") %>% group_by(tree) %>% summarise(top = max(height_m), .groups = "drop")
  g <- left_join(g, tops, by = "tree")
  stopifnot(setequal(unique(g$tree), unique(tops$tree)))
  g
}

# Frustum surface density pi*d*sqrt(1+(d'/2)^2). Composite Simpson quadrature
# splits at every observed flux/diameter height, so irregular chamber spacing is
# represented by physical interval length. Tops stop at the highest stem chamber;
# no unmeasured zero-diameter tip is invented. Missing end diameters use the nearest
# measured diameter (explicit boundary assumption, including missing basal values).
stem_grid <- function(g, lower, upper, flux_heights = numeric(), max_step = 0.05) {
  stopifnot(nrow(g) >= 2, upper > lower)
  knots <- sort(unique(c(lower, upper, g$h[g$h > lower & g$h < upper],
    flux_heights[flux_heights > lower & flux_heights < upper])))
  edges <- sort(unique(unlist(lapply(seq_len(length(knots)-1), function(i)
    seq(knots[i], knots[i+1], length.out = ceiling((knots[i+1]-knots[i])/max_step)+1)))))
  left <- head(edges,-1); right <- tail(edges,-1); mid <- (left+right)/2
  dl <- approx(g$h,g$d,xout=left,rule=2)$y; dr <- approx(g$h,g$d,xout=right,rule=2)$y
  fac <- pi * sqrt(1 + ((dr-dl)/(2*(right-left)))^2) * (right-left)/6
  data.frame(h = c(left,mid,right), area = c(fac*dl, 4*fac*(dl+dr)/2, fac*dr))
}
profile_weights <- function(heights, grid) {
  heights <- sort(unique(heights)); stopifnot(length(heights) >= 1)
  if (length(heights)==1) return(1)
  vapply(seq_along(heights), function(i) sum(grid$area *
    approx(heights, as.numeric(seq_along(heights)==i), xout=grid$h, rule=2)$y)/sum(grid$area), numeric(1))
}
height_means <- function(d) {
  a <- aggregate(flux ~ height_m, d, mean)
  a[order(a$height_m),]
}
pooled_slope <- function(low) {
  pos <- low %>% group_by(tree) %>% filter(all(flux > 0), n() >= 2) %>% ungroup()
  if (!nrow(pos)) return(NA_real_)
  if (n_distinct(pos$tree)==1) return(unname(coef(lm(log(flux) ~ height_m, pos))[2]))
  unname(coef(lm(log(flux) ~ height_m + tree, pos))[2])
}
extrapolation_forms <- function(lo, b_pool = NA_real_) {
  h<-lo$height_m; f<-lo$flux; stopifnot(length(unique(h))>=2)
  top<-f[which.max(h)]; lin<-coef(lm(f~h)); rising<-lin[2]>0
  ex<-if(all(f>0)) coef(lm(log(f)~h)) else c(NA,NA)
  list(const_mean=function(x) rep(mean(f),length(x)), const_top=function(x) rep(top,length(x)),
    exp_decay=function(x) if(anyNA(ex)) rep(NA_real_,length(x)) else if(ex[2]>0) rep(top,length(x)) else exp(ex[1]+ex[2]*x),
    exp_pooled=function(x) if(is.na(b_pool)||b_pool>0) rep(top,length(x)) else top*exp(b_pool*(x-max(h))),
    linear=function(x) if(rising) rep(top,length(x)) else lin[1]+lin[2]*x,
    linear_zero=function(x) pmax(if(rising) rep(top,length(x)) else lin[1]+lin[2]*x,0),
    zero=function(x) rep(0,length(x)))
}
# Resample entire trees first, then closures within each tree/component/height.
# Duplicate selected trees receive distinct IDs; no group_by operation can collapse them.
resample_trees <- function(d, trees = unique(d$tree)) {
  ids <- sample(trees, length(trees), replace=TRUE)
  bind_rows(lapply(seq_along(ids),function(i) {
    x<-d[d$tree==ids[i],,drop=FALSE]; x$source_tree<-ids[i]; x$tree<-paste0("draw_",i)
    groups<-interaction(x$component,x$height_m,drop=TRUE)
    ii<-unlist(lapply(split(seq_len(nrow(x)),groups),function(j) j[sample.int(length(j),length(j),replace=TRUE)]),use.names=FALSE)
    x[ii,,drop=FALSE]
  }))
}
# Equal tree weighting, equal sampled-height weighting within tree (nonintegrated summaries).
tree_mean <- function(d) {
  d %>% group_by(tree,height_m) %>% summarise(f=mean(flux),.groups="drop") %>%
    group_by(tree) %>% summarise(f=mean(f),.groups="drop")
}
ci95 <- function(x) {
  if (!any(is.finite(x))) return(c(lo=NA_real_,hi=NA_real_))
  setNames(as.numeric(quantile(x[is.finite(x)],c(.025,.975))),c("lo","hi"))
}
