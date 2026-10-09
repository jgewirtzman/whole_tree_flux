# Leaf areas are assigned to shoots, then shared by all closures of that shoot.
# Original UniqueIDs remain join keys; field columns define tissue and height.
leaf_area_assignments <- function(root = '.') {
  p <- file.path(root, 'data processing')
  f <- read.csv(file.path(p, 'Field Data Entry - Clean Canopy Lift Total.csv'), stringsAsFactors=FALSE)
  a <- read.csv(file.path(p, 'leaf_areas.csv'), fileEncoding='UTF-8-BOM', stringsAsFactors=FALSE)
  names(a) <- sub('^\\.+', '', trimws(gsub('[^A-Za-z0-9_]', '.', names(a))))
  a$Tree <- trimws(a$Tree)
  map <- read.csv(file.path(p, 'leaf_area_mapping.csv'), stringsAsFactors=FALSE)
  stopifnot(!anyDuplicated(f$UniqueID), !anyDuplicated(map$UniqueID),
            !anyDuplicated(paste(a$Tree,a$Flux.Rep)), all(is.finite(a$Area)), all(a$Area>0))
  l <- f[grepl('leaf', f$Type, ignore.case=TRUE), ]
  stopifnot(setequal(map$UniqueID,l$UniqueID), all(nzchar(map$ShootID)))
  map <- map[match(l$UniqueID,map$UniqueID),]
  sp <- c(bg='Blackgum', hem='Hemlock', rm='Maple', ro='Oak')
  species <- unname(sp[l$Species]); stopifnot(!anyNA(species), all(species==map$Tree))
  av <- tapply(a$Area, a$Tree, mean)
  area <- unname(av[species]); method <- rep('species_mean',nrow(l))
  for(i in which(!is.na(map$Flux.Rep))) {
    k <- which(a$Tree==map$Tree[i] & a$Flux.Rep==map$Flux.Rep[i])
    stopifnot(length(k)==1L)
    area[i] <- a$Area[k]; method[i] <- 'measured'
  }
  stopifnot(all(is.finite(area)), all(area>0))
  for (s in unique(map$ShootID)) {
    j <- which(map$ShootID==s)
    if(length(unique(paste(species[j],map$Flux.Rep[j],area[j],method[j])))!=1L)
      stop('Inconsistent area assignment within shoot: ',s)
  }
  data.frame(UniqueID=l$UniqueID, ShootID=map$ShootID,
    Leaf_area_cm2=area, Leaf_area_source=method, Light_condition=map$Light_condition)
}

override_leaf_area <- function(df, root='.') {
  a <- leaf_area_assignments(root)
  leaf <- grepl('leaf',df$Type,ignore.case=TRUE)
  i <- match(df$UniqueID,a$UniqueID)
  if(any(leaf & is.na(i)) || any(!leaf & !is.na(i)))
    stop('Leaf classification differs between timing keys and field metadata')
  df$Area[leaf] <- a$Leaf_area_cm2[i[leaf]]
  df
}

validate_canopy_metadata <- function(root='.') {
  p <- file.path(root,'data processing')
  f <- read.csv(file.path(p,'Field Data Entry - Clean Canopy Lift Total.csv'))
  tissue <- function(x) ifelse(grepl('leaf',x,ignore.case=TRUE),'leaf',tolower(trimws(x)))
  for(k in 1:3) {
    d <- read.csv(file.path(p,'input',paste0('times_key_tree - Canopy Lift_LGR',k,' (2).csv')))
    i <- match(d$UniqueID,f$UniqueID)
    stopifnot(!anyNA(i),!anyDuplicated(d$UniqueID))
    bad <- tissue(d$Type)!=tissue(f$Type[i]) | abs(d$Height-f$Height_m[i])>1e-8
    if(anyNA(bad) || any(bad)) stop('Timing/field tissue or height mismatch: ',paste(d$UniqueID[which(bad)],collapse=', '))
  }
  invisible(TRUE)
}

# Scaling-only runs must not silently use results predating a metadata correction.
validate_compiled_canopy <- function(root='.') {
  p <- file.path(root,'data processing')
  f <- read.csv(file.path(p,'Field Data Entry - Clean Canopy Lift Total.csv'))
  a <- leaf_area_assignments(root)
  for(nm in c('canopy_flux_goFlux_compiled.csv','canopy_flux_goFlux_compiled_with_mdf.csv')) {
    d <- read.csv(file.path(p,'goFlux_reprocessing/results',nm))
    i <- match(d$UniqueID,f$UniqueID)
    if(anyNA(i) || !setequal(d$UniqueID,f$UniqueID) ||
       !identical(as.character(d$Type),as.character(f$Type[i])) ||
       !isTRUE(all.equal(d$Height_m,f$Height_m[i],tolerance=1e-10)))
      stop('Compiled field metadata are stale; rerun flux calculations: ',nm)
    j <- match(a$UniqueID,d$UniqueID)
    for(k in setdiff(names(a),'UniqueID')) {
      if(!k %in% names(d) || !isTRUE(all.equal(d[[k]][j],a[[k]],tolerance=1e-10)))
        stop('Compiled leaf-area assignments are stale; rerun flux calculations: ',nm)
    }
  }
  invisible(TRUE)
}
