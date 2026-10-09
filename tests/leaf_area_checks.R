# Field-authorized metadata corrections and independent area/unit expectations.
source('data processing/goFlux_reprocessing/leaf_area_helpers.R')
validate_canopy_metadata()
validate_compiled_canopy()
a <- leaf_area_assignments()
i <- match(c('2 10.5 leaf 1','2 10.5 leaf 2','2 11.7 leaf 1',
             '2 11.7 leaf 2','2 11.7 leaf_s 1','2 9.1 leaf 1','2 7 leaf 1'),a$UniqueID)
stopifnot(!anyNA(i),max(abs(a$Leaf_area_cm2[i]-c(286.4,334.2,286.4,286.4,286.4,228.8,286.4)))<1e-8,
          a$ShootID[i[1]]!=a$ShootID[i[2]], a$ShootID[i[4]]==a$ShootID[i[5]])
for(k in 1:3) {
  e<-new.env();load(paste0('data processing/goFlux_reprocessing/RData/manID_LGR',k,'.RData'),e)
  m<-e[[paste0('manID.LGR',k)]];j<-match(m$UniqueID,a$UniqueID);ok<-!is.na(j)
  stopifnot(all(abs(m$Area[ok]-a$Leaf_area_cm2[j[ok]])<1e-8))
}
f<-read.csv('data processing/Field Data Entry - Clean Canopy Lift Total.csv')
stopifnot(f$Type[f$UniqueID=='2 11.7 stem 1']=='stem',
          f$Type[f$UniqueID=='2 11.7 leaf 1']=='leaf')
cat('PASS: corrected tissue/height metadata, within-species leaf areas and shared-shoot assignments\n')
