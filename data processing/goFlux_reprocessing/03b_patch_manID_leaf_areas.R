# Optional metadata refresh without changing selected concentration windows.
# Rebuild auxfiles with 03_build_auxfiles.R first. The standard flux step also
# performs this synchronization; no separate leaf-area lookup is maintained.
source("data processing/goFlux_reprocessing/00_setup.R")
for (inst in paste0("LGR",1:3)) {
  e <- new.env()
  load(file.path(rdata_dir,paste0("aux_",inst,".RData")),e)
  load(file.path(rdata_dir,paste0("manID_",inst,".RData")),e)
  nm <- paste0("manID.",inst)
  e[[nm]] <- sync_manID(e[[nm]],e[[paste0("aux.",inst)]])
  save(list=nm,envir=e,file=file.path(rdata_dir,paste0("manID_",inst,".RData")))
}
