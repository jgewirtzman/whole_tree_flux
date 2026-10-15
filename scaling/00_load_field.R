# 00_load_field.R — field CH4 fluxes (six Harvard Forest trees + the Yale Myers black oak) used by all scaling steps.
# Returns `F`: one row per measurement with site, tree, species, component, height_m, flux (nmol m-2 s-1), below_mdf.
suppressPackageStartupMessages({library(dplyr); library(stringr)})
ROOT <- "."   # run all scripts from the repository root
DATA <- ROOT
hf <- read.csv(file.path(DATA, "flux/out/harvard_forest_fluxes.csv"), stringsAsFactors = FALSE) %>%
  filter(!(Species == "bg" & Tree_Tag == 3)) %>%   # black gum tag 3: three basal stem closures only, not a profiled tree
  mutate(component = ifelse(Type == "leaf (shaded)", "leaf", Type)) %>%
  filter(!is.na(CH4_best.flux), !is.na(Height_m)) %>%
  transmute(uid = UniqueID, site = Site, tree = paste(Site, Species, Tree_Tag), species = Species, component, height_m = Height_m,
            flux = CH4_best.flux, below_mdf = CH4_below_MDF_wass95)
ymf <- read.csv(file.path(DATA, "flux/out/yale_myers_fluxes.csv"), stringsAsFactors = FALSE) %>%
  mutate(height_m = suppressWarnings(as.numeric(str_extract(as.character(Height_m), "^[0-9.]+")))) %>%
  filter(!is.na(height_m), !is.na(CH4_best.flux))
ymf <- ymf %>% transmute(uid = UniqueID, site = "Yale Myers", tree = "YMF black oak", species = "qv",
                         component = if ("Type" %in% names(ymf)) Type else "stem", height_m, flux = CH4_best.flux, below_mdf = CH4_below_MDF_wass95)
F <- bind_rows(hf, ymf)


# Refuse flux tables that predate the current closure metadata or leaf areas; rerun the flux steps if this fails.
source("flux/leaf_area_helpers.R")
local({
  for (site in c("harvard_forest", "yale_myers")) {
    m <- read.csv(file.path("data/metadata", paste0(site, "_closures.csv")), stringsAsFactors = FALSE)
    d <- read.csv(file.path("flux/out", paste0(site, "_fluxes.csv")), stringsAsFactors = FALSE)
    i <- match(m$UniqueID, d$UniqueID)
    for (k in c("Type", "Height_m", "Chamber_volume_L", "Chamber_area_cm2", "Window_start", "Window_end"))
      if (anyNA(i) || nrow(d) != nrow(m) || !isTRUE(all.equal(d[[k]][i], m[[k]])))
        stop("Flux table is stale for ", site, " (", k, "); rerun flux/ before scaling")
  }
  a <- leaf_area_assignments(); d <- read.csv("flux/out/harvard_forest_fluxes.csv", stringsAsFactors = FALSE)
  if (!isTRUE(all.equal(d$Leaf_area_cm2[match(a$UniqueID, d$UniqueID)], a$Leaf_area_cm2)))
    stop("Leaf areas changed since the flux run; rerun flux/ before scaling")
})
