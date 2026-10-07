# 00_load_field.R — GRL field CH4 fluxes (Harvard Forest 6 trees + Yale Myers black oak), as loaded for Figure 1.
# Returns `F`: one row per measurement with site, tree, species, component, height_m, flux (nmol m-2 s-1), below_mdf.
suppressPackageStartupMessages({library(dplyr); library(stringr)})
ROOT <- "."   # run all scripts from the repository root
# Flux outputs with the analyzer-volume convention 0.028 L (vtot_addition 0.057 L): on main since 2026-10-06 (cherry-pick of
# a512162 as 4d23b0a; verified byte-identical to the branch outputs).
DATA <- ROOT
hf <- read.csv(file.path(DATA, "data processing/goFlux_reprocessing/results/canopy_flux_goFlux_compiled_with_mdf.csv"), stringsAsFactors = FALSE) %>%
  filter(!(Species == "bg" & Tree_Tag == 3)) %>%            # same exclusion as figure1_composite.R
  mutate(component = ifelse(Type == "leaf (shaded)", "leaf", Type)) %>%
  filter(!is.na(CH4_best.flux), !is.na(Height_m)) %>%
  transmute(site = Site, tree = paste(Site, Species, Tree_Tag), species = Species, component, height_m = Height_m,
            flux = CH4_best.flux, below_mdf = CH4_below_MDF_wass95)
ymf <- read.csv(file.path(DATA, "data processing/goFlux_reprocessing/ymf_black_oak/results/ymf_black_oak_flux_compiled_with_mdf.csv"), stringsAsFactors = FALSE) %>%
  mutate(height_m = suppressWarnings(as.numeric(str_extract(as.character(Height_m), "^[0-9.]+")))) %>%
  filter(!is.na(height_m), !is.na(CH4_best.flux))
ymf <- ymf %>% transmute(site = "Yale Myers", tree = "YMF black oak", species = "qv",
                         component = if ("Type" %in% names(ymf)) Type else "stem", height_m, flux = CH4_best.flux, below_mdf = CH4_below_MDF_wass95)
F <- bind_rows(hf, ymf)
