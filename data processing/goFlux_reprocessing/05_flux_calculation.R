# =============================================================================
# 05_flux_calculation.R
# Calculate CO2 and CH4 fluxes using goFlux() and select best model
# using best.flux().
#
# For each measurement, goFlux calculates:
#   - Linear model (LM) flux
#   - Non-linear Hutchinson & Mosier (HM) flux
# Then best.flux selects the best estimate based on objective criteria.
#
# Units:
#   CO2 flux: umol m-2 s-1  (from CO2dry_ppm)
#   CH4 flux: nmol m-2 s-1  (from CH4dry_ppb)
# =============================================================================

# Source setup (works both when source()'d and run interactively in RStudio)
setup_path <- file.path(
  "/Users/jongewirtzman/My Drive/Research/whole_tree_flux",
  "data processing", "goFlux_reprocessing", "00_setup.R")
source(setup_path)

# --- Load manual ID results ---------------------------------------------------

load(file.path(rdata_dir, "manID_LGR1.RData"))
load(file.path(rdata_dir, "manID_LGR2.RData"))
load(file.path(rdata_dir, "manID_LGR3.RData"))

# =============================================================================
# LGR1
# =============================================================================

message("\n=== LGR1: Calculating fluxes ===")
CO2_flux.LGR1 <- goFlux(manID.LGR1, "CO2dry_ppm")
CH4_flux.LGR1 <- goFlux(manID.LGR1, "CH4dry_ppb")

message("LGR1: Selecting best flux estimates...")
CO2_best.LGR1 <- best.flux(CO2_flux.LGR1, flux_criteria)
CH4_best.LGR1 <- best.flux(CH4_flux.LGR1, flux_criteria)

message("LGR1: ", nrow(CO2_best.LGR1), " CO2 fluxes, ",
        nrow(CH4_best.LGR1), " CH4 fluxes")

# =============================================================================
# LGR2
# =============================================================================

message("\n=== LGR2: Calculating fluxes ===")
CO2_flux.LGR2 <- goFlux(manID.LGR2, "CO2dry_ppm")
CH4_flux.LGR2 <- goFlux(manID.LGR2, "CH4dry_ppb")

message("LGR2: Selecting best flux estimates...")
CO2_best.LGR2 <- best.flux(CO2_flux.LGR2, flux_criteria)
CH4_best.LGR2 <- best.flux(CH4_flux.LGR2, flux_criteria)

message("LGR2: ", nrow(CO2_best.LGR2), " CO2 fluxes, ",
        nrow(CH4_best.LGR2), " CH4 fluxes")

# =============================================================================
# LGR3
# =============================================================================

message("\n=== LGR3: Calculating fluxes ===")
CO2_flux.LGR3 <- goFlux(manID.LGR3, "CO2dry_ppm")
CH4_flux.LGR3 <- goFlux(manID.LGR3, "CH4dry_ppb")

message("LGR3: Selecting best flux estimates...")
CO2_best.LGR3 <- best.flux(CO2_flux.LGR3, flux_criteria)
CH4_best.LGR3 <- best.flux(CH4_flux.LGR3, flux_criteria)

message("LGR3: ", nrow(CO2_best.LGR3), " CO2 fluxes, ",
        nrow(CH4_best.LGR3), " CH4 fluxes")

# =============================================================================
# Combine across instruments
# =============================================================================

CO2_best_all <- rbind(CO2_best.LGR1, CO2_best.LGR2, CO2_best.LGR3)
CH4_best_all <- rbind(CH4_best.LGR1, CH4_best.LGR2, CH4_best.LGR3)

message("\n=== Combined results ===")
message("Total CO2 fluxes: ", nrow(CO2_best_all))
message("Total CH4 fluxes: ", nrow(CH4_best_all))

# Model selection summary
message("\nCO2 model selection: ",
        sum(CO2_best_all$model == "LM"), " LM, ",
        sum(CO2_best_all$model == "HM"), " HM")
message("CH4 model selection: ",
        sum(CH4_best_all$model == "LM"), " LM, ",
        sum(CH4_best_all$model == "HM"), " HM")

# =============================================================================
# Save results
# =============================================================================

# Per-instrument
save(CO2_best.LGR1, CH4_best.LGR1,
     file = file.path(rdata_dir, "flux_results_LGR1.RData"))
save(CO2_best.LGR2, CH4_best.LGR2,
     file = file.path(rdata_dir, "flux_results_LGR2.RData"))
save(CO2_best.LGR3, CH4_best.LGR3,
     file = file.path(rdata_dir, "flux_results_LGR3.RData"))

# Combined
save(CO2_best_all, CH4_best_all,
     file = file.path(rdata_dir, "flux_results_all.RData"))

# Excel exports
write.xlsx(CO2_best_all, file.path(results_dir, "CO2_best_all.xlsx"))
write.xlsx(CH4_best_all, file.path(results_dir, "CH4_best_all.xlsx"))

message("\nFlux results saved to: ", rdata_dir)
message("Excel exports saved to: ", results_dir)
message("Proceed to 06_compile_results.R and/or 07_quality_plots.R")
