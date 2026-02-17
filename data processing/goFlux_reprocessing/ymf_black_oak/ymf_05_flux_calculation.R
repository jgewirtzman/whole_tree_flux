# =============================================================================
# ymf_05_flux_calculation.R
# Calculate CO2 and CH4 fluxes using goFlux() and select best model
# using best.flux().
#
# Units:
#   CO2 flux: umol m-2 s-1  (from CO2dry_ppm)
#   CH4 flux: nmol m-2 s-1  (from CH4dry_ppb)
# =============================================================================

source(file.path(
  "/Users/jongewirtzman/My Drive/Research/whole_tree_flux",
  "data processing", "goFlux_reprocessing", "ymf_black_oak", "ymf_00_setup.R"))

# --- Load manual ID results ---------------------------------------------------

load(file.path(ymf_rdata_dir, "manID_YMF.RData"))

# =============================================================================
# Calculate fluxes
# =============================================================================

message("\n=== YMF: Calculating CO2 fluxes ===")
CO2_flux.YMF <- goFlux(manID.YMF, "CO2dry_ppm")

message("=== YMF: Calculating CH4 fluxes ===")
CH4_flux.YMF <- goFlux(manID.YMF, "CH4dry_ppb")

message("\nYMF: Selecting best flux estimates...")
CO2_best.YMF <- best.flux(CO2_flux.YMF, flux_criteria)
CH4_best.YMF <- best.flux(CH4_flux.YMF, flux_criteria)

message("YMF: ", nrow(CO2_best.YMF), " CO2 fluxes, ",
        nrow(CH4_best.YMF), " CH4 fluxes")

# Model selection summary
message("\nCO2 model selection: ",
        sum(CO2_best.YMF$model == "LM"), " LM, ",
        sum(CO2_best.YMF$model == "HM"), " HM")
message("CH4 model selection: ",
        sum(CH4_best.YMF$model == "LM"), " LM, ",
        sum(CH4_best.YMF$model == "HM"), " HM")

# =============================================================================
# Save results
# =============================================================================

save(CO2_best.YMF, CH4_best.YMF,
     file = file.path(ymf_rdata_dir, "flux_results_YMF.RData"))

write.xlsx(CO2_best.YMF, file.path(ymf_results_dir, "CO2_best_YMF.xlsx"))
write.xlsx(CH4_best.YMF, file.path(ymf_results_dir, "CH4_best_YMF.xlsx"))

message("\nFlux results saved to: ", ymf_rdata_dir)
message("Excel exports saved to: ", ymf_results_dir)
message("Proceed to ymf_06_compile_results.R")
