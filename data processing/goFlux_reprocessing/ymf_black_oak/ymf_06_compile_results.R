# =============================================================================
# ymf_06_compile_results.R
# Merge goFlux flux results with field data from the Excel sheet to create
# one final compiled file for the YMF black oak tree.
# =============================================================================

source(file.path(
  "/Users/jongewirtzman/My Drive/Research/whole_tree_flux",
  "data processing", "goFlux_reprocessing", "ymf_black_oak", "ymf_00_setup.R"))

# --- Load flux results and field data -----------------------------------------

load(file.path(ymf_rdata_dir, "flux_results_YMF.RData"))
load(file.path(ymf_rdata_dir, "ymf_field_data.RData"))

message("Field data: ", nrow(ymf_field_data), " rows")

# --- Extract key flux columns -------------------------------------------------

# CO2 results: prefix all columns with CO2_
co2_summary <- CO2_best.YMF %>%
  mutate(UniqueID = trimws(UniqueID)) %>%
  select(
    UniqueID,
    CO2_best.flux    = best.flux,
    CO2_model        = model,
    CO2_quality.check = quality.check,
    CO2_LM.flux      = LM.flux,
    CO2_LM.r2        = LM.r2,
    CO2_LM.p.val     = LM.p.val,
    CO2_LM.RMSE      = LM.RMSE,
    CO2_HM.flux      = HM.flux,
    CO2_HM.r2        = HM.r2,
    CO2_HM.RMSE      = HM.RMSE,
    CO2_LM.diagnose  = LM.diagnose,
    CO2_HM.diagnose  = HM.diagnose,
    CO2_nb.obs       = nb.obs,
    CO2_flux.term    = flux.term,
    CO2_MDF          = MDF,
    CO2_g.fact       = g.fact
  )

# CH4 results: prefix all columns with CH4_
ch4_summary <- CH4_best.YMF %>%
  mutate(UniqueID = trimws(UniqueID)) %>%
  select(
    UniqueID,
    CH4_best.flux    = best.flux,
    CH4_model        = model,
    CH4_quality.check = quality.check,
    CH4_LM.flux      = LM.flux,
    CH4_LM.r2        = LM.r2,
    CH4_LM.p.val     = LM.p.val,
    CH4_LM.RMSE      = LM.RMSE,
    CH4_HM.flux      = HM.flux,
    CH4_HM.r2        = HM.r2,
    CH4_HM.RMSE      = HM.RMSE,
    CH4_LM.diagnose  = LM.diagnose,
    CH4_HM.diagnose  = HM.diagnose,
    CH4_nb.obs       = nb.obs,
    CH4_flux.term    = flux.term,
    CH4_MDF          = MDF,
    CH4_g.fact       = g.fact
  )

# --- Merge CO2 and CH4 -------------------------------------------------------

flux_combined <- full_join(co2_summary, ch4_summary, by = "UniqueID")
message("Flux results: ", nrow(flux_combined), " unique measurements")

# --- Merge with field data ----------------------------------------------------

final_data <- left_join(ymf_field_data, flux_combined, by = "UniqueID")

# --- Check for unmatched records ----------------------------------------------

unmatched_field <- final_data %>%
  filter(is.na(CO2_best.flux) & is.na(CH4_best.flux))
if (nrow(unmatched_field) > 0) {
  message("\nWARNING: ", nrow(unmatched_field),
          " field records had NO matching flux data:")
  message("  ", paste(unmatched_field$UniqueID, collapse = "\n  "))
}

unmatched_flux <- flux_combined %>%
  filter(!UniqueID %in% ymf_field_data$UniqueID)
if (nrow(unmatched_flux) > 0) {
  message("\nWARNING: ", nrow(unmatched_flux),
          " flux records had NO matching field data:")
  message("  ", paste(unmatched_flux$UniqueID, collapse = "\n  "))
}

# --- Summary ------------------------------------------------------------------

message("\n=== Final compiled dataset ===")
message("Rows: ", nrow(final_data))
message("Columns: ", ncol(final_data))

n_with_co2 <- sum(!is.na(final_data$CO2_best.flux))
n_with_ch4 <- sum(!is.na(final_data$CH4_best.flux))
message("Rows with CO2 flux: ", n_with_co2, " / ", nrow(final_data))
message("Rows with CH4 flux: ", n_with_ch4, " / ", nrow(final_data))

# --- Save final output --------------------------------------------------------

write.csv(final_data,
          file.path(ymf_results_dir, "ymf_black_oak_flux_compiled.csv"),
          row.names = FALSE)

write.xlsx(final_data,
           file.path(ymf_results_dir, "ymf_black_oak_flux_compiled.xlsx"))

message("\nFinal compiled data saved to:")
message("  CSV:  ", file.path(ymf_results_dir, "ymf_black_oak_flux_compiled.csv"))
message("  XLSX: ", file.path(ymf_results_dir, "ymf_black_oak_flux_compiled.xlsx"))
