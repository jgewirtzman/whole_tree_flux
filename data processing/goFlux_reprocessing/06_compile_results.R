# =============================================================================
# 06_compile_results.R
# Merge goFlux flux results with field data to create one final compiled file.
#
# The final output contains:
#   - All columns from "Field Data Entry - Clean Canopy Lift Total.csv"
#   - CO2 and CH4 flux estimates (best.flux, model, diagnostics)
#   - Joined by UniqueID
# =============================================================================

# Source setup (works both when source()'d and run interactively in RStudio)
setup_path <- file.path(
  "/Users/jongewirtzman/My Drive/Research/whole_tree_flux",
  "data processing", "goFlux_reprocessing", "00_setup.R")
source(setup_path)

# --- Load flux results -------------------------------------------------------

load(file.path(rdata_dir, "flux_results_all.RData"))

# --- Read field data ----------------------------------------------------------

field_data <- read.csv(field_data_path, stringsAsFactors = FALSE)
field_data$UniqueID <- trimws(field_data$UniqueID)

message("Field data: ", nrow(field_data), " rows, ",
        ncol(field_data), " columns")

# --- Extract key flux columns -------------------------------------------------

# CO2 results: prefix all columns with CO2_
co2_cols <- CO2_best_all
co2_cols$UniqueID <- trimws(co2_cols$UniqueID)

co2_summary <- co2_cols %>%
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
ch4_cols <- CH4_best_all
ch4_cols$UniqueID <- trimws(ch4_cols$UniqueID)

ch4_summary <- ch4_cols %>%
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

final_data <- left_join(field_data, flux_combined, by = "UniqueID")

# --- Check for unmatched records ----------------------------------------------

unmatched_field <- final_data %>%
  filter(is.na(CO2_best.flux) & is.na(CH4_best.flux))
if (nrow(unmatched_field) > 0) {
  message("\nWARNING: ", nrow(unmatched_field),
          " field records had NO matching flux data:")
  message("  ", paste(unmatched_field$UniqueID, collapse = "\n  "))
}

unmatched_flux <- flux_combined %>%
  filter(!UniqueID %in% field_data$UniqueID)
if (nrow(unmatched_flux) > 0) {
  message("\nWARNING: ", nrow(unmatched_flux),
          " flux records had NO matching field data:")
  message("  ", paste(unmatched_flux$UniqueID, collapse = "\n  "))
}

# --- Summary ------------------------------------------------------------------

message("\n=== Final compiled dataset ===")
message("Rows: ", nrow(final_data))
message("Columns: ", ncol(final_data))
message("Field data columns: ", ncol(field_data))
message("CO2 flux columns added: ", sum(grepl("^CO2_", names(final_data))))
message("CH4 flux columns added: ", sum(grepl("^CH4_", names(final_data))))

n_with_co2 <- sum(!is.na(final_data$CO2_best.flux))
n_with_ch4 <- sum(!is.na(final_data$CH4_best.flux))
message("Rows with CO2 flux: ", n_with_co2, " / ", nrow(final_data))
message("Rows with CH4 flux: ", n_with_ch4, " / ", nrow(final_data))

# --- Save final output --------------------------------------------------------

write.csv(final_data,
          file.path(results_dir, "canopy_flux_goFlux_compiled.csv"),
          row.names = FALSE)

write.xlsx(final_data,
           file.path(results_dir, "canopy_flux_goFlux_compiled.xlsx"))

message("\nFinal compiled data saved to:")
message("  CSV:  ", file.path(results_dir, "canopy_flux_goFlux_compiled.csv"))
message("  XLSX: ", file.path(results_dir, "canopy_flux_goFlux_compiled.xlsx"))
