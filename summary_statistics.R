#!/usr/bin/env Rscript
# =============================================================================
# summary_statistics.R
# Calculates all summary statistics cited in the manuscript.
# Run from the project root (whole_tree_flux/).
# =============================================================================

library(readxl)
library(dplyr)
library(stringr)

base_dir <- here::here()  # or set manually to your project root
if (!file.exists(file.path(base_dir, "synthesis.R"))) {
  base_dir <- "/Users/jongewirtzman/My Drive/Research/whole_tree_flux"
}

# =============================================================================
# 1. LOAD FIELD DATA
# =============================================================================

# --- Harvard Forest ---
hf <- read.csv(file.path(base_dir, "data processing", "goFlux_reprocessing",
                          "results", "canopy_flux_goFlux_compiled.csv"),
               stringsAsFactors = FALSE)

# Remove Nyssa sylvatica tag 3
hf <- hf %>% filter(!(Species == "bg" & Tree_Tag == 3))

# Combine "leaf (shaded)" into "leaf"
hf <- hf %>% mutate(Component = ifelse(Type == "leaf (shaded)", "leaf", Type))

species_map <- c(bg = "Nyssa sylvatica", rm = "Acer rubrum",
                 ro = "Quercus rubra", hem = "Tsuga canadensis")
hf <- hf %>%
  mutate(Species_name = species_map[Species],
         Site = "Harvard Forest") %>%
  filter(!is.na(CH4_best.flux), !is.na(Height_m))

# --- Yale Myers Forest ---
ymf <- read.csv(file.path(base_dir, "data processing", "goFlux_reprocessing",
                           "ymf_black_oak", "results",
                           "ymf_black_oak_flux_compiled.csv"),
                stringsAsFactors = FALSE)

parse_ymf_height <- function(h) {
  s <- str_trim(as.character(h))
  if (grepl("\\(", s)) s <- str_trim(str_extract(s, "^[^(]+"))
  suppressWarnings(as.numeric(s))
}

ymf <- ymf %>%
  mutate(Height_m = sapply(Height_m, parse_ymf_height)) %>%
  filter(!is.na(Height_m), !is.na(CH4_best.flux)) %>%
  mutate(Component = "stem",
         Species_name = "Quercus velutina",
         Tree_Tag = "YMF_1",
         Site = "Yale Myers Forest")

# --- Combine ---
hf <- hf %>% mutate(Tree_Tag = as.character(Tree_Tag))
cols <- c("Site", "Tree_Tag", "Species_name", "Height_m", "Component", "CH4_best.flux")
field <- bind_rows(hf[, cols], ymf[, cols])

# =============================================================================
# 2. LOAD WU ET AL. (2024) SYNTHESIS
# =============================================================================

wu_file <- file.path(base_dir, "1-s2.0-S0168192324000911-mmc2.xlsx")
upland  <- read_xlsx(wu_file, sheet = "upland-stem")  %>% mutate(ecosystem = "Upland")
wetland <- read_xlsx(wu_file, sheet = "wetland-stem") %>% mutate(ecosystem = "Wetland")

colnames(upland)[1]  <- "Reference"
colnames(upland)[9]  <- "Height_raw"
colnames(upland)[10] <- "CH4_flux"
colnames(wetland)[1]  <- "Reference"
colnames(wetland)[9]  <- "Height_raw"
colnames(wetland)[10] <- "CH4_flux"

wu <- bind_rows(
  upland  %>% select(Reference, Height_raw, CH4_flux, ecosystem),
  wetland %>% select(Reference, Height_raw, CH4_flux, ecosystem)
)

parse_height <- function(x) {
  if (is.na(x) || x == "None") return(NA_real_)
  s <- str_trim(as.character(x))
  s <- str_replace(s, "^0\\.\\.6", "0.6")
  m <- str_match(s, "^(-?[0-9.]+)\\s*-\\s*([0-9.]+)$")
  if (!is.na(m[1, 1])) return(mean(c(as.numeric(m[1, 2]), as.numeric(m[1, 3]))))
  suppressWarnings(as.numeric(s))
}

wu <- wu %>%
  mutate(Height = sapply(Height_raw, parse_height)) %>%
  filter(!is.na(Height), !is.na(CH4_flux))

# =============================================================================
# 3. FIELD DATA STATISTICS
# =============================================================================

# Open a sink to also write all output to a text file
out_file <- file.path(base_dir, "summary_statistics_output.txt")
sink(out_file, split = TRUE)  # split = TRUE prints to console AND file

cat("========================================================================\n")
cat("FIELD DATA SUMMARY (Harvard Forest + Yale Myers Forest)\n")
cat("========================================================================\n\n")

cat(sprintf("Total measurements:  %d\n", nrow(field)))
cat(sprintf("Trees:               %d\n", n_distinct(field$Tree_Tag)))
cat(sprintf("Species:             %d  (%s)\n\n",
            n_distinct(field$Species_name),
            paste(sort(unique(field$Species_name)), collapse = ", ")))

cat("Measurements by component:\n")
print(table(field$Component))
cat("\n")

# Below vs above 2m
below2 <- field %>% filter(Height_m < 2)
above2 <- field %>% filter(Height_m >= 2)
cat(sprintf("Mean CH4 flux < 2 m:   %.4f nmol m-2 s-1  (n=%d)\n",
            mean(below2$CH4_best.flux), nrow(below2)))
cat(sprintf("Mean CH4 flux >= 2 m:  %.4f nmol m-2 s-1  (n=%d)\n",
            mean(above2$CH4_best.flux), nrow(above2)))
cat(sprintf("Ratio (below/above):   %.1fx\n\n",
            mean(below2$CH4_best.flux) / mean(above2$CH4_best.flux)))

# Stem only at >= 2m
stem_above2 <- field %>% filter(Component == "stem", Height_m >= 2)
n_neg <- sum(stem_above2$CH4_best.flux < 0)
cat(sprintf("Stem measurements >= 2 m:    %d\n", nrow(stem_above2)))
cat(sprintf("  Negative (uptake):          %d  (%.1f%%)\n\n", n_neg,
            100 * n_neg / nrow(stem_above2)))

# Per-tree at >= 2m
cat("Per-tree summary (all components, >= 2 m):\n")
field %>%
  filter(Height_m >= 2) %>%
  group_by(Tree_Tag, Species_name) %>%
  summarise(n = n(),
            neg = sum(CH4_best.flux < 0),
            pct_neg = 100 * neg / n,
            .groups = "drop") %>%
  arrange(Species_name) %>%
  { for (i in 1:nrow(.)) {
      cat(sprintf("  %-25s tag=%-8s  n=%3d  neg=%2d (%.1f%%)\n",
                  .$Species_name[i], as.character(.$Tree_Tag[i]),
                  .$n[i], .$neg[i], .$pct_neg[i]))
    }
  }
cat("\n")

# Leaf and branch medians
for (comp in c("leaf", "branch")) {
  sub <- field %>% filter(Component == comp)
  if (nrow(sub) > 0) {
    cat(sprintf("Median CH4 flux (%s): %.5f nmol m-2 s-1  (n=%d)\n",
                comp, median(sub$CH4_best.flux), nrow(sub)))
  }
}
cat("\n")

# Height threshold breakdown (stem only)
cat("Stem measurements by height threshold:\n")
for (h_thresh in c(2, 5, 10)) {
  stem_above <- field %>% filter(Component == "stem", Height_m >= h_thresh)
  n_trees    <- n_distinct(stem_above$Tree_Tag)
  n_meas     <- nrow(stem_above)
  n_neg_h    <- sum(stem_above$CH4_best.flux < 0)
  cat(sprintf("  >= %2d m:  %d trees  |  %3d measurements  |  %d negative (%.1f%%)\n",
              h_thresh, n_trees, n_meas, n_neg_h,
              ifelse(n_meas > 0, 100 * n_neg_h / n_meas, 0)))
}
cat("\n")

# =============================================================================
# 4. WU ET AL. STATISTICS
# =============================================================================

cat("========================================================================\n")
cat("WU ET AL. (2024) SYNTHESIS\n")
cat("========================================================================\n\n")

total_studies <- n_distinct(wu$Reference)
total_obs     <- nrow(wu)
cat(sprintf("Total studies in compilation:      %d\n", total_studies))
cat(sprintf("Total observations:                %d\n\n", total_obs))

# Studies with any measurement >= 2m
refs_2m <- wu %>% filter(Height >= 2) %>% distinct(Reference) %>% pull(Reference)
n_refs_2m <- length(refs_2m)
cat(sprintf("Studies with >= 2 m measurements:  %d  (%.1f%% of %d studies)\n\n",
            n_refs_2m, 100 * n_refs_2m / total_studies, total_studies))

# Observations at >= 2m from those studies
wu_above2 <- wu %>% filter(Reference %in% refs_2m, Height >= 2)
n_neg_wu  <- sum(wu_above2$CH4_flux < 0)
cat(sprintf("Observations at >= 2 m (from those studies):  %d\n", nrow(wu_above2)))
cat(sprintf("  Negative (uptake):  %d  (%.1f%%)\n\n", n_neg_wu,
            100 * n_neg_wu / nrow(wu_above2)))

cat("By study (>= 2 m):\n")
wu_above2 %>%
  group_by(Reference) %>%
  summarise(n = n(),
            neg = sum(CH4_flux < 0),
            pct_neg = 100 * neg / n,
            .groups = "drop") %>%
  mutate(label = str_replace(Reference, "_[A-Za-z]+$", "") %>%
           str_replace_all("_", " ")) %>%
  arrange(label) %>%
  { for (i in 1:nrow(.)) {
      cat(sprintf("  %-35s  n=%3d  neg=%2d (%5.1f%%)\n",
                  .$label[i], .$n[i], .$neg[i], .$pct_neg[i]))
    }
  }

majority_neg <- wu_above2 %>%
  group_by(Reference) %>%
  summarise(neg = sum(CH4_flux < 0), n = n(), .groups = "drop") %>%
  filter(neg > n / 2)
cat(sprintf("\nStudies with majority negative at >= 2 m: %d\n\n", nrow(majority_neg)))

# Height threshold breakdown
cat("Studies and observations by height threshold:\n")
for (h_thresh in c(2, 5, 10)) {
  refs_h  <- wu %>% filter(Height >= h_thresh) %>% distinct(Reference) %>% pull(Reference)
  obs_h   <- wu %>% filter(Height >= h_thresh)
  n_neg_h <- sum(obs_h$CH4_flux < 0)
  cat(sprintf("  >= %2d m:  %2d of %d studies (%4.1f%%)  |  %3d obs  |  %d negative (%4.1f%%)\n",
              h_thresh, length(refs_h), total_studies,
              100 * length(refs_h) / total_studies,
              nrow(obs_h), n_neg_h,
              ifelse(nrow(obs_h) > 0, 100 * n_neg_h / nrow(obs_h), 0)))
}
cat(sprintf("  Max height in compilation: %.1f m\n", max(wu$Height, na.rm = TRUE)))

cat("\n========================================================================\n")
cat("DONE\n")
cat("========================================================================\n")

# Close sink — output saved to file
sink()
cat(sprintf("\nSummary written to: %s\n", out_file))
