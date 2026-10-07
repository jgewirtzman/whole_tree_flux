# 09_tree_weighting.R — how much do the wet (transitional) trees drive the per-area rates? Pooled vs tree-weighted means,
# leave-one-tree-out, upland-only vs transitional-only. Output: scaling/out/tree_weighting.csv
source("scaling/00_load_field.R")
suppressPackageStartupMessages(library(readr))
WET <- c("EMS rm", "Swamp Rd")   # transitional: saturated pocket in upland stand + Swamp Rd margin
G <- F %>% mutate(comp = case_when(component == "stem" & height_m < 2 ~ "stem_lo", component == "stem" ~ "stem_up", TRUE ~ component),
                  group = ifelse(grepl("Swamp Rd", tree) | grepl("EMS.* rm|rm.*EMS", tree), "transitional", "upland"))
print(G %>% distinct(tree, group))
treemean <- G %>% group_by(tree, group, comp) %>% summarise(m = mean(flux), n = n(), .groups = "drop")
print(treemean %>% tidyr::pivot_wider(id_cols = c(tree, group), names_from = comp, values_from = m), width = 200)
summ <- function(d, lab) d %>% group_by(comp) %>% summarise(pooled = mean(flux), .groups = "drop") %>%
  left_join(d %>% group_by(tree, comp) %>% summarise(m = mean(flux), .groups = "drop") %>% group_by(comp) %>% summarise(tree_wt = mean(m), ntree = n()), by = "comp") %>% mutate(set = lab)
R <- bind_rows(summ(G, "all trees"), summ(filter(G, group == "upland"), "upland only"), summ(filter(G, group == "transitional"), "transitional only"),
  bind_rows(lapply(unique(G$tree), function(t) summ(filter(G, tree != t), paste("drop", t)))))
write_csv(R, file.path(ROOT, "scaling/out/tree_weighting.csv"))
print(R %>% tidyr::pivot_wider(id_cols = set, names_from = comp, values_from = tree_wt) %>% mutate(across(where(is.numeric), ~ round(.x, 3))), width = 200)
cat("\npooled (observation-weighted):\n"); print(R %>% tidyr::pivot_wider(id_cols = set, names_from = comp, values_from = pooled) %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>% head(3))
