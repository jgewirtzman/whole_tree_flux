#!/usr/bin/env python3
"""
Summary statistics for tree methane flux manuscript.
Calculates all numbers cited in the text from the raw data.

Datasets:
  1. Harvard Forest (6 trees, July-Aug 2023) + Yale Myers Forest (1 tree, Oct 2022)
  2. Wu et al. (2024) compiled synthesis data
"""

import pandas as pd
import numpy as np
import re
import os

BASE = os.path.dirname(os.path.abspath(__file__))

# =============================================================================
# 1. LOAD OUR FIELD DATA
# =============================================================================

# --- Harvard Forest ---
hf = pd.read_csv(os.path.join(BASE, "data processing", "goFlux_reprocessing",
                               "results", "canopy_flux_goFlux_compiled.csv"))

# Remove Nyssa sylvatica (bg) tag 3
hf = hf[~((hf["Species"] == "bg") & (hf["Tree_Tag"] == 3))].copy()

# Combine "leaf (shaded)" into "leaf"
hf["Component"] = hf["Type"].replace("leaf (shaded)", "leaf")

species_map = {"bg": "Nyssa sylvatica", "rm": "Acer rubrum",
               "ro": "Quercus rubra", "hem": "Tsuga canadensis"}
hf["Species_name"] = hf["Species"].map(species_map)

# Keep rows with valid CH4 flux and height
hf = hf.dropna(subset=["CH4_best.flux", "Height_m"])

# --- Yale Myers Forest ---
ymf = pd.read_csv(os.path.join(BASE, "data processing", "goFlux_reprocessing",
                                "ymf_black_oak", "results",
                                "ymf_black_oak_flux_compiled.csv"))

def parse_ymf_height(h):
    s = str(h).strip()
    if "(" in s:
        s = s.split("(")[0].strip()
    try:
        return float(s)
    except ValueError:
        return np.nan

ymf["Height_m"] = ymf["Height_m"].apply(parse_ymf_height)
ymf = ymf.dropna(subset=["Height_m", "CH4_best.flux"])
ymf["Component"] = "stem"
ymf["Species_name"] = "Quercus velutina"
ymf["Tree_Tag"] = "YMF_1"
ymf["Site"] = "Yale Myers Forest"

# Combine
hf["Site"] = "Harvard Forest"
cols = ["Site", "Tree_Tag", "Species_name", "Height_m", "Component", "CH4_best.flux"]
field = pd.concat([hf[cols], ymf[cols]], ignore_index=True)

# =============================================================================
# 2. LOAD WU ET AL. (2024) SYNTHESIS
# =============================================================================

wu_file = os.path.join(BASE, "1-s2.0-S0168192324000911-mmc2.xlsx")
upland  = pd.read_excel(wu_file, sheet_name="upland-stem")
wetland = pd.read_excel(wu_file, sheet_name="wetland-stem")

upland["ecosystem"]  = "Upland"
wetland["ecosystem"] = "Wetland"

# Standardise column names (positional: col 0=Reference, col 8=Height, col 9=CH4)
for df in [upland, wetland]:
    df.rename(columns={df.columns[0]: "Reference",
                       df.columns[8]: "Height_raw",
                       df.columns[9]: "CH4_flux"}, inplace=True)

wu = pd.concat([upland[["Reference", "Height_raw", "CH4_flux", "ecosystem"]],
                wetland[["Reference", "Height_raw", "CH4_flux", "ecosystem"]]],
               ignore_index=True)

def parse_height(x):
    if pd.isna(x) or str(x).strip().lower() in ("none", ""):
        return np.nan
    s = str(x).strip().replace("0..6", "0.6")
    m = re.match(r"^(-?[\d.]+)\s*-\s*([\d.]+)$", s)
    if m:
        return (float(m.group(1)) + float(m.group(2))) / 2
    try:
        return float(s)
    except ValueError:
        return np.nan

wu["Height"] = wu["Height_raw"].apply(parse_height)
wu = wu.dropna(subset=["Height", "CH4_flux"])

# =============================================================================
# 3. FIELD DATA STATISTICS
# =============================================================================

print("=" * 72)
print("FIELD DATA SUMMARY (Harvard Forest + Yale Myers Forest)")
print("=" * 72)

n_obs   = len(field)
n_trees = field["Tree_Tag"].nunique()
n_spp   = field["Species_name"].nunique()
print(f"Total measurements:  {n_obs}")
print(f"Trees:               {n_trees}")
print(f"Species:             {n_spp}  ({', '.join(sorted(field['Species_name'].unique()))})")
print()

# Height breakdown
print("Measurements by component:")
print(field["Component"].value_counts().to_string())
print()

# Below vs above 2 m (all components)
below2 = field[field["Height_m"] < 2]
above2 = field[field["Height_m"] >= 2]
print(f"Mean CH4 flux < 2 m:   {below2['CH4_best.flux'].mean():.4f} nmol m-2 s-1  (n={len(below2)})")
print(f"Mean CH4 flux >= 2 m:  {above2['CH4_best.flux'].mean():.4f} nmol m-2 s-1  (n={len(above2)})")
ratio = below2["CH4_best.flux"].mean() / above2["CH4_best.flux"].mean()
print(f"Ratio (below/above):   {ratio:.1f}x")
print()

# Stem only at >= 2m
stem_above2 = field[(field["Component"] == "stem") & (field["Height_m"] >= 2)]
n_neg_stem = (stem_above2["CH4_best.flux"] < 0).sum()
print(f"Stem measurements >= 2 m:    {len(stem_above2)}")
print(f"  Negative (uptake):          {n_neg_stem}  ({100*n_neg_stem/len(stem_above2):.1f}%)")
print()

# Per-tree: majority negative above 2m?
print("Per-tree summary (all components, >= 2 m):")
for (tag, spp), grp in field[field["Height_m"] >= 2].groupby(["Tree_Tag", "Species_name"]):
    n = len(grp)
    nn = (grp["CH4_best.flux"] < 0).sum()
    print(f"  {spp:<25s} tag={str(tag):<8s}  n={n:>3d}  neg={nn:>2d} ({100*nn/n:.1f}%)")
print()

# Leaf / branch medians
for comp in ["leaf", "branch"]:
    sub = field[field["Component"] == comp]
    if len(sub):
        print(f"Median CH4 flux ({comp}): {sub['CH4_best.flux'].median():.5f} nmol m-2 s-1  (n={len(sub)})")
print()

# =============================================================================
# 4. WU ET AL. STATISTICS
# =============================================================================

print("=" * 72)
print("WU ET AL. (2024) SYNTHESIS")
print("=" * 72)

total_studies = wu["Reference"].nunique()
total_obs     = len(wu)
print(f"Total studies in compilation:      {total_studies}")
print(f"Total observations:                {total_obs}")
print()

# Studies with any measurement >= 2m
refs_with_2m = wu[wu["Height"] >= 2]["Reference"].unique()
n_refs_2m    = len(refs_with_2m)
pct_refs_2m  = 100 * n_refs_2m / total_studies
print(f"Studies with >= 2 m measurements:  {n_refs_2m}  ({pct_refs_2m:.1f}% of {total_studies} studies)")
print()

# Subset to those studies, heights >= 2m
wu_above2 = wu[(wu["Reference"].isin(refs_with_2m)) & (wu["Height"] >= 2)]
n_wu2     = len(wu_above2)
n_neg_wu  = (wu_above2["CH4_flux"] < 0).sum()
print(f"Observations at >= 2 m (from those studies):  {n_wu2}")
print(f"  Negative (uptake):  {n_neg_wu}  ({100*n_neg_wu/n_wu2:.1f}%)")
print()

print("By study (>= 2 m):")
for ref, grp in wu_above2.groupby("Reference"):
    n  = len(grp)
    nn = (grp["CH4_flux"] < 0).sum()
    label = re.sub(r"_[A-Za-z]+$", "", ref).replace("_", " ")
    print(f"  {label:<35s}  n={n:>3d}  neg={nn:>2d} ({100*nn/n:5.1f}%)")

majority_neg = [ref for ref, grp in wu_above2.groupby("Reference")
                if (grp["CH4_flux"] < 0).sum() > len(grp) / 2]
print(f"\nStudies with majority negative at >= 2 m: {len(majority_neg)}")
print()

print("=" * 72)
print("DONE")
print("=" * 72)
