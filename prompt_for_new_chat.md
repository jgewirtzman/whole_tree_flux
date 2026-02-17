# Task: Write a GRL short report on tree methane fluxes

## Title

"Better evidence is needed for widespread methane uptake by tree surfaces"

Author: J. Gewirtzman, Yale School of the Environment. Sole author for now; co-authors may be added later.

## Core argument

Recent work has suggested that tree surfaces switch from CH₄ emission to net uptake at canopy heights. This paper argues that the evidence base for this claim is insufficient — on two fronts:

1. **The measurement gap.** The vast majority of stem CH₄ flux studies measure only near the base of the tree. Only 20% of studies in the Wu et al. (2024) synthesis measured above 2m. And even at 10m, you've captured only ~4% of a tree's total surface area (stem + branches + leaves). We are extrapolating from almost nothing.

2. **What the measurements actually show.** Among the few studies and measurements that do exist above 2m — including our own new data from 7 trees — uptake is rare. Only 6–7% of above-2m measurements are negative, and no tree or study shows a consistent pattern of net uptake.

**The key framing: absence of evidence is not evidence of absence.** We are NOT claiming that canopy uptake doesn't happen. We ARE saying that (a) it has barely been measured, and (b) where it has been measured, the data mostly show continued emission. The field needs targeted canopy-height measurements before drawing conclusions.

## Outline

**Introduction** — The claim: recent work suggests trees may switch from CH₄ emission to uptake at canopy heights. Why it matters (scaling implications for global CH₄ budgets). But how strong is the evidence?

**The measurement gap** — Wu et al. reanalysis: only 10 of 50 studies (20%) measure above 2m. The truncation figure: a representative 26m tree modeled as a cone — at 2m you've sampled 14.8% of stem SA but only 1.0% of total tree SA (stem + branches + leaves). At 10m: 62.1% of stem, but only 4.2% of total. We are drawing sweeping conclusions from a tiny fraction of the tree.

**What the measurements actually show** — Our 7 trees at Harvard Forest and Yale Myers Forest: persistent emission from ground to canopy (up to 22m). Of 83 stem measurements at ≥2m, only 5 are negative (6.0%). Mean flux declines with height (3.4× higher below 2m) but does not cross zero. Wu et al. above-2m data: 9 of 127 measurements negative (7.1%). Zero studies show a majority of negative fluxes at height.

**Discussion** — Absence of evidence ≠ evidence of absence. Methanotrophy on bark is plausible; reduced emission ≠ net uptake. Hotspots above the base complicate simple monotonic models. Limitations: 7 trees, two temperate sites, snapshots in time. What's needed: continuous, multi-height measurements extending into the canopy.

**Conclusions** — Current evidence does not support the claim. Better measurements are needed before we can conclude anything about canopy-height uptake.

## Key statistics

All from `summary_statistics_output.txt`:

**Our field data (Harvard Forest + Yale Myers Forest):**
- 141 total measurements, 7 trees, 5 species (Acer rubrum, Nyssa sylvatica, Quercus rubra, Quercus velutina, Tsuga canadensis)
- 109 stem, 11 branch, 21 leaf measurements
- Mean CH₄ flux below 2m: 0.88 nmol m⁻² s⁻¹ (n=26); at ≥2m: 0.26 nmol m⁻² s⁻¹ (n=115). Ratio: 3.4×
- Stem measurements by height threshold: ≥2m: 83 meas, 5 negative (6.0%); ≥5m: 61 meas, 5 negative (8.2%); ≥10m: 31 meas, 4 negative (12.9%). All 7 trees represented at every threshold.
- Median leaf flux: 0.007 nmol m⁻² s⁻¹; median branch flux: 0.070 nmol m⁻² s⁻¹

**Wu et al. (2024) reanalysis:**
- 50 studies, 1010 total observations
- Studies by height threshold: ≥2m: 10 of 50 (20%); ≥5m: 3 of 50 (6%); ≥10m: 0 of 50 (0%). Max height in compilation: 8.0m.
- 127 observations at ≥2m; 9 negative (7.1%). At ≥5m: 12 obs, 0 negative.
- 0 studies with majority negative fluxes at ≥2m

**Truncation analysis (representative 26m tree, 40cm DBH, cone geometry):**
- Whittaker & Woodwell (1967) surface area ratios: stem bark 0.45, branch bark 1.70, LAI 4.5 (m² per m² ground)
- At 2m: 14.8% of stem SA captured, 3.1% of stem+branch, 1.0% of stem+branch+leaf
- At 10m: 62.1% of stem, 13.0% of stem+branch, 4.2% of stem+branch+leaf

## Figures

**Figure 1 (Harvard Forest flux profiles):** 6 panels, one per individual tree. Height on y-axis, CH₄ flux on x-axis. Points colored by component (stem, branch, leaf). LOESS smooth on stem only. Dashed line at zero flux. Shows persistent positive fluxes at all heights, with occasional hotspots.

**Figure 2 (Yale Myers Forest flux profile):** Single panel for Quercus velutina. Same orientation. Stem only, heights 0.5–10m. Same pattern — positive throughout.

**Figure 3 (Wu et al. reanalysis):** Faceted by study (10 studies that measured ≥2m). Points colored by ecosystem (upland/wetland). LOESS smooths. Shows declining but mostly positive fluxes across the literature.

**Figure 4 (Truncation figure):** 3-panel layout. (a) Cone schematic of a tree with gradient fill (blue at base → red at top), dashed ellipses at 2m and 10m. (b) Capture-fraction curves: three lines showing % of surface area captured as a function of measurement height, for three denominators (stem only, stem+branches, stem+branches+leaves). (c) Bar chart: same three categories at 2m and 10m. This figure makes the measurement-gap argument visually.

## Format

Write the manuscript as a Markdown file (.md). Start from scratch — do not reference any existing draft. Target GRL length: ~12 publication units max (roughly 4800 words + figures). Include 3 key points (≤140 chars each) at the top.
