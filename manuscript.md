# Better evidence is needed for widespread methane uptake by tree surfaces

**J. Gewirtzman**

Yale School of the Environment, New Haven, CT, USA

---

## Key Points

1. Only 20% of stem CH₄ flux studies measure above 2 m, capturing as little as 1% of total tree surface area.
2. Among the few above-2 m measurements that exist, uptake is rare: 6–7% of fluxes are negative.
3. Claims of canopy-height CH₄ uptake by trees outpace the evidence; targeted measurements are needed.

---

## Abstract

Recent studies have proposed that tree surfaces transition from methane (CH₄) emission near the base to net uptake at canopy heights, with potentially large implications for global CH₄ budgets. Here we evaluate the evidence for this claim on two fronts. First, we show that the existing measurement base is truncated: only 20% of studies in a recent global synthesis measured stem fluxes above 2 m, and measurements at 10 m height capture only 4.2% of a representative tree's total surface area (stem, branches, and leaves combined). Second, we present new field data from seven trees at Harvard Forest and Yale Myers Forest measured up to 22 m and reanalyze above-2 m data from the literature. In both datasets, CH₄ uptake is rare: only 6.0% of our above-2 m stem measurements and 7.1% of literature above-2 m observations are negative, and no study shows a consistent pattern of net uptake at height. While mean fluxes decline 3.4-fold above 2 m, they do not cross zero. We conclude that current data do not support widespread canopy-height CH₄ uptake by tree surfaces. Absence of evidence is not evidence of absence, but neither is it grounds for assuming uptake occurs. Continuous, multi-height flux measurements extending into the upper canopy are needed before conclusions can be drawn about this process.

## 1. Introduction

Trees transport methane (CH₄) from the subsurface to the atmosphere. Over the past two decades, studies have documented CH₄ emissions from tree stems across ecosystems ranging from tropical wetlands to temperate upland forests (Pangala et al., 2017; Covey et al., 2012; Barba et al., 2019; Gauci, 2025). These emissions are thought to originate primarily from methanogenesis in anoxic soil microsites, with CH₄ transported upward through the xylem or bark tissues and released through the stem surface (Covey & Megonigal, 2019), though internal methanogenesis within heartwood can also contribute (Wang et al., 2016; Gewirtzman et al., 2025).

A consistent observation across this literature is that stem CH₄ fluxes decline with measurement height. Fluxes measured near the base of the tree — typically below 1–2 m — tend to be substantially larger than those measured higher on the stem. This pattern is generally attributed to the proximity of lower stem surfaces to the soil CH₄ source and to the potential for diffusive loss and oxidation during upward transport.

Recently, several authors have extended this declining-flux pattern to suggest that tree surfaces may switch from net emission to net uptake at canopy heights (Welch et al., 2022). Gauci et al. (2024) reported net CH₄ uptake on woody surfaces above approximately 2 m in tropical, temperate, and boreal forests and estimated a global woody-surface CH₄ sink of 24.6–49.9 Tg yr⁻¹. Jeffrey et al. (2021a) demonstrated that methane-oxidizing bacteria residing in tree bark can reduce stem CH₄ emissions by 36 ± 5%. Other studies have reported net CH₄ uptake by stems or foliage under certain conditions (Machacova et al., 2021; Gorgolewski et al., 2023a). If canopy-height uptake is widespread, it would affect how tree CH₄ fluxes are scaled to landscape and global levels, since a tree that emits at the base but absorbs CH₄ across its canopy surface area would have a different net flux than one estimated from basal measurements alone. However, most tree surface area — comprising stems, branches, and leaves — resides above the heights at which measurements are typically made (Gewirtzman, 2026).

Here we evaluate the evidence for this transition. We show that the empirical basis for canopy-height CH₄ uptake is limited, and that the measurements which do exist mostly do not support it.

## 2. Methods

### 2.1 Field measurements

We measured CH₄ fluxes on tree stems, branches, and leaves at two sites: Harvard Forest (Petersham, MA; 42.53°N, 72.19°W) and Yale Myers Forest (Union, CT). At Harvard Forest, six trees of four species (*Acer rubrum*, *Nyssa sylvatica*, *Quercus rubra*, and *Tsuga canadensis*) were measured on 18–19 July and 16–17 August 2023. Heights ranged from ground level to 22 m. The upper canopy was accessed using a canopy aerial lift (DinoLift). At Yale Myers Forest, a single mature *Quercus velutina* was measured on 4 October 2022 at seven heights along the stem (0.5, 1.25, 2, 4, 6, 8, and 10 m), with heights above 2 m accessed via arborist climbing equipment.

### 2.2 Chamber design and gas analysis

Stem CH₄ fluxes were measured using rigid static chambers constructed from transparent plastic containers (Rubbermaid) with arcs cut to conform to the stem surface. Chamber volumes (0.5–2.0 L) were determined by gas standard dilution in the laboratory (Jeffrey et al., 2021b; Siegenthaler et al., 2016). Enclosed stem surface area was calculated as the product of the chamber's planar area and arc length. Chambers were sealed to stems using potting clay (Amaco) and secured with lashing straps. Seal integrity was verified by visual inspection, by blowing on the chamber perimeter and monitoring for CO₂ spikes, and by checking concentration linearity during measurement.

Chambers were connected via 5 mm internal-diameter PVC tubing (Bev-a-Line) to a portable off-axis integrated cavity output spectroscopy (OA-ICOS) analyzer (GLA131-GGA, Los Gatos Research) measuring CO₂, CH₄, and H₂O concentrations at 1 Hz with CH₄ precision <0.9 ppb. The analyzer's 5 L min⁻¹ flow rate provided complete chamber volume turnover within approximately 30 s. Measurements lasted 3–10 min per location, with the initial 30 s excluded to allow equilibration.

### 2.3 Meteorological data

Air temperature and barometric pressure used in flux calculations were obtained from the Fisher Meteorological Station at Harvard Forest (Boose & VanScoy, 2026; HF001 v.35). The Fisher station records 15-minute averages of air temperature (Vaisala HMP45C, 2.2 m height), barometric pressure (Vaisala CS105), and other variables on a Campbell Scientific CR10X datalogger (Boose & VanScoy, 2026). Each flux measurement was matched to the nearest 15-minute meteorological observation by timestamp. Pressure was converted from millibar to kPa.

### 2.4 Flux calculation

CH₄ and CO₂ fluxes were calculated using the goFlux R package (Rheault et al., 2024), which fits both linear (LM) and non-linear Hutchinson-Mosier (HM) models to the time series of chamber headspace concentration, with automatic correction for water vapor dilution (Hutchinson & Mosier, 1981; Hüppi et al., 2018). Flux (F, in nmol m⁻² s⁻¹ for CH₄) is computed as:

F = (dC/dt) × (V_c / A_c) × (P / RT) × (1 / (1 − X_H₂O))

where dC/dt is the rate of concentration change, V_c is chamber volume (L), A_c is enclosed surface area (m²), P is atmospheric pressure (kPa), R is the gas constant (8.314 L kPa K⁻¹ mol⁻¹), T is air temperature (K), and X_H₂O is the water vapor mole fraction.

Best-fit model selection was performed using the goFlux `best.flux()` function with criteria including goodness-of-fit metrics (MAE, RMSE relative to instrument precision; AICc), physical constraints (g-factor < 2.0, where g = HM flux / LM flux; κ/κ_max ≤ 1.0), and quality thresholds (P < 0.05; flux > minimal detectable flux; n ≥ 60 observations). CO₂ flux R² was used as a chamber seal quality metric, since tree stem respiration produces consistently positive CO₂ fluxes.

### 2.5 Wu et al. (2024) reanalysis

We reanalyzed the compiled dataset from Wu et al. (2024), accessed from the supplementary material of the original publication (https://doi.org/10.1016/j.agrformet.2024.109976). The dataset contains observations from 50 studies across upland and wetland ecosystems. We used the "upland-stem" and "wetland-stem" sheets, which together contain 1,010 observations with associated measurement heights and CH₄ flux values.

Measurement heights reported as ranges (e.g., "2–4") were converted to their midpoint. Observations with missing height or flux values were excluded. We classified each study by its maximum reported measurement height and tallied observations at or above height thresholds of 2, 5, and 10 m. No deduplication across studies was performed, as individual observations in the compilation represent independent measurements. All reanalysis code is available in the project repository.

### 2.6 Surface area truncation model

To quantify what fraction of a tree's total surface area is captured by stem measurements up to a given height, we modeled a representative temperate canopy tree as a cone (height = 26 m, DBH = 40 cm). The stem radius was assumed to taper linearly from its base value to zero at the tree top. Lateral surface area of the cone below height h was calculated analytically. Total tree surface area was partitioned using the surface area indices of Whittaker and Woodwell (1967): stem bark = 0.45 m² m⁻² ground area, branch bark = 1.70 m² m⁻² ground area, and leaf area index (LAI) = 4.5 m² m⁻² ground area. The capture fraction at height h is the stem surface area below h divided by the relevant total (stem only, stem + branches, or stem + branches + leaves).

### 2.7 Software

All analyses were conducted in R. Flux calculations used goFlux (Rheault et al., 2024). Data processing used dplyr, readr, tidyr, stringr, lubridate, readxl, and openxlsx. Figures were produced with ggplot2 and patchwork. The truncation model was computed analytically in R.

## 3. Results

### 3.1 Field measurements

The combined dataset from Harvard Forest and Yale Myers Forest (Section 2) comprises 141 measurements on seven trees of five species, including 109 stem, 11 branch, and 21 leaf measurements spanning heights from ground level to 22 m.

Consistent with prior work, mean CH₄ flux was substantially higher below 2 m (0.88 nmol m⁻² s⁻¹, n = 26) than at or above 2 m (0.26 nmol m⁻² s⁻¹, n = 115) — a ratio of 3.4×. Fluxes generally declined with height (Figure 1), and occasional hotspots (localized zones of elevated emission) were observed at various heights along some stems, complicating a simple monotonic decline model.

> **Figure 1.** Vertical profiles of CH₄ flux versus height. (a) Six individual trees at Harvard Forest, with measurement height on the y-axis and CH₄ flux (nmol m⁻² s⁻¹) on the x-axis. Points are colored by component (stem, branch, leaf). LOESS smooths are fitted to stem measurements only. The dashed vertical line indicates zero flux. (b) *Quercus velutina* at Yale Myers Forest; heights range from 0.5 to 10 m. (c) Canopy access via aerial lift (DinoLift) at Harvard Forest. (d) Arborist climbing at Yale Myers Forest. Fluxes are predominantly positive at all heights in both datasets.

However, the declining fluxes did not cross zero. Among 83 stem measurements at or above 2 m, only 5 were negative (6.0%). At higher thresholds, the results are similar: at ≥5 m, 5 of 61 measurements (8.2%) were negative; at ≥10 m, 4 of 31 (12.9%). All seven trees are represented at every height threshold, so this pattern is not driven by sampling bias toward particular individuals.

Branch and leaf fluxes were small in magnitude. Median leaf CH₄ flux was 0.007 nmol m⁻² s⁻¹ (n = 21) and median branch flux was 0.070 nmol m⁻² s⁻¹ (n = 11). While these values are much lower than basal stem fluxes, they are predominantly positive — not negative. This contrasts with the daytime foliar uptake (−0.54 nmol m⁻² s⁻¹) reported by Gorgolewski et al. (2023a) in upland forests in Ontario, where uptake was attributed to endophytic methanotrophic bacteria and was linked to transpiration. Our leaf measurements were made during the day, suggesting that foliar uptake is not a universal feature of upland tree canopies.

On a per-tree basis, no individual showed a consistent pattern of net uptake at height. The tree with the highest proportion of negative fluxes above 2 m was a *Tsuga canadensis* individual with 4 of 13 measurements negative (30.8%), but even this tree had a majority of positive fluxes. Most trees had few or no negative measurements at height.

### 3.2 Reanalysis of above-2 m literature data

Using the Wu et al. (2024) synthesis (Section 2.5; 1,010 observations from 50 studies), we classified each study by its maximum measurement height. The measurement record is heavily truncated: only 10 of 50 studies (20%) include any measurements at or above 2 m: Churkina et al. (2018), Epron et al. (2022), Gorgolewski et al. (2023b), Iddris et al. (2021), Machacova et al. (2023), Moldaschl et al. (2021), Sjögersten et al. (2020), Vainio et al. (2022), van Haren et al. (2021), and Wang et al. (2016). Coverage thins further at greater heights: 3 studies (6%) measured at or above 5 m, and none reached 10 m. The maximum measurement height in the entire compilation is 8.0 m.

We examined the 127 above-2 m observations from these 10 studies (Figure 2). Of these, 9 were negative (7.1%), similar to the 6.0% in our own data. At ≥5 m, there were 12 observations across 3 studies, and none were negative.

> **Figure 2.** Above-2 m CH₄ flux observations from the Wu et al. (2024) synthesis, faceted by study (10 studies with measurements at ≥2 m). Points are colored by ecosystem type (upland/wetland). LOESS smooths are shown per study. Fluxes decline with height but remain predominantly positive across the literature.

No study in the compilation had a majority of negative fluxes at or above 2 m. The study with the highest proportion of negative fluxes was Moldaschl et al. (2021), with 1 of 3 measurements negative (33.3%), but this is based on a very small sample. The largest above-2 m dataset is from Iddris et al. (2021), with 44 observations and 5 negative (11.4%). Across the remaining studies — spanning boreal (Vainio et al., 2022), temperate (Epron et al., 2022; Machacova et al., 2023; Churkina et al., 2018; Gorgolewski et al., 2023b), and tropical ecosystems (Sjögersten et al., 2020; van Haren et al., 2021; Wang et al., 2016) — the proportion of negative fluxes at height ranges from 0 to 20%.

Across both our new data and the existing literature, the above-2 m measurement record shows declining but predominantly positive CH₄ fluxes. Negative fluxes, where observed, are sporadic and represent a small minority of measurements.

## 4. Discussion

### 4.1 The surface area truncation problem

To put these measurement heights in perspective, we modeled a representative temperate deciduous tree as a cone with a height of 26 m and a DBH of 40 cm (Section 2.6). Using the surface area partitioning of Whittaker and Woodwell (1967) — which gives ratios of 0.45 m² m⁻² for stem bark, 1.70 m² m⁻² for branch bark, and 4.5 m² m⁻² (leaf area index) for leaves — we calculated cumulative capture fractions for three nested surface area denominators: stem only, stem plus branches, and stem plus branches plus leaves (Figure 3).

At 2 m — the height below which most measurements are concentrated — one has sampled 14.8% of the stem surface area, 3.1% of stem-plus-branch area, and 1.0% of total tree surface area (stem, branches, and leaves). At 10 m, these fractions reach 62.1%, 13.0%, and 4.2%, respectively (Figure 3b,c). Most of a tree's surface area is in the branches and leaves of the upper canopy, which has rarely been measured.

> **Figure 3.** The measurement truncation problem. (a) Schematic of a representative 26 m tree modeled as a cone, with gradient fill (blue at base to red at top) and dashed ellipses at 2 m and 10 m measurement heights. (b) Capture-fraction curves showing the cumulative percentage of surface area sampled as a function of measurement height, for three surface area denominators: stem only, stem plus branches, and stem plus branches plus leaves. (c) Bar chart showing the same capture fractions at 2 m and 10 m for the three surface area categories.

The claim that tree surfaces transition to CH₄ uptake at canopy heights is therefore being made about a portion of the tree for which almost no data exist. Current conclusions about whole-tree fluxes rest on measurements of a small fraction of the stem.

The scaling implication of this truncation is straightforward. Most of a tree's surface area lies above the heights at which fluxes are typically measured. As a result, the whole-tree methane balance depends disproportionately on assumptions about canopy fluxes. If canopy surfaces are assumed to take up methane, the net tree flux can shift toward neutrality or even become a sink. If canopy fluxes remain small but positive, the tree remains a source. Because this large canopy-area term is largely unmeasured, different extrapolation choices can produce different whole-tree flux directions. Current observations show declining flux with height, but they do not demonstrate a sign reversal across the canopy.

Importantly, even fluxes that are small in magnitude at height can strongly influence whole-tree estimates when multiplied by the large surface area of branches and leaves. The scaling question is therefore not only whether uptake occurs, but whether it is sufficiently large and spatially consistent to outweigh basal emissions when integrated over the canopy.

### 4.2 Absence of evidence is not evidence of absence

We are not claiming that methanotrophic CH₄ uptake on tree surfaces does not occur. Methanotrophy on bark is biologically plausible: methane-oxidizing bacteria have been detected on bark surfaces and shown to reduce stem emissions by more than a third in laboratory and field experiments (Jeffrey et al., 2021a). Machacova et al. (2021) reported net CH₄ uptake by tree stems in a tropical rain forest on volcanic Réunion Island, indicating that tree-level uptake can occur under specific conditions. Reduced emission at height is consistent with partial methanotrophic consumption of CH₄ diffusing through bark, even if the net flux remains positive.

Our point is that the existing evidence does not support the broader claim that tree surfaces generally exhibit net CH₄ uptake at canopy heights — as proposed, for example, by Gauci et al. (2024), who estimated a global woody-surface CH₄ sink of 24.6–49.9 Tg yr⁻¹. The measurement record above 2 m is thin, and the data that exist show mostly positive fluxes. Extrapolating a declining trend to assume it crosses zero is a modeling inference, not an empirical observation.

The distinction between reduced emission and net uptake matters for scaling. A tree whose canopy surfaces emit small but positive CH₄ fluxes contributes to landscape emissions. A tree whose canopy surfaces take up CH₄ partially offsets its own basal emissions. The sign of the flux, not just its magnitude, determines the direction of the scaling correction.

### 4.3 Complications to simple height-flux models

Our data also show complications for smooth monotonic decline models. Several trees had hotspots — localized zones of elevated emission at heights well above the base — likely associated with wounds, branch junctions, or internal decay columns. Mochidome et al. (2025) recently demonstrated that local methanogenesis within upper trunk tissues can drive emissions at heights of 4–12 m, with 44–89% of total trunk CH₄ emissions originating above 3 m. Mochidome and Epron (2024) documented up to 15-fold intra-individual spatial variation in trunk CH₄ fluxes, driven largely by variation in sapwood CH₄ production rates. Wounds and branches can also be significant emission sources (Gorgolewski et al., 2023b). These features produce fluxes at mid-stem and upper-stem heights that rival or exceed basal values and are not predictable from height alone. Parameterizing tree CH₄ fluxes as a simple function of height will need to account for this heterogeneity.

### 4.4 Limitations

Our field data are limited to seven trees at two temperate forest sites in the northeastern United States, measured during the growing season. We cannot rule out the possibility that uptake is more prevalent in other biomes, seasons, or environmental conditions. The measurements are also snapshots in time rather than continuous records, and fluxes may vary diurnally and seasonally.

These same limitations apply to the studies that have been cited as evidence for canopy-height uptake. The literature record above 2 m is thinner and less systematic than what we present here.

### 4.5 What is needed

**Table 1.** Qualitative uncertainty in scaling whole-tree CH₄ fluxes.

| Component | Surface area certainty | Flux certainty | Primary limitation |
|---|---|---|---|
| Lower stem (≤2 m) | Relatively high | Moderate | Frequently measured; simple geometry |
| Upper stem (>2 m) | Relatively high | Low | Sparse vertical measurements |
| Branches | Low | Low | Architecture complex; limited area and flux data |
| Leaves | High (LAI-based) | Low | Leaf area index known; direct CH₄ flux data rare |

Uncertainty in whole-tree scaling is dominated by canopy components. Stem surface area can be estimated geometrically with relatively high confidence, and lower-stem flux has been measured extensively. In contrast, upper-stem, branch, and leaf fluxes remain sparsely quantified despite comprising most of total tree surface area. As a result, the least constrained terms in the scaling equation correspond to the largest surface areas.

Resolving the question of canopy-height CH₄ fluxes will require a different measurement approach: continuous or repeated measurements at multiple heights on the same trees, extending from the base into the upper canopy and including branches and leaves. As Gauci (2025) emphasizes, tree CH₄ exchange responses to environmental changes remain poorly quantified, and reducing uncertainty in the global CH₄ budget requires better spatial and temporal coverage. Automated chamber systems on canopy towers or accessed via tree-climbing techniques could provide the resolution needed to detect whether and when uptake occurs. Pairing flux measurements with characterization of bark-associated microbial communities (Jeffrey et al., 2021a; Gewirtzman et al., 2025) would help connect observed flux patterns to their biological drivers.

Independent constraints on canopy-height exchange could also be obtained by integrating bottom-up and top-down approaches. Eddy covariance (EC) measurements over upland forests capture net ecosystem CH₄ exchange. In systems where soil CH₄ flux is independently measured, the residual flux could, in principle, constrain the aggregate tree component. However, such partitioning compounds uncertainties from multiple terms (soil, litter, stem, foliage) and requires careful co-location and temporal alignment of measurements. Coordinated EC–chamber studies would therefore provide an important test of canopy-uptake hypotheses at ecosystem scale.

Surface area itself is an under-constrained term in many scaling exercises. Stem area can be estimated geometrically or from diameter–height allometries with relatively low uncertainty. In contrast, branch surface area — often exceeding stem area — is poorly quantified in most forests. Terrestrial LiDAR and structure-from-motion photogrammetry now permit direct three-dimensional reconstruction of tree architecture and offer a path toward empirically constraining woody surface area. Incorporating such measurements would reduce geometric uncertainty in whole-tree CH₄ budgets and clarify the magnitude of the canopy-area term.

Until such measurements are made, canopy-height CH₄ uptake by trees should be treated as a hypothesis, not an established finding.

## 5. Conclusions

The idea that tree surfaces switch from CH₄ emission to net uptake at canopy heights has implications for how tree CH₄ fluxes are scaled to landscape and global levels. We have shown that this claim has a limited empirical basis: only 20% of studies measure above 2 m, capturing as little as 1% of a tree's total surface area. Among the above-2 m measurements that do exist — including our new data from seven trees measured up to 22 m — uptake is rare, representing 6–7% of observations. Mean fluxes decline with height but do not cross zero.

Absence of evidence is not evidence of absence, but current observations do not yet broadly demonstrate a sign reversal with height. Targeted, multi-height measurements extending into the upper canopy are needed before conclusions about this process can be drawn.

## Acknowledgments

Field measurements at Harvard Forest were conducted at the Harvard Forest Long Term Ecological Research site. Yale Myers Forest measurements were facilitated by the Yale School of the Environment.

## References

Barba, J., Bradford, M. A., Brewer, P. E., Bruhn, D., Covey, K., van Haren, J., ... & Vargas, R. (2019). Methane emissions from tree stems: a new frontier in the global carbon cycle. *New Phytologist*, 222(1), 18–28.

Boose, E., & VanScoy, M. (2026). Fisher Meteorological Station at Harvard Forest since 2001. Harvard Forest Data Archive: HF001 (v.35). Environmental Data Initiative. https://doi.org/10.6073/pasta/00eae8d70316abe9b27a18760093d3fb.

Churkina, A. I., Mochenov, S. Yu., Sabrekov, S. F., Glagolev, M. V., Il'yasov, D. V., Terentieva, I. E., & Maksyutov, S. S. (2018). Soils in seasonally flooded forests as methane sources: A case study of West Siberian South taiga. *IOP Conference Series: Earth and Environmental Science*, 138, 012012.

Covey, K. R., & Megonigal, J. P. (2019). Methane production and emissions in trees and forests. *New Phytologist*, 222(1), 35–51.

Covey, K. R., Wood, S. A., Warren, R. J., Lee, X., & Bradford, M. A. (2012). Elevated methane concentrations in trees of an upland forest. *Geophysical Research Letters*, 39, L15705.

Epron, D., Mochidome, T., Tanabe, T., Dannoura, M., & Sakabe, A. (2022). Variability in stem methane emissions and wood methane production of different tree species in a cold temperate mountain forest. *Ecosystems*, 26, 784–799.

Gauci, V. (2025). Tree methane exchange in a changing world. *Nature Reviews Earth & Environment*, 6, 471–483.

Gauci, V., Pangala, S. R., Shenkin, A., Barba, J., Bastviken, D., Figueiredo, V., ... & Malhi, Y. (2024). Global atmospheric methane uptake by upland tree woody surfaces. *Nature*, 631, 796–800.

Gewirtzman, J. (2026). The Global Woody Surface: A planetary interface for biodiversity, ecosystem function, and climate. *Global Change Biology*, 32, e70699.

Gewirtzman, J., Arnold, W., Taylor, M., Burrows, H., Merenstein, C., Woodbury, D., ... & Bradford, M. A. (2025). Tree microbiomes and methane emissions in upland forests. *bioRxiv* preprint. https://doi.org/10.1101/2025.09.30.679632.

Gorgolewski, A. S., Caspersen, J. P., Vantellingen, J., & Thomas, S. C. (2023a). Tree foliage is a methane sink in upland temperate forests. *Ecosystems*, 26, 174–186.

Gorgolewski, A. S., Vantellingen, J., Caspersen, J. P., & Thomas, S. C. (2023b). Overlooked sources of methane emissions from trees: branches and wounds. *Canadian Journal of Forest Research*, 52, 1165–1175.

Hutchinson, G. L., & Mosier, A. R. (1981). Improved soil cover method for field measurement of nitrous oxide fluxes. *Soil Science Society of America Journal*, 45(2), 311–316.

Hüppi, R., Felber, R., Krauss, M., Six, J., Leifeld, J., & Fuß, R. (2018). Restricting the nonlinearity parameter in soil greenhouse gas flux calculation for more reliable flux estimates. *PLOS ONE*, 13(7), e0200876.

Iddris, N. A., Corre, M. D., Yemefack, M., van Straaten, O., & Veldkamp, E. (2021). Stem and soil nitrous oxide fluxes from rainforest and cacao agroforest on highly weathered soils in the Congo Basin. *Journal of Geophysical Research: Biogeosciences*, 126, e2020JG006014.

Jeffrey, L. C., Maher, D. T., Chiri, E., Leung, P. M., Nauer, P. A., Arndt, S. K., Tait, D. R., Greening, C., & Johnston, S. G. (2021a). Bark-dwelling methanotrophic bacteria decrease methane emissions from trees. *Nature Communications*, 12, 2127.

Jeffrey, L. C., Maher, D. T., Tait, D. R., & Johnston, S. G. (2021b). A small nimble *in situ* fine-scale flux method for measuring tree stem greenhouse gas emissions and processes. *Ecosystems*, 24, 2021–2037.

Machacova, K., Borak, L., Agyei, T., Schindler, T., Soosaar, K., Mander, Ü., & Ah-Peng, C. (2021). Trees as net sinks for methane (CH₄) and nitrous oxide (N₂O) in the lowland tropical rain forest on volcanic Réunion Island. *New Phytologist*, 229, 1559–1574.

Machacova, K., Warlo, H., Svobodová, K., Agyei, T., Uchytilová, T., Horáček, P., & Lang, F. (2023). Methane emission from stems of European beech (*Fagus sylvatica*) offsets as much as half of methane oxidation in soil. *New Phytologist*, 238, 584–597.

Mochidome, T., & Epron, D. (2024). Drivers of intra-individual spatial variability in methane emissions from tree trunks in upland forest. *Trees*, 38, 625–636.

Mochidome, T., Hölttä, T., Asakawa, S., Watanabe, T., Dannoura, M., & Epron, D. (2025). Local methanogenesis drives significant methane emissions from upper tree trunks in a cool-temperate upland forest. *New Phytologist*. https://doi.org/10.1111/nph.70331.

Moldaschl, E., Kitzler, B., Machacova, K., Schindler, T., & Schindlbacher, A. (2021). Stem CH₄ and N₂O fluxes of *Fraxinus excelsior* and *Populus alba* trees along a flooding gradient. *Plant and Soil*, 461, 407–420.

Pangala, S. R., Enrich-Prast, A., Basso, L. S., Peixoto, R. B., Bastviken, D., Hornibrook, E. R., ... & Gauci, V. (2017). Large emissions from floodplain trees close the Amazon methane budget. *Nature*, 552, 230–234.

Rheault, K., Bhatt, M., Bhiry, N., & Garneau, M. (2024). goFlux: An R package for calculating greenhouse gas fluxes from static chambers. *Journal of Open Source Software*, 9(99), 6547.

Siegenthaler, A., Welch, B., Pangala, S. R., Peacock, M., & Gauci, V. (2016). Semi-rigid chambers for methane gas flux measurements on tree stems. *Wetlands*, 36, 47–53.

Sjögersten, S., Siegenthaler, A., Lopez, O. R., Aplin, P., Turner, B., & Gauci, V. (2020). Methane emissions from tree stems in neotropical peatlands. *New Phytologist*, 225, 769–781.

Vainio, E., Haikarainen, I. P., Machacova, K., Putkinen, A., Pihlatie, M., & Pumpanen, J. (2022). Soil-tree-atmosphere CH₄ flux dynamics of boreal birch and spruce trees during spring leaf-out. *Plant and Soil*, 478, 391–407.

van Haren, J., Brewer, P. E., Kurtzberg, L., Wehr, R. N., Springer, V. L., Tello Espinoza, R., ... & Cadillo-Quiroz, H. (2021). A versatile gas flux chamber reveals high tree stem CH₄ emissions in Amazonian peatland. *Agricultural and Forest Meteorology*, 297, 108244.

Wang, Z. P., Gu, Q., Deng, F. D., Huang, J. H., Megonigal, J. P., Yu, Q., ... & Han, X. G. (2016). Methane emissions from the trunks of living trees on upland soils. *New Phytologist*, 211, 429–439.

Welch, B., Gauci, V., & Sayer, E. J. (2022). Tree stem methane fluxes: an emerging frontier in the global carbon cycle. *Frontiers in Forests and Global Change*, 5, 867112.

Whittaker, R. H., & Woodwell, G. M. (1967). Surface area relations of woody plants and forest communities. *American Journal of Botany*, 54(8), 931–939.

Wu, J., Zhang, H., Cheng, X., & Liu, G. (2024). Tree stem methane emissions: Global patterns and controlling factors. *Agricultural and Forest Meteorology*, 350, 109976.

