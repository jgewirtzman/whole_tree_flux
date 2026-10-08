# Tree illustrations and figure layouts

This folder contains the selected Figures 1 and S5, presentation graphics, and alternative layouts. The main analysis rebuilds both selected figures through `scaling/12_tree_illustrations.R`. Figure 1 uses option E: faint genus silhouettes to the left of the profiles, followed by photographs and the pooled height summary.

## Outputs

- `presentations/T_upland_row.*`: seven upland trees, with stem and branch flux colors.
- `presentations/T_all_two_rows.*`: eight trees including the wetland reference, with flux colors.
- Corresponding `*_gray_branches.*` files retain stem flux colors and use neutral branches.
- `presentations/T_structure.*`: neutral woody geometry, without flux colors.
- `presentations/genus_cluster_*.*`: overlapping genus illustrations for slide decoration.
- **Selected Figure S5:** `si_candidates/FigS5_tree_flux.*` shows all eight trees in the original two-row layout, with one shared nonlinear color scale, gray branches, and labels below. No measurement ticks or Yale Myers subtitle are drawn.
- `si_candidates/upland_stem_flux_sampling_two_rows.*` and `all_stem_flux_sampling_two_rows.*` preserve the alternative seven-tree and shared-scale eight-tree layouts.
- `si_candidates/upland_stem_flux_sampling.*`: upland row with slightly darker gray branches and small notches marking measured stem heights.
- `figure1_options/F1_A_data_only.*`: measured profiles, habitat brackets, photographs and pooled height summary.
- `figure1_options/F1_B_faint_aligned.*`: the same layout with faint simple crown silhouettes.
- `figure1_options/F1_C_data_plus_tree_row.*`: A with the upland illustration as a third row.
- `figure1_options/F1_D_faint_plus_tree_row.*`: B with the same third row.
- `figure1_options/F1_E_faint_side_left.*` and `F1_F_faint_side_right.*`: data-only profiles with a faint simple genus cluster beside panel a. Tree bases align exactly with y=0; all heights use the same axis transformation and expansion as the profiles.

- `figure1_options/F1_G_summary_first_left.*` and `F1_H_summary_first_right.*`: E/F with the second row reordered to height summary, lift photograph, climbing photograph (relabeled b–d).

PNG and vector PDF versions are provided. `*_panel_a.*` files show only the profile panels. Option E is the selected Figure 1; the other layouts remain alternatives.

## Rebuild

From the repository root, using the analysis environment:

```sh
Rscript scaling/tree_illustrations/scripts/render_final_si.R
Rscript scaling/tree_illustrations/scripts/render_figure1_options.R
Rscript scaling/tree_illustrations/scripts/render_clusters.R
Rscript scaling/tree_illustrations/scripts/render_side_silhouettes.R
```

The first command rebuilds the selected Figure S5 and companion presentation graphics. The second rebuilds Figure 1 alternatives and the earlier illustration variants without replacing Figure S5. Scripts accept a project root and optional illustration-bundle directory as the first and second arguments. They read the existing analysis outputs; rerun the analysis first if observations or processing change. They do not regenerate the numerical analyses. The fourth command rebuilds selected Figure 1 (E) and alternatives F–H. Dependencies are the existing plotting packages in `renv.lock`; color conversion uses base R.

The model geometry is a frozen illustration asset stored as transparent CSV inputs, not an additional fitted scientific result. Rendering from these inputs is reproducible; changes to crown architecture require updating the illustration inputs deliberately. The renderer checks the branch:stem surface-area ratios and flux interpolation. The Figure 1 renderer additionally verifies that silhouettes leave every observation, curve and horizontal flux scale unchanged.

## Interpretation

- Standalone tree layouts share the same horizontal and vertical spatial scale. Total heights are Sullivan species proxies, increased to the highest sample where necessary; they are not measured total heights of these individuals.
- Stem diameters use the measured profiles for the seven upland trees. The wetland stem diameter is illustrative.
- Crown widths use the saved Bechtold allometric inputs. Crown depth combines reference proportions with a shared within-species allowance for the lowest observed branch/leaf measurement. Branch topology, spacing and recursive branching rules are illustrative, not individual crown reconstructions.
- Every woody model has branch:stem lateral surface area `1.55 / 0.55 = 2.81818`. This is the common stand-scenario ratio imposed on the drawings, not a measured individual or species ratio. Segment frusta omit junction overlaps and end caps. Leaves are not drawn.
- Stem color interpolates measured height means linearly and holds endpoint rates constant below/above the sampled profile. It does not represent a fitted whole-tree budget or uncertainty. In the sampling-notch version, each notch marks a measured stem height, regardless of its detection flag.
- Colored branches use their own tree's measured branch height means where available. Asterisks mark trees whose branches instead use the measured upland tree-weighted mean. Gray branches carry no flux meaning; they do not indicate zero flux or the absence of branch measurements.
- The upland-only color range is −0.04 to 3 nmol m⁻² s⁻¹. The presentation layout sharing one scale with the wetland reference spans −0.05 to 160. Selected Figure S5 uses the same shared −0.05 to 160 scale as the eight-tree presentation layout. White is zero; equal absolute rates have matched color intensity. Color spacing uses an asinh transformation. Separate-scale versions remain exploratory alternatives.
- Simple Figure 1 silhouettes use the same height proxies and lower crown extents as the detailed drawings. Crown widths and stem widths are adapted to the flux-panel layout. Their horizontal positions have no spatial or flux meaning. They share the profiles' height-axis compression above 15 m.

Input parameter and surface-area tables in `inputs/` distinguish published crown parameters from assumed drawing rules. Output assignment and color-check tables accompany the presentation exports.

## Figure S5 color options

The compact two-row options remove all sampling-height ticks and the Yale Myers subtitle. Model geometry and flux interpolation remain unchanged. A thin neutral stem edge makes near-zero white colors visible. These are archived alternatives; the selected figure uses the earlier shared asinh scale.

- `S5_A_linear_red`: separate linear scales, red for positive rates in both groups.
- `S5_B_linear_purple`: separate linear scales, with the wetland reference in purple.
- `S5_C_gentle_purple`: upland signed-power scaling (exponent 0.65); wetland linear purple.
- `S5_D_shared_scale`: one common blue–red asinh(flux/0.5) scale for all eight trees.

Upland limits remain −0.04 to 3; the wetland scale is 0 to 160. The shared alternative spans −0.04 to 160. Zero is white; negative values are pale because their absolute magnitudes are small. Transformations affect display only. The narrower common spatial bounds (±5.4 m, 0–23 m) retain every tree polygon.

## Simplified shared-scale layout review

`si_candidates/S5_F_shared_gray_branches_labels_below.*` retains exactly the earlier `T_all_two_rows_gray_branches` stem palette and gray branch color. All eight trees use one shared nonlinear scale (−0.05 to 160 nmol m⁻² s⁻¹). Species labels sit below the stems; only the wetland reference has a subtitle. No measurement-height ticks are drawn. A single coordinate system preserves relative dimensions while reducing row/column padding. This is a review candidate, not a replacement for the selected manuscript figure.

Rebuild with `Rscript scaling/tree_illustrations/scripts/render_simple_si.R`.

The final selection restores the original equal-width panels and spacing. The compact labels-below candidate remains available for comparison; only `render_final_si.R` writes the selected `FigS5_tree_flux` files.
