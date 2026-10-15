# Tree illustrations (Figures 1 and S3)

The main analysis rebuilds both figures through `scaling/12_tree_illustrations.R`, writing `figures/Figure1.*` and
`figures/FigureS3.*`. The drawings are illustrations built from explicit geometry inputs, not fitted results.

## Contents

- `inputs/`: frozen geometry. Published crown parameters (`published_crown_parameters.csv`) are kept apart from
  assumed drawing rules (`branching_parameters_ASSUMED.csv`). Also included: per-tree crown inputs, stem and branch
  frusta, shared lower-branch allowances and a surface-area audit.
- `scripts/render_trees.R`: shared geometry, flux painting and color scale.
- `scripts/render_final_si.R`: Figure S3, with all eight trees in two rows, one shared color scale and gray branches.
- `scripts/render_side_silhouettes.R`: Figure 1, with faint genus silhouettes to the left of the measured profiles.
  It reuses plotting functions from `scaling/08_figures.R` and `scripts/render_figure1_options.R`.
- `checks/`: audit tables written on each run. They cover flux assignment, color balance, color-scale ranges and the
  Figure 1 silhouette alignment.

Set `EXPORT_ILLUSTRATION_ALTERNATIVES=1` to also write presentation layouts and alternative SI and Figure 1 layouts
to `presentations/`, `si_candidates/` and `figure1_options/`. These folders are not tracked.

## Rebuild

From the repository root, after the analysis outputs exist:

```sh
Rscript scaling/12_tree_illustrations.R
```

## Interpretation

- **Heights and scale:** all trees share one spatial scale. Total heights are species proxies from Sullivan et al.
  (2017), raised to the highest sample where necessary. They are not measured heights of these individuals.
- **Diameters:** stem diameters use the measured profiles of the seven upland trees. The wetland stem diameter is
  illustrative.
- **Crowns:** crown widths use Bechtold allometric inputs. Crown depth combines reference proportions with a shared
  within-species allowance for the lowest observed branch or leaf measurement. Branch topology is illustrative.
- **Surface area:** each woody model has branch:stem lateral surface area 1.55 / 0.55. This is the stand-scenario
  ratio imposed on the drawings, not a measured individual or species ratio. Leaves are not drawn.
- **Stem color:** linear interpolation between measured height means, held constant beyond the sampled range. This is
  a display rule, not a fitted whole-tree budget.
- **Color scale:** one shared asinh scale, −0.05 to 160 nmol m⁻² s⁻¹. White is zero, and equal absolute rates have
  matched color intensity.
- **Branches:** gray branches carry no flux meaning.
- **Figure 1 silhouettes:** these share the profiles' height axis (compressed above 15 m). Their horizontal positions
  carry no meaning.
