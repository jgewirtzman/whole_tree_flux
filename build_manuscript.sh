#!/usr/bin/env bash
# manuscript.md is the private source. Generated Word files remain gitignored.
set -euo pipefail
cd "$(dirname "$0")"
case "${1:-}" in
  --analysis) Rscript run_analysis.R ;;
  --figures)
    Rscript scaling/03_tree_component.R
    Rscript scaling/04_main_figures.R
    Rscript scaling/08_figures_v3c.R
    Rscript scaling/10_basal_by_site_and_blackgum.R
    Rscript scaling/11_manuscript_statistics.R
    ;;
  "") ;;
  *) echo "Usage: ./build_manuscript.sh [--analysis|--figures]" >&2; exit 2 ;;
esac
command -v pandoc >/dev/null || { echo "pandoc is required" >&2; exit 1; }
test -f manuscript.md || { echo "Private manuscript.md is required" >&2; exit 1; }
# Generated statistics and manuscript prose should be reviewed together after a rerun.
Rscript tests/analysis_checks.R
for part in combined main si; do
  suffix=""
  if [[ "$part" != combined ]]; then suffix="_$part"; fi
  output="DRAFT_v3d_Ground-to-Canopy_2026-10-07${suffix}.docx"
  MANUSCRIPT_PART="$part" pandoc manuscript.md --from markdown-implicit_figures --to docx --standalone \
    --resource-path=. --reference-doc=manuscript_reference.docx \
    --lua-filter=manuscript_layout.lua --output="$output"
  printf 'Built %s from manuscript.md.\n' "$output"
done
