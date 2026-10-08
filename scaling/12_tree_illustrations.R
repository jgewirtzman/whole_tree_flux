# Rebuild Figure S5 and companion presentation graphics from frozen drawing inputs.
# A subprocess isolates the illustration renderer from the numerical analysis.
status <- system2(file.path(R.home("bin"), "Rscript"),
                  "scaling/tree_illustrations/scripts/render_final_si.R")
if (status != 0L) stop("Tree illustration rendering failed")
