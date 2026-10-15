# Rebuild Figures 1 and S3 (tree illustrations).
# Subprocesses isolate illustration renderers from the numerical analysis.
for (script in c("render_final_si.R", "render_side_silhouettes.R")) {
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    file.path("scaling/tree_illustrations/scripts", script))
  if (status != 0L) stop("Figure rendering failed: ", script)
}
