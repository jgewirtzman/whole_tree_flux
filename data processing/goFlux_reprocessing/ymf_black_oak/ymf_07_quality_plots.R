# =============================================================================
# ymf_07_quality_plots.R
# Generate quality-check flux plots and export as PDF for YMF black oak.
# =============================================================================

source(file.path(
  "/Users/jongewirtzman/My Drive/Research/whole_tree_flux",
  "data processing", "goFlux_reprocessing", "ymf_black_oak", "ymf_00_setup.R"))

# --- Load results -------------------------------------------------------------

load(file.path(ymf_rdata_dir, "manID_YMF.RData"))
load(file.path(ymf_rdata_dir, "flux_results_YMF.RData"))

# --- Patch flux2pdf class check bug -------------------------------------------
# goFlux v0.2.0's flux2pdf() validation uses sapply(plot.list, class), which
# returns a matrix when ggplot objects have multiple classes (ggplot2 >= 3.5
# uses S7 and reports 5 classes per object). The grep-based length check then
# fails. We patch the function to use a corrected check.

flux2pdf_patched <- function(plot.list, outfile = NULL, width = 11.6, height = 8.2) {
  if (missing(plot.list)) stop("'plot.list' is required")
  if (!is.list(plot.list)) stop("'plot.list' must be of class list")
  is_gg <- sapply(plot.list, function(x) inherits(x, "ggplot"))
  if (!all(is_gg)) {
    stop("all elements of 'plot.list' must be of class 'gg, ggplot'")
  }
  if (!is.null(outfile) && !is.character(outfile))
    stop("'outfile' must be of class character")
  if (!is.numeric(width)) stop("'width' must be of class numeric")
  if (!is.numeric(height)) stop("'height' must be of class numeric")

  . <- NULL
  group_plot.list <- rlist::list.group(plot.list, .$plot_env$UniqueID)
  outplot <- pbapply::pblapply(group_plot.list, function(p) {
    title <- grid::textGrob(paste("Unique ID:", p[[1]]$plot_env$UniqueID),
                            gp = grid::gpar(fontsize = 16))
    footnote <- grid::textGrob(
      paste("page", which(names(group_plot.list) == p[[1]]$plot_env$UniqueID),
            "of", length(group_plot.list)),
      gp = grid::gpar(fontface = 3, fontsize = 12), hjust = 1)
    n.plot <- length(p)
    nrow <- ifelse(n.plot <= 2, 1, ifelse(n.plot <= 4, 2, 4))
    ncol <- ifelse(n.plot <= 1, 1, 2)
    gridExtra::marrangeGrob(grobs = p, ncol = ncol, nrow = nrow,
                            top = title, bottom = footnote)
  })
  if (is.null(outfile)) {
    outfile <- paste0(getwd(), "/", deparse(substitute(plot.list)), ".pdf")
  }
  if (length(outplot) > 10000)
    stop("The outfile contains more than 10,000 pages.")
  pdf(file = outfile, width = width, height = height)
  max.print <- getOption("max.print")
  options(max.print = 10000)
  suppressMessages(print(outplot))
  options(max.print = max.print)
  dev.off()
}

# --- Plot settings ------------------------------------------------------------

plot.legend  <- c("MAE", "RMSE", "AICc", "k.ratio", "g.factor")
plot.display <- c("MDF", "prec", "nb.obs", "flux.term")
quality.check <- TRUE

# =============================================================================
# Generate plots
# =============================================================================

message("=== Generating YMF plots ===")

CO2_plots.YMF <- flux.plot(
  CO2_best.YMF, manID.YMF, "CO2dry_ppm",
  shoulder = 20, plot.legend = plot.legend,
  plot.display = plot.display, quality.check = quality.check
)

CH4_plots.YMF <- flux.plot(
  CH4_best.YMF, manID.YMF, "CH4dry_ppb",
  shoulder = 20, plot.legend = plot.legend,
  plot.display = plot.display, quality.check = quality.check
)

flux2pdf_patched(c(CO2_plots.YMF, CH4_plots.YMF),
         outfile = file.path(ymf_plots_dir, "YMF_flux_results.pdf"))
message("YMF plots saved")

# =============================================================================
# Summary
# =============================================================================

message("\n=== Quality plots complete ===")
message("PDF saved to: ", file.path(ymf_plots_dir, "YMF_flux_results.pdf"))
message("\nReview these plots to verify flux estimates and identify")
message("any measurements that may need re-processing.")
