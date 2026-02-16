# =============================================================================
# 07_quality_plots.R
# Generate quality-check flux plots and export as PDF.
#
# Each measurement gets a 2-panel plot:
#   Left:  Time series with flagged data points
#   Right: Linear and non-linear model fits overlaid
# Annotations show the selected flux, model statistics, and quality flags.
#
# One PDF is generated per instrument, containing both CO2 and CH4 plots.
# =============================================================================

# Source setup (works both when source()'d and run interactively in RStudio)
setup_path <- file.path(
  "/Users/jongewirtzman/My Drive/Research/whole_tree_flux",
  "data processing", "goFlux_reprocessing", "00_setup.R")
source(setup_path)

# --- Load results -------------------------------------------------------------

load(file.path(rdata_dir, "manID_LGR1.RData"))
load(file.path(rdata_dir, "manID_LGR2.RData"))
load(file.path(rdata_dir, "manID_LGR3.RData"))
load(file.path(rdata_dir, "flux_results_LGR1.RData"))
load(file.path(rdata_dir, "flux_results_LGR2.RData"))
load(file.path(rdata_dir, "flux_results_LGR3.RData"))

# --- Patch flux2pdf class check bug -------------------------------------------
# goFlux v0.2.0's flux2pdf() validation uses sapply(plot.list, class), which
# returns a matrix when ggplot objects have multiple classes (ggplot2 >= 3.5
# uses S7 and reports 5 classes per object). The grep-based length check then
# fails. We patch the function to use a corrected check.

flux2pdf_patched <- function(plot.list, outfile = NULL, width = 11.6, height = 8.2) {
  if (missing(plot.list)) stop("'plot.list' is required")
  if (!is.list(plot.list)) stop("'plot.list' must be of class list")
  # Fixed validation: check each element individually
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
# LGR1
# =============================================================================

message("=== Generating LGR1 plots ===")
CO2_plots.LGR1 <- flux.plot(
  CO2_best.LGR1, manID.LGR1, "CO2dry_ppm",
  shoulder = 20, plot.legend = plot.legend,
  plot.display = plot.display, quality.check = quality.check
)

CH4_plots.LGR1 <- flux.plot(
  CH4_best.LGR1, manID.LGR1, "CH4dry_ppb",
  shoulder = 20, plot.legend = plot.legend,
  plot.display = plot.display, quality.check = quality.check
)

flux2pdf_patched(c(CO2_plots.LGR1, CH4_plots.LGR1),
         outfile = file.path(plots_dir, "LGR1_flux_results.pdf"))
message("LGR1 plots saved")

# =============================================================================
# LGR2
# =============================================================================

message("=== Generating LGR2 plots ===")
CO2_plots.LGR2 <- flux.plot(
  CO2_best.LGR2, manID.LGR2, "CO2dry_ppm",
  shoulder = 20, plot.legend = plot.legend,
  plot.display = plot.display, quality.check = quality.check
)

CH4_plots.LGR2 <- flux.plot(
  CH4_best.LGR2, manID.LGR2, "CH4dry_ppb",
  shoulder = 20, plot.legend = plot.legend,
  plot.display = plot.display, quality.check = quality.check
)

flux2pdf_patched(c(CO2_plots.LGR2, CH4_plots.LGR2),
         outfile = file.path(plots_dir, "LGR2_flux_results.pdf"))
message("LGR2 plots saved")

# =============================================================================
# LGR3
# =============================================================================

message("=== Generating LGR3 plots ===")
CO2_plots.LGR3 <- flux.plot(
  CO2_best.LGR3, manID.LGR3, "CO2dry_ppm",
  shoulder = 20, plot.legend = plot.legend,
  plot.display = plot.display, quality.check = quality.check
)

CH4_plots.LGR3 <- flux.plot(
  CH4_best.LGR3, manID.LGR3, "CH4dry_ppb",
  shoulder = 20, plot.legend = plot.legend,
  plot.display = plot.display, quality.check = quality.check
)

flux2pdf_patched(c(CO2_plots.LGR3, CH4_plots.LGR3),
         outfile = file.path(plots_dir, "LGR3_flux_results.pdf"))
message("LGR3 plots saved")

# =============================================================================
# Summary
# =============================================================================

message("\n=== Quality plots complete ===")
message("PDFs saved to: ", plots_dir)
message("  LGR1_flux_results.pdf")
message("  LGR2_flux_results.pdf")
message("  LGR3_flux_results.pdf")
message("\nReview these plots to verify flux estimates and identify")
message("any measurements that may need re-processing.")
