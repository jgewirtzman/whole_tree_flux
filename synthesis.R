library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)
library(scales)

# ── Read data ──
file_path <- "/Users/jongewirtzman/My Drive/Research/whole_tree_flux/1-s2.0-S0168192324000911-mmc2.xlsx"

upland <- read_xlsx(file_path, sheet = "upland-stem") %>%
  mutate(ecosystem = "Upland")
wetland <- read_xlsx(file_path, sheet = "wetland-stem") %>%
  mutate(ecosystem = "Wetland")

# Harmonize column names
colnames(upland)[1]  <- "Reference"
colnames(upland)[9]  <- "Height_raw"
colnames(upland)[10] <- "CH4_flux"
colnames(wetland)[1]  <- "Reference"
colnames(wetland)[9]  <- "Height_raw"
colnames(wetland)[10] <- "CH4_flux"

# Combine
dat <- bind_rows(
  upland  %>% select(Reference, Height_raw, CH4_flux, ecosystem),
  wetland %>% select(Reference, Height_raw, CH4_flux, ecosystem)
)

# ── Parse height ──
parse_height <- function(x) {
  if (is.na(x) || x == "None") return(NA_real_)
  x <- as.character(x) %>% str_trim()
  x <- str_replace(x, "^0\\.\\.6", "0.6")
  
  range_match <- str_match(x, "^(-?[0-9.]+)\\s*-\\s*([0-9.]+)$")
  if (!is.na(range_match[1, 1])) {
    a <- as.numeric(range_match[1, 2])
    b <- as.numeric(range_match[1, 3])
    if (!is.na(a) && !is.na(b)) return(mean(c(a, b)))
  }
  
  val <- suppressWarnings(as.numeric(x))
  return(val)
}

dat <- dat %>%
  mutate(Height = sapply(Height_raw, parse_height))

# ── Filter to projects with at least one measurement >= 2m ──
projects_2m <- dat %>%
  filter(!is.na(Height), Height >= 2) %>%
  distinct(Reference, ecosystem)

dat_filt <- dat %>%
  semi_join(projects_2m, by = c("Reference", "ecosystem")) %>%
  filter(!is.na(Height))

# Facet label
dat_filt <- dat_filt %>%
  mutate(facet_label = paste0(Reference, "\n(", ecosystem, ")"))

cat("Projects with measurements >= 2m:", nrow(projects_2m), "\n")
cat("Observations in filtered set:", nrow(dat_filt), "\n")

# ── Custom arcsinh transform ──
asinh_trans <- function() {
  trans_new(
    name = "asinh",
    transform = asinh,
    inverse = sinh,
    breaks = function(x) {
      rng <- range(x, na.rm = TRUE)
      raw_breaks <- pretty(asinh(rng), n = 5)
      sinh(raw_breaks)
    },
    format = function(x) {
      signif(x, 2)
    }
  )
}

# ── Shared plot skeleton ──
# x = Height, y = CH4_flux, then coord_flip
# After coord_flip: height displays on y-axis (vertical), flux on x-axis (horizontal)
# In ggplot with coord_flip:
#   scale_x_ controls Height (displayed vertically after flip)
#   scale_y_ controls CH4_flux (displayed horizontally after flip)
#   "free_y" in facet_wrap frees the flux axis (originally y, displayed as x after flip)
#   "free_x" in facet_wrap frees the height axis (originally x, displayed as y after flip)
base_plot <- function(data) {
  ggplot(data, aes(x = Height, y = CH4_flux, color = ecosystem)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
    geom_point(alpha = 0.5, size = 1.3) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 0.8, na.rm = TRUE) +
    scale_color_manual(values = c("Upland" = "steelblue", "Wetland" = "darkorange")) +
    labs(
      x = "Sampling height (m)",
      color = "Ecosystem"
    ) +
    coord_flip() +
    theme_bw(base_size = 10) +
    theme(
      strip.text = element_text(size = 7),
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
}

ncol_facet <- 4
n_facets <- length(unique(dat_filt$facet_label))
plot_h <- ceiling(n_facets / ncol_facet) * 3.5 + 1.5
plot_w <- 16

# ── Version 1: Free flux axis ──
# free_y frees the y aesthetic (CH4_flux), which displays horizontally after coord_flip
p1 <- base_plot(dat_filt) +
  facet_wrap(~ facet_label, scales = "free_x", ncol = ncol_facet) +
  labs(
    title = "Stem CH\u2084 flux vs. height (projects with measurements \u2265 2 m)",
    subtitle = "Free flux axis",
    y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1}))
  )

ggsave("flux_by_height_free.pdf", p1, width = plot_w, height = plot_h, limitsize = FALSE)

# ── Version 2: Free arcsinh-transformed axis ──
# scale_y_continuous transforms the y aesthetic (CH4_flux)
p2 <- base_plot(dat_filt) +
  scale_y_continuous(trans = asinh_trans()) +
  facet_wrap(~ facet_label, scales = "free_x", ncol = ncol_facet) +
  labs(
    title = "Stem CH\u2084 flux vs. height (projects with measurements \u2265 2 m)",
    subtitle = "Free arcsinh-transformed flux axis",
    y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})~"[arcsinh scale]")
  )

ggsave("flux_by_height_asinh_free.pdf", p2, width = plot_w, height = plot_h, limitsize = FALSE)

# ── Version 3: Fixed arcsinh-transformed axis ──
p3 <- base_plot(dat_filt) +
  scale_y_continuous(trans = asinh_trans()) +
  facet_wrap(~ facet_label, scales = "fixed", ncol = ncol_facet) +
  labs(
    title = "Stem CH\u2084 flux vs. height (projects with measurements \u2265 2 m)",
    subtitle = "Fixed arcsinh-transformed flux axis",
    y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})~"[arcsinh scale]")
  )

ggsave("flux_by_height_asinh_fixed.pdf", p3, width = plot_w, height = plot_h, limitsize = FALSE)

cat("\nSaved:\n")
cat("  flux_by_height_free.pdf\n")
cat("  flux_by_height_asinh_free.pdf\n")
cat("  flux_by_height_asinh_fixed.pdf\n")