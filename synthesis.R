library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)

# ── Read data ──
file_path <- "/Users/jongewirtzman/My Drive/Research/whole_tree_flux/1-s2.0-S0168192324000911-mmc2.xlsx"
upland <- read_xlsx(file_path, sheet = "upland-stem") %>%
  mutate(ecosystem = "Upland")
wetland <- read_xlsx(file_path, sheet = "wetland-stem") %>%
  mutate(ecosystem = "Wetland")
colnames(upland)[1]  <- "Reference"
colnames(upland)[9]  <- "Height_raw"
colnames(upland)[10] <- "CH4_flux"
colnames(wetland)[1]  <- "Reference"
colnames(wetland)[9]  <- "Height_raw"
colnames(wetland)[10] <- "CH4_flux"
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
  suppressWarnings(as.numeric(x))
}
dat <- dat %>% mutate(Height = sapply(Height_raw, parse_height))

# ── Filter to projects with >= 2m measurements ──
projects_2m <- dat %>%
  filter(!is.na(Height), Height >= 2) %>%
  distinct(Reference, ecosystem)
dat_filt <- dat %>%
  semi_join(projects_2m, by = c("Reference", "ecosystem")) %>%
  filter(!is.na(Height))

# Clean reference labels
dat_filt <- dat_filt %>%
  mutate(
    ref_short = str_replace(Reference, "_[A-Za-z]+$", ""),
    ref_short = str_replace_all(ref_short, "_", " "),
    facet_label = paste0(ref_short, " (", ecosystem, ")")
  )

# ── Shared theme ──
shared_theme <- theme_bw(base_size = 11) +
  theme(
    legend.position = c(0.85, 0.04),
    legend.background = element_rect(fill = "white", color = "grey70", linewidth = 0.3),
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 9),
    strip.text = element_text(size = 8, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3, color = "grey90"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    plot.margin = margin(10, 10, 10, 10)
  )

shared_colors <- scale_color_manual(
  values = c("Upland" = "#d7191c", "Wetland" = "#2c7bb6"), name = NULL
)
shared_fill <- scale_fill_manual(
  values = c("Upland" = "#d7191c", "Wetland" = "#2c7bb6"), name = NULL
)

n_facets <- length(unique(dat_filt$facet_label))
plot_h <- ceiling(n_facets / 4) * 3.2 + 0.8

# ── Plot 1: LOESS ──
p_loess <- ggplot(dat_filt, aes(x = Height, y = CH4_flux, color = ecosystem, fill = ecosystem)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_point(alpha = 0.45, size = 1.4, shape = 16) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 0.7,
              alpha = 0.12, na.rm = TRUE) +
  shared_colors + shared_fill +
  facet_wrap(~ facet_label, scales = "free_x", ncol = 4) +
  coord_flip() +
  labs(
    x = "Sampling height (m)",
    y = expression(Stem~CH[4]~flux~(nmol~m^{-2}~s^{-1}))
  ) +
  shared_theme

ggsave("flux_by_height_loess.pdf", p_loess, width = 14, height = plot_h, limitsize = FALSE)
cat("Saved: flux_by_height_loess.pdf\n")

# ── Plot 2: Exponential fit ──
# Fit CH4_flux ~ a * exp(b * Height) + c per facet
# Use nls with start values estimated from the data

# Generate exponential predictions per group
exp_preds <- dat_filt %>%
  group_by(facet_label, ecosystem) %>%
  group_modify(~ {
    d <- .x
    if (nrow(d) < 4) return(tibble())
    
    height_seq <- seq(min(d$Height), max(d$Height), length.out = 100)
    
    # Try nls with exponential: flux = a * exp(b * Height) + c
    # Use robust starting values
    tryCatch({
      # Estimate starting values
      flux_range <- range(d$CH4_flux)
      a_start <- diff(flux_range)
      c_start <- min(d$CH4_flux)
      b_start <- -0.5  # expect decay with height
      
      fit <- nls(CH4_flux ~ a * exp(b * Height) + c, data = d,
                 start = list(a = a_start, b = b_start, c = c_start),
                 control = nls.control(maxiter = 200, warnOnly = TRUE))
      
      pred <- predict(fit, newdata = data.frame(Height = height_seq))
      tibble(Height = height_seq, CH4_flux_pred = pred)
    }, error = function(e) {
      # Fallback: simpler 2-parameter exponential flux = a * exp(b * Height)
      tryCatch({
        a_start <- max(abs(d$CH4_flux))
        b_start <- -0.3
        fit <- nls(CH4_flux ~ a * exp(b * Height), data = d,
                   start = list(a = a_start, b = b_start),
                   control = nls.control(maxiter = 200, warnOnly = TRUE))
        pred <- predict(fit, newdata = data.frame(Height = height_seq))
        tibble(Height = height_seq, CH4_flux_pred = pred)
      }, error = function(e2) {
        tibble()
      })
    })
  }) %>%
  ungroup()

p_exp <- ggplot(dat_filt, aes(x = Height, y = CH4_flux, color = ecosystem)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_point(alpha = 0.45, size = 1.4, shape = 16) +
  shared_colors +
  facet_wrap(~ facet_label, scales = "free_x", ncol = 4) +
  coord_flip() +
  labs(
    x = "Sampling height (m)",
    y = expression(Stem~CH[4]~flux~(nmol~m^{-2}~s^{-1}))
  ) +
  shared_theme

# Add exponential curves if predictions exist
if (nrow(exp_preds) > 0) {
  p_exp <- p_exp +
    geom_line(data = exp_preds,
              aes(x = Height, y = CH4_flux_pred, color = ecosystem),
              linewidth = 0.7, alpha = 0.85, na.rm = TRUE)
}

ggsave("flux_by_height_exp.pdf", p_exp, width = 14, height = plot_h, limitsize = FALSE)
cat("Saved: flux_by_height_exp.pdf\n")
cat("Projects:", n_facets, "| Observations:", nrow(dat_filt), "\n")
cat("Exponential fits generated for",
    length(unique(exp_preds$facet_label)), "of", n_facets, "panels\n")