# Load necessary libraries
library(tidyverse)
library(ggplot2)
library(broom)

# Load the data files
O2_standards <- read.csv("/Users/jongewirtzman/Google Drive/Research/YMF Tree Microbiomes & Methane/Tree Methane Lab/YMF Black Oak/gas/o2_standards.csv")
GHG_standards <- read.csv("/Users/jongewirtzman/Google Drive/Research/YMF Tree Microbiomes & Methane/Tree Methane Lab/YMF Black Oak/gas/ghg_standards.csv")
GC_data <- read.csv("/Users/jongewirtzman/Google Drive/Research/YMF Tree Microbiomes & Methane/Tree Methane Lab/YMF Black Oak/gas/quve_gc_data.csv")


GC_data$O2.Area[which(GC_data$Sample.ID=="N2")]<-NA

# Rename columns for easier processing (remove spaces)
colnames(GC_data) <- make.names(colnames(GC_data))

# Merge the O2 and GHG standards with the GC data based on Sample/Sample ID
# O2 Standards Merge
O2_data <- O2_standards %>%
  left_join(GC_data, by = c("Sample" = "Sample.Type"))

# GHG Standards Merge
GHG_data <- GHG_standards %>%
  left_join(GC_data, by = c("Sample" = "Sample.Type"))

# Convert relevant columns to numeric to avoid coercion issues
# Remove commas from concentration columns and convert to numeric
O2_data$O2.Area <- as.numeric(O2_data$O2.Area)
O2_data$O2_ppm <- as.numeric(gsub(",", "", O2_data$X.O2...ppm.))  # Convert O2 concentration to numeric

GHG_data$N2O.Area <- as.numeric(GHG_data$N2O.Area)
GHG_data$N2O_ppm <- as.numeric(gsub(",", "", GHG_data$X.N2O...ppm.))  # Convert N2O concentration to numeric

GHG_data$CH4.Area <- as.numeric(GHG_data$CH4.Area)
GHG_data$CH4_ppm <- as.numeric(gsub(",", "", GHG_data$X.CH4...ppm.))  # Convert CH4 concentration to numeric

GHG_data$CO2.Area <- as.numeric(GHG_data$CO2.Area)
GHG_data$CO2_ppm <- as.numeric(gsub(",", "", GHG_data$X.CO2...ppm.))  # Convert CO2 concentration to numeric

# Filter out any NAs or missing values
GHG_data <- GHG_data %>%
  filter(!is.na(N2O.Area) & !is.na(N2O_ppm)) %>%
  filter(!is.na(CH4.Area) & !is.na(CH4_ppm)) %>%
  filter(!is.na(CO2.Area) & !is.na(CO2_ppm))

# Define the low-range and high-range standards for O2
low_range_o2_standards <- low_range_o2_standards <- c("Outdoor Air 2", "Outdoor Air 3", "Outdoor Air 4",
                            "Oxygen Standard 1", "Oxygen Standard 2", "Oxygen Standard 3", "Oxygen Standard 4")

high_range_o2_standards <- c("Outdoor Air 5", "Outdoor Air 6", "Outdoor Air 7", 
                             "Oxygen Standard 5")

# Filter O2 data to separate low-range and high-range standards
low_range_o2_data <- O2_data %>%
  filter(Sample %in% low_range_o2_standards)

high_range_o2_data <- O2_data %>%
  filter(Sample %in% high_range_o2_standards)

# Combine both low and high range standards for high-range curve
all_high_range_o2_data <- bind_rows(low_range_o2_data, high_range_o2_data)

# Create low-range and high-range standard curves for O2
O2_low_curve <- lm(O2_ppm ~ O2.Area, data = low_range_o2_data)
O2_high_curve <- lm(O2_ppm ~ O2.Area, data = all_high_range_o2_data)

# Get the maximum peak area for the low-range O2 curve
max_O2_low_area <- max(low_range_o2_data$O2.Area, na.rm = TRUE)

# Define the low-range and high-range standards for GHGs
low_range_standards <- c("N2", "SB1", "SB2", "SB3")
high_range_standards <- c("SB4", "SB5", "SB6")

# Filter GHG data to separate low-range and high-range standards
low_range_data <- GHG_data %>%
  filter(Sample %in% low_range_standards)

high_range_data <- GHG_data %>%
  filter(Sample %in% high_range_standards)

# Create low-range and high-range standard curves for low-range and all data for high-range
N2O_low_curve <- lm(N2O_ppm ~ N2O.Area, data = low_range_data)
N2O_high_curve <- lm(N2O_ppm ~ N2O.Area, data = GHG_data)

CH4_low_curve <- lm(CH4_ppm ~ CH4.Area, data = low_range_data)
CH4_high_curve <- lm(CH4_ppm ~ CH4.Area, data = GHG_data)

CO2_low_curve <- lm(CO2_ppm ~ CO2.Area, data = low_range_data)
CO2_high_curve <- lm(CO2_ppm ~ CO2.Area, data = GHG_data)

# Get the maximum peak areas for the low-range curves
#max_N2O_low_area <- max(low_range_data$N2O.Area, na.rm = TRUE)
max_N2O_low_area <- 1500
#max_CH4_low_area <- max(low_range_data$CH4.Area, na.rm = TRUE)
max_CH4_low_area <- 1000
#max_CO2_low_area <- max(low_range_data$CO2.Area, na.rm = TRUE)
max_CO2_low_area <- 2000

# Function to calculate concentration based on peak area and curve selection
calculate_concentration <- function(peak_area, max_low_area, low_curve, high_curve) {
  if (peak_area <= max_low_area) {
    return(predict(low_curve, newdata = data.frame(peak_area = peak_area)))
  } else {
    return(predict(high_curve, newdata = data.frame(peak_area = peak_area)))
  }
}

# Apply the calculation without rowwise(), using if_else and vectorized operations
GC_data <- GC_data %>%
  mutate(
    O2_concentration = if_else(
      O2.Area <= max_O2_low_area,
      predict(O2_low_curve, newdata = data.frame(O2.Area = O2.Area)),
      predict(O2_high_curve, newdata = data.frame(O2.Area = O2.Area))
    ),
    
    N2O_concentration = if_else(
      N2O.Area <= max_N2O_low_area,
      predict(N2O_low_curve, newdata = data.frame(N2O.Area = N2O.Area)),
      predict(N2O_high_curve, newdata = data.frame(N2O.Area = N2O.Area))
    ),
    CH4_concentration = if_else(
      CH4.Area <= max_CH4_low_area,
      predict(CH4_low_curve, newdata = data.frame(CH4.Area = CH4.Area)),
      predict(CH4_high_curve, newdata = data.frame(CH4.Area = CH4.Area))
    ),
    CO2_concentration = if_else(
      CO2.Area <= max_CO2_low_area,
      predict(CO2_low_curve, newdata = data.frame(CO2.Area = CO2.Area)),
      predict(CO2_high_curve, newdata = data.frame(CO2.Area = CO2.Area))
    )
  )


# Function to plot standard curve with equation and R^2
plot_standard_curve <- function(low_data, high_data, low_curve, high_curve, analyte, x, y, range_type) {
  low_summary <- summary(low_curve)
  high_summary <- summary(high_curve)
  
  low_equation <- paste("Low:", round(coef(low_curve)[2], 3), "*x +", round(coef(low_curve)[1], 3))
  high_equation <- paste("High:", round(coef(high_curve)[2], 3), "*x +", round(coef(high_curve)[1], 3))
  
  low_r2 <- paste("R² (low) =", round(low_summary$r.squared, 3))
  high_r2 <- paste("R² (high) =", round(high_summary$r.squared, 3))
  
  ggplot() +
    geom_point(data = high_data, aes_string(x = x, y = y), color = "red") +
    geom_smooth(data = high_data, aes_string(x = x, y = y), method = "lm", se = FALSE, color = "red") +
    geom_point(data = low_data, aes_string(x = x, y = y), color = "blue") +
    geom_smooth(data = low_data, aes_string(x = x, y = y), method = "lm", se = FALSE, color = "blue") +
    annotate("text", x = Inf, y = Inf, label = paste(low_equation, low_r2, sep = "\n"), 
             hjust = 1.2, vjust = 2, size = 3.5, color = "blue", parse = FALSE) +
    annotate("text", x = Inf, y = Inf, label = paste(high_equation, high_r2, sep = "\n"), 
             hjust = 1.2, vjust = 1, size = 3.5, color = "red", parse = FALSE) +
    labs(title = paste(analyte, range_type, "Standard Curve"),
         x = "Peak Area",
         y = "Concentration (ppm)") +
    theme_minimal()
}

# Create separate plots for low-range and high-range standards for each analyte
# O2 Low-Range and High-Range Curves
O2_low_plot <- plot_standard_curve(low_range_o2_data, O2_data, O2_low_curve, O2_high_curve, "O2", "O2.Area", "O2_ppm", "Low-Range")
O2_high_plot <- plot_standard_curve(all_high_range_o2_data, O2_data, O2_low_curve, O2_high_curve, "O2", "O2.Area", "O2_ppm", "High-Range")

# Low-Range Curves for GHGs
N2O_low_plot <- plot_standard_curve(low_range_data, GHG_data, N2O_low_curve, N2O_high_curve, "N2O", "N2O.Area", "N2O_ppm", "Low-Range")
CH4_low_plot <- plot_standard_curve(low_range_data, GHG_data, CH4_low_curve, CH4_high_curve, "CH4", "CH4.Area", "CH4_ppm", "Low-Range")
CO2_low_plot <- plot_standard_curve(low_range_data, GHG_data, CO2_low_curve, CO2_high_curve, "CO2", "CO2.Area", "CO2_ppm", "Low-Range")

# High-Range Curves for GHGs (including all standards)
N2O_high_plot <- plot_standard_curve(GHG_data, GHG_data, N2O_low_curve, N2O_high_curve, "N2O", "N2O.Area", "N2O_ppm", "High-Range")
CH4_high_plot <- plot_standard_curve(GHG_data, GHG_data, CH4_low_curve, CH4_high_curve, "CH4", "CH4.Area", "CH4_ppm", "High-Range")
CO2_high_plot <- plot_standard_curve(GHG_data, GHG_data, CO2_low_curve, CO2_high_curve, "CO2", "CO2.Area", "CO2_ppm", "High-Range")

# Display the plots in RStudio
print(O2_low_plot)
print(O2_high_plot)

print(N2O_low_plot)
print(N2O_high_plot)

print(CH4_low_plot)
print(CH4_high_plot)

print(CO2_low_plot)
print(CO2_high_plot)


# Function to plot standard curve with equation, R^2, and axis limits for low-range
plot_low_range_curve <- function(low_data, low_curve, analyte, x, y) {
  low_summary <- summary(low_curve)
  low_equation <- paste("Low:", round(coef(low_curve)[2], 3), "*x +", round(coef(low_curve)[1], 3))
  low_r2 <- paste("R² (low) =", round(low_summary$r.squared, 3))
  
  # Define the x-axis limits based on the low-range data
  x_min <- min(low_data[[x]], na.rm = TRUE)
  x_max <- max(low_data[[x]], na.rm = TRUE)
  
  ggplot() +
    geom_point(data = low_data, aes_string(x = x, y = y), color = "blue") +
    geom_smooth(data = low_data, aes_string(x = x, y = y), method = "lm", se = FALSE, color = "blue") +
    annotate("text", x = Inf, y = Inf, label = paste(low_equation, low_r2, sep = "\n"), 
             hjust = 1.2, vjust = 2, size = 3.5, color = "blue", parse = FALSE) +
    scale_x_continuous(limits = c(x_min, x_max)) +  # Set x-axis limits for low-range
    labs(title = paste(analyte, "Standard Curve (Low-Range)"),
         x = "Peak Area",
         y = "Concentration (ppm)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels
}

# Low-Range Plots for O2, N2O, CH4, and CO2
O2_low_plot <- plot_low_range_curve(low_range_o2_data, O2_low_curve, "O2", "O2.Area", "O2_ppm")
N2O_low_plot <- plot_low_range_curve(low_range_data, N2O_low_curve, "N2O", "N2O.Area", "N2O_ppm")
CH4_low_plot <- plot_low_range_curve(low_range_data, CH4_low_curve, "CH4", "CH4.Area", "CH4_ppm")
CO2_low_plot <- plot_low_range_curve(low_range_data, CO2_low_curve, "CO2", "CO2.Area", "CO2_ppm")

O2_low_plot
N2O_low_plot
CH4_low_plot
CO2_low_plot
