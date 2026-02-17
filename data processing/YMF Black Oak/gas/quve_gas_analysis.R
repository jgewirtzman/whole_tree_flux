quve_gas<-GC_data %>% filter(is.na(Lab.ID)==F)

int_gas<-quve_gas%>%filter(Tree.Tissue=="Trunk Gas")
int_gas$Tree.Height<-as.numeric(int_gas$Tree.Height)


incubation<-quve_gas%>%filter(Tree.Tissue!="Trunk Gas")
aerobic<-quve_gas%>%filter(Sample.Type=="Aerobic (UZA)")
anaerobic<-quve_gas%>%filter(Sample.Type=="Anaerobic (N2)")


mcra<-read.csv("/Users/jongewirtzman/Downloads/black_oak_aligned_mcrA_FAM_Oak1_Ben_Meso_20221114_171135_766.csv")
mcra_clean<- mcra %>% filter(Component %in% c("Heartwood", "Sapwood"))

flux <- c(0.023, 0.066, 0.024, 0.276, 0.030, 0.031, 0.002)
flux_df<-data.frame(height=int_gas$Tree.Height, flux=rev(flux), ch4=int_gas$CH4_concentration)

# Define common breaks for y-axis (Tree.Height)
height_breaks <- sort(unique(int_gas$Tree.Height))

# Plot 1: O2 concentration scaled from red (high) to blue (low)
p1a <- ggplot(int_gas, aes(x = Tree.Height, y = CH4_concentration)) +
  geom_smooth(se = FALSE, color = "black") +
  geom_point(aes(fill = O2_concentration / 1e6 * 100), size = 3, shape = 21, color = "black", stroke = 1, alpha=0.8) +
  theme_classic() +
  scale_fill_viridis_c(option="inferno", direction=-1, name = "O2 (%)") +  # Blue to red scale
  coord_flip() +
  ylab("CH4 (ppm)") +
  xlab("Height (m)") +
  scale_x_continuous(breaks = height_breaks, minor_breaks = NULL) +  # Consistent breaks
  theme(legend.position = "right")

# Plot 1: O2 concentration scaled from red (high) to blue (low)
p1b <- ggplot(int_gas, aes(x = Tree.Height, y = O2_concentration)) +
  geom_smooth(se = FALSE, color = "black") +
  geom_point(aes(fill = CH4_concentration / 1e6 * 100), size = 3, shape = 21, color = "black", stroke = 1, alpha=0.8) +
  theme_classic() +
  scale_fill_viridis_c(option="inferno", direction=-1, name = "CH4 (%)") +  # Blue to red scale
  coord_flip() +
  ylab("O2 (ppm)") +
  xlab("Height (m)") +
  scale_x_continuous(breaks = height_breaks, minor_breaks = NULL) +  # Consistent breaks
  theme(legend.position = "right")

# Plot 2: Heartwood and Sapwood colors updated
p2 <- ggplot(mcra_clean, aes(y = Conc.copies.µL., x = Height..cm./100)) +
  geom_smooth(se = FALSE, aes(color=Component)) +
  geom_jitter(size=3, shape = 21, aes(fill = Component), color = "black", stroke = 1, alpha=0.8) +
  theme_classic() +
  scale_color_manual(values = c("Heartwood" = "#a6611a", "Sapwood" = "#1f78b4")) +
  scale_fill_manual(values = c("Heartwood" = "#a6611a", "Sapwood" = "#1f78b4")) +
  coord_flip() +
  ylab("mcrA (copies/µL)") +
  xlab("") +
  scale_x_continuous(breaks = height_breaks, minor_breaks = NULL) +  # Consistent breaks
  theme(legend.position = "right")

# Plot 3: Flux with CH4 color scale
p3 <- ggplot(flux_df, aes(y = flux, x = height, fill=ch4)) +
  geom_smooth(se = FALSE, color="black") +
  geom_point(size=3, shape = 21, color = "black", stroke = 1, alpha=0.8) +
  theme_classic() +
  scale_fill_viridis_c(option="E", direction=-1, name = "CH4 (ppm)") +  # Blue to red scale
  coord_flip() +
  xlab("") +
  ylab("CH4 Flux") +
  scale_x_continuous(breaks = height_breaks, minor_breaks = NULL) +  # Consistent breaks
  theme(legend.position = "right")

library(patchwork)
# Combine the three plots with patchwork and consistent grid breaks
p1a + p2 + p3 + plot_layout(nrow = 1, guides = "collect") &  theme(legend.position = "bottom")



###


ggplot(aerobic, aes(y=CH4_concentration, x=Tree.Height))+
  geom_point()+
  geom_line()+
  facet_wrap(~Sample.Type)+
  scale_color_viridis_d(option="E")+
  coord_flip()+
  facet_wrap(~Tree.Tissue)

ggplot(anaerobic, aes(y=CH4_concentration, x=Tree.Height))+
  geom_point()+
  geom_line()+
  facet_wrap(~Sample.Type)+
  scale_color_viridis_d(option="E")+
  coord_flip()+
  facet_wrap(~Tree.Tissue)


ggplot(int_gas, aes(y=(CH4_concentration), x=(O2_concentration), color=N2O_concentration))+
  geom_point()+
  geom_smooth(method="lm")

mcra_clean$height<-mcra_clean$Height..cm./100
all_data<-left_join(mcra_clean, flux_df, by="height")

library(ggpmisc)
ggplot(all_data, aes(y=ch4, x=Copies.20µLWell, color=flux))+
  geom_point()+
  geom_smooth(method="lm") +
  stat_poly_eq(aes(label = paste(..rr.label..)), 
               formula = y ~ x, 
               parse = TRUE,
               label.x = "right",
               label.y = "top")+
  theme_classic()

