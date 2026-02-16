canopy1<-read.csv("/Users/jongewirtzman/Downloads/fluxdata_canopy_lgr1.csv")

canopy2<-read.csv("/Users/jongewirtzman/Downloads/fluxdata_canopy_lgr2.csv")

canopy3<-read.csv("/Users/jongewirtzman/Downloads/fluxdata_canopy_lgr3.csv")

canopy_all<-rbind(canopy1, canopy2, canopy3)

field<-read.csv("/Users/jongewirtzman/Downloads/Field Data Entry - Clean Canopy Lift Total.csv")

library(tidyverse)

canopy_data<-left_join(field, canopy_all, by="UniqueID")

canopy_data<-canopy_data[-which(canopy_data$Tree_Tag==3),]
canopy_data$Type[which(canopy_data$Type=="leaf (shaded)")]<-"leaf"

canopy_data$CH4_flux<-canopy_data$CH4_flux*1000

canopy_data <- canopy_data %>%
  mutate(name = case_when(
    Tree_Tag == 2 ~ "Swamp Blackgum",
    Tree_Tag == 4 ~ "Swamp Hemlock",
    Tree_Tag == 5 ~ "Swamp Red Maple",
    Tree_Tag == 300607 ~ "Upland Red Oak",
    Tree_Tag == 321071 ~ "Upland Red Maple",
    Tree_Tag == 321902 ~ "Upland Hemlock",
    TRUE ~ NA_character_  # Default case if no match is found
  ))

# Define the desired order
desired_order <- c("Swamp Blackgum", "Swamp Hemlock",  "Swamp Red Maple",
                   "Upland Red Oak",  "Upland Hemlock", "Upland Red Maple")

# Ensure all desired levels are present (even if some might not appear in your data)
canopy_data$name <- factor(canopy_data$name, levels = desired_order)

canopy_data <- canopy_data %>% rename(Tissue = Type)


###

ggplot(canopy_data,
       aes(y=Height_m, x=log(CH4_flux+1), color=Tissue))+
  #geom_smooth(method="lm", formula= (y~log(x)))+
  geom_point()+
  facet_wrap(~name,  scales="free_x")+
  xlab(expression("Log of" ~ CH[4] ~ "Flux" ~ (nmol ~ m^-2 ~ s^-1)))+ylab("Height Above Ground (m)")

ggplot(canopy_data,
       aes(y=Height_m, x=CH4_flux, color=Tissue))+
  #geom_smooth(method="lm", formula= (y~log(x)))+
  geom_point()+
  facet_wrap(~name)+
  ylim(0,25)+ 
  xlab(expression(CH[4] ~ "Flux" ~ (nmol ~ m^-2 ~ s^-1)))+ylab("Height Above Ground (m)")


ggplot(canopy_data,
       aes(y=Height_m, x=CH4_flux, color=Tissue))+
  #geom_smooth(method="lm", formula= (y~log(x)))+
  geom_vline(xintercept = 2, linetype = "dashed")+
  geom_point()+
  facet_wrap(~name,  scales="free_x")+
  xlab(expression(CH[4] ~ "Flux" ~ (nmol ~ m^-2 ~ s^-1)))+ylab("Height Above Ground (m)")


ggplot(canopy_data,
       aes(y=Height_m, x=CH4_flux, color=Tissue))+
  #geom_smooth(method="lm", formula= (y~log(x)))+
  geom_vline(xintercept = 2, linetype = "dashed")+
  geom_point()+
  facet_wrap(~name)+
  xlab(expression(CH[4] ~ "Flux" ~ (nmol ~ m^-2 ~ s^-1)))+ylab("Height Above Ground (m)")+
  theme_minimal()


ggplot(canopy_data,
       aes(y=Height_m, x=log(CH4_flux), color=Tissue))+
  #geom_smooth(method="lm", formula= (y~log(x)))+
  geom_vline(xintercept = log(2), linetype = "dashed")+
  geom_point()+
  facet_wrap(~name,  scales="free_x")+
  xlab(expression("Log of" ~ CH[4] ~ "Flux" ~ (nmol ~ m^-2 ~ s^-1)))+ylab("Height Above Ground (m)")+
  theme_minimal()
