#histogram by site, box plot/histogram by species, correlation matrix, 
library(dbplyr)
library(ggplot2)
#Color Codes for Consistency-----------------
#tree spec, bg, hem, rm, ro: "#8674ed", "#3c9e4c","#dba553", "#56B4E9"
#Site, bg, ems: "#4aa6dd", "#eaae59"

#Merging, DF, and Cleaning Data-----------------------------------
clean_data <- read.csv("input/clean_data.csv", header=TRUE)
flux_data <- read.csv("output/fluxdatasummer.csv")

flux_envfull <- merge(flux_data, clean_data, by="UniqueID") 

#filter, meanrow thingy, throw it back in

flux_env <- flux_envfull[-c(89),] %>% #18-2 was not sealed
            filter(between(CH4_r2, .1, 1))
names(flux_env)[27] <- "soilsat" 
            
flux_env[flux_env == "na"] <- NA

flux_env$VWC_1 <- str_replace_all(flux_env$VWC_1,
                                  startsWith(flux_env$VWC_1, "2."), 
                                  NA)

flux_env['VWC_1'][flux_env['VWC_1'] startsWith(flux_env$VWC_1, 2.)] <- NA


dry <- filter(flux_env, soilsat == "n")

       
dry$VWC_1 <- as.numeric(dry$VWC_1)
dry$VWC_2 <- as.numeric(dry$VWC_2)
dry$VWC_3 <- as.numeric(dry$VWC_3)

         
dry_averaged <- dry %>% 
  mutate(VWC_average=select(.,c("VWC_1","VWC_2", "VWC_3")) %>% 
  rowMeans())

flux_env$SiteSpecies <-  interaction(flux_env$Site, flux_env$Species)

flux_EMS <- filter(flux_env, Site== "EMS")

flux_bg <- filter(flux_env, Site == "BG")

rm<- filter(flux_env, Species == "rm")

bg <- filter(flux_env, Species == "bg")

ro <- filter(flux_env, Species == "ro")

hem <- filter(flux_env, Species == "hem")


hemrm <- filter(flux_env, Species %in% c("hem", "rm"))

hemrm$SiteSpecies <- interaction(hemrm$Site, hemrm$Species)

#Summary--------------
flux_summary <- flux_env %>%
  group_by(Species)%>%
  summarize(count=length(CH4_flux),
            mean=mean(CH4_flux),
            min=min(CH4_flux),
            q1=quantile(CH4_flux, prob=.25),
            median=median(CH4_flux),
            q3=quantile(CH4_flux, prob= .75),
            max=max(CH4_flux))
site_summary <- flux_env %>%
  group_by(Site)%>%
  summarize(count=length(CH4_flux),
            mean=mean(CH4_flux),
            min=min(CH4_flux),
            q1=quantile(CH4_flux, prob=.25),
            median=median(CH4_flux),
            q3=quantile(CH4_flux, prob= .75),
            max=max(CH4_flux))

#By Site--------------------------------------------------------

ggplot(flux_env, aes(x=CH4_flux, fill=Site))+
  geom_histogram()

ggplot(flux_env, aes(x=Site, y=CH4_flux, fill=Site))+
  geom_boxplot()

ggplot(justrm, aes(x=Plot, y=CH4_flux, fill=Plot))+
  geom_boxplot()

#By Species--------------------------------
ggplot(flux_env, aes(x=log(CH4_flux*1000), fill=Site,  colour=Site))+
  geom_histogram(alpha=.3, position="identity")+
  scale_fill_manual(values=c("#4aa6dd", "#eaae59"))+
  scale_color_manual(values=c("#4aa6dd", "#eaae59"))+
  labs(
    title= "Methane Flux vs. Site and Species",
    x= "expression(log(CH[4]~Flux)~(nmol/m^2~s))"
  )

ggplot(flux_env, aes(x=log(CH4_flux*1000), fill=Species,  colour=Species))+
  geom_histogram(alpha=.3, position="identity")+
  ylim(0,6)+
  scale_fill_manual(values=c("#8674ed", "#000000","#000000", "#000000"))+
  scale_color_manual(values=c("#8674ed", "#000000","#000000", "#000000"))+
  labs(
    title= "Methane Flux vs. Site and Species",
    x= "log(CH4 Flux) (nmol m^-2 s^-1)"
  )
  
  

ggplot(flux_env, aes(x=Species, y=CH4_flux, fill=Species))+
  geom_boxplot()
#This graph!
ggplot(flux_env, aes(x=Site, y=log(CH4_flux*1000), fill=Species))+
  geom_boxplot() +
  scale_fill_manual(values=c("#8674ed", "#3c9e4c","#dba553", "#56B4E9"))+
  labs(
    title= "Methane Flux vs. Site and Species",
    y= expression(log(CH[4]~Flux)~(nmol/m^2~s))
  )
#not logged, also this one
ggplot(flux_env, aes(x=Site, y=CH4_flux*1000, fill=Species))+
  geom_boxplot() +
  scale_fill_manual(values=c("#8674ed", "#3c9e4c","#dba553", "#56B4E9"))+
  labs(
    title= "Methane Flux vs. Site and Species",
    y= expression(CH[4]~Flux~(nmol/m^2~s))
  )

ggplot(hemrm, aes(x=SiteSpecies, y=CH4_flux, fill=Species))+
  geom_boxplot()

#Moisture--------------
flux_env$Moisture <-  interaction(flux_env$Site, flux_env$soilsat)

ggplot(flux_env, aes(x=soilsat, y=log(CH4_flux), fill=Site))+
  geom_boxplot()

ggplot(flux_env, aes(x=soilsat, y=log(CH4_flux), fill=soilsat))+
  geom_boxplot()

ggplot(flux_EMS, aes(x=soilsat, y=CH4_flux, fill=soilsat))+
  geom_boxplot()

#this graph!
ggplot(dry_averaged, aes(x=VWC_average, y=log(CH4_flux*1000), color=Site))+
  geom_point() +
  geom_smooth(method=lm, se = TRUE)+
  labs(title="Methane Flux vs. Soil Moisture",
       x="Average Moisture, VWC%", y=expression(log(CH[4]~Flux)~(nmol/m^2~s)))+
  scale_color_manual(values=c("#4aa6dd", "#eaae59"))
  
bgvwc <- filter(dry_averaged, Site == "BG")
emsvwc <- filter(dry_averaged, Site == "EMS")

summary(lm(VWC_average ~ CH4_flux, data=bgvwc)) #significant
summary(lm(VWC_average ~ CH4_flux, data=emsvwc)) #notsignificant!

ggplot(bg, aes(x=water_table_height, y=CH4_flux, color=Species))+
  geom_point() +
  geom_smooth(method=lm, se = TRUE)
summary(lm(water_table_height ~ CH4_flux, data=bg))

#By Size------------------
ggplot(hemrm, aes(x=DBH_cm, y=CH4_flux, color=Species))+
  geom_point() +
  geom_smooth(method=lm, se = TRUE)

#This graph!
ggplot(flux_env,aes(x=DBH_cm, y=log(CH4_flux*1000), color=Species))+
  geom_point() +  
  geom_smooth(method=lm, se=TRUE)+
  labs(title="Methane Flux vs. DBH by Species", 
       x="Diameter at Breast Height, DBH cm", y= "log CH4 flux, nmol m^-2 s^-1")+
  scale_color_manual(values=c("#8674ed", "#3c9e4c","#dba553", "#56B4E9"))

bg.size.lm <- lm( DBH_cm ~ CH4_flux, data=bg)
summary(bg.size.lm)

summary(lm(DBH_cm ~ CH4_flux, data=rm))
summary(lm(DBH_cm ~ CH4_flux, data=hem))
summary(lm(DBH_cm ~ CH4_flux, data=ro))

#Time-------------------
weather <- data.frame(jday=180:206,
                     precipitation_mm= c(.3, .3, 0,31, 16, 8.9, 0, 0, 0, 5.6, 1, 42.2, .5, 3.8, 24.6, 4.1,0,27.7,.3, 10.2, 0, 0, 39.9, 0, 0, 4.8, 1.8))
#Include? Send to Jon and Jackie
ggplot()+
  geom_point(data=flux_env, aes(jday,CH4_flux*1000, color=Site)) +
  geom_bar(data=weather, aes(jday, precipitation_mm/2,alpha=.05), stat="identity")+
  scale_y_continuous(sec.axis= sec_axis(~.*2, name= "Precipitation (mm)"))+
  labs(title="Day vs Precipitation and Methane Flux", 
       x="Day of Year", y=expression(CH[4]~Flux~(nmol/m^2~s)), alpha="Precipitation")+
  scale_color_manual(values=c("#4aa6dd", "#eaae59"))
  
summary(lm(jday ~ CH4_flux, data = flux_bg))
summary(lm(jday ~ CH4_flux, data = flux_EMS))


ggplot(bg) +
  geom_smooth(aes(jday, water_table_height), method="lm")

#Statistic Tings --------------------------
#run methane against site species
by.species <- aov(CH4_flux ~ Species, data = flux_env)
summary(by.species)
plot(by.species)

summary(lm(CH4_flux ~ Site, data = flux_env))
summary(lm(CH4_flux ~ SiteSpecies, data = flux_env))
summary(lm(CH4_flux ~ Moisture, data = flux_env))

summary(lm(CH4_flux ~ SiteSpecies:DBH_cm, data = flux_env))

by.site <- aov(CH4_flux ~ Site, data = flux_env)
summary(by.site)

by.sitespecies <- aov(CH4_flux ~ SiteSpecies, data = flux_env)
summary(by.sitespecies)

by.soilsat <- aov(CH4_flux ~ Moisture, data = flux_bg)
summary(by.soilsat)

by.rm <- aov(CH4_flux ~ Site, data = rm)
summary(by.rm)

summary(lm(CH4_flux ~ VWC_average, data=dry_averaged))

library(lme4)
library(nlme)
lm(CH4_flux ~ moisture + species + DBH)
m1 <- lmer(CH4_flux ~ Species * Site * jday+ 
     (1|Tree_Tag), data=flux_env)
summary(m1, digit=4)
m2 <- lmer(CH4_flux ~ Species * DBH_cm  * Site * jday+ 
            (1|Tree_Tag), data=flux_env)
summary(m2)
anova(m2)
#correlation matrix------------------------------------------------
library(GGally)
library(corrr)
library(ggcorrplot)
library(psych)

corinput <- flux_env[, -c(1,2,3,5,6,8,9,10, 12,13,14, 15,16,17, 18, 20, 24,25,26, 28, 30,31,32,33,34)] %>%
            na.omit(flux_env)
  #add vwc average data
corinput['soilsat'][corinput['soilsat'] == "n"] <- 0
corinput['soilsat'][corinput['soilsat'] == "y"] <- 1
corinput['Site'][corinput['Site'] == "EMS"] <- 0
corinput['Site'][corinput['Site'] == "BG"] <- 1
corinput <- sapply(corinput,as.numeric)
  
normal <- scale(corinput)
#correlation matrix
cor_matx <- cor(normal)
ggcorrplot(cor_matx, hc.order = TRUE, type = "lower", lab=TRUE) +
  scale_fill_gradient2(low ="pink3", high = "blue3",name = "") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(fill = guide_colourbar(barwidth = 10, barheight = .8)) +
  theme(legend.direction="horizontal")
#p.mat is significance , p.mat = corr.test(cor_matx)$p, lab=True is numbers

#Canopy!----------------
clean_canopy<- read.csv('input/times_key_canopy.csv')
flux_c <-merge(flux_canopy, clean_canopy, by="UniqueID")
stemheight <- filter(flux_c, Type == "stem")
tree2 <- filter(flux_c, Tree_tag=="2", Type %in% c("branch", "stem"))
leaves <- filter(flux_c, Type == "leaf")

#Do these if we get more data points
ggplot(flux_c[-c(10),], aes(x=Type, y=CH4_flux, fill=Type))+
  geom_boxplot()

stemheight <- stemheight %>% 
  mutate(location = ifelse(Height > 5, "Above", "Below"))

library(plotrix)

stem_grouped <- stemheight %>% 
  group_by(location) %>%
  summarize(mean = mean(CH4_flux*1000),
            se = std.error(CH4_flux*1000))
#Box plot grouped by height
ggplot(stem_grouped, aes(x = location, y = mean, fill = location)) +
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))+
  labs(title="Average Methane Flux Above and Below 5 meters",
       x= "Height relative to 5 meters",
       y= expression(CH[4]~Flux~(nmol/m^2~s)),
       fill="Relative Position")


ggplot(stemheight, aes(x=CH4_flux, y=Height))+
  geom_bar(stat="identity")

ggplot(tree2, aes(x=Height, y=CH4_flux*1000, color=Type))+
  geom_point() +
  geom_smooth(method="lm", color="blue3", se=FALSE, formula= (y ~ exp(-x)), linetype = 2)+
  ylim(0,1.5)+
  labs(
    title="Methane Flux of a Black Gum vs. Height",
    x="Height (m)",
    y=expression(CH[4]~Flux~(nmol/m^2~s))
  )

mean(leaves$CH4_flux)

t.test(leaves$CH4_flux, y = NULL,
       alternative = "less",
       mu = 0, paired = FALSE, var.equal = FALSE,
       conf.level = 0.95)

