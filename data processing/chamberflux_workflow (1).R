# Workflow for processing raw LGR/Picarro data and calculating
# chamber flux on different dates at different soil chambers

library(tidyverse)
library(lubridate)
library(data.table)

# Load local functions
file_sources <- list.files("/Users/jongewirtzman/Downloads/functions/", pattern="*.R", full.names = TRUE)
sapply(file_sources, source, .GlobalEnv)

# File names for raw data from the LGR/Picarro analyzers
data_path <- "/Users/jongewirtzman/Downloads/input/LGR3" # local path to raw data
raw_files <- list.files(data_path, full.names=TRUE)

# # Load file with sampling dates and start times for flux measurements
# FM test time file - Picarro # might not need to rename if rep_vol_L in orig

# JG test time file - LGR
date_time <- read_csv("/Users/jongewirtzman/Downloads/times_key_tree - Canopy Lift_LGR3 (2) copy.csv") %>%
  mutate(dates = lubridate::mdy(Date),
         UniqueID = paste(UniqueID),
         start_time = lubridate::ymd_hms(paste0(dates,comp_start_time)),
         end_time = lubridate::ymd_hms(paste0(dates,comp_end_time)),
         rep_vol_L = Rep_vol_L) %>%
  filter(!is.na(comp_start_time)) 
  #UniqueID = paste(UniqueID) #trees
#UniqueID = paste(Tree_tag, Height, Type, repnumber),

# Flux processing settings - change these for your application
init <- list()
init$analyzer <- "lgr" # can be "picarro" or "lgr"
init$data_path <- data_path # path to analyzer files
init$startdelay <- 0 # 20s delay for Picarro
init$fluxend   <- 3 # minutes to include data after start (will ignore if end times are in date_time)
#init$surfarea  <- pi*(4*2.54/2)^2 / 100^2 #m^2, 4-inch pvc collars
init$vol_system <- .2 # interior volume of LGR = 1.969 L included in data sheet ??? who wrote this, actually LGR V = .2L
init$plotslope <- 1 # make a plot with the slope: 0 = off, save images?? pdf (looP) dev.off
init$outputfile <- 1 # write an output file: 0 = off
init$outfilename <- "fluxdata_lgr3.csv"

# Calculate soil CO2 & CH4 flux for each measurement date & replicate
#pdf("Fluxplots.pdf")
flux_data_monthly <- calculate_chamber_flux(raw_files, date_time, init)   
#dev.off()



# AFTER THIS YOU CAN EDIT TO MERGE WHATEVER OTHER DATA YOU WANT
# by the UniqueID column
########
# Merge temperature & moisture dataset with clean flux data
temp_moist  <- read.csv("input/TF_insttempmoisture.csv", 
                        header=TRUE, stringsAsFactors = FALSE) 

temp_moist$date <- as.Date(temp_moist$date, "%m/%d/%y") %>% 
  format(.,"%Y-%m") %>% 
  as.character()

flux_env <- merge(flux_clean, temp_moist, by.x = c("fmonth", "fid"),
                     by.y = c("date","collar"), all=TRUE) 

