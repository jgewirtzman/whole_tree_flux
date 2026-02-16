# This function reads a raw text output file from a Picarro GasScouter Analyzer and exports the 
# time, CO2, and CH4 concentration data. 
#
# The format_Picarro_output() function requires two inputs:
#   1. raw_files = a vector of raw file names to be processed
# Output of this function is a data frame with measurement times and CO2 & CH4 concentration data

format_Picarro_output <- function(raw_files, date_time){
  
  # Aggregate the files if there is more than one file per day.
  if(length(raw_files) > 1){
    for(f in 1:length(raw_files)){ 
      
      # Grab date from file name
      date = lubridate::ymd(substr(gsub("\\D", "", raw_files[f]),2,9))
      
      # Read in raw text Picarro data for each file
      dat  <- read.table(raw_files[f], header=TRUE) %>% 
        mutate(times2 = lubridate::ymd_hms(paste0(date, TIME)))
      
      # Format time column #times(dat$TIME)
      Pic_tmp_times <- dat %>% 
        mutate(times2 = lubridate::ymd_hms(paste0(date, TIME)))
        
      # Aggregate each file into a big daily file
      if(f == 1){
        Pic_data  <- dat
        #Pic_times <- Pic_tmp_times
      } else {
        Pic_data  <- rbind(Pic_data,dat)
        #Pic_times <- c(Pic_times,Pic_tmp_times)
      }
    }
  } else {
    Pic_data  <- read.table(raw_files, header=TRUE)
    #Pic_times <- times(Pic_data$TIME)
  }
  
  # Format output for times, CO2, and CH4
  Pic_df <- data.frame(times = Pic_data$times2,
                        CO2 = Pic_data$CO2_dry, CH4 = Pic_data$CH4_dry)
  return(Pic_df)
}

