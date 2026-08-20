# MSc Project: A little sedge goes a long way: emergent vegetation as a key driver of littoral methane flux

# Function for merging and cleaning relevant datasets (LICOR single inject,
# Environmental Variables, BCWS Weather)

# 07/15/2025 (Last Edit: 08/17/2026)
# Author: Kelsey McGuire
# kmcgu@mail.ubc.ca; kmcguire.9@outlook.com

# LOAD RELEVANT PACKAGES/LOCAL FUNCTIONS
library(dplyr)
library(lubridate)
library(zoo)
library(purrr)
library(here)

source(here::here("code_DataCleaning", "Cawley_Dissolved_Gases_Function_2018.R"))
source(here::here("code_DataCleaning", "EVCH4_SlopeCalc.R"))

# reduce package conflicts which use the same function to perform different tasks
library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

CleanLICOR <- function(csv_filepath) {
  # ------------------------------------------------------------------------ #
  # Arguments:
  # - csv [str]: .csv file path as a string, differs depending on if the
  #              working directory is set or not - csv should contain the
  #             data from LICOR single-inject processing
  # - colnames_list [List[str]]: List of desired column names for the processed
  #                              dataframe.
  #
  # Purpose:
  # - Cleans up the column names and splits apart ID names for easier visualization
  #   and processing of LICOR data.
  #
  # Output:
  # - dataframe with licor data that now has site data split apart for easier
  #   understanding.
  # ------------------------------------------------------------------------ #

  # read in the licor csv as specified in the argument section and rename columns
  # based on their placement and if they contain CH4 or CO2 in the name
  licor_data <- read.csv(csv_filepath) %>%
    select(Sample.Date, Processed.Date, Sample.ID, LICOR.Interval, Time, CH4Gas_ppm, CO2.conc...ppm., CH4_corrected_ppm, CO2_corrected_ppm, NOTES) %>%
    rename(
      sample_date = Sample.Date,
      processed_date = Processed.Date,
      sample_id = Sample.ID,
      chamber_int = LICOR.Interval,
      time = Time,
      raw_ch4_ppm = CH4Gas_ppm,
      raw_co2_ppm = CO2.conc...ppm.,
      ch4_ppm = CH4_corrected_ppm,
      co2_ppm = CO2_corrected_ppm,
      notes = NOTES
    )

  # get specific site (i.e., 1CO, 3AV)
  licor_data$site <- sapply(strsplit(licor_data[,3], split='-', fixed=T), function(x) (x[1]))
  licor_data$zone <- substr(licor_data$site, 1, 1)

  # classify if the sites are vegetated or open
  licor_data$veg_class <- NA
  for (i in 1:length(licor_data$sample_id)) {
    if (grepl('V', licor_data$sample_id[i], fixed = T) == TRUE) {
      licor_data$veg_class[i] <- 'Vegetated'
    }
    if (grepl('O', licor_data$sample_id[i], fixed = T) == TRUE) {
      licor_data$veg_class[i] <- 'Open'
    }
  }

  # classify which subsample (i.e., A, B, C) plot it is
  licor_data$subsample <- NA
  licor_data$subsample <- ifelse(grepl("A", licor_data$sample_id, fixed = TRUE), "A",
                                 ifelse(grepl("B", licor_data$sample_id, fixed = TRUE), "B",
                                        ifelse(grepl("C", licor_data$sample_id, fixed = TRUE), "C", NA)))

  # classify which analysis type it is
  licor_data$sample_type <- NA
  licor_data$sample_type <- ifelse(grepl("FC", licor_data$sample_id, fixed = TRUE), "Chamber",
                                   ifelse(grepl("DISS", licor_data$sample_id, fixed = TRUE), "Dissolved",
                                          ifelse(grepl("BUB", licor_data$sample_id, fixed = TRUE), "Bubble", NA)))

  # classify which chamber time interval it is
  licor_data$chamber_time <- NA
  licor_data$chamber_time <- ifelse(grepl("T0", licor_data$sample_id, fixed = TRUE), "t0",
                                    ifelse(grepl("T1", licor_data$sample_id, fixed = TRUE), "t1",
                                           ifelse(grepl("T2", licor_data$sample_id, fixed = TRUE), "t2", NA)))

  # multiply the last methane time series value by the dilution factor (5/2) for a more accurate concentration
  for (i in 1:length(licor_data$chamber_time)) {
   if (!is.na(licor_data$chamber_time[i]) && licor_data$chamber_time[i] == "t2") {
     licor_data$ch4_ppm[i] <- licor_data$ch4_ppm[i] * (5/2)
   }
  }

  # add a sample_date column
  licor_data$sample_date <- mdy(licor_data$sample_date)
  
  # filter out high ppm concentrations
  licor_data <- licor_data %>%
    filter(ch4_ppm < 40)

  return(licor_data) # return the processed dataset
}

CleanBubble <- function(bubcsv_file) {
  # ------------------------------------------------------------------------ #
  # Arguments:
  # - csv_file [str]: .csv file path as a string, differs depending on if the
  #              working directory is set or not - should be the csv containing
  #              bubble data from sampling days.
  #
  # Purpose:
  # - Cleans the bubble data so that it may be more easily combined with
  #   the other datasets (i.e., LICOR data, Weather Data ...)
  #
  # Output:
  # - dataframe with bubble data that is now cleaned and ready to be
  #   bound with other datasets
  # ------------------------------------------------------------------------ #
  bub <- read.csv(bubcsv_file)
  
  cleaned <- bub %>%
    dplyr::select(-c(LICOR.Remark, Process.Date, Push.Time, Measure.Time, 
                     CO2.Con...ppm....Baseline, CO2.Con...ppm....Measurement, CH4.Con...ppb....Baseline, 
                     Delt_CH4, Delt_CO2, Vol_Inj, Veff, 
                     CO2.Con...Vial.Dilution., CH4.Con...Vial.Dilution., N.Syringe.Dilution, 
                     CH4.Con...N.Dil.., CO2.Con...N.dil., Notes)) %>%
    rename(zone = Zone, subsample = Subsample, veg_class = Veg.Presence,
           set_datetime = Set.Date.Time, sample_datetime = Sample.Date.Time,
           days_deploy_t = Days.Deployed..t., sample_taken = Sample.Taken.,
           tot_vol_ml = Total.Volume..mL., vol_vial_ml = Volume.in.Vial..mL.,
           funnel_area_m2 = Funnel.Area..m.2., depth_bt = Depth, reset_datetime = Reset.Time,
           co2_ppm = CO2_corrected_ppm, ch4_ppb = CH4.Con...ppb....Measurement, ch4_ppm = CH4_corrected_ppm,
           ch4_ppm_postdil = CH4.Con...N...vial.Dil.., co2_ppm_postdil = CO2.Con...N...vial.Dil., 
           perc_ch4 = X..Methane, perc_co2 = X..Carbon.Dioxide, bub_flux_lm2d1 = Bubble.Flux..L.m.2.d.1., 
           ebull_flux_ch4 = Flux..mg.C.CH4.m.2.d.1.,
           ebull_flux_co2 = Flux..mg.C.CO2.m.2.d.1.) %>%
    mutate(site = paste0(zone, subsample, veg_class),
           veg_class = factor(veg_class,
                              levels = c("O", "V"),
                              labels = c("Open", "Vegetated")),
           set_datetime = mdy_hms(set_datetime),
           sample_datetime = mdy_hms(sample_datetime),
           sample_date = as.Date(sample_datetime),
           reset_datetime = mdy_hms(reset_datetime),
           ebull_flux_ch4 = case_when(ebull_flux_ch4 == "#VALUE!" ~ NA_real_,
                                      TRUE ~ as.numeric(ebull_flux_ch4)),
           ebull_flux_co2 = case_when(ebull_flux_co2 == "#VALUE!" ~ NA_real_,
                                      TRUE ~ as.numeric(ebull_flux_co2))) %>%
    dplyr::select(set_datetime, sample_datetime, reset_datetime, site, everything())
  
  return(cleaned)
}

CleanENV <- function(csv_filepath) {
  # ------------------------------------------------------------------------ #
  # Arguments:
  # - csv [str]: .csv file path as a string, differs depending on if the
  #              working directory is set or not - should be the csv containing
  #              environmental data from sampling days.
  #
  # Purpose:
  # - Cleans the environmental data so that it may be more easily combined with
  #   the other datasets (i.e., LICOR data, Weather Data ...)
  #
  # Output:
  # - dataframe with environmental data that is now cleaned and ready to be
  #   bound with other datasets
  # ------------------------------------------------------------------------ #

  # read in field environmental site data and clean it up
  env_data <- read.csv(csv_filepath)

  # create a site variable to match other datasets
  env_data$site <- paste0(env_data$Zone, env_data$Subsample, env_data$Vegetation.Presence)

  # remove now unnecessary columns and rename them
  env_data <- env_data %>%
    select(!c(Zone, Subsample, Vegetation.Presence, Notes)) %>%
    rename(sample_date = 1, sample_time = 2, depth = 3, ph = 4,
           cond_mScm = 5, tds_mgL = 6, wat_temp_C = 7, air_temp_C = 8,
           press_kPa = 9, DO = 10, site = 11)

  # convert the date and time columns to more useable formats
  env_data$sample_date <- mdy(env_data$sample_date)
  env_data$sample_time <- hms(env_data$sample_time)

  return(env_data) # return slightly cleaned env. data
}

CleanBCWS <- function(folder_path, stationcode) {
  # ------------------------------------------------------------------------ #
  # Arguments:
  # - folder_path [str]: folder path to where all the weather station data from
  #   BCWS is contained (https://www.for.gov.bc.ca/ftp/HPR/external/!publish/BCWS_DATA_MART/2025/)
  # - stationcode [int]: integer number for relevant weather station code.
  #   BCWS files contain all data from their weather stations so best to narrow
  #   it to the relevant stations (i.e., Wheeler Lake uses Boya Station data = 445)
  #
  # Purpose:
  # - Combines the csv's for relevant sampling days from BCWS weather station,
  #   narrows the data to only the relevant station, and cleans the data to
  #   prep it for easier merging between other datasets.
  #
  # Output:
  # - dataframe with weather data from the chosen station which is in proper
  #   format for post-processing and visualization.
  # ------------------------------------------------------------------------ #

  # List all CSV files with full paths
  bcws_files <- list.files(folder_path, full.names = TRUE)

  # Initialize empty dataframe to accumulate all results
  tachila_weather <- data.frame()

  # Loop through each file
  for (csv_path in bcws_files) {
    # Read in one CSV
    data <- read.csv(csv_path)

    tachila_daily_weather <- data %>%
      dplyr::filter(STATION_CODE == stationcode) %>%
      select(c(STATION_CODE, STATION_NAME, DATE_TIME, HOURLY_PRECIPITATION, HOURLY_TEMPERATURE,
                      HOURLY_RELATIVE_HUMIDITY, HOURLY_WIND_SPEED, HOURLY_WIND_DIRECTION, PRECIP_RGT))

    # clean date and time data
    tachila_daily_weather$DATE_TIME <- as.character(tachila_daily_weather$DATE_TIME)
    tachila_daily_weather$DATE <- ymd(str_sub(tachila_daily_weather$DATE_TIME, 1, 8))
    tachila_daily_weather$HOUR <- as.numeric(str_sub(tachila_daily_weather$DATE_TIME, 9, 10))

    # bind all weather csv's together
    tachila_weather <- bind_rows(tachila_weather, tachila_daily_weather)
  }

  return(tachila_weather)
}

CleanWaterSamples <- function(folderpath_watercsvs) {
	# ------------------------------------------------------------------------ #
	# Arguments:
	# - folderpath_watercsvs [str]: folder path to where all relevant water chemistry
	#		csv's (i.e., Anions, Nutrients) are contained. This is based on the csv structure
	#		used by U of Alberta NRAL Lab. 
	#
	# Purpose:
	# - Combines the different water chemistry columns togehter to form one full dataset
	#
	# Output:
	# - Combined water chemistry dataset from all relevant samples.
	# ------------------------------------------------------------------------ #
  water_csvs <- list.files(folderpath_watercsvs, full.names = T) # get full folder path names of csv's in a list
  ## water_csvs will list files alphabetically, so processing follows that order

  anions <- read.csv(water_csvs[1], skip = 7, nrows = 18) # only keeping necessary data rows
  anions <- anions[, -c(8:13)]
  names(anions) <- c('sample_id', 'lab_id', 'F_mgL', 'Cl_mgL', 'Br_mgL', 'NO3-N_mgL', 'S04-S_mgL')

  cn <- read.csv(water_csvs[2], skip = 8, nrows = 18)
  names(cn) <- c('sample_id', 'DOC_mgL', 'TOC_mgL', 'TC_mgL', 'TIC_mgL', 'TN_mgL')

  icpoes <- read.csv(water_csvs[3], skip = 11, nrows = 18) # read in where there is a chemical row and unit row so need to resolve that
  columns_icpoes <- names(icpoes) # get the column names
  columns_icpoes[1] <- 'sample_id' # replace first and second column names with those in second row
  columns_icpoes[2] <- 'lab_id'
  names(icpoes) <- columns_icpoes # rename the columns so that all names are in one row
  icpoes <- icpoes[-1,] # remove additional row

  pn <- read.csv(water_csvs[5], skip = 7, nrows = 18) # skips a number because of there being two ICPOES csv's
  
  names(pn) <- c('sample_id', 'lab_id', 'NH4-N_ugL', 'NO2-N_ugL', 'PO4-P_ugL')

  # water meta data is the google sheet that was sent to the lab that contains the details of the samples
  water_meta <- read.csv('~/Desktop/masters/data/msc-aquatic-ch4/surfwater-samples/ubc-watersurf-aug2025.csv', nrows = 18)
  water_meta <- water_meta %>%
    select(c(Sample.ID.on.Vial, Sample.Date))
  names(water_meta) <- c('sample_id', 'sample_date')
  water_meta$zone <- substr(water_meta$sample_id, 9, 9)
  water_meta$veg_presence <- substr(water_meta$sample_id, 10, 10)

  water_chem <- plyr::join(water_meta, anions, by = 'sample_id')
  water_chem <- plyr::join(water_chem, cn, by = 'sample_id')
  water_chem <- plyr::join(water_chem, pn, by = c('sample_id', 'lab_id'))
  water_chem <- plyr::join(water_chem, icpoes, by = c('sample_id', 'lab_id'))

  water_chem[water_chem == '<LOD'] <- NA # replace any undetected samples with NA for easier processing

  return(water_chem)
}

CleanVeg <- function(field_veg_data, surface_area_csv, rgb_folderpath) {
	# ------------------------------------------------------------------------ #
	# Arguments:
	# - field_veg_data [str]: folder path to the csv containing the field vegetation
	#		data, such as number of shoots, biomass, etc. 
	# - surface_area_csv [str]: folder path of the csv from ImageJ that contains the 
	#		data on percent cover and area measures from field photos.
	# - rgb_folderpath [str]: folder path of the csv from ImageJ that contains all RGB
	# 	band information.
	#
	# Purpose:
	# - To clean and combine all relevant data relating to vegetation metrics specifically
	#		and have it in a useable format
	#
	# Output:
	# - Dataframe that contains any relevant vegetation data (without CH4 fluxes)
	# ------------------------------------------------------------------------ #
  field_veg_data <- read.csv(field_veg_data)
  # print(names(field_veg_data))
  veg_surfarea <- read.csv(surface_area_csv)

  rgb_filesfn <- list.files(rgb_folderpath,
                            full.names = T)
  rgb_files <- list.files(rgb_folderpath)

  rgb_list <- list()
  for (file in 1:length(rgb_filesfn)) {
    file_name <- tools::file_path_sans_ext(rgb_files[file]) # split the filename to get only the date, time, and site info
    # print(paste("Processing:", file_name))
    rgb_data <- read.csv(rgb_filesfn[file])

    split_name <- strsplit(file_name, "-", fixed = TRUE)[[1]]

    # correct datetime into more useable formatting
    datetime_str <- paste0(split_name[1], "-", split_name[2])  # "20250817-133030"
    datetime_val <- as.POSIXct(datetime_str, format = "%Y%m%d-%H%M%S",
                               tz = "America/Los_Angeles")

    # edit site name
    site_raw <- split_name[3]  # "1CTOP"
    site <- substr(site_raw, 1, 2) # extract "1C"

    rgb_data <- rgb_data %>%
      select(Red, Green, Blue) %>%
      summarise(across(where(is.numeric), list(mean = mean,
                                               sd = sd))) %>%
      mutate(datetime = datetime_val,
             site = site,
             GCC = (Green_mean/(Red_mean + Green_mean + Blue_mean)), # Greenness Index/Greenness Chromatic Coordinate
             ExG = (2*Green_mean - Red_mean + Blue_mean)) # Excess Greenness

    # print(head(rgb_data))

    rgb_list[[file_name]] <- rgb_data
  }

  all_rgbs <- bind_rows(rgb_list)

  field_veg_data <- field_veg_data %>%
    select(-c(Photo.Downloaded)) %>%
    filter(!is.na(Time)) %>%
    mutate(datetime = as.POSIXct(paste0(Date, Time), format = "%m/%d/%Y%H:%M:%S",
                                 tz = "America/Los_Angeles"),
           site = paste0(as.character(Zone),Subsample)) %>%
    rename(num_shoots = X..of.Shoots,
           biomass_g = Biomass..g.DW.m.2.,
           abovewater_vol_cm3 = Approx..Vol.cm3..cylinder.volume.AW.,
           total_vol_cm3 = Approx..Vol.cm3..cylinder.volume.Total.) %>%
    select(datetime, site, num_shoots, biomass_g, abovewater_vol_cm3, total_vol_cm3)

  veg_surfarea <- veg_surfarea %>%
    select(-c(X, Area, Mean, StdDev, Min, Max)) %>%
    separate(Label, c("date", "time", "site"), extra = "drop") %>%
    mutate(datetime = as.POSIXct(paste0(date, time), format = "%Y%m%d%H%M%S",
                                 tz = "America/Los_Angeles"),
           site = substr(site, 1, 2),
           perc_cover = case_when(X.Area < 50 ~ X.Area,
                                  X.Area > 50 ~ 100 - X.Area)) %>%
    select(-c(date, time, X.Area))

  veg_data <- left_join(all_rgbs, veg_surfarea, by = c("datetime", "site")) %>%
    right_join(field_veg_data, by = c("datetime", "site")) %>%
    filter(!is.na(num_shoots)) %>%
    mutate(site = paste0(site, "V"),
           date = as.Date(datetime, tz = "America/Los_Angeles")) %>%
    select(site, datetime, everything())

  return(veg_data)
}

CleanHOBOs <- function(hobo_folderpath, hobo_rangemetafile) {
	# ------------------------------------------------------------------------ #
	# Arguments:
	# - hobo_folderpath [str]: folder where relevant HOBO *csv* files are found. 
	#		May need to convert .hobo files to csv before running function.
	# - hobo_rangemetafile [str]: the file that contains the information on the 
	#		date and times in which hobo loggers are launched and removed, to filter 
	#		data that falls beyond that
	#
	# Purpose:
	# - Clean HOBO logger data to make it readable, and then calculate daily, 
	#		3-day and 7-day averages of temperature from the 2-hour interval data. 
	#
	# Output:
	# - List of dataframes with various specified uses including (1) raw dataset
	#		with all 2-hour data; (2) dataset with temperature data based on the
	#		placement (i.e., Air, Sediment); (3) dataset with temperature based on 
	#		specific hobo (i.e., Zone 1 Air Open); (4) dataset with placement as columns.
	# ------------------------------------------------------------------------ #
  # list out the hobo file names including and excluding the path
  hobo_csvs_fn <- list.files(hobo_folderpath, full.names = T)
  hobo_csvs <- list.files(hobo_folderpath)

  # read in the meta data that includes when the hobos were set and taken out
  hobo_meta <- read.csv(hobo_rangemetafile)

  # create a loop for processing individual hobo files
  hobo_list <- list()
  for (file in 1:length(hobo_csvs)) {
    hobo_name <- sapply(strsplit(hobo_csvs[file], split='.', fixed=T), function(x) (x[1])) # split the filename to get only the HOBO name
    hobo_data <- read.csv(hobo_csvs_fn[file], skip = 1, sep = ",", quote = "\"") # read in the data with special case for the separator being quotes and commas (i.e., ",")

    if (grepl("Intensity", names(hobo_data)[4], fixed = TRUE)) { # process based on if it's a light + temp HOBO or just temp
      hobo_data <- hobo_data[,-c(1,5:length(names(hobo_data)))] # remove unnecessary columns - keeping only datetime, temperature, and light
      names(hobo_data) <- c('datetime', 'temp_c', 'light_lux')
    } else {
      hobo_data <- hobo_data[,-c(1, 4:length(names(hobo_data)))] # remove unnecessary columns - keeping only datetime and temperature
      names(hobo_data) <- c('datetime', 'temp_c')
    }

    # correct datetime into more useable formatting, and split into Date and Time data
    hobo_data$datetime <- as.POSIXct(hobo_data$datetime,
                                     format = "%m/%d/%y %I:%M:%S %p")

    # use the hobo meta data to find ranges of deployment time and filter based on those conditions
    hobo_range <- hobo_meta %>%
      filter(placement == hobo_name)
    paste(hobo_range$date_in, hobo_range$time_in)

    datetime_in <- as.POSIXct(paste(hobo_range$date_in, hobo_range$time_in), format = "%m/%d/%Y %H:%M:%S")
    datetime_out <- as.POSIXct(paste(hobo_range$date_out, hobo_range$time_out), format = "%m/%d/%Y %H:%M:%S")

    hobo_list[[hobo_name]] <- hobo_data %>%
      filter(datetime_in <= datetime & datetime <= datetime_out) %>%
      mutate(hobo_name = hobo_name,
             zone = substr(hobo_name, 6, 6),
             placement = substr(hobo_name, 8, (length(hobo_name) - 4)))
  }

  # create a dataframe of the hobos that were processed above into a list
  filtered_combined_hobos <- bind_rows(hobo_list)

  # find the daily average values for temp and lux based on the specific HOBO ...
  daily_average_hobo <- filtered_combined_hobos %>%
    mutate(date = as.Date(datetime)) %>%
    group_by(date, hobo_name) %>%
    summarise(across(where(is.numeric), list(min = min,
                                             max = max,
                                             mean = mean,
                                             sd = sd)))

  # ... or the specific placement (i.e., air, bottom, etc.)
  daily_average_placement <- filtered_combined_hobos %>%
    mutate(date = as.Date(datetime)) %>%
    group_by(date, placement) %>%
    summarise(across(where(is.numeric), list(min = min,
                                             max = max,
                                             mean = mean,
                                             sd = sd)))

  # do the same for the three day ...
  threeday_avg_hobo <- daily_average_hobo %>%
    group_by(hobo_name) %>%
    arrange(date) %>%
    mutate(threeday_avg = map_dbl(seq_along(temp_c_mean), function(i) {
      mean(temp_c_mean[max(1, i-2):i], na.rm = T)
    })) %>%
    ungroup() %>%
    select(date, hobo_name, threeday_avg)

  threeday_avg_placement <- daily_average_placement %>%
    group_by(placement) %>%
    arrange(date) %>%
    mutate(threeday_avg = map_dbl(seq_along(temp_c_mean), function(i) {
      mean(temp_c_mean[max(1, i-2):i], na.rm = T)
    })) %>%
    ungroup() %>%
    select(date, placement, threeday_avg)

  # and seven day averages, making it an average from the past seven days rather than centering the day
  sevenday_avg_hobo <- daily_average_hobo %>%
    group_by(hobo_name) %>%
    arrange(date) %>%
    mutate(sevenday_avg = map_dbl(seq_along(temp_c_mean), function(i) {
      mean(temp_c_mean[max(1, i-6):i], na.rm = T)
    })) %>%
    ungroup() %>%
    select(date, hobo_name, sevenday_avg)

  sevenday_avg_placement <- daily_average_placement %>%
    group_by(placement) %>%
    arrange(date) %>%
    mutate(sevenday_avg = map_dbl(seq_along(temp_c_mean), function(i) {
      mean(temp_c_mean[max(1, i-6):i], na.rm = T)
    })) %>%
    ungroup() %>%
    select(date, placement, sevenday_avg)

  # join together the average daily, 3-day and 7-day temperatures into one dataframe for hobo ...
  average_temps_hobo <- left_join(
    daily_average_hobo,
    threeday_avg_hobo,
    by = c('date', 'hobo_name')
  )

  average_temps_hobo <- left_join(
    average_temps_hobo,
    sevenday_avg_hobo,
    by = c('date', 'hobo_name')
  )

  # ... and for the placement
  average_temps_placement <- left_join(
    daily_average_placement,
    threeday_avg_placement,
    by = c('date', 'placement')
  )

  average_temps_placement <- left_join(
    average_temps_placement,
    sevenday_avg_placement,
    by = c('date', 'placement')
  )

  # classify by zone in the hobo for easier binding
  average_temps_hobo$zone <- substr(average_temps_hobo$hobo_name, 6, 6)

  # create a wider tibble that seperates placements into their own columns for easier processing
  average_temps <- average_temps_placement %>%
    select(c(date, placement, temp_c_mean, light_lux_mean, threeday_avg, sevenday_avg)) %>%
    pivot_wider(names_from = placement,
                values_from = c(temp_c_mean, light_lux_mean, threeday_avg, sevenday_avg),
                names_glue = "{.value}_{placement}") %>%
    select(-c(light_lux_mean_AIR, light_lux_mean_SED_OPEN, light_lux_mean_SED_VEG, light_lux_mean_SUR_VEG))

  # create final list of multiple datasets to be used when needed
  hobo_data <- list(raw_hobo = filtered_combined_hobos,
                    byhobo = average_temps_hobo,
                    byplacement = average_temps_placement,
                    widerplacement = average_temps)

  return(hobo_data)
}

CleanData <- function(licor_csvpath,
                      env_csvpath,
                      bcws_folderpath, stationcode,
                      hobo_folderpath, hobo_rangemetafile,
                      field_veg_data, surface_area_csv, rgb_folderpath,
                      bubcsv_file) {
  # ------------------------------------------------------------------------ #
  # Arguments:
  # - licor_csvpath [str]: .csv file path as a string, differs depending on if
  #   the working directory is set or not - csv should contain the data from
  #   LICOR single-inject processing
  # - env_csvpath [str]: .csv file path as a string, differs depending on if the
  #   working directory is set or not - should be the csv containing
  #   environmental data from sampling days.
  # - bcws_folderpath [str]: folder path to where all the weather station data from
  #   BCWS is contained (https://www.for.gov.bc.ca/ftp/HPR/external/!publish/BCWS_DATA_MART/2025/)
  # - stationcode [int]: integer number for relevant weather station code.
  #   BCWS files contain all data from their weather stations so best to narrow
  #   it to the relevant stations (i.e., Wheeler Lake uses Boya Station data = 445)
	# - hobo_folderpath [str]: folder where relevant HOBO *csv* files are found. 
	#		May need to convert .hobo files to csv before running function.
	# - hobo_rangemetafile [str]: the file that contains the information on the 
	#		date and times in which hobo loggers are launched and removed, to filter 
	#		data that falls beyond that.
	# - field_veg_data [str]: folder path to the csv containing the field vegetation
	#		data, such as number of shoots, biomass, etc. 
	# - surface_area_csv [str]: folder path of the csv from ImageJ that contains the 
	#		data on percent cover and area measures from field photos.
	# - rgb_folderpath [str]: folder path of the csv from ImageJ that contains all RGB
	# 	band information.
	# - bubcsv_file [str]: .csv file path as a string, differs depending on if the
	#		working directory is set or not - should be the csv containing bubble data 
	#		from sampling days.
  #
  # Purpose:
  # - Takes in all relevant CSV's to form one cleaned dataset with relevant data.
	#   For specific data cleaning processes and requirements, look at the relevant 
	#		function descriptions above. 
  #
  # Output:
  # - list containing 5 datasets. (1) raw data is the combined data before any daily
	#		averages or cleaning. (2) clean data is the data that has been averaged on a
	#		daily time scale, has been filtered for outliers, and keeps data from the relevant
	#		sampling days. (3) dissolved only data with environmental variables over the 
	#   entire sampling period. (4) Chamber-based data with environmental variables across
	# 	the sampling period. (5) Bubble only data from the entire sampling period.
  # ------------------------------------------------------------------------ #

  # load in each of the independently cleaned datasets (go to their respective
  # functions for more information)
  licor_data <- CleanLICOR(licor_csvpath)
  env_data <- CleanENV(env_csvpath)
  bcws_weather <- CleanBCWS(bcws_folderpath, stationcode)
  hobo_data <- CleanHOBOs(hobo_folderpath, hobo_rangemetafile)
  veg_data <- CleanVeg(field_veg_data, surface_area_csv, rgb_folderpath)
  bubbledata <- CleanBubble(bubcsv_file)

  # pull out the pivoted hobo data tibble
  hobo_wider <- hobo_data$widerplacement

  # join together the licor and environmental data
  lake_data <- plyr::join(licor_data,
                          env_data,
                          by = c('site', 'sample_date'))
  
  # print("test 1:")
  # print(head(lake_data))

  # Extract hour data
  lake_data <- lake_data %>%
    mutate(sample_hour = hour(sample_time))  # sample_hour will be numeric

  # find the hour from the bcws data - uncomment if you want to include tachila data
  bcws_weather <- bcws_weather %>%
    mutate(weather_hour = as.numeric(str_extract(HOUR, "^\\d+")))

  # join licor, environmental, and BCWS data together
  lake_data <- left_join(lake_data, bcws_weather,
                         by = c("sample_date" = "DATE", "sample_hour" = "weather_hour"))

  # join the hobo and lake data together
  lake_data <- left_join(lake_data,
                         hobo_wider,
                         by = c("sample_date" = "date"))

  lake_data <- left_join(lake_data,
                         veg_data,
                         by = c("sample_date" = "date",
                                "site"))

  # rename and process the temperature columns so that they are named more clearly, and they are sorted by
  # vegetation and open areas
  lake_data <- lake_data %>%
    mutate(
      # Daily mean temperatures
      temp_air_daily = temp_c_mean_AIR,
      temp_bottom_daily = temp_c_mean_BOT,
      temp_surface_daily = ifelse((veg_class == "Vegetated"), temp_c_mean_SUR_VEG, temp_c_mean_SUR), # if the vegetation class is vegetated, apply the surface_veg hobo, if not the sur open hobo
      temp_sed_daily = ifelse((veg_class == "Vegetated"), temp_c_mean_SED_VEG, temp_c_mean_SED_OPEN), # same condition as above but with sediments - this carries on with the 3-day and weekly temperatures

      # 3-day averages
      temp_air_3day = threeday_avg_AIR,
      temp_bottom_3day = threeday_avg_BOT,
      temp_surface_3day = ifelse((veg_class == "Vegetated"), threeday_avg_SUR_VEG, threeday_avg_SUR),
      temp_sed_3day = ifelse((veg_class == "Vegetated"), threeday_avg_SED_VEG, threeday_avg_SED_OPEN),

      # 7-day averages
      temp_air_7day = sevenday_avg_AIR,
      temp_bottom_7day = sevenday_avg_BOT,
      temp_surface_7day = ifelse((veg_class == "Vegetated"), sevenday_avg_SUR_VEG, sevenday_avg_SUR),
      temp_sed_7day = ifelse((veg_class == "Vegetated"), sevenday_avg_SED_VEG, sevenday_avg_SED_OPEN)
    ) %>%
    select(-c(temp_c_mean_AIR, temp_c_mean_BOT, temp_c_mean_SED_OPEN, temp_c_mean_SED_VEG, temp_c_mean_SUR, temp_c_mean_SUR_VEG,
              threeday_avg_AIR, threeday_avg_BOT, threeday_avg_SED_OPEN, threeday_avg_SED_VEG, threeday_avg_SUR, threeday_avg_SUR_VEG,
              sevenday_avg_AIR, sevenday_avg_BOT, sevenday_avg_SED_OPEN, sevenday_avg_SED_VEG, sevenday_avg_SUR, sevenday_avg_SUR_VEG)) # remove any repetitive columns

  # create a dissolved concentration specific dataset to put into the CAWLEY function
  diss_data <- lake_data %>%
    filter(sample_type == 'Dissolved') %>%
    mutate(gasVolume = if_else(sample_date >= as.Date("2025-07-15"), 5, 15), # add in the specific columns requested - NEED TO CHANGE THE VALUES BASED ON DATES WHEN SWITCHED BETWEEN THE 30 TO THE 10 ML
           waterVolume = if_else(sample_date >= as.Date("2025-07-15"), 5, 15),
           barometricPressure = if_else(is.na(press_kPa), 101.33, press_kPa), # NEED TO CONFIRM NUMBERS FOR THIS - MAY USE HOBO FROM WHEELER
           waterTemp = wat_temp_C,
           headspaceTemp = temp_air_daily, # CAN USE HOBO TEMPS WHEN SORTED OUT HOW TO BIND LOL
           concentrationCO2Gas = if_else(chamber_int == 'Measurement', co2_ppm, NA),
           concentrationCO2Air = if_else(chamber_int == 'Baseline', raw_co2_ppm, NA),
           concentrationCH4Gas = if_else(chamber_int == 'Measurement', ch4_ppm, NA),
           concentrationCH4Air = if_else(chamber_int == 'Baseline', raw_ch4_ppm, NA),
           concentrationN2OGas = NA,
           concentrationN2OAir = NA)

  # keep necessary columns for the processing of dissolved concentrations
  diss_data_fxn <- diss_data %>%
    select(c(sample_date, sample_id, time, gasVolume, waterVolume, barometricPressure, waterTemp, headspaceTemp,
    concentrationCO2Gas, concentrationCO2Air, concentrationCH4Gas, concentrationCH4Air, concentrationN2OGas, concentrationN2OAir))

  # use the cawley function to find out the dissolved concentrations of CO2 and CH4 in the water column
  diss_data_calc <- def.calc.sdg.conc(diss_data_fxn) %>%
    select(-c("concentrationN2OGas", "concentrationN2OAir", "dissolvedN2O")) %>% # remove N2O columns since this wasn't measured
    mutate(dissolvedCH4 = (dissolvedCH4*(1000^3)),
           dissolvedCO2 = (dissolvedCO2*(1000^3))) # convert the dissolved values into nmol L-1

  # join the dissolved data to the lake data (licor, environmental, bcws) to the dissolved data that was just calculated
  lake_data <- left_join(
    lake_data,
    diss_data_calc,
    by = c("sample_id", "sample_date", "time", "wat_temp_C" = "waterTemp")
  )

  # clean up the lake data slightly to make it easier to communicate
  cleaned_lake_data <- lake_data %>%
    select(-c(notes, time, sample_hour, STATION_CODE, STATION_NAME, DATE_TIME, HOUR)) %>%
    group_by(processed_date, sample_date, sample_id, site, zone, veg_class,
             subsample, sample_type, chamber_time, chamber_int) %>%
    summarise(across(everything(), ~ first(na.omit(.x)), .names = "{.col}"), .groups = "drop")
  
  # create a dissolved data independent dataset
  dissolved_data <- cleaned_lake_data %>%
    filter(sample_type == "Dissolved")
  
  write.csv(dissolved_data, file = "~/Desktop/masters/data/msc-aquatic-ch4/cleaned-data/wheeler_dissolved_2025.csv")

  # create a chamber independent dataset (diffusive measurements)
  chamber_data <- cleaned_lake_data %>%
    filter(sample_type == "Chamber") %>%
    select(-c(gasVolume, waterVolume, barometricPressure, headspaceTemp, concentrationCO2Gas, concentrationCO2Air,
              concentrationCH4Gas, concentrationCH4Air, dissolvedCO2, dissolvedCH4))

  # calculate the slopes from the chamber measurements
  chamber_data <- ChamberSlopeCalc(chamber_data)

  chamber_data <- chamber_data %>%
    group_by(sample_date, site) %>%
    select(c(sample_date, site, ch4_slope, co2_slope,
             ch4_flux_umolm2s1, ch4_flux_mgm2d1, cch4_flux_mgm2d1, plant_ch4_flux_mgm2d1, plant_cch4_flux_mgm2d1,
             co2_flux_umolm2s1, co2_flux_mgm2d1, cco2_flux_mgm2d1, plant_co2_flux_mgm2d1, plant_cco2_flux_mgm2d1)) %>%
    summarise(across(where(is.numeric), mean))

  # get daily, summarised data
  cleaned_lake_data <- cleaned_lake_data %>%
    filter(chamber_int == "Measurement") %>% # only keep measurements
    group_by(sample_date, zone, subsample, site, veg_class) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
    left_join(chamber_data,
              by = c("sample_date", "site")) %>%
    filter(site != "REF") %>%
    mutate(EffB = cch4_flux_mgm2d1/biomass_g) %>%
    select(-c(gasVolume, waterVolume, barometricPressure, headspaceTemp, concentrationCO2Gas,
              concentrationCO2Air, concentrationCH4Gas, concentrationCH4Air))
  
  # print(head(cleaned_lake_data$dissolvedCH4))

  bub <- bubbledata %>%
    mutate(zone = as.factor(zone),
           ebull_flux_ch4 = case_when(
             bub_flux_lm2d1 == 0 ~ 0,
             TRUE ~ ebull_flux_ch4),
           log_ebull_flux_ch4 = log(ebull_flux_ch4 + 1),
           log_ebull_flux_co2 = log(ebull_flux_ch4 + 1)) %>%
    filter(!is.na(sample_date)) %>%
    select(sample_date, site, zone, subsample, veg_class, days_deploy_t, sample_taken,
           tot_vol_ml, depth_bt, 
           perc_ch4, perc_co2, co2_ppm, ch4_ppb, ch4_ppm, ch4_ppm_postdil, co2_ppm_postdil,
           bub_flux_lm2d1, ebull_flux_ch4, ebull_flux_co2, log_ebull_flux_ch4, log_ebull_flux_co2)
  
  # fill in missing bub values
  bub <- bub %>%
    group_by(zone, veg_class) %>%
    arrange(sample_date) %>%
    mutate(
      avg_ch4_ppm = mean(ch4_ppm_postdil, na.rm = T), 
      avg3day_ch4_ppm = map_dbl(sample_date, function(d) {
        mean(
          ch4_ppm_postdil[abs(sample_date - d) <= 1],
          na.rm = TRUE
        )
      }),
      avg7day_ch4_ppm = map_dbl(sample_date, function(d) {
        mean(
          ch4_ppm_postdil[abs(sample_date - d) <= 3],
          na.rm = TRUE
        )
      }),
      ch4_ppm_postdil = case_when(
        !is.na(ch4_ppm_postdil) ~ ch4_ppm_postdil,
        !is.na(avg3day_ch4_ppm) ~ avg3day_ch4_ppm,
        TRUE ~ avg_ch4_ppm
      ),
      ebull_flux_ch4 = case_when(
        is.na(ebull_flux_ch4) ~ ((ch4_ppm_postdil * bub_flux_lm2d1 * 16.04)/22.4/1000),
        TRUE ~ ebull_flux_ch4
      ),
      perc_ch4 = case_when(
        is.na(perc_ch4) ~ (ch4_ppm_postdil / 10000),
        TRUE ~ perc_ch4
      )
    ) %>%
    ungroup()
  
  write.csv(bub, file = "~/Desktop/masters/data/msc-aquatic-ch4/cleaned-data/wheeler_bubble_2025.csv")
  
  cleaned_lake_data <- full_join(cleaned_lake_data, bub,
                                  by = c("sample_date", "site",
                                         "zone", "subsample", "veg_class")) %>%
    mutate(veg_class = factor(veg_class, levels = c("Open", "Vegetated"))) %>%
    filter(cch4_flux_mgm2d1 > 0 & cch4_flux_mgm2d1 <= 300)

  # write a clean csv for distribution
  write.csv(cleaned_lake_data, file = "~/Desktop/masters/data/msc-aquatic-ch4/cleaned-data/daily_wheeler_envflux_2025.csv")
  
  # create a list of the processed (and raw) datasets that may be useful during visualization
  datasets <- list(raw_data = lake_data,
                   clean_data = cleaned_lake_data,
                   dissolved_data = dissolved_data,
                   chamber_data = chamber_data,
                   bubble_data = bub)

  return(datasets)
}
