# MSc Project: A little sedge goes a long way: emergent vegetation as a key driver of littoral methane flux

# Function for calculating GHG slopes from interval chamber data.

# 08/15/2025 (Last Edit: 08/17/2026)
# Author: Kelsey McGuire
# kmcgu@mail.ubc.ca; kmcguire.9@outlook.com

# LOAD RELEVANT PACKAGES ----
library(dplyr)
library(lubridate)
library(broom)

# FUNCTION DESCRIPTION ----
ChamberSlopeCalc <- function(chamber_data) {
  # ------------------------------------------------------------------------ #
  # Arguments:
  # - merged_dataframe [df]: dataframe from the CleanData function.
  #
  # Purpose:
  # - Calculate the slopes of provided chamber flux data. Set up for a three
  #   time point floating chamber measurement.
  #
  # Output:
  # - dataframe with chamber slopes and accompanying environmental data
  # ------------------------------------------------------------------------ #
  chamber_slopes <- data.frame(site = character(),
                               sample_date = as.Date(character()),
                               veg_class = character(),
                               ph = numeric(),
                               wat_temp_C = numeric(),
                               do = numeric(),
                               ch4_slope = numeric(),
                               co2_slope = numeric(),
                               stringsAsFactors = FALSE) # create an empty data frame for computed values to go into

  sites <- unique(chamber_data$site) # get a list of the sites
  for (i in 1:length(sites)) { # create a loop to go over each site
    # subset data only from that specific site
    site_data <- chamber_data %>%
      filter(site == sites[i] & chamber_int == "Measurement") %>%
      dplyr::select(sample_date, sample_id, chamber_time, ch4_ppm, co2_ppm) %>% # keep only the relevant variables for computing flux
      mutate(sample_time = case_when(chamber_time == "t0" ~ 0,
                                     chamber_time == "t1" ~ 15,
                                     chamber_time == "t2" ~ 30)) # give times to the intervals

    # print(site_data)
    dates <- unique(site_data$sample_date) # get a list of the dates that sites were sampled on
    
    for (date in 1:length(dates)) {
      daily_site_data <- site_data %>% # pull data from the specific sampling day
        filter(sample_date == dates[date],
               chamber_time %in% c("t0", "t1", "t2")) %>% # make sure the chamber time corresponds to one of the time intervals
        distinct(chamber_time, .keep_all = TRUE)  %>% # only keep 3 - no repeats!
        arrange(chamber_time)
      
      # print(daily_site_data)
      
      if (length(daily_site_data$chamber_time) >= 2) {
        ch4_lm <- lm(ch4_ppm ~ sample_time, data = daily_site_data) # find the slope of the three time points to get an idea of conc over time
        co2_lm <- lm(co2_ppm ~ sample_time, data = daily_site_data) # find the slope of the three time points to get an idea of conc over time
        
        ch4_slope <- coef(ch4_lm)[["sample_time"]] / 60 # pull the slope from the linear regression above - convert to ppm over seconds
        co2_slope <- coef(co2_lm)[["sample_time"]] / 60 # pull the slope from the linear regression above - convert to seconds 
      }

      chamber_slopes <- rbind(chamber_slopes, data.frame( # create a data frame for just the slopes of from the chambers
        site = sites[i],
        sample_date = dates[date],
        ch4_slope = ch4_slope,
        co2_slope = co2_slope
      ))
      
    }
  }

  # choose to write a slope only csv if you'd like
  # write_csv(chamber_slopes, "~/Desktop/masters/data/msc-aquatic-ch4/flux-calculations/chamber_calculations/chamber_slopes.csv") # write a csv for slopes to manually see if there's anything wonky

  # set up constants and other values needed for flux calculations
  IGL <- 0.08206
  chamber_data$RT_light <- (chamber_data$temp_air_daily+273.5)*IGL
  chamber_data <- dplyr::left_join(chamber_data, chamber_slopes, by = c('site', 'sample_date'))

  chamber_data$co2_flux_umolm2s1 <- NA # create empty columns for fluxes
  chamber_data$ch4_flux_umolm2s1 <- NA
  chamber_data$co2_flux_mgm2d1 <- NA
  chamber_data$ch4_flux_mgm2d1 <- NA
  chamber_data$cco2_flux_mgm2d1 <- NA
  chamber_data$cch4_flux_mgm2d1 <- NA

  # create for loops over sites
  for (i in 1:length(chamber_data$site)) {
    if (grepl('V', chamber_data$site[i], fixed = T) == TRUE) { # if it's a vegetated site ...
      # ... use the measurements from the vegetated chamber (i.e., tall boy)
      veg_chamber_hgt <- 28.8/100 #m
      veg_chamber_rad <- (mean(25.7, 27, 28.6)/2)/100 #m

      volume <- (pi * veg_chamber_rad^2 * veg_chamber_hgt * 1000) - (chamber_data$abovewater_vol_cm3[i]/1000000) # account for the amount of volume taken up by the plants within the chamber
        area <- pi * veg_chamber_rad^2
    }

    if (grepl('O', chamber_data$site[i], fixed = T) == TRUE) { # if it's an open site ...
      # use either the measurements from the first short chamber (thing 1) ...
      if (chamber_data$sample_date[i] < as.Date("08/14/2025", format = "%m/%d/%Y")) {
        open_chamber_hgt <- 18.7/100 #m
        open_chamber_rad <- (mean(26.3, 28.4, 27.1)/2)/100 #m

        volume <- pi * open_chamber_rad^2 * open_chamber_hgt * 1000
        area <- pi * open_chamber_rad^2
      } else {
        # or the replacement short chamber (thing 2) ...
        open_chamber_hgt <- 20.4/100 #m
        open_chamber_rad <- (mean(27.8, 26, 26)/2)/100 #m

        volume <- pi * open_chamber_rad^2 * open_chamber_hgt * 1000
        area <- pi * open_chamber_rad^2
      }
    }
    if (is.na(chamber_data$press_kPa[i])) {
      atm_press <- 0.96
    } else {
      atm_press <- chamber_data$press_kPa[i]/101.3
    }

    # calculate fluxes using information from above in ...
    # umol ...
    chamber_data$co2_flux_umolm2s1[i] <- ((chamber_data$co2_slope[i]*atm_press)/
                                            chamber_data$RT_light[i]) * (volume/area)

    chamber_data$ch4_flux_umolm2s1[i] <- ((chamber_data$ch4_slope[i]*atm_press)/
                                            chamber_data$RT_light[i]) * (volume/area)

    # mg of co2/ch4 or ...
    chamber_data$co2_flux_mgm2d1[i] <- ((chamber_data$co2_flux_umolm2s1[i]*86400)*44.01)/1000

    chamber_data$ch4_flux_mgm2d1[i] <- ((chamber_data$ch4_flux_umolm2s1[i]*86400)*16.04)/1000

    # mg of c in the co2/ch4
    chamber_data$cco2_flux_mgm2d1[i] <- ((chamber_data$co2_flux_umolm2s1[i]*86400)*12.011)/1000

    chamber_data$cch4_flux_mgm2d1[i] <- ((chamber_data$ch4_flux_umolm2s1[i]*86400)*12.011)/1000
  }
  
  chamber_data <- chamber_data %>%
    mutate(
      zone = str_extract(site, "^\\d+"),
      site_group = str_extract(site, "^\\d+[A-Z]"),
      type = ifelse(str_detect(site, "V$"), "V", "O")
    ) 
    
  # find outlier boundaries
  bounds <- chamber_data %>%
    group_by(zone, type) %>%
    summarise(
      q1 = quantile(cch4_flux_mgm2d1, 0.25, na.rm = TRUE),
      q3 = quantile(cch4_flux_mgm2d1, 0.75, na.rm = TRUE),
      iqr = IQR(cch4_flux_mgm2d1, na.rm = TRUE),
      lower_bound = q1 - 3 * iqr,
      upper_bound = q3 + 3 * iqr,
      .groups = "drop"
    )
  
  print(bounds)
  
  # Filter out outliers (keep only normal data)
  chamber_data_clean <- chamber_data %>%
    left_join(bounds, by = c("zone", "type")) %>%
    filter(
      cch4_flux_mgm2d1 >= lower_bound,
      cch4_flux_mgm2d1 <= upper_bound
    ) %>%
    select(-q1, -q3, -iqr, -lower_bound, -upper_bound)
  
  # print(head(chamber_data_clean))
  
  plant_med <- chamber_data_clean %>% # get plant-mediated specific measurements
    
    group_by(sample_date, site_group, zone, type) %>% # find the average flux over the zone
    summarise(
      ch4_flux_mgm2d1  = mean(ch4_flux_mgm2d1, na.rm = TRUE),
      co2_flux_mgm2d1  = mean(co2_flux_mgm2d1, na.rm = TRUE),
      cch4_flux_mgm2d1 = mean(cch4_flux_mgm2d1, na.rm = TRUE),
      cco2_flux_mgm2d1 = mean(cco2_flux_mgm2d1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    
    pivot_wider( 
      names_from = type, # pivot values longer so vegetated and open measurements are side by side based on site
      values_from = c(
        ch4_flux_mgm2d1,
        co2_flux_mgm2d1,
        cch4_flux_mgm2d1,
        cco2_flux_mgm2d1
      )
    ) %>%
    
    group_by(zone) %>% # subtract open values from vegetated to get the final plant-med value
    mutate(
      plant_ch4_flux_mgm2d1  = ch4_flux_mgm2d1_V - ch4_flux_mgm2d1_O,
      plant_co2_flux_mgm2d1  = co2_flux_mgm2d1_V - co2_flux_mgm2d1_O,
      plant_cch4_flux_mgm2d1 = cch4_flux_mgm2d1_V - cch4_flux_mgm2d1_O,
      plant_cco2_flux_mgm2d1 = cco2_flux_mgm2d1_V - cco2_flux_mgm2d1_O,
      ) %>%
    mutate(
      plant_cch4_flux_mgm2d1 = case_when( # make sure the erreonous data is accounted for ...
        plant_cch4_flux_mgm2d1 < 0 ~ cch4_flux_mgm2d1_V,# if any plant-mediated are below zero ... 
        is.na(plant_cch4_flux_mgm2d1) ~ cch4_flux_mgm2d1_V, # or missing, replace the value with the vegetated zone average
      TRUE ~ plant_cch4_flux_mgm2d1 # if that's not the case keep computed value
      )) %>%
    ungroup() %>%
    mutate(site = paste0(site_group, "V")) %>%
    select(
      sample_date, zone, site,
      plant_ch4_flux_mgm2d1,
      plant_cch4_flux_mgm2d1,
      plant_co2_flux_mgm2d1,
      plant_cco2_flux_mgm2d1,
    ) # keep only relevant columns

  chamber_data <- chamber_data %>% # bind the plant-mediated data with other chamber data
    left_join(plant_med %>%
                select(sample_date, site, zone,
                       plant_ch4_flux_mgm2d1, plant_co2_flux_mgm2d1,
                       plant_cch4_flux_mgm2d1, plant_cco2_flux_mgm2d1), 
              by = c("sample_date", "site", "zone")) 

  return(chamber_data)
}


