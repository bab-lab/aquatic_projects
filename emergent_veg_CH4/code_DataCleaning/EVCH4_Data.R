# MSc Project: A little sedge goes a long way: emergent vegetation as a key driver of littoral methane flux

# Script containing any and all subsets of data and post-processing after intitial data cleaning.
# Incubation and Isotope Data cleaning is contained within this script.

# 02/07/2026 (Last Update: 08/17/2026)
# Author: Kelsey McGuire
# kmcgu@mail.ubc.ca; kmcguire.9@outlook.com

# INSTALL LOCAL FUNCTIONS/LIBRARIES
source(here::here('code_DataCleaning', 'CleanData.R'))
source(here::here('code_DataCleaning', 'SlopeCalc.R'))
source(here::here('code_PaperFiguresStats', 'ModelRunning.R'))

# INSTALL/LOAD GLOBAL LIBRARIES
# install.packages("datetime")
# remotes::install_github("wilkelab/cowplot")
# install.packages("colorspace", repos = "http://R-Forge.R-project.org")
# remotes::install_github("clauswilke/colorblindr")

## relevant packages
library(ggplot2)
library(tidyverse)
library(nlme)
library(ggsignif)
library(ggpubr)
# install.packages("multcomp")
library(multcomp)
library(lme4)
library(lmerTest)
conflicted::conflicts_prefer(lmerTest::lmer)
library(broom)
library(janitor)

# LOAD IN DATA ----
# replace csv paths with your own ...
datasets <- CleanData(licor_csvpath = '~/Desktop/masters/data/msc-aquatic-ch4/licor-data-methane2025/WHE-methane-LICOR-jan19-clean.csv',
                      env_csvpath = '~/Desktop/masters/data/msc-aquatic-ch4/flux-calculations/chamber_calculations/environmental-variables-aug18.csv',
                      bcws_folderpath = '~/Desktop/masters/data/msc-aquatic-ch4/tachila-weather', 445,
                      hobo_folderpath = '~/Desktop/masters/data/msc-aquatic-ch4/Wheeler_HOBO_data/hobo_csvs/',
                      hobo_rangemetafile = '~/Desktop/masters/data/msc-aquatic-ch4/Wheeler_HOBO_data/hobo_whe_2025.csv',
                      field_veg_data = '~/Desktop/masters/data/msc-aquatic-ch4/veg-data/whe-veg-data.csv',
                      surface_area_csv = '~/Desktop/masters/data/msc-aquatic-ch4/veg-data/veg-photos/top-view/imageJ/area/surf-area.csv',
                      rgb_folderpath = '~/Desktop/masters/data/msc-aquatic-ch4/veg-data/veg-photos/top-view/imageJ/RGB',
                      bubcsv_file = '~/Desktop/masters/data/msc-aquatic-ch4/flux-calculations/bubble_calculations/whe-bubble-data-feb02.csv')

# get cleaned (chamber) ----
clean_data <- datasets$clean %>%
  filter(!is.na(cch4_flux_mgm2d1)) %>% # filter any NA values for flux since this will be our key variable
  mutate(log_cch4_flux_mgm2d1 = log(cch4_flux_mgm2d1 + 1),
         log_plant_cch4_flux_mgm2d1 = log(plant_cch4_flux_mgm2d1 + 1))

# remove outliers from whole dataset ----
# find outlier boundaries
bounds <- clean_data %>%
  group_by(zone, veg_class) %>%
  summarise(
    q1 = quantile(cch4_flux_mgm2d1, 0.25, na.rm = TRUE),
    q3 = quantile(cch4_flux_mgm2d1, 0.75, na.rm = TRUE),
    iqr = IQR(cch4_flux_mgm2d1, na.rm = TRUE),
    lower_bound = q1 - 3 * iqr,
    upper_bound = q3 + 3 * iqr,
    .groups = "drop"
  )

no_clean_data <- clean_data %>%
  left_join(bounds, by = c("zone", "veg_class")) %>%
  group_by(zone, veg_class) %>%
  filter(cch4_flux_mgm2d1 <= upper_bound & # IQR 75
           cch4_flux_mgm2d1 >= lower_bound) %>% # IQR 25
  ungroup() %>%
  select(-c(q1, q3, iqr, lower_bound, upper_bound))
  
## data cleaning specifically for modelling ----
### create vegetation only model data ----
veg_ch4 <- clean_data %>%
  filter(!is.na(veg_class) & !is.na(cch4_flux_mgm2d1)) %>%
  filter(veg_class == "Vegetated")

no_veg_ch4 <- no_clean_data %>%
  filter(veg_class == "Vegetated")

### create open only model data ----
open_ch4 <- clean_data %>%
  filter(!is.na(veg_class) & !is.na(cch4_flux_mgm2d1)) %>%
  filter(veg_class == "Open")

no_open_ch4 <- no_clean_data %>%
  filter(veg_class == "Open") 

### make perfectly paired data ----
no_perfect_pairs <- no_clean_data %>%
  filter(!is.na(cch4_flux_mgm2d1)) %>%
  group_by(sample_date, zone, subsample) %>%
  filter(setequal(veg_class, c("Open", "Vegetated"))) %>%
  ungroup()

perfect_pairs <- clean_data %>%
  filter(!is.na(cch4_flux_mgm2d1)) %>%
  group_by(sample_date, zone, subsample) %>%
  filter(setequal(veg_class, c("Open", "Vegetated"))) %>%
  ungroup()

perf_veg_ch4 <- perfect_pairs %>%
  filter(veg_class == "Vegetated")

perf_open_ch4 <- perfect_pairs %>%
  filter(veg_class == "Open")

no_perf_veg_ch4 <- no_perfect_pairs %>%
  filter(veg_class == "Vegetated")

no_perf_open_ch4 <- no_perfect_pairs %>%
  filter(veg_class == "Open")

### find zone averages ----
zone_avg <- clean_data %>%
  filter(!is.na(cch4_flux_mgm2d1)) %>%
  group_by(sample_date, zone, veg_class) %>%
  summarise(across(where(is.numeric), list(n = ~sum(!is.na(.x)),
                                           mean = ~mean(.x, na.rm = TRUE),
                                           se = ~ sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))))) 
no_zone_avg <- no_clean_data %>%
  filter(!is.na(cch4_flux_mgm2d1)) %>%
  group_by(sample_date, zone, veg_class) %>%
  summarise(across(where(is.numeric), list(n = ~sum(!is.na(.x)),
                                           mean = ~mean(.x, na.rm = TRUE),
                                           se = ~ sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))))) 

perf_zone_avg <- perfect_pairs %>%
  filter(!is.na(cch4_flux_mgm2d1)) %>%
  group_by(sample_date, zone, veg_class) %>%
  summarise(across(where(is.numeric), list(n = ~sum(!is.na(.x)),
                                           mean = ~mean(.x, na.rm = TRUE),
                                           se = ~ sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))))) 

veg_za <- zone_avg %>%
  filter(!is.na(veg_class) & !is.na(cch4_flux_mgm2d1_mean)) %>%
  filter(veg_class == "Vegetated")

no_veg_za <- no_zone_avg %>%
  filter(!is.na(veg_class) & !is.na(cch4_flux_mgm2d1_mean)) %>%
  filter(veg_class == "Vegetated")

perf_veg_za <- perf_zone_avg %>%
  filter(!is.na(veg_class) & !is.na(cch4_flux_mgm2d1_mean)) %>%
  filter(veg_class == "Vegetated")

open_za <- zone_avg %>%
  filter(!is.na(veg_class) & !is.na(cch4_flux_mgm2d1_mean)) %>%
  filter(veg_class == "Open")

no_open_za <- no_zone_avg %>%
  filter(!is.na(veg_class) & !is.na(cch4_flux_mgm2d1_mean)) %>%
  filter(veg_class == "Open")

perf_open_za <- perf_zone_avg %>%
  filter(!is.na(veg_class) & !is.na(cch4_flux_mgm2d1_mean)) %>%
  filter(veg_class == "Open")

depth_modelling <- clean_data %>%
  filter(
    !(
      (site == "2AO" & sample_date == "2025-07-18") | # remove the super low 5cm depth
        (site == "2AV" & sample_date %in% c("2025-07-18")) # remove it's pair (07/18), and the extremely high open measurement (08/17)
    )
  )

# bub data ----
bubdata <- datasets$bubble_data %>% 
  filter(!is.na(bub_flux_lm2d1)) %>%
  mutate(log_ebull_flux_ch4 = log(ebull_flux_ch4 + 1),
         log_ebull_flux_co2 = log(ebull_flux_ch4 + 1))

perf_bubdata <- bubdata %>% 
  filter(!is.na(bub_flux_lm2d1)) %>%
  group_by(sample_date, zone, subsample) %>%
  filter(setequal(veg_class, c("Open", "Vegetated"))) %>%
  ungroup()

perf_bubdata %>%
  count(sample_date, zone, subsample) %>%
  count(n)

perf_veg_bub <- perf_bubdata %>%
  filter(veg_class == "Vegetated")

perf_open_bub <- perf_bubdata %>%
  filter(veg_class == "Open")

bubdata_za <- bubdata %>%
  group_by(sample_date, zone, veg_class) %>%
  summarise(across(where(is.numeric), list(n = ~sum(!is.na(.x)),
                                           mean = ~mean(.x, na.rm = TRUE),
                                           se = ~ sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))))))

perf_bubdata_za <- bubdata_za %>% 
  filter(!is.na(bub_flux_lm2d1_mean)) %>%
  group_by(sample_date, zone) %>%
  filter(setequal(veg_class, c("Open", "Vegetated"))) %>%
  ungroup()

perf_veg_bub_za <- perf_bubdata_za %>%
  filter(veg_class == "Vegetated")

perf_open_bub_za <- perf_bubdata_za %>%
  filter(veg_class == "Open")

# ebull data ----
ebulldata <- bubdata %>% 
  filter(!is.na(ebull_flux_ch4))

perf_ebulldata <- ebulldata %>% 
  group_by(sample_date, zone, subsample) %>%
  filter(setequal(veg_class, c("Open", "Vegetated"))) %>%
  ungroup()

veg_ebull <- perf_ebulldata %>%
  filter(veg_class == "Vegetated")

open_ebull <- perf_ebulldata %>%
  filter(veg_class == "Open")

ebulldata_za <- ebulldata %>%
  group_by(sample_date, zone, veg_class) %>%
  summarise(across(where(is.numeric), list(n = ~sum(!is.na(.x)),
                                           mean = ~mean(.x, na.rm = TRUE),
                                           se = ~ sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))))))

perf_ebulldata_za <- ebulldata_za %>% 
  group_by(sample_date, zone) %>%
  filter(setequal(veg_class, c("Open", "Vegetated"))) %>%
  ungroup()

veg_ebull_za <- perf_ebulldata_za %>%
  filter(veg_class == "Vegetated")

open_ebull_za <- perf_ebulldata_za %>%
  filter(veg_class == "Open")

## model bubble ----
model_bub <- datasets$clean_data %>%
  filter(!is.na(ebull_flux_ch4))

model_bub_za <- datasets$clean_data %>%
  filter(!is.na(ebull_flux_ch4)) %>%
  group_by(sample_date, zone, veg_class) %>%
  summarise(across(where(is.numeric), list(n = ~sum(!is.na(.x)),
                                           mean = ~mean(.x, na.rm = TRUE),
                                           se = ~ sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))))))

model_bub_open <- model_bub %>%
  filter(veg_class == "Open")

model_bub_veg <- model_bub %>%
  filter(veg_class == "Vegetated")

model_bub_za_open <- model_bub_za %>%
  filter(veg_class == "Open")

model_bub_za_veg <- model_bub_za %>%
  filter(veg_class == "Vegetated")

# diss data ----
diss_clean_data <- datasets$clean_data %>%
  filter(!is.na(dissolvedCH4)) %>%
  group_by(sample_date, zone, subsample) %>%
  filter(setequal(veg_class, c("Open", "Vegetated"))) %>%
  ungroup()

## data cleaning specifically for modelling ----
### create vegetation only model data ----
veg_diss_ch4 <- diss_clean_data %>%
  filter(veg_class == "Vegetated")

### create open only model data ----
open_diss_ch4 <- diss_clean_data %>%
  filter(veg_class == "Open")

### find zone averages ----
diss_zone_avg <- diss_clean_data %>%
  filter(!is.na(dissolvedCH4)) %>%
  group_by(sample_date, zone, veg_class) %>%
  summarise(across(where(is.numeric), list(n = ~sum(!is.na(.x)),
                                           mean = ~mean(.x, na.rm = TRUE),
                                           se = ~ sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))))))      

diss_veg_za <- diss_zone_avg %>%
  filter(!is.na(veg_class) & !is.na(dissolvedCH4_mean)) %>%
  filter(veg_class == "Vegetated")

diss_open_za <- diss_zone_avg %>%
  filter(!is.na(veg_class) & !is.na(dissolvedCH4_mean)) %>%
  filter(veg_class == "Open")

## relative contributions ----
rel_summary <- no_clean_data %>%
  select(sample_date, zone, site, veg_class,
         diffusive_cch4 = cch4_flux_mgm2d1, diffusive_cco2 = cco2_flux_mgm2d1,
         plantmed_cch4 = plant_cch4_flux_mgm2d1, plantmed_cco2 = plant_cco2_flux_mgm2d1,
         ebullitive_cch4 = ebull_flux_ch4, ebullitive_cco2 = ebull_flux_co2) %>%
  ungroup() %>%
  group_by(zone, veg_class) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = T)))

## get relative summaries ----
ch4_rel_summary <- rel_summary %>%
  pivot_longer(cols = c(diffusive_cch4, plantmed_cch4, ebullitive_cch4),
               names_to = "process",
               values_to = "flux") %>%
  group_by(zone) %>%
  mutate(
    open_diffusive = flux[veg_class == "Open" & process == "diffusive_cch4"],
    
    flux_adj = case_when(
      # OPEN
      veg_class == "Open" & process == "plantmed_cch4" ~ 0,
      veg_class == "Open" & process == "diffusive_cch4" ~ flux,
      veg_class == "Open" & process == "ebullitive_cch4" ~ flux,
      
      # VEGETATED
      veg_class == "Vegetated" & process == "diffusive_cch4" ~ open_diffusive,
      veg_class == "Vegetated" & process != "diffusive_cch4" ~ flux
    )
  ) %>%
  ungroup() %>%
  select(-c(diffusive_cco2, plantmed_cco2, ebullitive_cco2, open_diffusive)) %>%
  mutate(process = factor(process,
                          levels = c("diffusive_cch4", "plantmed_cch4", "ebullitive_cch4"),
                          labels = c("Diffusive", "Plant-Mediated", "Ebullitive")))

relcont_ch4_sum <- ch4_rel_summary %>%
  group_by(process, veg_class) %>%
  summarise(across(where(is.numeric),
                   list(
                     mean = ~mean(.x, na.rm = TRUE),
                     sd = ~sd(.x, na.rm = TRUE),
                     se = ~sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))
                   ),
                   .names = "{.col}_{.fn}"
  ))

total_ch4_cont <- relcont_ch4_sum %>%
  group_by(veg_class) %>%
  summarise(across(where(is.numeric),
                   list(
                     sum = ~sum(.x, na.rm = TRUE),
                     sd = ~sd(.x, na.rm = TRUE),
                     se   = ~sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))
                   ),
                   .names = "{.col}_{.fn}"
  ))

# INCUBATION ----
IncLICOR <- read.csv("~/Desktop/masters/data/msc-aquatic-ch4/incubations/data/licor-measurements.csv",
                     skip = 6)
IncProduction <- read.csv("~/Desktop/masters/data/msc-aquatic-ch4/incubations/data/production_calcs_apr3.csv")

IncLOI <- read.csv("~/Desktop/masters/data/msc-aquatic-ch4/incubations/data/inc-loi.csv")[1:12,]

IncCN <- read.csv("~/Desktop/masters/data/msc-aquatic-ch4/incubations/data/inc-cn.csv")

## clean up LICOR columns ----
CleanIncLICOR <- IncLICOR %>%
  clean_names() %>% # make column names better
  filter(!grepl("-2", sample_name)) %>% # remove replicate samples
  rename(licor_time_base = licor_time,
         licor_time_measure = licor_time_1) %>% # rename licor times for better practice
  filter(sample_name != "ref-air") %>% # remove ref air samples
  separate_wider_delim(col = sample_name, 
                       delim = "-",
                       names = c("treatment_short", "temperature", "replicate"),
                       cols_remove = FALSE) %>%
  filter(!is.na(measurement_ch4_ppb) & treatment_day != "t14") %>%
  mutate(sample_date = as.Date(sample_date, format = "%m/%d/%Y"),
         temperature = recode(temperature,
                              "t20" = "Average Field Temperature (20°C)",
                              "t30" = "Potential Warming Temperature (30°C)"),
         treatment_day = fct_recode(treatment_day,
                                    "0"  = "t0", "1"  = "t1", "2"  = "t2", "3"  = "t3", 
                                    "7"  = "t7", "10" = "t10"),
         treatment_day = fct_relevel(treatment_day, "10", after = Inf),
         licor_time_base = hms(licor_time_base),
         licor_time_measure = hms(licor_time_measure),
         delt_ch4 = as.numeric(delt_ch4),
         delt_co2 = as.numeric(delt_co2),
         ch4_corrected_ppm = as.numeric(ch4_corrected_ppm),
         co2_corrected_ppm = as.numeric(co2_corrected_ppm),
         delt_ch4_post_veff = as.numeric(delt_ch4_post_veff),
         delt_co2_post_veff = as.numeric(delt_co2_post_veff),
         ch4_corrected_ppm_post_veff = as.numeric(ch4_corrected_ppm_post_veff),
         co2_corrected_ppm_post_veff = as.numeric(co2_corrected_ppm_post_veff))

## clean up PRODUCTION columns ----
CleanIncProduction <- IncProduction %>%
  clean_names() %>% # make column names cleaner
  filter(!is.na(ch4_ppb) & time_step != "t14") %>% 
  dplyr::select(-c(x, x_1, x_2)) %>%
  separate_wider_delim(col = sample_id, 
                       delim = "-",
                       names = c("treatment_short", "t_c", "replicate"),
                       cols_remove = FALSE) %>%
  mutate(treatment = fct_recode(treatment_short,
                                "Control"  = "C", "Water + Benthic"  = "WB", 
                                "Water + Sediment + Benthic"  = "WSB", "Water + Sediment"  = "WS"),
         treatment = fct_relevel(treatment,
                                 "Water + Sediment",
                                 "Water + Sediment + Benthic",
                                 "Water + Benthic"),
         time_step = fct_recode(time_step,
                                "1"  = "t1", "2"  = "t2", "3"  = "t3", "7"  = "t7", "10" = "t10"),
         time_step = fct_relevel(time_step, "10", after = Inf),
         time_between_steps_hrs = as.numeric(time_between_steps_hrs),
         t_c = fct_recode(t_c,
                          "Average Field Temperature (20°C)" = "t20",
                          "Potential Warming Temperature (30°C)" = "t30"),
         log_nmol_ch4_dry_gw_1_hr_1 = log(as.numeric(nmol_ch4_dry_gw_1_hr_1)),
         log_umol_ch4_dry_gw_1_hr_1 = log(as.numeric(umol_ch4_dry_gw_1_hr_1)),
         plog_nmol_ch4_dry_gw_1_hr_1 = log(as.numeric(nmol_ch4_dry_gw_1_hr_1) + 1),
         plog_umol_ch4_dry_gw_1_hr_1 = log(as.numeric(umol_ch4_dry_gw_1_hr_1) + 1),
         log_nmol_ch4_dry_gw_1_d_1 = log(as.numeric(nmol_ch4_dry_gw_1_d_1)),
         log_umol_ch4_dry_gw_1_d_1 = log(as.numeric(umol_ch4_dry_gw_1_d_1)),
         plog_nmol_ch4_dry_gw_1_d_1 = log(as.numeric(nmol_ch4_dry_gw_1_d_1) + 1),
         plog_umol_ch4_dry_gw_1_d_1 = log(as.numeric(umol_ch4_dry_gw_1_d_1) + 1),
         inc_split = case_when(time_step == "1" ~ "Early",
                               time_step == "2" ~ "Early",
                               time_step == "3" ~ "Early",
                               TRUE ~ "Late")) %>%
  select(treatment, treatment_short, t_c, t_k, replicate, sample_id, time_step, inc_split,
         ch4_ppb, nmol_ch4_dry_gw_1_hr_1, umol_ch4_dry_gw_1_hr_1, nmol_ch4_dry_gw_1_d_1, umol_ch4_dry_gw_1_d_1,
         log_nmol_ch4_dry_gw_1_hr_1, log_umol_ch4_dry_gw_1_hr_1, log_nmol_ch4_dry_gw_1_d_1, log_umol_ch4_dry_gw_1_d_1,
         plog_nmol_ch4_dry_gw_1_hr_1, plog_umol_ch4_dry_gw_1_hr_1, plog_nmol_ch4_dry_gw_1_d_1, plog_umol_ch4_dry_gw_1_d_1)

## remove high and low values ----
NoOutIncProd <- data.frame()

unique_temps <- unique(CleanIncProduction$t_c)
for (z in 1:length(unique_temps)) {
  temperature_data <- CleanIncProduction %>%
    filter(t_c == unique_temps[z])
  
  unique_time_steps <- unique(CleanIncProduction$time_step)
  for (i in 1:length(unique_time_steps)) {
    daily_temp_data <- temperature_data %>%
      filter(time_step == unique_time_steps[i])
    
    unique_treatments <- unique(daily_temp_data$treatment)
    for (x in 1:length(unique_treatments)) {
      treatment_daily <- daily_temp_data %>%
        filter(treatment == unique_treatments[x])
      
      negatives <- sum(treatment_daily$nmol_ch4_dry_gw_1_hr_1 < 0)
      treatment <- unique(treatment_daily$treatment_short)
      
      if (negatives > 1 & treatment %in% c("WS", "WSB", "WB")) {
        print(treatment_daily)
        treatment_daily <- treatment_daily %>%
          filter(nmol_ch4_dry_gw_1_hr_1 > 0)
      } else {
        treatment_daily <- treatment_daily %>%
          filter(!nmol_ch4_dry_gw_1_hr_1 %in% range(nmol_ch4_dry_gw_1_hr_1))
      }
      
      NoOutIncProd <- rbind(NoOutIncProd, treatment_daily)
    }
  }
}

write.csv(NoOutIncProd, file = "~/Desktop/masters/paper_submission/data/NOIncubationProduction.csv")


## Q10 ----
Q10IncProd <- NoOutIncProd %>%
  dplyr::select(treatment, time_step, inc_split, t_k, nmol_ch4_dry_gw_1_hr_1)

Q10IncProd_20deg <- Q10IncProd %>%
  filter(t_k == 293.15) %>%
  group_by(treatment, time_step) %>%
  summarise(
    mean_dryprod = mean(nmol_ch4_dry_gw_1_hr_1, na.rm = TRUE),
    se_dryprod = sd(nmol_ch4_dry_gw_1_hr_1, na.rm = TRUE) / sqrt(sum(!is.na(nmol_ch4_dry_gw_1_hr_1))),
  )

Q10IncProd_30deg <- Q10IncProd %>%
  filter(t_k == 303.15) %>%
  group_by(treatment, time_step) %>%
  summarise(
    mean_dryprod = mean(nmol_ch4_dry_gw_1_hr_1, na.rm = TRUE),
    se_dryprod = sd(nmol_ch4_dry_gw_1_hr_1, na.rm = TRUE) / sqrt(sum(!is.na(nmol_ch4_dry_gw_1_hr_1))),
  )

Q10IncProd <- Q10IncProd_20deg %>%
  full_join(Q10IncProd_30deg, 
            by = c("treatment", "time_step"),
            suffix = c("_20", "_30")) %>%
  mutate(
    q10_dry = (mean_dryprod_30/mean_dryprod_20)*exp(10/(303.15-293.15)),
  ) %>%
  filter(q10_dry > 0 & q10_dry < 30) %>%
  dplyr::select(treatment, time_step, q10_dry) 

### Q10 IN EARLY V. LATE
SplitQ10 <- Q10IncProd %>%
  mutate(inc_split = case_when(time_step == "1" ~ "Early",
                               time_step == "2" ~ "Early",
                               time_step == "3" ~ "Early",
                               TRUE ~ "Late")) %>%
  group_by(inc_split, treatment) %>%
  summarise(
    mean_q10_dry = mean(q10_dry, na.rm = TRUE),
    se_q10_dry = sd(q10_dry, na.rm = TRUE) / sqrt(sum(!is.na(q10_dry)))
  )

AvgQ10 <- Q10IncProd %>%
  group_by(treatment) %>%
  summarise(
    mean_q10_dry = mean(q10_dry, na.rm = TRUE),
    se_q10_dry = sd(q10_dry, na.rm = TRUE) / sqrt(sum(!is.na(q10_dry)))
  )

## calculate cumulative production, clean up LOI, and bind ----
TotalIncProduction <- data.frame()

unique_temps <- unique(CleanIncProduction$t_c)
for (x in 1:length(unique_temps)) {
  temperature_data <- CleanIncProduction %>%
    filter(t_c == unique_temps[x])
  
  unique_treatments <- unique(temperature_data$treatment)
  for (y in 1:length(unique_treatments)) {
    treatment_by_temp <- temperature_data %>%
      filter(treatment == unique_treatments[y])
    
    treatment_by_temp <- treatment_by_temp %>%
      arrange(replicate, time_step) %>%
      group_by(replicate) %>%
      mutate(
        day_diff = as.numeric(time_step) - dplyr::lag(as.numeric(time_step)),
        nmol_ch4 = if_else(
          day_diff == 1,
          nmol_ch4_dry_gw_1_d_1,
          (nmol_ch4_dry_gw_1_d_1 + dplyr::lag(nmol_ch4_dry_gw_1_d_1)) / day_diff
        ),
        cumulative_prod = sum(nmol_ch4, na.rm = T)
      )
    
    print(treatment_by_temp$nmol_ch4)
    
    TotalIncProduction <- rbind(TotalIncProduction, treatment_by_temp)
  }
}
## find avg. cumulative ----
SummaryCumulativeProduction <- TotalIncProduction %>%
  select(treatment, t_c, cumulative_prod) %>%
  group_by(treatment, t_c) %>%
  summarise(avg_cumulative_prod = mean(cumulative_prod, na.rm = TRUE),
            sd_cumulative_prod = sd(cumulative_prod, na.rm = TRUE))

## find cumulative production totals for Nutrient ----
### CLEAN LOI ----
CleanIncLOI <- IncLOI %>%
  clean_names() %>%
  separate_wider_delim(col = sample_name, 
                       delim = "-",
                       names = c("treatment_short", "t_c", "replicate"),
                       cols_remove = FALSE) %>%
  mutate(sed_om = soil_organic_matter_content,
         t_c = fct_recode(t_c,
                          "Average Field Temperature (20°C)" = "t20",
                          "Potential Warming Temperature (30°C)" = "t30")) %>%
  select(treatment_short, t_c, sample_name, replicate, sed_om) 

loi_reps <- unique(CleanIncLOI$sample_name)

### CLEAN CN ----
CleanIncCN <- IncCN %>%
  clean_names() %>%
  mutate(sample_id = sub("-[^-]*$", "", name)) %>%
  select(-c("name")) %>%
  separate_wider_delim(col = sample_id, 
                       delim = "-",
                       names = c("treatment_short", "t_c", "replicate"),
                       cols_remove = FALSE) %>%
  mutate(sed_c = x_total_c,
         sed_n = x_total_n,
         sed_cn = c_n,
         t_c = fct_recode(t_c,
                          "Average Field Temperature (20°C)" = "t20",
                          "Potential Warming Temperature (30°C)" = "t30")) %>%
  
  select(treatment_short, t_c, sample_id, replicate, sed_c, sed_n, sed_cn) 

cn_reps <- unique(CleanIncCN$sample_id)

TotalIncProduction <- TotalIncProduction %>%
  filter(sample_id %in% c(loi_reps, cn_reps)) %>% 
  select(treatment_short, t_c, cumulative_prod) %>%
  group_by(treatment_short, t_c) %>%
  summarise(avg_cumulative_prod = mean(cumulative_prod, na.rm = TRUE),
            se_cumulative_prod = sd(cumulative_prod, na.rm = TRUE) / sqrt(sum(!is.na(cumulative_prod))))

CleanIncLOISum <- CleanIncLOI %>%
  select(treatment_short, t_c, sed_om) %>%
  group_by(treatment_short, t_c) %>%
  summarise(avg_sed_om = mean(sed_om, na.rm = TRUE),
            se_sed_om = sd(sed_om, na.rm = TRUE) / sqrt(sum(!is.na(sed_om))))

CleanIncCNSum <- CleanIncCN %>%
  select(treatment_short, t_c, sed_c, sed_n, sed_cn) %>%
  group_by(treatment_short, t_c) %>%
  summarise(avg_sed_c = mean(sed_c, na.rm = TRUE),
            se_sed_c = sd(sed_c, na.rm = TRUE) / sqrt(sum(!is.na(sed_c))),
            
            avg_sed_n = mean(sed_n, na.rm = TRUE),
            se_sed_n = sd(sed_n, na.rm = TRUE) / sqrt(sum(!is.na(sed_n))),
            
            avg_sed_cn = mean(sed_cn, na.rm = TRUE),
            se_sed_cn = sd(sed_cn, na.rm = TRUE) / sqrt(sum(!is.na(sed_cn))))

CleanIncNutrient <- TotalIncProduction %>%
  left_join(CleanIncLOISum, by = c("treatment_short", "t_c")) %>%
  left_join(CleanIncCNSum, by = c("treatment_short", "t_c"))

## get daily averages ----
### LICOR ----
DailyIncLICOR <- CleanIncLICOR %>%
  group_by(treatment, temperature, treatment_day) %>%
  summarise(avg_conc = mean(ch4_corrected_ppm, na.rm = TRUE),
            se = sd(ch4_corrected_ppm, na.rm = TRUE) / sqrt(sum(!is.na(ch4_corrected_ppm))),
            avg_conc_pVeff = mean(ch4_corrected_ppm_post_veff, na.rm = TRUE),
            se_pVeff = sd(ch4_corrected_ppm_post_veff, na.rm = TRUE) / sqrt(sum(!is.na(ch4_corrected_ppm_post_veff))))

### PRODUCTION ----
DailyIncProduction <- CleanIncProduction %>%
  group_by(treatment, t_c, time_step) %>%
  summarise(avg_production_dry = mean(nmol_ch4_dry_gw_1_d_1, na.rm = TRUE),
            sd_dry = sd(nmol_ch4_dry_gw_1_d_1, na.rm = TRUE),
            se_dry = sd(nmol_ch4_dry_gw_1_d_1, na.rm = TRUE) / sqrt(sum(!is.na(nmol_ch4_dry_gw_1_hr_1))),
            log_avg_production_dry = mean(log_nmol_ch4_dry_gw_1_d_1, na.rm = TRUE),
            log_se_dry = sd(log_nmol_ch4_dry_gw_1_d_1, na.rm = TRUE) / sqrt(sum(!is.na(log_nmol_ch4_dry_gw_1_hr_1))),
            avg_production_dry_umol = mean(umol_ch4_dry_gw_1_d_1, na.rm = TRUE),
            se_dry_umol = sd(umol_ch4_dry_gw_1_d_1, na.rm = TRUE) / sqrt(sum(!is.na(umol_ch4_dry_gw_1_hr_1))),
            log_avg_production_dry_umol = mean(log_umol_ch4_dry_gw_1_d_1, na.rm = TRUE),
            log_se_dry_umol = sd(log_umol_ch4_dry_gw_1_d_1, na.rm = TRUE) / sqrt(sum(!is.na(log_umol_ch4_dry_gw_1_hr_1))),
            plog_avg_production_dry = mean(plog_nmol_ch4_dry_gw_1_d_1, na.rm = TRUE),
            plog_se_dry = sd(plog_nmol_ch4_dry_gw_1_d_1, na.rm = TRUE) / sqrt(sum(!is.na(plog_nmol_ch4_dry_gw_1_hr_1))),
            plog_avg_production_dry_umol = mean(plog_umol_ch4_dry_gw_1_d_1, na.rm = TRUE),
            plog_se_dry_umol = sd(plog_umol_ch4_dry_gw_1_d_1, na.rm = TRUE) / sqrt(sum(!is.na(plog_umol_ch4_dry_gw_1_hr_1))))

DailyNoOutIncProd <- NoOutIncProd %>%
  group_by(treatment, t_c, time_step) %>%
  summarise(avg_production_dry = mean(nmol_ch4_dry_gw_1_d_1, na.rm = TRUE),
            sd_dry = sd(nmol_ch4_dry_gw_1_d_1, na.rm = TRUE),
            se_dry = sd(nmol_ch4_dry_gw_1_d_1, na.rm = TRUE) / sqrt(sum(!is.na(nmol_ch4_dry_gw_1_hr_1))),
            log_avg_production_dry = mean(log_nmol_ch4_dry_gw_1_hr_1, na.rm = TRUE),
            log_sd_dry = sd(log_nmol_ch4_dry_gw_1_hr_1, na.rm = TRUE),
            # avg_production_dry_umol = mean(umol_ch4_dry_gw_1_hr_1, na.rm = TRUE),
            # se_dry_umol = sd(umol_ch4_dry_gw_1_hr_1, na.rm = TRUE) / sqrt(sum(!is.na(umol_ch4_dry_gw_1_hr_1))),
            # log_avg_production_dry_umol = mean(log_umol_ch4_dry_gw_1_hr_1, na.rm = TRUE),
            # log_se_dry_umol = sd(log_umol_ch4_dry_gw_1_hr_1, na.rm = TRUE) / sqrt(sum(!is.na(log_umol_ch4_dry_gw_1_hr_1))),
            plog_avg_production_dry = mean(plog_nmol_ch4_dry_gw_1_hr_1, na.rm = TRUE),
            plog_sd_dry = sd(plog_nmol_ch4_dry_gw_1_hr_1, na.rm = TRUE))
            # plog_avg_production_dry_umol = mean(plog_umol_ch4_dry_gw_1_hr_1, na.rm = TRUE),
            # plog_se_dry_umol = sd(plog_umol_ch4_dry_gw_1_hr_1, na.rm = TRUE) / sqrt(sum(!is.na(plog_umol_ch4_dry_gw_1_hr_1))))

AvgNoOutIncProd <- NoOutIncProd %>%
  group_by(treatment, t_c) %>% 
  summarise(avg_production_dry = mean(nmol_ch4_dry_gw_1_hr_1, na.rm = TRUE),
            se_dry = sd(nmol_ch4_dry_gw_1_hr_1, na.rm = TRUE) / sqrt(sum(!is.na(nmol_ch4_dry_gw_1_hr_1))),
            log_avg_production_dry = mean(log_nmol_ch4_dry_gw_1_hr_1, na.rm = TRUE),
            log_se_dry = sd(log_nmol_ch4_dry_gw_1_hr_1, na.rm = TRUE) / sqrt(sum(!is.na(log_nmol_ch4_dry_gw_1_hr_1))),
            avg_production_dry_umol = mean(umol_ch4_dry_gw_1_hr_1, na.rm = TRUE),
            se_dry_umol = sd(umol_ch4_dry_gw_1_hr_1, na.rm = TRUE) / sqrt(sum(!is.na(umol_ch4_dry_gw_1_hr_1))),
            log_avg_production_dry_umol = mean(log_umol_ch4_dry_gw_1_hr_1, na.rm = TRUE),
            log_se_dry_umol = sd(log_umol_ch4_dry_gw_1_hr_1, na.rm = TRUE) / sqrt(sum(!is.na(log_umol_ch4_dry_gw_1_hr_1))),
            plog_avg_production_dry = mean(plog_nmol_ch4_dry_gw_1_hr_1, na.rm = TRUE),
            plog_se_dry = sd(plog_nmol_ch4_dry_gw_1_hr_1, na.rm = TRUE) / sqrt(sum(!is.na(plog_nmol_ch4_dry_gw_1_hr_1))),
            plog_avg_production_dry_umol = mean(plog_umol_ch4_dry_gw_1_hr_1, na.rm = TRUE),
            plog_se_dry_umol = sd(plog_umol_ch4_dry_gw_1_hr_1, na.rm = TRUE) / sqrt(sum(!is.na(plog_umol_ch4_dry_gw_1_hr_1))))

SplitNoOutIncProd <- NoOutIncProd %>%
  group_by(treatment, t_c, inc_split) %>%
  summarise(avg_production_dry = mean(nmol_ch4_dry_gw_1_hr_1, na.rm = TRUE),
            se_dry = sd(nmol_ch4_dry_gw_1_hr_1, na.rm = TRUE) / sqrt(sum(!is.na(nmol_ch4_dry_gw_1_hr_1))),
            log_avg_production_dry = mean(log_nmol_ch4_dry_gw_1_hr_1, na.rm = TRUE),
            log_se_dry = sd(log_nmol_ch4_dry_gw_1_hr_1, na.rm = TRUE) / sqrt(sum(!is.na(log_nmol_ch4_dry_gw_1_hr_1))),
            avg_production_dry_umol = mean(umol_ch4_dry_gw_1_hr_1, na.rm = TRUE),
            se_dry_umol = sd(umol_ch4_dry_gw_1_hr_1, na.rm = TRUE) / sqrt(sum(!is.na(umol_ch4_dry_gw_1_hr_1))),
            log_avg_production_dry_umol = mean(log_umol_ch4_dry_gw_1_hr_1, na.rm = TRUE),
            log_se_dry_umol = sd(log_umol_ch4_dry_gw_1_hr_1, na.rm = TRUE) / sqrt(sum(!is.na(log_umol_ch4_dry_gw_1_hr_1))))

# ISOTOPES ----
isotopes <- read.csv("~/Desktop/masters/data/msc-aquatic-ch4/isotopes/results/McGuire_Lakes_d13C-CH4.csv")

## clean up data for plotting ----
clean_iso <- isotopes %>%
  clean_names() %>%
  dplyr::select(-c(pathway, x, x_1)) %>%
  separate_wider_delim(col = sample, 
                       delim = "-",
                       names = c("date", "site", "pathway")) %>%
  mutate(
    type = case_when(
      grepl("O", type) ~ "Open",
      TRUE ~ "Vegetated"),
    zone = substr(site, 1, 1),
    subzone = substr(site, 2, 2),
    season = case_when(
      date %in% c("JUL14", "JULY14", "JUL16") ~ "Early",
      date %in% c("JULY31", "AUG02") ~ "Middle",
      TRUE ~ "Late"),
    date = case_when(
      date == "JUL14" ~ as.Date("2025-07-14"),
      date == "JULY14" ~ as.Date("2025-07-14"),
      date == "JUL16" ~ as.Date("2025-07-16"),
      date == "JULY16" ~ as.Date("2025-07-14"),
      date == "JULY31" ~ as.Date("2025-07-31"),
      date == "AUG02" ~ as.Date("2025-08-02"),
      date == "AUG14" ~ as.Date("2025-08-14"),
      date == "AUG17" ~ as.Date("2025-08-17"),
      date == "AUG18" ~ as.Date("2025-08-18")),
    season = fct_relevel(season, "Early", "Middle", "Late"),
    pathway = fct_recode(pathway,
                         "Chamber-Based" = "FC", "Ebullitive" = "BUB"),
    pathway = fct_relevel(pathway, "Ebullitive", after = Inf),
  ) %>%
  filter(date != "2025-08-02" | zone != "2" | subzone != "B") # remove low veg value and the pair

d13FC_veg <- clean_iso %>%
  filter(type == "Vegetated" & pathway == "Chamber-Based")

d13FC_open <- clean_iso %>%
  filter(type == "Open" & pathway == "Chamber-Based")

d13BUB_veg <- clean_iso %>%
  filter(type == "Vegetated" & pathway == "Ebullitive")

d13BUB_open <- clean_iso %>%
  filter(type == "Open" & pathway == "Ebullitive") 

d13BUB_open <- d13BUB_open[!(d13BUB_open$date == "2025-07-31" &
                               d13BUB_open$site == "1AO"),]

iso_stats <- clean_iso %>%
  group_by(pathway, type) %>%
  summarise(min_d13c = min(d13c_ch4, na.rm = TRUE),
            max_d13c = max(d13c_ch4, na.rm = TRUE),
            avg_d13c = mean(d13c_ch4, na.rm = TRUE),
            se_d13c = sd(d13c_ch4, na.rm = TRUE) / sqrt(sum(!is.na(d13c_ch4))))

## filter to chamber isos ----
FCIso <- clean_iso %>%
  filter(pathway == "Chamber-Based")

FCIsoVeg <- FCIso %>%
  group_by(type) %>%
  summarise(avg_d13 = mean(d13c_ch4, na.rm = TRUE),
            se = sd(d13c_ch4, na.rm = TRUE) / sqrt(sum(!is.na(d13c_ch4))))

FCIsoDaily <- FCIso %>%
  group_by(date, type) %>%
  summarise(avg_d13 = mean(d13c_ch4, na.rm = TRUE),
            se = sd(d13c_ch4, na.rm = TRUE) / sqrt(sum(!is.na(d13c_ch4))))

FCIsoSeason <- FCIso %>%
  group_by(season, type) %>%
  summarise(avg_d13 = mean(d13c_ch4, na.rm = TRUE),
            se = sd(d13c_ch4, na.rm = TRUE) / sqrt(sum(!is.na(d13c_ch4))))

## filter to bubble isos ----
BUBIso <- clean_iso %>%
  filter(pathway == "Ebullitive")

BUBIsoVeg <- BUBIso %>%
  group_by(type) %>%
  summarise(avg_d13 = mean(d13c_ch4, na.rm = TRUE),
            se = sd(d13c_ch4, na.rm = TRUE) / sqrt(sum(!is.na(d13c_ch4))))

BUBIsoDaily <- BUBIso %>%
  group_by(date, type) %>%
  summarise(avg_d13 = mean(d13c_ch4, na.rm = TRUE),
            se = sd(d13c_ch4, na.rm = TRUE) / sqrt(sum(!is.na(d13c_ch4))))

BUBIsoSeason <- BUBIso %>%
  group_by(season, type) %>%
  summarise(avg_d13 = mean(d13c_ch4, na.rm = TRUE),
            se = sd(d13c_ch4, na.rm = TRUE) / sqrt(sum(!is.na(d13c_ch4))))

# MULTI-STUDY COMPS ----
## wheeler (this study)
whe_summary <- no_clean_data %>%
  ungroup() %>%
  select(veg_class, 
         # plantmed_cch4_mgm2d1,
         cch4_flux_mgm2d1) %>%
  mutate(dataset = "This Study (2025)") %>%
  filter(!is.na(cch4_flux_mgm2d1))

## kyzivat (2022) 
kyzivat <- read.csv("~/Desktop/masters/data/mock-data/kyzivat-mock-data/Field_flux_data.csv")

kyzivat_clean <- kyzivat %>%
  select(CH4Flux..mg.m2.d., Vegetation) %>%
  rename(ch4_flux_mgm2d1 = CH4Flux..mg.m2.d.,
         veg_class = Vegetation) %>%
  mutate(dataset = "Kyzivat (2022)",
         veg_class = factor(veg_class, 
                            levels = c("Open water", "Vegetation"),
                            labels = c("Open", "Vegetated")),
         cch4_flux_mgm2d1 = as.numeric(ch4_flux_mgm2d1) / 16.04 * 12.01) %>%
  select(-ch4_flux_mgm2d1) %>%
  filter(!is.na(veg_class) &  cch4_flux_mgm2d1 < 2000) # remove outliers

## desrosiers 2022
desrosiers <- read.csv("~/Desktop/masters/data/mock-data/desrosiers-data/surface_fluxes_desrosiers.csv")

desrosiers_clean <- desrosiers %>%
  select(Habitat, Diffusive.CH4.flux..mmol.m.2.d.1., Plant.mediated.CH4.flux..mmol.m.2.d.1.) %>%
  mutate(
    dataset = "Desrosiers (2022)",
    cch4_flux_mgm2d1 = (Diffusive.CH4.flux..mmol.m.2.d.1. + Plant.mediated.CH4.flux..mmol.m.2.d.1.) * 12.011,
    veg_class = case_when(
      Habitat == "Typha" ~ "Vegetated",
      Habitat == "Brasenia" ~ "Vegetated - FL",
      Habitat == "Unvegetated Littoral" ~ "Open"
    ),
    veg_class = factor(veg_class, levels = c("Open", "Vegetated", "Vegetated - FL"))) %>%
  select(-c(Habitat, Diffusive.CH4.flux..mmol.m.2.d.1., Plant.mediated.CH4.flux..mmol.m.2.d.1.)) %>%
  filter(!is.na(veg_class) & veg_class != "Vegetated - FL")

## BAWLD 2021
bawld <- read.csv("~/Desktop/masters/data/mock-data/kuhn-bawld-data/Kuhn_BAWLD_CH4_Aquatic.csv")

bawld_clean <- bawld %>%
  filter(LAKE.TYPE == "GP" & D.METHOD == "CH" & SIZE == "M") %>%
  select(CH4.D.FLUX) %>%
  mutate(dataset = "BAWLD-CH4 (2021)",
         veg_class = factor("Open"),
         cch4_flux_mgm2d1 = as.numeric(CH4.D.FLUX)/16.04 * 12.01) %>%
  select(-CH4.D.FLUX) %>%
  filter(!is.na(cch4_flux_mgm2d1))

bawld_summary <- bawld %>%
  filter(PATHWAY == "D" | PATHWAY == "E" &
           LAKE.TYPE == "GP") %>%
  select(LAKE.TYPE, SIZE, PATHWAY, BOTTOM, CH4.D.FLUX, CH4.E.FLUX) %>%
  mutate(cch4_flux_mgm2d1 = as.numeric(CH4.D.FLUX) * 12.011/16.04,
         ebull_cch4_mgm2d1 = as.numeric(CH4.E.FLUX) * 12.011/16.04) %>%
  select(-c(CH4.D.FLUX, CH4.E.FLUX)) %>%
  filter(BOTTOM != "U") %>%
  group_by(LAKE.TYPE, SIZE, PATHWAY, BOTTOM) %>%
  summarise(across(everything(), 
                   list(n = ~sum(!is.na(.x)),
                        min = ~min(.x, na.rm = T),
                        max = ~max(.x, na.rm = T),
                        mean = ~mean(.x, na.rm = TRUE),
                        se = ~ sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))))))

bawld_summary_cont <- bawld %>%
  filter((PATHWAY == "D" | PATHWAY == "E") &
           !is.na(CH4.D.FLUX) & !is.na(CH4.E.FLUX)) %>%
  select(LAKE.TYPE, SIZE, PATHWAY, BOTTOM, CH4.D.FLUX, CH4.E.FLUX) %>%
  mutate(cch4_flux_mgm2d1 = as.numeric(CH4.D.FLUX) * 12.011/16.04,
         ebull_cch4_mgm2d1 = as.numeric(CH4.E.FLUX) * 12.011/16.04) %>%
  select(-c(CH4.D.FLUX, CH4.E.FLUX, BOTTOM, PATHWAY)) %>%
  group_by(LAKE.TYPE, SIZE) %>%
  summarise(across(everything(), 
                   list(n = ~sum(!is.na(.x)),
                        min = ~min(.x, na.rm = T),
                        max = ~max(.x, na.rm = T),
                        mean = ~mean(.x, na.rm = TRUE),
                        se = ~ sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))))) %>%
  mutate(diff_contribution = case_when((!is.na(cch4_flux_mgm2d1_mean) & !is.na(ebull_cch4_mgm2d1_mean) ~ 
                                          cch4_flux_mgm2d1_mean / (cch4_flux_mgm2d1_mean + ebull_cch4_mgm2d1_mean))),
         ebull_contribution = case_when((!is.na(cch4_flux_mgm2d1_mean) & !is.na(ebull_cch4_mgm2d1_mean) ~ 
                                           ebull_cch4_mgm2d1_mean / (cch4_flux_mgm2d1_mean + ebull_cch4_mgm2d1_mean))))

mean(bawld_summary_cont$diff_contribution)

## combine the different datasets
multistudy <- rbind(whe_summary, desrosiers_clean, bawld_clean, kyzivat_clean)

study_summary <- multistudy %>%
  select(dataset, veg_class, cch4_flux_mgm2d1) %>%
  group_by(dataset, veg_class) %>%
  summarise(across(everything(), 
                   list(n = ~sum(!is.na(.x)),
                        min = ~min(.x, na.rm = TRUE),
                        max = ~max(.x, na.rm = TRUE),
                        mean = ~mean(.x, na.rm = TRUE),
                        median = ~median(.x, na.rm = T),
                        sd = ~ sd(.x, na.rm = TRUE),
                        se = ~ sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))))))

open_veg_rat <- study_summary %>%
  select(dataset, veg_class, 
         cch4_flux_mgm2d1_median, 
         cch4_flux_mgm2d1_sd) %>%
  pivot_wider(names_from = veg_class,
              values_from = c(cch4_flux_mgm2d1_median, 
                              cch4_flux_mgm2d1_sd)) %>%
  mutate(ov_ratio = cch4_flux_mgm2d1_median_Vegetated/cch4_flux_mgm2d1_median_Open) %>%
  select(dataset, ov_ratio) %>%
  mutate(dataset = factor(dataset,
                          levels = c("Desrosiers (2022)",
                                     "BAWLD-CH4 (2021)",
                                     "Kyzivat (2022)",
                                     "This Study (2025)")),
         ecosystem = case_when(
           dataset %in% c("BAWLD-CH4 (2021)", "Kyzivat (2022)") ~ "Arctic-boreal",
           dataset == "This Study (2025)" ~ "Northern Boreal",
           dataset == "Desrosiers (2022)" ~ "Temperate"
         ))

observations <- multistudy %>%
  group_by(dataset, veg_class) %>%
  summarise(n = n())

multistudy <- multistudy %>%
  mutate(dataset = factor(dataset,
                          levels = c("Desrosiers (2022)",
                                     "BAWLD-CH4 (2021)",
                                     "Kyzivat (2022)",
                                     "This Study (2025)")),
         ecosystem = case_when(
           dataset %in% c("BAWLD-CH4 (2021)", "Kyzivat (2022)") ~ "Arctic-boreal",
           dataset == "This Study (2025)" ~ "Northern Boreal",
           dataset == "Desrosiers (2022)" ~ "Temperate"
         ))   # order you want