# MSc Project: A little sedge goes a long way: emergent vegetation as a key driver of littoral methane flux

# Script containing the stats contained within thesis and in prep. paper.

# 03/14/2026 (Last Update: 08/17/2026)
# Author: Kelsey McGuire
# kmcgu@mail.ubc.ca; kmcguire.9@outlook.com

# LOAD LOCAL PACKAGES ----
source(here::here('code_DataCleaning', 'EVCH4_Data.R'))
source(here::here('code_PaperFiguresStats', 'EVCH4_ModelRunning.R'))

# install.packages("glmmTMB")
# install.packages("statmod")
# install.packages("tweedie")
library(statmod)
library(glmmTMB)
library(emmeans)

# GENERALIZED LINEAR MIXED MODELS ----
## DIFF X VEG X ZONE ----
DCH4VegGLM <- glmmTMB(log_cch4_flux_mgm2d1 ~ veg_class + zone + (1|site),
                      data = no_clean_data,
                      family = Gamma(link = "log"))

aov(DCH4VegGLM) # this is not a true anova, it is looking at fixed effects and their impact
summary(DCH4VegGLM) # effect sizes

# Residuals
plot(fitted(DCH4VegGLM), resid(DCH4VegGLM), 
      ylab="Residuals", xlab="CH4 Flux", 
      main="DCH4 Residuals - Veg + Zone") 
abline(0, 0)                  # the horizon

# Pairs
emmeans(DCH4VegGLM, pairwise ~ veg_class) # shows pairwise comps for veg
emmeans(DCH4VegGLM, pairwise ~ zone) # for zone, but not veg type specific
emmeans(DCH4VegGLM, pairwise ~ veg_class * zone) # for zone, veg type specific

## EBULL X VEG X ZONE ----
ECH4VegGLMzi <- glmmTMB(ebull_flux_ch4 ~ veg_class + zone + (1|site),
                        data = bubdata,
                        family = ziGamma(link = "log"), 
                        ziformula = ~ 1)
aov(ECH4VegGLMzi)
summary(ECH4VegGLMzi)

emmeans(ECH4VegGLMzi, pairwise ~ veg_class) # shows pairwise comps for veg
emmeans(ECH4VegGLMzi, pairwise ~ zone) # for zone, but not veg type specific
emmeans(ECH4VegGLMzi, pairwise ~ veg_class * zone) # for zone, veg type specific

## BUB X VEG X ZONE ----
BCH4VegGLMzi <- glmmTMB(log(bub_flux_lm2d1 + 1) ~ veg_class + zone + (1|site),
                        data = bubdata,
                        family = ziGamma(link = "log"), 
                        ziformula = ~ 1)
aov(BCH4VegGLMzi)
summary(BCH4VegGLMzi)

# Residuals
plot(fitted(BCH4VegGLMzi), resid(BCH4VegGLMzi), 
     ylab="Residuals", xlab="BUB Flux", 
     main="BCH4 Residuals - Veg + Zone") 
abline(0, 0)                  # the horizon

emmeans(BCH4VegGLMzi, pairwise ~ veg_class) # shows pairwise comps for veg
emmeans(BCH4VegGLMzi, pairwise ~ zone) # for zone, but not veg type specific
emmeans(BCH4VegGLMzi, pairwise ~ veg_class * zone) # for zone, veg type specific

names(no_clean_data)

## TOTAL FLUX X VEG ----
hist(ch4_rel_summary$flux_adj)

tt_rel_cont <- ch4_rel_summary %>%
  group_by(process, veg_class) %>%
  summarise(across(where(is.numeric),
                   list(
                     mean = ~mean(.x, na.rm = TRUE),
                     sd = ~sd(.x, na.rm = TRUE),
                     se = ~sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))
                   ),
                   .names = "{.col}_{.fn}"
  )) %>%
  select(process, veg_class, flux_adj_mean) %>%
  pivot_wider(names_from = veg_class,
              values_from = flux_adj_mean)

TCH4VegWR <- wilcox.test(tt_rel_cont$Vegetated, tt_rel_cont$Open,
                         paired = TRUE)
TCH4VegWR # diff

# ENV PARA. X VEG X ZONE ----
## pH - NS ----
pHVegGLM <- glmmTMB(ph ~ veg_class + (1|site),
                       data = no_clean_data,
                       family = Gamma(link = "log"))

aov(pHVegGLM)
summary(pHVegGLM) # NS

## Cond - NS----
CondVegGLM <- glmmTMB(cond_mScm ~ veg_class + (1|site),
                       data = no_clean_data,
                       family = Gamma(link = "log"))

aov(CondVegGLM)
summary(CondVegGLM) # NS

## Depth - SIG (veg + zone differences) ----
DepthVegGLM <- glmmTMB(depth ~ veg_class + zone + (1|site),
                    data = no_clean_data,
                    family = Gamma(link = "log"))

aov(DepthVegGLM)
summary(DepthVegGLM) # SIG - VEG

# Residuals
plot(fitted(DepthVegGLM), resid(DepthVegGLM), 
     ylab="Residuals", xlab="DEPTH", 
     main="Depth Residuals - Veg + Zone") 
abline(0, 0)                  # the horizon

emmeans(DepthVegGLM, pairwise ~ veg_class) # shows pairwise comps for veg
emmeans(DepthVegGLM, pairwise ~ zone) # for zone, but not veg type specific
emmeans(DepthVegGLM, pairwise ~ veg_class * zone) # for zone, veg type specific

## DO - NS ----
DOVegGLM <- glmmTMB(DO ~ veg_class + zone + (1|site),
                    data = no_clean_data,
                    family = Gamma(link = "log"))

aov(DOVegGLM)
summary(DOVegGLM) # NS

emmeans(DOVegGLM, pairwise ~ veg_class) # shows pairwise comps for veg
emmeans(DOVegGLM, pairwise ~ zone) # for zone, but not veg type specific
emmeans(DOVegGLM, pairwise ~ veg_class * zone) # for zone, veg type specific

## TDS - NS ----
TDSVegGLM <- glmmTMB(tds_mgL ~ veg_class + zone + (1|site),
                    data = no_clean_data,
                    family = Gamma(link = "log"))

aov(TDSVegGLM)
summary(TDSVegGLM) # NS

## sed temps; NS ----
### daily 
DSedTempVegGLM <- glmmTMB(temp_sed_daily ~ veg_class + zone + (1|site),
                     data = no_clean_data,
                     family = Gamma(link = "log"))

aov(DSedTempVegGLM)
summary(DSedTempVegGLM) # NS

### 3-day 
D3SedTempVegGLM <- glmmTMB(temp_sed_3day ~ veg_class + zone + (1|site),
                          data = no_clean_data,
                          family = Gamma(link = "log"))

aov(D3SedTempVegGLM)
summary(D3SedTempVegGLM) # NS

### 7-day 
D7SedTempVegGLM <- glmmTMB(temp_sed_7day ~ veg_class + zone + (1|site),
                           data = no_clean_data,
                           family = Gamma(link = "log"))

aov(D7SedTempVegGLM)
summary(D7SedTempVegGLM) # NS

## bot temps; NS ----
### daily 
DBotTempVegGLM <- glmmTMB(temp_bottom_daily ~ veg_class + zone + (1|site),
                          data = no_clean_data,
                          family = Gamma(link = "log"))

aov(DBotTempVegGLM)
summary(DBotTempVegGLM) # NS

### 3-day 
D3BotTempVegGLM <- glmmTMB(temp_bottom_3day ~ veg_class + zone + (1|site),
                           data = no_clean_data,
                           family = Gamma(link = "log"))

aov(D3BotTempVegGLM)
summary(D3BotTempVegGLM) # NS

### 7-day 
D7BotTempVegGLM <- glmmTMB(temp_bottom_7day ~ veg_class + zone + (1|site),
                           data = no_clean_data,
                           family = Gamma(link = "log"))

aov(D7BotTempVegGLM)
summary(D7BotTempVegGLM) # NS

## surf temps; NS ----
### daily 
DSurfTempVegGLM <- glmmTMB(temp_surface_daily ~ veg_class + zone + (1|site),
                          data = no_clean_data,
                          family = Gamma(link = "log"))

aov(DSurfTempVegGLM)
summary(DSurfTempVegGLM) # NS

### 3-day 
D3SurfTempVegGLM <- glmmTMB(temp_surface_3day ~ veg_class + zone + (1|site),
                           data = no_clean_data,
                           family = Gamma(link = "log"))

aov(D3SurfTempVegGLM)
summary(D3SurfTempVegGLM) # NS

### 7-day 
D7SurfTempVegGLM <- glmmTMB(temp_surface_7day ~ veg_class + zone + (1|site),
                           data = no_clean_data,
                           family = Gamma(link = "log"))

aov(D7SurfTempVegGLM)
summary(D7SurfTempVegGLM) # NS

# DI,NUTRIENTS X VEG X ZONE ----
water_chem <- CleanWaterSamples('~/Desktop/masters/data/msc-aquatic-ch4/surfwater-samples/surfwater-csv') # function from the CleanData file

## Fluoride - NS ----
FVegGLM <- glmmTMB(F_mgL ~ veg_presence + zone + (1|sample_date),
                    data = water_chem,
                    family = Gamma(link = "log"))

aov(FVegGLM)
summary(FVegGLM) # NS

## Chlorine - NS ----
ClVegGLM <- glmmTMB(Cl_mgL ~ veg_presence + zone + (1|sample_date),
                   data = water_chem,
                   family = Gamma(link = "log"))

aov(ClVegGLM)
summary(ClVegGLM) 

emmeans(ClVegGLM, pairwise ~ veg_presence) # shows pairwise comps for veg
emmeans(ClVegGLM, pairwise ~ zone) # shows pairwise comps for zone

## TN - SIG VEG, LESS IN OPEN ----
TNVegGLM <- glmmTMB(`TN_mgL` ~ veg_presence + zone + (1|sample_date),
                    data = water_chem,
                    family = Gamma(link = "log"))

aov(TNVegGLM)
summary(TNVegGLM)

emmeans(TNVegGLM, pairwise ~ veg_presence) # shows pairwise comps for zone

## NH4 - SIG, ZONE 3 + VEG ----
NH4VegGLM <- glmmTMB(`NH4-N_ugL` ~ veg_presence + zone + (1|zone),
                    data = water_chem,
                    family = Gamma(link = "log"))

aov(NH4VegGLM)
summary(NH4VegGLM)

emmeans(NH4VegGLM, pairwise ~ veg_presence) # shows pairwise comps for veg
emmeans(NH4VegGLM, pairwise ~ zone) # shows pairwise comps for zone

## NO3 - VEG SIG ----
NO3VegGLM <- glmmTMB(`NO3-N_mgL` ~ veg_presence * zone + (1|zone),
                    data = water_chem,
                    family = Gamma(link = "log"))

aov(NO3VegGLM)
summary(NO3VegGLM) 

emmeans(NO3VegGLM, pairwise ~ zone) # shows pairwise comps for zone
emmeans(NO3VegGLM, pairwise ~ veg_presence) # shows pairwise comps for veg

## SO4 - NS ----
SO4VegGLM <- glmmTMB(`S04-S_mgL` ~ veg_presence + zone + (1|zone),
                     data = water_chem,
                     family = Gamma(link = "log"))

aov(SO4VegGLM)
summary(SO4VegGLM) 

## PO4 - SIG Z2-Z3 ----
PO4VegGLM <- glmmTMB(`PO4-P_ugL` ~ veg_presence + zone + (1|zone),
                     data = water_chem,
                     family = Gamma(link = "log"))

aov(PO4VegGLM)
summary(PO4VegGLM) 

emmeans(PO4VegGLM, pairwise ~ zone) # shows pairwise comps for zone

## DOC - NS ----
DOCVegGLM <- glmmTMB(`DOC_mgL` ~ veg_presence + zone + (1|zone),
                     data = water_chem,
                     family = Gamma(link = "log"))

aov(DOCVegGLM)
summary(DOCVegGLM) 

## TOC - NS ----
TOCVegGLM <- glmmTMB(`TOC_mgL` ~ veg_presence + zone + (1|zone),
                     data = water_chem,
                     family = Gamma(link = "log"))

aov(TOCVegGLM)
summary(TOCVegGLM)

## TC - SIG Z2 ----
TCVegGLM <- glmmTMB(`TC_mgL` ~ veg_presence + zone + (1|zone),
                     data = water_chem,
                     family = Gamma(link = "log"))

aov(TCVegGLM)
summary(TCVegGLM)

emmeans(TCVegGLM, pairwise ~ zone) # shows pairwise comps for zone

# DISS X VEG X ZONE ----
## DISS CH4 - NS
DissCH4VegGLM <- glmmTMB(dissolvedCH4 ~ veg_class * zone + (1|site),
                   data = no_clean_data,
                   family = Gamma(link = "log"))

aov(DissCH4VegGLM)
summary(DissCH4VegGLM) # NS

emmeans(DissCH4VegGLM, pairwise ~ veg_class) # shows pairwise comps for veg
emmeans(DissCH4VegGLM, pairwise ~ zone) # for zone, but not veg type specific
emmeans(DissCH4VegGLM, pairwise ~ veg_class * zone) # for zone, veg type specific



# PLANT METRICS X ZONE ----
## Biomass - SIG
BioZoneGLM <- glmmTMB(biomass_g ~ zone + (1|site),
                      data = no_veg_ch4,
                      family = Gamma(link = "log"))

aov(BioZoneGLM)
summary(BioZoneGLM) # SIG 

# Residuals
plot(fitted(BioZoneGLM), resid(BioZoneGLM), 
     ylab="Residuals", xlab="Biomass", 
     main="Biomass Residuals - Zone") 
abline(0, 0)                  # the horizon

# Pairs
emmeans(BioZoneGLM, pairwise ~ zone) # for zone, but not veg type specific

## Num. Shoots - SIG (same as biomass)
NShootsZoneGLM <- glmmTMB(num_shoots ~ zone + (1|site),
                      data = no_veg_ch4,
                      family = Gamma(link = "log"))

aov(NShootsZoneGLM)
summary(NShootsZoneGLM) # SIG 

# Residuals
plot(fitted(NShootsZoneGLM), resid(NShootsZoneGLM), 
     ylab="Residuals", xlab="Num. Shoots", 
     main="Num Shoots Residuals - Zone") 
abline(0, 0)                  # the horizon

# Pairs
emmeans(NShootsZoneGLM, pairwise ~ zone) # for zone, but not veg type specific

## Abovewater Volume - SIG
ABWVZoneGLM <- glmmTMB(abovewater_vol_cm3 ~ zone + (1|site),
                      data = no_veg_ch4,
                      family = Gamma(link = "log"))

aov(ABWVZoneGLM)
summary(ABWVZoneGLM) # SIG 

# Residuals
plot(fitted(ABWVZoneGLM), resid(ABWVZoneGLM), 
     ylab="Residuals", xlab="ABW Volume", 
     main="ABW Volume Residuals - Zone") 
abline(0, 0)                  # the horizon

# Pairs
emmeans(ABWVZoneGLM, pairwise ~ zone) # for zone, but not veg type specific

## Total Volume - SIG (same as ABWV)
TVZoneGLM <- glmmTMB(total_vol_cm3 ~ zone + (1|site),
                       data = no_veg_ch4,
                       family = Gamma(link = "log"))

aov(TVZoneGLM)
summary(TVZoneGLM) # SIG 

# Residuals
plot(fitted(TVZoneGLM), resid(TVZoneGLM), 
     ylab="Residuals", xlab="Total Volume", 
     main="Total Volume Residuals - Zone") 
abline(0, 0)                  # the horizon

# Pairs
emmeans(TVZoneGLM, pairwise ~ zone) # for zone, but not veg type specific

## GCC - SIG 
GCCZoneGLM <- glmmTMB(GCC ~ zone + (1|site),
                     data = no_veg_ch4,
                     family = Gamma(link = "log"))

aov(GCCZoneGLM)
summary(GCCZoneGLM) # SIG 

# Residuals
plot(fitted(GCCZoneGLM), resid(GCCZoneGLM), 
     ylab="Residuals", xlab="GCC", 
     main="GCC Residuals - Zone") 
abline(0, 0)                  # the horizon

# Pairs
emmeans(GCCZoneGLM, pairwise ~ zone) # for zone, but not veg type specific


# ENV. VARIABLES ----
## DIFFUSIVE (ENTIRE) ----
DCH4Envglmm <- runglmm(no_clean_data, log_cch4_flux_mgm2d1, 
                       c("depth", "ph", "cond_mScm", "tds_mgL", "wat_temp_C", "air_temp_C", "press_kPa", "DO", 
                         "light_lux_mean_BOT", "light_lux_mean_SUR", 
                         "temp_air_daily", "temp_bottom_daily", "temp_surface_daily", "temp_sed_daily", 
                         "temp_air_3day", "temp_bottom_3day", "temp_surface_3day", "temp_sed_3day", 
                         "temp_air_7day", "temp_bottom_7day", "temp_surface_7day", "temp_sed_7day"),
                       random_effect = "site")
### SIG MODELS
DCH4Envglmm$summary$depth

#### TDS
DCH4Envglmm$summary$tds_mgL
performance::r2(DCH4Envglmm$model$tds_mgL) # CR2: 0.326; MR2: 0.054; p-val: 0.02096

# Residuals
plot(fitted(DCH4Envglmm$model$tds_mgL), resid(DCH4Envglmm$model$tds_mgL), 
     ylab="Residuals", xlab="TDS", 
     main="CH4 Flux x TDS Residuals") 
abline(0, 0) # the horizon

## ZA DIFFUSIVE (ENTIRE) ----
DCH4ZAEnvglmm <- runzaglmm(no_zone_avg, log_cch4_flux_mgm2d1_mean, 
                       c("depth_mean", "ph_mean", "cond_mScm_mean", "tds_mgL_mean", "press_kPa_mean", "DO_mean", 
                         "light_lux_mean_BOT_mean", "light_lux_mean_SUR_mean", "HOURLY_WIND_SPEED_mean",
                         "temp_air_daily_mean", "temp_bottom_daily_mean", "temp_surface_daily_mean", "temp_sed_daily_mean", 
                         "temp_air_3day_mean", "temp_bottom_3day_mean", "temp_surface_3day_mean", "temp_sed_3day_mean", 
                         "temp_air_7day_mean", "temp_bottom_7day_mean", "temp_surface_7day_mean", "temp_sed_7day_mean"),
                       random_effect = "zone")
### SIG MODELS
DCH4ZAEnvglmm$summary

#### DEPTH
DCH4ZAEnvglmm$summary$depth_mean # p - 0.00136
performance::r2(DCH4ZAEnvglmm$model$depth_mean) # Conditional R2: 0.255; Marginal R2: 0.237

# Residuals
plot(fitted(DCH4ZAEnvglmm$model$depth_mean), resid(DCH4ZAEnvglmm$model$depth_mean), 
     ylab="Residuals", xlab="Depth", 
     main="CH4 Flux x Depth Residuals") 
abline(0, 0) # the horizon

## DIFFUSIVE (VEG) ----
DCH4EnvVegglmm <- runglmm(no_veg_ch4, log_cch4_flux_mgm2d1, 
                       c("depth", "ph", "cond_mScm", "tds_mgL", "wat_temp_C", "air_temp_C", "press_kPa", "DO", 
                         "light_lux_mean_BOT", "light_lux_mean_SUR", "HOURLY_WIND_SPEED",
                         "temp_air_daily", "temp_bottom_daily", "temp_surface_daily", "temp_sed_daily", 
                         "temp_air_3day", "temp_bottom_3day", "temp_surface_3day", "temp_sed_3day", 
                         "temp_air_7day", "temp_bottom_7day", "temp_surface_7day", "temp_sed_7day",
                         "GCC", "ExG", "biomass_g", "abovewater_vol_cm3", "total_vol_cm3", "perc_cover"),
                       random_effect = "site")
### SIG MODELS - NADA
DCH4EnvVegglmm$summary$temp_sed_3day

## ZA DIFFUSIVE (VEG) ----
DCH4ZAEnvVegglmm <- runzaglmm(no_veg_za, log_cch4_flux_mgm2d1_mean, 
                         c("depth_mean", "ph_mean", "cond_mScm_mean", "tds_mgL_mean", "press_kPa_mean", "DO_mean", 
                           "light_lux_mean_BOT_mean", "light_lux_mean_SUR_mean", "HOURLY_WIND_SPEED_mean",
                           "temp_air_daily_mean", "temp_bottom_daily_mean", "temp_surface_daily_mean", "temp_sed_daily_mean", 
                           "temp_air_3day_mean", "temp_bottom_3day_mean", "temp_surface_3day_mean", "temp_sed_3day_mean", 
                           "temp_air_7day_mean", "temp_bottom_7day_mean", "temp_surface_7day_mean", "temp_sed_7day_mean",
                           "GCC_mean", "ExG_mean", "biomass_g_mean", "abovewater_vol_cm3_mean", "total_vol_cm3_mean", "perc_cover_mean"),
                         random_effect = "zone")
### SIG MODELS - NADA
DCH4ZAEnvVegglmm$summary

## DIFFUSIVE (OPEN) ----
DCH4EnvOpenglmm <- runglmm(no_open_ch4, log_cch4_flux_mgm2d1, 
                          c("depth", "ph", "cond_mScm", "tds_mgL", "wat_temp_C", "air_temp_C", "press_kPa", "DO", 
                            "light_lux_mean_BOT", "light_lux_mean_SUR", "HOURLY_WIND_SPEED",
                            "temp_air_daily", "temp_bottom_daily", "temp_surface_daily", "temp_sed_daily", 
                            "temp_air_3day", "temp_bottom_3day", "temp_surface_3day", "temp_sed_3day", 
                            "temp_air_7day", "temp_bottom_7day", "temp_surface_7day", "temp_sed_7day"),
                          random_effect = "site")
### SIG MODELS - NADA
DCH4EnvOpenglmm$summary$depth

## ZA DIFFUSIVE (OPEN) ----
DCH4ZAEnvOpenglmm <- runzaglmm(no_open_za, log_cch4_flux_mgm2d1_mean, 
                            c("depth_mean", "ph_mean", "cond_mScm_mean", "tds_mgL_mean", "press_kPa_mean", "DO_mean", 
                              "light_lux_mean_BOT_mean", "light_lux_mean_SUR_mean", "HOURLY_WIND_SPEED_mean",
                              "temp_air_daily_mean", "temp_bottom_daily_mean", "temp_surface_daily_mean", "temp_sed_daily_mean", 
                              "temp_air_3day_mean", "temp_bottom_3day_mean", "temp_surface_3day_mean", "temp_sed_3day_mean", 
                              "temp_air_7day_mean", "temp_bottom_7day_mean", "temp_surface_7day_mean", "temp_sed_7day_mean"),
                            random_effect = "zone")
### SIG MODELS - NADA
DCH4ZAEnvOpenglmm$summary

#### conductivity 
DCH4ZAEnvOpenglmm$summary$cond_mScm_mean
performance::r2(DCH4ZAEnvglmm$model$cond_mScm_mean) # Marginal R2: 0.073

# Residuals
plot(fitted(DCH4ZAEnvglmm$model$cond_mScm_mean), resid(DCH4ZAEnvglmm$model$cond_mScm_mean), 
     ylab="Residuals", xlab="Conductivity", 
     main="CH4 Flux x Conductivity Residuals") 
abline(0, 0) # the horizon

## PM ----
PMCH4Envglmm <- runglmm(no_veg_ch4, log_plant_cch4_flux_mgm2d1, 
                        c("depth", "ph", "cond_mScm", "tds_mgL", "press_kPa", "DO", 
                          "light_lux_mean_BOT", "light_lux_mean_SUR", "HOURLY_WIND_SPEED",
                          "temp_air_daily", "temp_bottom_daily", "temp_surface_daily", "temp_sed_daily", 
                          "temp_air_3day", "temp_bottom_3day", "temp_surface_3day", "temp_sed_3day", 
                          "temp_air_7day", "temp_bottom_7day", "temp_surface_7day", "temp_sed_7day",
                          "GCC", "ExG", "perc_cover", "num_shoots", "biomass_g", "abovewater_vol_cm3", "total_vol_cm3"),
                        random_effect = "site")
### SIG MODELS - NADA
PMCH4Envglmm$summary$GCC

## PM (ZA) ----
PMZAEnvglmm <- runzaglmm(no_veg_za, log_plant_cch4_flux_mgm2d1_mean, 
                            c("depth_mean", "ph_mean", "cond_mScm_mean", "tds_mgL_mean", "press_kPa_mean", "DO_mean", 
                              "light_lux_mean_BOT_mean", "light_lux_mean_SUR_mean", "HOURLY_WIND_SPEED_mean",
                              "temp_air_daily_mean", "temp_bottom_daily_mean", "temp_surface_daily_mean", "temp_sed_daily_mean", 
                              "temp_air_3day_mean", "temp_bottom_3day_mean", "temp_surface_3day_mean", "temp_sed_3day_mean", 
                              "temp_air_7day_mean", "temp_bottom_7day_mean", "temp_surface_7day_mean", "temp_sed_7day_mean",
                              "GCC_mean", "ExG_mean", "biomass_g_mean", "abovewater_vol_cm3_mean", "total_vol_cm3_mean", "perc_cover_mean"),
                            random_effect = "zone")
### SIG MODELS - NADA
PMZAEnvglmm$summary

## PM (no extreme biomass) ----
veg_ch4_clean <- no_veg_ch4 %>%
  filter(abovewater_vol_cm3 < 400) # remove extreme biomass

PMCH4CleanEnvglmm <- runglmm(veg_ch4_clean, log_plant_cch4_flux_mgm2d1, 
                          c("depth", "ph", "cond_mScm", "tds_mgL", "press_kPa", "DO", 
                            "light_lux_mean_BOT", "light_lux_mean_SUR", 
                            "temp_air_daily", "temp_bottom_daily", "temp_surface_daily", "temp_sed_daily", 
                            "temp_air_3day", "temp_bottom_3day", "temp_surface_3day", "temp_sed_3day", 
                            "temp_air_7day", "temp_bottom_7day", "temp_surface_7day", "temp_sed_7day",
                            "GCC", "ExG", "perc_cover", "num_shoots", "biomass_g", "abovewater_vol_cm3", "total_vol_cm3"),
                          random_effect = "site")
### SIG MODELS - NADA
PMCH4CleanEnvglmm$summary$GCC

## PM (ZA - no extreme bio) ----
veg_za_clean <- no_veg_za %>%
  filter(abovewater_vol_cm3_mean < 300) 

PMZAEnvglmm <- runzaglmm(veg_za_clean, log_plant_cch4_flux_mgm2d1_mean, 
                       c("depth_mean", "ph_mean", "cond_mScm_mean", "tds_mgL_mean", "press_kPa_mean", "DO_mean", 
                         "light_lux_mean_BOT_mean", "light_lux_mean_SUR_mean", 
                         "temp_air_daily_mean", "temp_bottom_daily_mean", "temp_surface_daily_mean", "temp_sed_daily_mean", 
                         "temp_air_3day_mean", "temp_bottom_3day_mean", "temp_surface_3day_mean", "temp_sed_3day_mean", 
                         "temp_air_7day_mean", "temp_bottom_7day_mean", "temp_surface_7day_mean", "temp_sed_7day_mean",
                         "GCC_mean", "ExG_mean", "biomass_g_mean", "abovewater_vol_cm3_mean", "total_vol_cm3_mean", "perc_cover_mean"),
                       random_effect = "zone")
### SIG MODELS - NADA
PMZAEnvglmm$summary

## EBULL (ENTIRE) ----
ECH4Envglmm <- runziglm(no_clean_data, log_ebull_flux_ch4, 
                             c("depth_bt", "ph", "cond_mScm", "tds_mgL", "wat_temp_C", "air_temp_C", "press_kPa", "DO", 
                               "light_lux_mean_BOT", "light_lux_mean_SUR", 
                               "temp_air_daily", "temp_bottom_daily", "temp_surface_daily", "temp_sed_daily", 
                               "temp_air_3day", "temp_bottom_3day", "temp_surface_3day", "temp_sed_3day", 
                               "temp_air_7day", "temp_bottom_7day", "temp_surface_7day", "temp_sed_7day"),
                             random_effect = "site")
## SIG MODELS - NADA
ECH4Envglmm$summary

## EBULL (ENTIRE - ZA) ----
ECH4ZAEnvglmm <- runzigzalm(no_zone_avg, log_ebull_flux_ch4_mean, 
                        c("depth_bt_mean", "depth_mean", "ph_mean", "cond_mScm_mean", "tds_mgL_mean", "press_kPa_mean", "DO_mean", 
                          "light_lux_mean_BOT_mean", "light_lux_mean_SUR_mean", "HOURLY_WIND_SPEED_mean",
                          "temp_air_daily_mean", "temp_bottom_daily_mean", "temp_surface_daily_mean", "temp_sed_daily_mean", 
                          "temp_air_3day_mean", "temp_bottom_3day_mean", "temp_surface_3day_mean", "temp_sed_3day_mean", 
                          "temp_air_7day_mean", "temp_bottom_7day_mean", "temp_surface_7day_mean", "temp_sed_7day_mean"),
                        random_effect = "zone")
## SIG MODELS - NADA
ECH4ZAEnvglmm$summary

## EBULL (Veg) ----
ECH4EnvVegglmm <- runziglm(veg_ch4, log_ebull_flux_ch4, 
                             c("depth_bt", "ph", "cond_mScm", "tds_mgL", "wat_temp_C", "air_temp_C", "press_kPa", "DO", 
                               "light_lux_mean_BOT", "light_lux_mean_SUR", "HOURLY_WIND_SPEED",
                               "temp_air_daily", "temp_bottom_daily", "temp_surface_daily", "temp_sed_daily", 
                               "temp_air_3day", "temp_bottom_3day", "temp_surface_3day", "temp_sed_3day", 
                               "temp_air_7day", "temp_bottom_7day", "temp_surface_7day", "temp_sed_7day"),
                          random_effect = "site")
## SIG MODELS - NADA
ECH4EnvVegglmm$summary$temp_sed_3day

ECH4VegDepthGLMM <- glmECH4VegDepthGLMM <- glmECH4VegDepthGLMM <- glmmTMB(log_ebull_flux_ch4 ~ depth_bt + zone + (1|site),
                         data = veg_ebull,
                         family = ziGamma(link = "log"), 
                         ziformula = ~ 1)
summary(ECH4VegDepthGLMM)

## EBULL (VEG - ZA) ----
ECH4ZAVegEnvglmm <- runzigzalm(no_veg_za, log_ebull_flux_ch4_mean, 
                          c("depth_bt_mean", "ph_mean", "cond_mScm_mean", "tds_mgL_mean", "press_kPa_mean", "DO_mean", 
                            "light_lux_mean_BOT_mean", "light_lux_mean_SUR_mean", "HOURLY_WIND_SPEED_mean",
                            "temp_air_daily_mean", "temp_bottom_daily_mean", "temp_surface_daily_mean", "temp_sed_daily_mean", 
                            "temp_air_3day_mean", "temp_bottom_3day_mean", "temp_surface_3day_mean", "temp_sed_3day_mean", 
                            "temp_air_7day_mean", "temp_bottom_7day_mean", "temp_surface_7day_mean", "temp_sed_7day_mean"),
                          random_effect = "zone")
## SIG MODELS - NADA
ECH4ZAVegEnvglmm$summary

## EBULL (Open) ----
ECH4EnvOpenglmm <- runzigzalm(no_open_ch4, log_ebull_flux_ch4, 
                            c("depth_bt", "ph", "cond_mScm", "tds_mgL", "wat_temp_C", "air_temp_C", "press_kPa", "DO", 
                              "light_lux_mean_BOT", "light_lux_mean_SUR", "HOURLY_WIND_SPEED",
                              "temp_air_daily", "temp_bottom_daily", "temp_surface_daily", "temp_sed_daily", 
                              "temp_air_3day", "temp_bottom_3day", "temp_surface_3day", "temp_sed_3day", 
                              "temp_air_7day", "temp_bottom_7day", "temp_surface_7day", "temp_sed_7day"),
                             random_effect = "site")

## SIG MODELS
ECH4EnvOpenglmm$summary$temp_sed_3day

ECH4OpenDepthGLMM <- glmmTMB(log_ebull_flux_ch4 ~ depth_bt + (1|subsample),
                             data = open_ebull,
                             family = ziGamma(link = "log"),
                             ziformula = ~ 1)
summary(ECH4OpenDepthGLMM)
performance::r2(ECH4OpenDepthGLMM)

# Residuals
plot(fitted(ECH4OpenDepthGLMM), resid(ECH4OpenDepthGLMM), 
     ylab="Residuals", xlab="Depth", 
     main="ECH4 Flux x Depth Residuals") 
abline(0, 0) # the horizon

## EBULL (OPEN - ZA) ----
ECH4EnvZAOpenglmm <- runzigzalm(no_open_za, log_ebull_flux_ch4_mean, 
                            c("depth_bt_mean", "ph_mean", "cond_mScm_mean", "tds_mgL_mean", "press_kPa_mean", "DO_mean", 
                              "light_lux_mean_BOT_mean", "light_lux_mean_SUR_mean",  "HOURLY_WIND_SPEED_mean",
                              "temp_air_daily_mean", "temp_bottom_daily_mean", "temp_surface_daily_mean", "temp_sed_daily_mean", 
                              "temp_air_3day_mean", "temp_bottom_3day_mean", "temp_surface_3day_mean", "temp_sed_3day_mean", 
                              "temp_air_7day_mean", "temp_bottom_7day_mean", "temp_surface_7day_mean", "temp_sed_7day_mean"),
                            random_effect = "zone")

## SIG MODELS
ECH4EnvZAOpenglmm$summary

# TIME SERIES ----
## DCH4 ----
DCH4Time <- glmmTMB(cch4_flux_mgm2d1_mean ~ sample_date,
                   data = no_zone_avg,
                   family = Gamma(link = "log"))
summary(DCH4Time)

DCH4OpenTime <- glmmTMB(log_cch4_flux_mgm2d1 ~ sample_date,
                        data = no_open_ch4,
                        family = Gamma(link = "log"))
summary(DCH4OpenTime)

DCH4VegTime <- glmmTMB(log_cch4_flux_mgm2d1 ~ sample_date,
                    data = no_veg_ch4,
                    family = Gamma(link = "log"))
summary(DCH4VegTime)

## PMCH4 ----
PMCH4Time <- glmmTMB(log_plant_cch4_flux_mgm2d1 ~ sample_date,
                    data = no_veg_ch4,
                    family = Gamma(link = "log"))
summary(PMCH4Time)

## ECH4 ----
ECH4Time <- glmmTMB(ebull_flux_ch4 ~ sample_date,
                     data = bubdata,
                    family = ziGamma(link = "log"), 
                    ziformula = ~ 1)
summary(ECH4Time)

open_bub_za <- bubdata_za %>%
  filter(veg_class == "Open")

ECH4OpenTime <- glmmTMB(log_ebull_flux_ch4 ~ sample_date,
                    data = open_ebull,
                    family = ziGamma(link = "log"), 
                    ziformula = ~ 1)
summary(ECH4OpenTime)

veg_bub_za <- bubdata_za %>%
  filter(veg_class == "Vegetated")

ECH4VegTime <- glmmTMB(log_ebull_flux_ch4 ~ sample_date,
                        data = veg_ebull,
                        family = ziGamma(link = "log"), 
                        ziformula = ~ 1)
summary(ECH4VegTime)

## GCC ----
GCCTime <- glmmTMB(GCC ~ sample_date + (1|site),
                   data = no_clean_data,
                   family = Gamma(link = "log"))
summary(GCCTime)
performance::r2(GCCTime)

# Residuals
plot(fitted(GCCTime), resid(GCCTime), 
     ylab="Residuals", xlab="GCC", 
     main="GCC x Sample Date Residuals") 
abline(0, 0) # the horizon

## FINAL MODELS ----

# ISOTOPES ----
hist(clean_iso$d13c_ch4) # normal, non-independent!

## floating chamber ----
d13CH4VegTTFC <- t.test(x = d13FC_veg$d13c_ch4, y = d13FC_open$d13c_ch4,
                        paired = T)
d13CH4VegTTFC # t = -0.62855, df = 10, p-value = 0.5437

## bubbles  ----
d13CH4VegTTBub <- t.test(x = d13BUB_veg$d13c_ch4, y = d13BUB_open$d13c_ch4,
                         paired = T)
d13CH4VegTTBub # t = 1.0923, df = 7, p-value = 0.3109

## LMER
## all iso ----
d13CH4VegLMER <- glmmTMB(d13c_ch4 ~ pathway * type + (1|site), 
                         data = clean_iso)

summary(d13CH4VegLMER) # pathway -> 0.000224 ***

emmeans(d13CH4VegLMER, pairwise ~ pathway) # shows pairwise comps for veg
emmeans(d13CH4VegLMER, pairwise ~ type) # shows pairwise comps for veg
emmeans(d13CH4VegLMER, pairwise ~ pathway * type) # shows pairwise comps for veg

# INCUBATION (NP) ----
# ## test for differences ----
# DailyIncProductionWide <- DailyIncProduction %>%
#   select(block, treatment, t_c, avg_production_dry) %>%
#   pivot_wider(names_from = treatment,
#               values_from = avg_production_dry)
# 
# DailyIncProductionW20 <- DailyIncProductionWide %>%
#   filter(t_c == "Average Field Temperature (20°C)")
# 
# DailyIncProductionW30 <- DailyIncProductionWide %>%
#   filter(t_c == "Maximum Field Temperature (30°C)")
# 
# ## test treatment pairs ----
# # CONTROL
# c_temp <- wilcox.test(DailyIncProductionW20$Control, 
#                       DailyIncProductionW30$Control, 
#                       paired = TRUE)
# c_temp
# 
# # SED
# s_temp <- wilcox.test(DailyIncProductionW20$`Water + Sediment`, 
#                       DailyIncProductionW30$`Water + Sediment`, 
#                       paired = TRUE)
# s_temp
# 
# # BENTHIC
# b_temp <- wilcox.test(DailyIncProductionW20$`Water + Benthic`, 
#                       DailyIncProductionW30$`Water + Benthic`, 
#                       paired = TRUE)
# b_temp
# 
# # SED + BENTHIC
# bs_temp <- wilcox.test(DailyIncProductionW20$`Water + Sediment + Benthic`, 
#                        DailyIncProductionW30$`Water + Sediment + Benthic`, 
#                        paired = TRUE)
# bs_temp
# 
# ## test treatment differences in temp controls ----
# DailyIncProduction20 <- DailyIncProduction %>%
#   filter(t_c == "Average Field Temperature (20°C)")
# 
# DailyIncProduction30 <- DailyIncProduction %>%
#   filter(t_c == "Maximum Field Temperature (30°C)")
# 
# ## 20 deg.
# treatments_20 <- friedman.test(avg_production_dry ~ treatment | time_step,
#                                data = DailyIncProduction20)
# treatments_20 # chi-squared = 9, df = 3, p-value = 0.02929
# 
# ### CONTROLS (NS)
# cs_20 <- wilcox.test(DailyIncProductionW20$`Control`, DailyIncProductionW20$`Water + Sediment`, paired = TRUE)
# cb_20 <- wilcox.test(DailyIncProductionW20$`Control`, DailyIncProductionW20$`Water + Benthic`, paired = TRUE)
# csb_20 <- wilcox.test(DailyIncProductionW20$`Control`, DailyIncProductionW20$`Water + Sediment + Benthic`, paired = TRUE)
# 
# ### SEDIMENT (NS)
# sb_20 <- wilcox.test(DailyIncProductionW20$`Water + Sediment`, DailyIncProductionW20$`Water + Benthic`, paired = TRUE)
# ssb_20 <- wilcox.test(DailyIncProductionW20$`Water + Sediment`, DailyIncProductionW20$`Water + Sediment + Benthic`, paired = TRUE)
# 
# ### BENTHIC (NS)
# bsb_20 <- wilcox.test(DailyIncProductionW20$`Water + Benthic`, DailyIncProductionW20$`Water + Sediment + Benthic`, paired = TRUE)
# 
# ## 30 deg.
# treatments_30 <- friedman.test(avg_production_dry ~ treatment | time_step,
#                                data = DailyIncProduction30)
# treatments_30 # chi-squared = 9, df = 3, p-value = 0.02929
# 
# ### CONTROLS (NS)
# cs_30 <- wilcox.test(DailyIncProductionW30$`Control`, DailyIncProductionW30$`Water + Sediment`, paired = TRUE)
# cb_30 <- wilcox.test(DailyIncProductionW30$`Control`, DailyIncProductionW30$`Water + Benthic`, paired = TRUE)
# csb_30 <- wilcox.test(DailyIncProductionW30$`Control`, DailyIncProductionW30$`Water + Sediment + Benthic`, paired = TRUE)
# 
# ### SEDIMENT (NS)
# sb_30 <- wilcox.test(DailyIncProductionW30$`Water + Sediment`, DailyIncProductionW30$`Water + Benthic`, paired = TRUE)
# ssb_30 <- wilcox.test(DailyIncProductionW30$`Water + Sediment`, DailyIncProductionW30$`Water + Sediment + Benthic`, paired = TRUE)
# 
# ### BENTHIC (NS)
# bsb_30 <- wilcox.test(DailyIncProductionW30$`Water + Benthic`, DailyIncProductionW30$`Water + Sediment + Benthic`, paired = TRUE)


# INCUBATION (PARA) ----
NoOutIncProd %>%
  ggplot() +
  geom_histogram(aes(x = nmol_ch4_dry_gw_1_d_1)) +
  facet_wrap(~t_c)

NoOutIncProd <- NoOutIncProd %>%
  filter(nmol_ch4_dry_gw_1_d_1 > 0)

NoOutIncProduction20 <- NoOutIncProd %>%
  filter(t_c == "Average Field Temperature (20°C)" &
           nmol_ch4_dry_gw_1_d_1 > 0)

NoOutIncProduction30 <- NoOutIncProd %>%
  filter(t_c == "Potential Warming Temperature (30°C)" &
           nmol_ch4_dry_gw_1_d_1 > 0)

## test differences on specific days ----
ProdGLM <- glmmTMB(nmol_ch4_dry_gw_1_d_1 ~ treatment + time_step + t_c + (1|replicate),
                   data = NoOutIncProd,
                   family = Gamma(link = "log"))
summary(ProdGLM)

emmeans::emmeans(ProdGLM, pairwise ~ t_c | treatment) # specific days 

## test differences on specific days (20 deg) ----
ProdGLM20 <- glmmTMB(nmol_ch4_dry_gw_1_d_1 ~ treatment * time_step + (1|replicate),
                     data = NoOutIncProduction20,
                     family = Gamma(link = "log"))
summary(ProdGLM20)

# find specific differences on days
emmeans::emmeans(ProdGLM20, pairwise ~ treatment | time_step)

emmeans::emmeans(ProdGLM20, pairwise ~ time_step | treatment)

## test differences on specific days (30 deg.)----
ProdGLM30 <- glmmTMB(nmol_ch4_dry_gw_1_d_1 ~ treatment * time_step + (1|replicate),
                     data = NoOutIncProduction30,
                     family = Gamma(link = "log"))
summary(ProdGLM30)

# find specific differences on days
emmeans::emmeans(ProdGLM30, pairwise ~ treatment | time_step)

emmeans::emmeans(ProdGLM30, pairwise ~ time_step | treatment)

## Q10 ----
Q10glmm <- glmmTMB(q10_dry ~ treatment,
                     data = Q10IncProd,
                     family = Gamma(link = "log"))
summary(Q10glmm)

## OM ----
OMglmm <- glmmTMB(avg_sed_om ~ treatment_short + t_c,
                   data = CleanIncNutrient,
                   family = Gamma(link = "log"))
summary(OMglmm)

emmeans::emmeans(OMglmm, pairwise ~ t_c | treatment_short) # within treatments, does it differ between temperatures?
emmeans::emmeans(OMglmm, pairwise ~ treatment_short | t_c) # within temperatures, are there differences in treatments?

## CN ----
CNglmm <- glmmTMB(avg_sed_cn ~ treatment_short + t_c,
                  data = CleanIncNutrient,
                  family = Gamma(link = "log"))
summary(CNglmm)

emmeans::emmeans(CNglmm, pairwise ~ t_c | treatment_short) # within treatments, does it differ between temperatures?
emmeans::emmeans(CNglmm, pairwise ~ treatment_short | t_c) # within temperatures, are there differences in treatments?


## PROD X OM ----
ProdOMglmm <- glmmTMB(avg_cumulative_prod ~ avg_sed_om + treatment_short + t_c,
                   data = CleanIncNutrient,
                   family = Gamma(link = "log"))
summary(ProdOMglmm)

emmeans::emmeans(ProdOMglmm, pairwise ~ t_c | treatment_short) # within treatments, does it differ between temperatures?
emmeans::emmeans(ProdOMglmm, pairwise ~ treatment_short | t_c) # within temperatures, are there differences in treatments?

## PROD X CN ----
ProdCNglmm <- glmmTMB(avg_cumulative_prod ~ avg_sed_cn + treatment_short + t_c,
                  data = CleanIncNutrient,
                  family = Gamma(link = "log"))
summary(ProdCNglmm)

emmeans::emmeans(ProdCNglmm, pairwise ~ t_c | treatment_short) # within treatments, does it differ between temperatures?
emmeans::emmeans(ProdCNglmm, pairwise ~ treatment_short | t_c) # within temperatures, are there differences in treatments?


