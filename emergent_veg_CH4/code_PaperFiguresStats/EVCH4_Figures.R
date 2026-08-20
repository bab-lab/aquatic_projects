# MSc Project: A little sedge goes a long way: emergent vegetation as a key driver of littoral methane flux

# Script containing the figure outputs for thesis and in prep. paper.

# 03/14/2026 (Last Update: 08/17/2026)
# Author: Kelsey McGuire
# kmcgu@mail.ubc.ca; kmcguire.9@outlook.com

# LOAD DATA AND LIBRARIES ----
source(here::here('code_DataCleaning', 'EVCH4_Data.R'))
source(here::here('code_PaperFiguresStats', 'EVCH4_ModelRunning.R'))
library(janitor)
library(patchwork)
library(ggnewscale)
library(scales)

# FIGURE 3: EBULL + DIFF ----
## pull observation data ----
diff_observations <- no_clean_data %>%
  group_by(veg_class) %>%
  summarise(n = n())

bub_observations <- bubdata %>%
  select(veg_class, bub_flux_lm2d1) %>%
  group_by(veg_class) %>%
  summarise(n = n())

ebull_observations <- perf_bubdata %>%
  select(veg_class, ebull_flux_ch4) %>%
  group_by(veg_class) %>%
  summarise(n = n())

d_means <- no_clean_data %>%
  group_by(veg_class) %>%
  summarise(mean = mean(cch4_flux_mgm2d1))

## DIFF ----
DCH4Veg_boxplot <- no_clean_data %>%
  ggplot(aes(x = veg_class, y = cch4_flux_mgm2d1)) +
  geom_jitter(aes(fill = veg_class), alpha = 0.5,
              width = 0.35, shape = 21, size = 2.5) +
  geom_boxplot(aes(fill = veg_class), alpha = 0.7,
               outlier.shape = NA) +
  geom_text(data = diff_observations, aes(x = veg_class, y = Inf, label = n),
            inherit.aes = FALSE,
            vjust = 6, size = 4) +
  geom_point(data = d_means, aes(x = veg_class, y = mean, fill = veg_class),
             shape = 4, size = 3) +
  geom_signif(comparisons = list(c("Open", "Vegetated")),
              map_signif_level=TRUE,
              y_position = 270) +
  scale_colour_manual(values = c("Open" = "black",
                                 "Vegetated" = "black")) +
  scale_fill_manual(values = c("Open" = "dodgerblue3",
                               "Vegetated" = "olivedrab3")) +
  theme_classic(base_size = 15)+
  coord_cartesian(ylim = c(0,300)) +
  ylab(expression(Chamber-Based~CH[4]~Flux~"("*mg~C~m^{-2}~day^{-1}*")")) + # changed from mg CH4 to mg C 
  theme(legend.position = "None",
        axis.title.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
# DCH4Veg_boxplot

## EBULL ----
e_means <- bubdata %>%
  group_by(veg_class) %>%
  summarise(mean = mean(ebull_flux_ch4))

ECH4Veg_boxplot <- bubdata %>%
  ggplot(aes(x = veg_class, y = ebull_flux_ch4)) +
  geom_jitter(aes(fill = veg_class), alpha = 0.5,
              width = 0.35, shape = 21, size = 2.5) +
  geom_boxplot(aes(fill = veg_class), alpha = 0.7,
               outlier.shape = NA) +
  geom_text(data = bub_observations, aes(x = veg_class, y = Inf, label = n),
            inherit.aes = FALSE,
            vjust = 5, size = 4) +
  geom_point(data = e_means, aes(x = veg_class, y = mean, fill = veg_class),
             shape = 4, size = 3) +
  #geom_signif(comparisons = list(c("Open", "Vegetated")),
  #            map_signif_level=TRUE,
  #            y_position = c(280),
  #            tip_length = 0.06) +
  scale_colour_manual(values = c("Open" = "black",
                                 "Vegetated" = "black")) +
  scale_fill_manual(values = c("Open" = "dodgerblue3",
                               "Vegetated" = "olivedrab3")) +
  theme_classic(base_size = 15)+
  ylab(expression(Ebullitive~CH[4]~Flux~"("*mg~C~m^{-2}~day^{-1}*")")) + # changed from mg CH4 to mg C 
  theme(legend.position = "None",
        axis.title.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
# ECH4Veg_boxplot

## PUT E+D TOGETHER ----
# DCH4Veg_boxplot | ECH4Veg_boxplot

## RELATIVE CONTRIBUTIONS ----
relativeCH4 <- ggplot(relcont_ch4_sum, 
       aes(x = veg_class, y = flux_adj_mean, 
           fill = fct_rev(process))) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.85) +
  scale_fill_manual(name = "Flux Pathway",
                    values = c("Diffusive" = 'dodgerblue',
                               "Plant-Mediated" = 'olivedrab3',
                               "Ebullitive" = "orangered3")) +
  ylab(expression(Average~CH[4]~Flux~"("*mg~C~m^{-2}~day^{-1}*")")) +
  #  ggtitle(expression("Wheeler Lake Total"~CH[4]~Flux)) +
  theme_classic(base_size = 15) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        axis.title.x = element_blank())
# relativeCH4

## put final plot together ----
fig3_fluxes <- DCH4Veg_boxplot / ECH4Veg_boxplot | relativeCH4
fig3_fluxes + plot_annotation(tag_levels = 'a',
                              tag_prefix = '(',
                              tag_suffix = ')')

# FIGURE 4: ISOTOPES ----
obs <- clean_iso %>%
  group_by(pathway, type) %>%
  summarise(n = n(),
            mean = mean(d13c_ch4))

fig4_isos <- clean_iso %>%
  ggplot(aes(x = pathway, y = d13c_ch4)) +
  geom_jitter(aes(fill = type), 
              position = position_jitterdodge(
                jitter.width = 0.4,
                dodge.width = 0.75
              ),
              alpha = 0.5, shape = 21, size = 2.5) +
  geom_boxplot(aes(fill = type),
               alpha = 0.7,
               outlier.shape = NA) +
  geom_text(data = obs, aes(x = pathway, y = Inf, label = n, group = type),
            position = position_dodge(width = 0.75),
            inherit.aes = FALSE,
            vjust = 5, size = 4) +
  geom_point(data = obs, aes(x = pathway, y = mean, group = type),
             position = position_dodge(width = 0.75),
             shape = 4, size = 3) +
  geom_signif(comparisons = list(c("Open", "Vegetated")),
              map_signif_level=c("***" = 0.001, 
                                 "**" = 0.01,
                                 "*" = 0.05),
              y_position = c(338),
              tip_length = 0.06) +
  scale_fill_manual(name = NULL,
                    values = c("Open" = "dodgerblue3",
                               "Vegetated" = "olivedrab3")) +
  scale_y_continuous(
    limits = c(min(clean_iso$d13c_ch4, na.rm = TRUE),
               max(clean_iso$d13c_ch4, na.rm = TRUE) + 10)
  ) + 
  theme_classic(base_size = 15)+
  ylab(expression(delta^{13}*C*"-"*CH[4]~"("*"\u2030"*")")) + 
  theme(axis.title.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
fig4_isos

# FIGURE 5: MIXED LINEAR EFFECTS ----
## depth ----
### PM DEPTH
pm_depth_sig_plot <- no_veg_ch4 %>%
  filter(depth > 25) %>%
  ggplot(aes(x=depth, 
             y=plant_cch4_flux_mgm2d1)) + 
  # geom_errorbar(aes(ymin=plant_cch4_flux_mgm2d1_mean-plant_cch4_flux_mgm2d1_se,
  #                   ymax=plant_cch4_flux_mgm2d1_mean+plant_cch4_flux_mgm2d1_se),
  #               width = 0, alpha = 0.6, colour = "grey70") +
  # geom_errorbarh(aes(xmin=depth_mean-depth_se,
  #                    xmax=depth_mean+depth_se),
  #                width = 0, alpha = 0.6, colour = "grey70") +
  # geom_jitter(data = depth_modelling,
  #            aes(x = depth, y = cch4_flux_mgm2d1, fill = veg_class),
  #            alpha = 0.3, shape = 21) +
  geom_hline(yintercept = 0, alpha = 0.6, linetype = "dashed") +
  geom_smooth(aes(x = depth,
                  y = (exp(log_cch4_flux_mgm2d1) - 1)), 
              fill="grey80", colour = "grey70",
              method = glm, alpha=0.3, show.legend = F)+
  geom_point(fill="darkolivegreen", colour="black", 
             shape=23, size = 2.5)+
  xlab(expression("Depth (cm)")) +
  ylab(expression(CH[4]~Flux~"("*mg~C~m^{-2}~day^{-1}*")")) +
  labs(fill = NULL) +
  theme_classic(base_size = 15) +
  theme(legend.position = "bottom", 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
# pm_depth_sig_plot

### DIFF/PM DEPTH (ZA)
diff_depth_sig_plot <- no_clean_data %>%
  filter(depth > 25) %>%
  ggplot(aes(x=depth, y=cch4_flux_mgm2d1, 
             colour=veg_class, fill = veg_class)) + 
  # geom_errorbar(aes(ymin=cch4_flux_mgm2d1_mean-cch4_flux_mgm2d1_se, 
  #                   ymax=cch4_flux_mgm2d1_mean+cch4_flux_mgm2d1_se),
  #               alpha = 0.6) +
  # geom_errorbarh(aes(xmin=depth_mean-depth_se, 
  #                    xmax=depth_mean+depth_se), width = 0,
  #                alpha = 0.6) +
  # geom_jitter(data = depth_modelling,
  #            aes(x = depth, y = cch4_flux_mgm2d1, fill = veg_class),
  #            alpha = 0.3, shape = 21) +
  geom_hline(yintercept = 0, alpha = 0.6, linetype = "dashed") +
  geom_smooth(aes(x = depth,
                  y = (exp(log_cch4_flux_mgm2d1) - 1),
                  fill=veg_class), method = glm, alpha=0.3, show.legend = F) +
  scale_fill_manual(values = c("Open"="grey80",
                               "Vegetated" = "grey80") )+
  new_scale_fill() + 
  geom_point(aes(fill=veg_class), colour="black", shape=21, size = 2.5)+
  scale_fill_manual(values = c("Open"="dodgerblue2",
                               "Vegetated"="olivedrab3"))+
  scale_color_manual(values = c("Open"="grey70",
                                "Vegetated"="grey70"))+
  xlab(expression("Depth (cm)")) +
  ylab(expression(atop(Chamber-Based~CH[4]~Flux,
                       "("*mg~C~m^{-2}~day^{-1}*")"))) +
  labs(fill = NULL) +
  theme_classic(base_size = 15) +
  theme(legend.position = "bottom", 
        axis.title.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
# diff_depth_sig_plot

### EBULL DEPTH (ZA)
ebull_depth_sig_plot <- bubdata %>%
  ggplot(aes(x=depth_bt, 
             y=ebull_flux_ch4, 
             colour=veg_class, fill = veg_class)) + 
  # geom_errorbar(aes(ymin=ebull_flux_ch4_mean-ebull_flux_ch4_se, 
  #                   ymax=ebull_flux_ch4_mean+ebull_flux_ch4_se),
  #               alpha = 0.7) +
  # geom_errorbarh(aes(xmin=depth_bt_mean-depth_bt_se, 
  #                    xmax=depth_bt_mean+depth_bt_se), width = 0,
  #                alpha = 0.7) +
  #geom_jitter(data = bubdata,
  #            aes(x = depth_bt, y = ebull_flux_ch4, fill = veg_class),
  #            alpha = 0.3, shape = 21, size = 1.5) +
  geom_hline(yintercept = 0, alpha = 0.6, linetype = "dashed") +
  geom_smooth(aes(x = depth_bt,
                  y = (exp(log_ebull_flux_ch4) - 1),
                  fill=veg_class), method = glm, alpha=0.3, show.legend = F) +
  scale_fill_manual(values = c("Open"="dodgerblue3",
                               "Vegetated" = "grey70") )+
  new_scale_fill() + 
  geom_point(aes(fill=veg_class), 
             colour="black", shape=21, size = 2.5)+
  scale_fill_manual(values = c("Open"="dodgerblue2",
                               "Vegetated"="olivedrab3"))+
  scale_color_manual(values = c("Open"="dodgerblue3",
                                "Vegetated"="grey70")) +
  xlab(expression("Depth (cm)")) +
  ylab(expression(CH[4]~Flux~"("*mg~C~m^{-2}~day^{-1}*")")) +
  labs(fill = NULL) +
  theme_classic(base_size = 15) +
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
# ebull_depth_sig_plot

depth <- pm_depth_sig_plot | diff_depth_sig_plot | ebull_depth_sig_plot

## sed temp ----
pm_sed_sig_plot <- no_veg_ch4 %>%
  ggplot(aes(x=temp_sed_3day, 
             y=plant_cch4_flux_mgm2d1, 
             colour=veg_class, fill = veg_class)) + 
  # geom_errorbar(aes(ymin=plant_cch4_flux_mgm2d1_mean-plant_cch4_flux_mgm2d1_se, 
  #                   ymax=plant_cch4_flux_mgm2d1_mean+plant_cch4_flux_mgm2d1_se),
  #               width = 0, alpha = 0.6) +
  # geom_errorbarh(aes(xmin=temp_sed_3day_mean-temp_sed_3day_se, 
  #                    xmax=temp_sed_3day_mean+temp_sed_3day_se), 
  #                width = 0, alpha = 0.6) +
  # geom_jitter(data = depth_modelling,
  #            aes(x = depth, y = cch4_flux_mgm2d1, fill = veg_class),
  #            alpha = 0.3, shape = 21) +
  geom_hline(yintercept = 0, alpha = 0.6, linetype = "dashed") +
  geom_smooth(aes(x = temp_sed_3day,
                  y = (exp(log_plant_cch4_flux_mgm2d1) - 1),
                  fill=veg_class), 
              method = glm, alpha=0.3, show.legend = F)+
  scale_fill_manual(values = c("Open"="grey80",
                               "Vegetated" = "grey80") )+
  new_scale_fill() + 
  geom_point(aes(fill=veg_class), colour="black", shape=23, size = 2.5)+
  scale_fill_manual(values = c("Open"="dodgerblue",
                               "Vegetated"="darkolivegreen"))+
  scale_color_manual(values = c("Open"="grey70",
                                "Vegetated"="grey70"))+
  xlab("Sediment Temperature (°C)") +
  ylab(expression(CH[4]~Flux~"("*mg~C~m^{-2}~day^{-1}*")")) +
  labs(fill = NULL) +
  theme_classic(base_size = 15) +
  theme(legend.position = "bottom", 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
# pm_sed_sig_plot

diff_sed_sig_plot <- no_clean_data %>%
  ggplot(aes(x=temp_sed_3day, 
             y=cch4_flux_mgm2d1, 
             colour=veg_class, fill = veg_class)) + 
  # geom_errorbar(aes(ymin=cch4_flux_mgm2d1_mean-cch4_flux_mgm2d1_se, 
  #                   ymax=cch4_flux_mgm2d1_mean+cch4_flux_mgm2d1_se),
  #               alpha = 0.6) +
  # geom_errorbarh(aes(xmin=temp_sed_3day_mean-temp_sed_3day_se, 
  #                    xmax=temp_sed_3day_mean+temp_sed_3day_se), width = 0,
  #                alpha = 0.6) +
  # geom_jitter(data = depth_modelling,
  #            aes(x = depth, y = cch4_flux_mgm2d1, fill = veg_class),
  #            alpha = 0.3, shape = 21) +
  geom_hline(yintercept = 0, alpha = 0.6, linetype = "dashed") +
  geom_smooth(aes(x = temp_sed_3day,
                  y = exp(log_cch4_flux_mgm2d1) - 1,
                  fill=veg_class), method = glm, alpha=0.3, show.legend = F)+
  scale_fill_manual(values = c("Open"="grey80",
                               "Vegetated" = "grey80") )+
  new_scale_fill() + 
  geom_point(aes(fill=veg_class), colour="black", shape=21, size = 2.5)+
  scale_fill_manual(values = c("Open"="dodgerblue2",
                               "Vegetated"="olivedrab3"))+
  scale_color_manual(values = c("Open"="grey70",
                                "Vegetated"="grey70"))+
  xlab("Sediment Temperature (°C)") +
  ylab(expression(CH[4]~Flux~"("*mg~C~m^{-2}~day^{-1}*")")) +
  labs(fill = NULL) +
  theme_classic(base_size = 15) +
  theme(legend.position = "bottom", 
        axis.title.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
# diff_sed_sig_plot

ebull_sed_sig_plot <- no_clean_data %>%
  ggplot(aes(x=temp_sed_3day,
             y=ebull_flux_ch4, 
             colour=veg_class, fill = veg_class)) + 
  # geom_errorbar(aes(ymin=ebull_flux_ch4_mean-ebull_flux_ch4_se, 
  #                   ymax=ebull_flux_ch4_mean+ebull_flux_ch4_se),
  #               width = 0, alpha = 0.6) +
  # geom_errorbarh(aes(xmin=temp_sed_3day_mean-temp_sed_3day_se, 
  #                    xmax=temp_sed_3day_mean+temp_sed_3day_se), 
  #                width = 0, alpha = 0.6) +
  # geom_jitter(data = depth_modelling,
  #            aes(x = depth, y = cch4_flux_mgm2d1, fill = veg_class),
  #            alpha = 0.3, shape = 21) +
  geom_hline(yintercept = 0, alpha = 0.6, linetype = "dashed") +
  geom_smooth(aes(x = temp_sed_3day,
                  y = exp(log_ebull_flux_ch4) - 1,
                  fill=veg_class), method = glm, alpha=0.3, show.legend = F)+
  scale_fill_manual(values = c("Open"="grey80",
                               "Vegetated" = "grey80") )+
  new_scale_fill() + 
  geom_point(aes(fill=veg_class), colour="black", shape=21, size = 2.5)+
  scale_fill_manual(values = c("Open"="dodgerblue2",
                               "Vegetated"="olivedrab3"))+
  scale_color_manual(values = c("Open"="grey70",
                                "Vegetated"="grey70"))+
  xlab("Sediment Temperature (°C)") +
  ylab(expression(CH[4]~Flux~"("*mg~C~m^{-2}~day^{-1}*")")) +
  labs(fill = NULL) +
  theme_classic(base_size = 15) +
  theme(legend.position = "bottom", 
        axis.title.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
# ebull_sed_sig_plot

sed <- pm_sed_sig_plot | diff_sed_sig_plot | ebull_sed_sig_plot

## biomass ----
pm_bio <- no_veg_ch4 %>%
  filter(biomass_g < 200) %>%
  ggplot(aes(x=biomass_g, 
             y=plant_cch4_flux_mgm2d1, 
             colour=veg_class, fill = veg_class)) + 
  # geom_errorbar(aes(ymin=plant_cch4_flux_mgm2d1_mean-plant_cch4_flux_mgm2d1_se, 
  #                   ymax=plant_cch4_flux_mgm2d1_mean+plant_cch4_flux_mgm2d1_se),
  #               width = 0, alpha = 0.6) +
  # geom_errorbarh(aes(xmin=biomass_g_mean-biomass_g_se, 
  #                    xmax=biomass_g_mean+biomass_g_se), 
  #                width = 0, alpha = 0.6) +
  # geom_jitter(data = depth_modelling,
  #            aes(x = depth, y = cch4_flux_mgm2d1, fill = veg_class),
  #            alpha = 0.3, shape = 21) +
  geom_hline(yintercept = 0, alpha = 0.6, linetype = "dashed") +
  geom_smooth(aes(x = biomass_g,
                  y = exp(log_plant_cch4_flux_mgm2d1) - 1,
                  fill=veg_class), method = glm, alpha=0.3, show.legend = F)+
  scale_fill_manual(values = c("Open"="grey80",
                               "Vegetated" = "grey80") )+
  new_scale_fill() + 
  geom_point(aes(fill=veg_class), colour="black", shape=23, size = 2.5)+
  scale_fill_manual(values = c("Open"="dodgerblue2",
                               "Vegetated"="darkolivegreen"))+
  scale_color_manual(values = c("Open"="grey70",
                                "Vegetated"="grey70"))+
  xlab(expression(Biomass~"("*g*")")) +
  ylab(expression(CH[4]~Flux~"("*mg~C~m^{-2}~day^{-1}*")")) +
  labs(fill = NULL) +
  theme_classic(base_size = 15) +
  theme(legend.position = "None", 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
# pm_bio

## aw vol ----
pm_abwvol <- no_veg_ch4 %>%
  filter(biomass_g < 200) %>%
  ggplot(aes(x=abovewater_vol_cm3, 
             y=plant_cch4_flux_mgm2d1, 
             colour=veg_class, fill = veg_class)) + 
  # geom_errorbar(aes(ymin=plant_cch4_flux_mgm2d1_mean-plant_cch4_flux_mgm2d1_se,
  #                   ymax=plant_cch4_flux_mgm2d1_mean+plant_cch4_flux_mgm2d1_se),
  #               width = 0, alpha = 0.6) +
  # geom_errorbarh(aes(xmin=biomass_g_mean-biomass_g_se,
  #                    xmax=biomass_g_mean+biomass_g_se),
  #                width = 0, alpha = 0.6) +
  # geom_jitter(data = depth_modelling,
  #            aes(x = depth, y = cch4_flux_mgm2d1, fill = veg_class),
  #            alpha = 0.3, shape = 21) +
  geom_hline(yintercept = 0, alpha = 0.6, linetype = "dashed") +
  geom_smooth(aes(x = abovewater_vol_cm3,
                  y = exp(log_plant_cch4_flux_mgm2d1) - 1,
                  fill=veg_class), method = glm, alpha=0.3, show.legend = F)+
  scale_fill_manual(values = c("Open"="grey80",
                               "Vegetated" = "grey70") )+
  new_scale_fill() + 
  geom_point(aes(fill=veg_class), colour="black", shape=23, size = 2.5)+
  scale_fill_manual(values = c("Open"="dodgerblue2",
                               "Vegetated"="darkolivegreen"))+
  scale_color_manual(values = c("Open"="grey70",
                                "Vegetated"="darkolivegreen"))+
  xlab(expression(Above~Water~Plant~Volume~"("*cm^{3}*")")) +
  ylab(expression(CH[4]~Flux~"("*mg~C~m^{-2}~day^{-1}*")")) +
  labs(fill = NULL) +
  theme_classic(base_size = 15) +
  theme(legend.position = "None", 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
# pm_abwvol

## GCC ----
pm_gcc_plot <- no_veg_ch4 %>%
  filter(biomass_g < 200) %>%
  ggplot(aes(x=GCC,
             y=plant_cch4_flux_mgm2d1, 
             colour=veg_class, fill = veg_class)) + 
  # geom_errorbar(aes(ymin=plant_cch4_flux_mgm2d1_mean-plant_cch4_flux_mgm2d1_se, 
  #                   ymax=plant_cch4_flux_mgm2d1_mean+plant_cch4_flux_mgm2d1_se),
  #               width = 0, alpha = 0.6) +
  # geom_errorbarh(aes(xmin=GCC_mean-GCC_se, 
  #                    xmax=GCC_mean+GCC_se), 
  #                width = 0, alpha = 0.6) +
  # geom_jitter(data = depth_modelling,
  #            aes(x = depth, y = cch4_flux_mgm2d1, fill = veg_class),
  #            alpha = 0.3, shape = 21) +
  geom_hline(yintercept = 0, alpha = 0.6, linetype = "dashed") +
  geom_smooth(aes(x = GCC,
                  y = exp(log_plant_cch4_flux_mgm2d1) - 1,
                  fill=veg_class), method = lm, alpha=0.3, show.legend = F)+
  scale_fill_manual(values = c("Open"="grey80",
                               "Vegetated" = "grey80") )+
  new_scale_fill() + 
  geom_point(aes(fill=veg_class), colour="black", shape=23, size = 2.5)+
  scale_fill_manual(values = c("Open"="dodgerblue2",
                               "Vegetated"="darkolivegreen"))+
  scale_color_manual(values = c("Open"="grey70",
                                "Vegetated"="grey70"))+
  xlab("Greenness Chromatic Coordinate") +
  ylab(expression(CH[4]~Flux~"("*mg~C~m^{-2}~day^{-1}*")")) +
  labs(fill = NULL) +
  theme_classic(base_size = 15) +
  theme(legend.position = "None", 
        axis.title.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
# pm_gcc_plot

## perc cover ----
pm_pc_plot <- no_veg_ch4 %>%
  filter(biomass_g < 200) %>%
  ggplot(aes(x=perc_cover, 
             y=plant_cch4_flux_mgm2d1, 
             colour=veg_class, fill = veg_class)) + 
  # geom_errorbar(aes(ymin=plant_cch4_flux_mgm2d1_mean-plant_cch4_flux_mgm2d1_se, 
  #                   ymax=plant_cch4_flux_mgm2d1_mean+plant_cch4_flux_mgm2d1_se),
  #               width = 0, alpha = 0.6) +
  # geom_errorbarh(aes(xmin=perc_cover_mean-perc_cover_se, 
  #                    xmax=perc_cover_mean+perc_cover_se), 
  #                width = 0, alpha = 0.6) +
  # geom_jitter(data = depth_modelling,
  #            aes(x = depth, y = cch4_flux_mgm2d1, fill = veg_class),
  #            alpha = 0.3, shape = 21) +
  geom_hline(yintercept = 0, alpha = 0.6, linetype = "dashed") +
  geom_smooth(aes(x = perc_cover,
                  y = exp(log_plant_cch4_flux_mgm2d1) - 1,
                  fill=veg_class), method = lm, alpha=0.3, show.legend = F)+
  scale_fill_manual(values = c("Open"="grey80",
                               "Vegetated" = "grey80") )+
  new_scale_fill() + 
  geom_point(aes(fill=veg_class), colour="black", shape=23, size = 2.5)+
  scale_fill_manual(values = c("Open"="dodgerblue2",
                               "Vegetated"="darkolivegreen"),
                    labels = c("Open" = "Open",
                               "Vegetated" = "Plant-Mediated"))+
  scale_color_manual(values = c("Open"="grey70",
                                "Vegetated"="grey70")) +
  xlab("Percent Cover (%)") +
  ylab(expression(CH[4]~Flux~"("*mg~C~m^{-2}~day^{-1}*")")) +
  labs(fill = NULL) +
  theme_classic(base_size = 15) +
  theme(legend.position = "None", 
        axis.title.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
# pm_pc_plot


## put final plot together ----
fig5_glmms <- depth / sed / (pm_bio | pm_pc_plot | pm_gcc_plot)
fig5_glmms + plot_annotation(tag_levels = 'a',
                             tag_prefix = '(',
                             tag_suffix = ')') +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# FIGURE 5 v2: MIXED LINEAR EFFECTS - ONLY SIG ----
## depth ----
### PM DEPTH
pm_depth_sig_plotv2 <- no_veg_ch4 %>%
  filter(depth > 25) %>%
  ggplot(aes(x=depth, 
             y=plant_cch4_flux_mgm2d1)) + 
  # geom_errorbar(aes(ymin=plant_cch4_flux_mgm2d1_mean-plant_cch4_flux_mgm2d1_se,
  #                   ymax=plant_cch4_flux_mgm2d1_mean+plant_cch4_flux_mgm2d1_se),
  #               width = 0, alpha = 0.6, colour = "grey70") +
  # geom_errorbarh(aes(xmin=depth_mean-depth_se,
  #                    xmax=depth_mean+depth_se),
  #                width = 0, alpha = 0.6, colour = "grey70") +
  # geom_jitter(data = depth_modelling,
  #            aes(x = depth, y = cch4_flux_mgm2d1, fill = veg_class),
  #            alpha = 0.3, shape = 21) +
  geom_hline(yintercept = 0, alpha = 0.6, linetype = "dashed") +
  #geom_smooth(aes(x = depth,
  #                y = (exp(log_cch4_flux_mgm2d1) - 1)), 
  #            fill="grey80", colour = "grey70",
  #            method = glm, alpha=0.3, show.legend = F)+
  geom_point(fill="darkolivegreen", colour="black", 
             shape=23, size = 2.5)+
  xlab(expression("Depth (cm)")) +
  ylab(expression(CH[4]~Flux~"("*mg~C~m^{-2}~day^{-1}*")")) +
  labs(fill = NULL) +
  theme_classic(base_size = 15) +
  theme(legend.position = "bottom", 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
# pm_depth_sig_plot

### DIFF/PM DEPTH (ZA)
diff_depth_sig_plotv2 <- no_clean_data %>%
  filter(depth > 25) %>%
  ggplot(aes(x=depth, y=cch4_flux_mgm2d1, 
             colour=veg_class, fill = veg_class)) + 
  # geom_errorbar(aes(ymin=cch4_flux_mgm2d1_mean-cch4_flux_mgm2d1_se, 
  #                   ymax=cch4_flux_mgm2d1_mean+cch4_flux_mgm2d1_se),
  #               alpha = 0.6) +
  # geom_errorbarh(aes(xmin=depth_mean-depth_se, 
  #                    xmax=depth_mean+depth_se), width = 0,
  #                alpha = 0.6) +
  # geom_jitter(data = depth_modelling,
  #            aes(x = depth, y = cch4_flux_mgm2d1, fill = veg_class),
  #            alpha = 0.3, shape = 21) +
  geom_hline(yintercept = 0, alpha = 0.6, linetype = "dashed") +
  #geom_smooth(aes(x = depth,
  #                y = (exp(log_cch4_flux_mgm2d1) - 1),
  #                fill=veg_class), method = glm, alpha=0.3, show.legend = F) +
  scale_fill_manual(values = c("Open"="grey80",
                               "Vegetated" = "grey80") )+
  new_scale_fill() + 
  geom_point(aes(fill=veg_class), colour="black", shape=21, size = 2.5)+
  scale_fill_manual(values = c("Open"="dodgerblue2",
                               "Vegetated"="olivedrab3"))+
  scale_color_manual(values = c("Open"="grey70",
                                "Vegetated"="grey70"))+
  xlab(expression("Depth (cm)")) +
  ylab(expression(atop(Chamber-Based~CH[4]~Flux,
                       "("*mg~C~m^{-2}~day^{-1}*")"))) +
  labs(fill = NULL) +
  theme_classic(base_size = 15) +
  theme(legend.position = "bottom", 
        axis.title.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
# diff_depth_sig_plot

### EBULL DEPTH (ZA)
ebull_depth_sig_plotv2 <- bubdata %>%
  ggplot(aes(x=depth_bt, 
             y=ebull_flux_ch4, 
             colour=veg_class, fill = veg_class)) + 
  # geom_errorbar(aes(ymin=ebull_flux_ch4_mean-ebull_flux_ch4_se, 
  #                   ymax=ebull_flux_ch4_mean+ebull_flux_ch4_se),
  #               alpha = 0.7) +
  # geom_errorbarh(aes(xmin=depth_bt_mean-depth_bt_se, 
  #                    xmax=depth_bt_mean+depth_bt_se), width = 0,
  #                alpha = 0.7) +
  #geom_jitter(data = bubdata,
  #            aes(x = depth_bt, y = ebull_flux_ch4, fill = veg_class),
  #            alpha = 0.3, shape = 21, size = 1.5) +
  geom_hline(yintercept = 0, alpha = 0.6, linetype = "dashed") +
  geom_smooth(aes(x = depth_bt,
                  y = (exp(log_ebull_flux_ch4) - 1),
                  fill=veg_class), method = glm, alpha=0.3, show.legend = F) +
  scale_fill_manual(values = c("Open"="dodgerblue3",
                               "Vegetated" = "transparent") )+
  new_scale_fill() + 
  geom_point(aes(fill=veg_class), 
             colour="black", shape=21, size = 2.5)+
  scale_fill_manual(values = c("Open"="dodgerblue2",
                               "Vegetated"="olivedrab3"))+
  scale_color_manual(values = c("Open"="dodgerblue3",
                                "Vegetated"="transparent")) +
  xlab(expression("Depth (cm)")) +
  ylab(expression(CH[4]~Flux~"("*mg~C~m^{-2}~day^{-1}*")")) +
  labs(fill = NULL) +
  theme_classic(base_size = 15) +
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
# ebull_depth_sig_plot

depthv2 <- pm_depth_sig_plotv2 | diff_depth_sig_plotv2 | ebull_depth_sig_plotv2

## sed temp ----
pm_sed_sig_plotv2 <- no_veg_ch4 %>%
  ggplot(aes(x=temp_sed_3day, 
             y=plant_cch4_flux_mgm2d1, 
             colour=veg_class, fill = veg_class)) + 
  # geom_errorbar(aes(ymin=plant_cch4_flux_mgm2d1_mean-plant_cch4_flux_mgm2d1_se, 
  #                   ymax=plant_cch4_flux_mgm2d1_mean+plant_cch4_flux_mgm2d1_se),
  #               width = 0, alpha = 0.6) +
  # geom_errorbarh(aes(xmin=temp_sed_3day_mean-temp_sed_3day_se, 
  #                    xmax=temp_sed_3day_mean+temp_sed_3day_se), 
  #                width = 0, alpha = 0.6) +
  # geom_jitter(data = depth_modelling,
  #            aes(x = depth, y = cch4_flux_mgm2d1, fill = veg_class),
  #            alpha = 0.3, shape = 21) +
  geom_hline(yintercept = 0, alpha = 0.6, linetype = "dashed") +
  #geom_smooth(aes(x = temp_sed_3day,
  #                y = (exp(log_plant_cch4_flux_mgm2d1) - 1),
  #                fill=veg_class), 
  #            method = glm, alpha=0.3, show.legend = F)+
  scale_fill_manual(values = c("Open"="grey80",
                               "Vegetated" = "grey80") )+
  new_scale_fill() + 
  geom_point(aes(fill=veg_class), colour="black", shape=23, size = 2.5)+
  scale_fill_manual(values = c("Open"="dodgerblue",
                               "Vegetated"="darkolivegreen"))+
  scale_color_manual(values = c("Open"="grey70",
                                "Vegetated"="grey70"))+
  xlab("Sediment Temperature (°C)") +
  ylab(expression(CH[4]~Flux~"("*mg~C~m^{-2}~day^{-1}*")")) +
  labs(fill = NULL) +
  theme_classic(base_size = 15) +
  theme(legend.position = "bottom", 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
# pm_sed_sig_plot

diff_sed_sig_plotv2 <- no_clean_data %>%
  ggplot(aes(x=temp_sed_3day, 
             y=cch4_flux_mgm2d1, 
             colour=veg_class, fill = veg_class)) + 
  # geom_errorbar(aes(ymin=cch4_flux_mgm2d1_mean-cch4_flux_mgm2d1_se, 
  #                   ymax=cch4_flux_mgm2d1_mean+cch4_flux_mgm2d1_se),
  #               alpha = 0.6) +
  # geom_errorbarh(aes(xmin=temp_sed_3day_mean-temp_sed_3day_se, 
  #                    xmax=temp_sed_3day_mean+temp_sed_3day_se), width = 0,
  #                alpha = 0.6) +
  # geom_jitter(data = depth_modelling,
  #            aes(x = depth, y = cch4_flux_mgm2d1, fill = veg_class),
  #            alpha = 0.3, shape = 21) +
  geom_hline(yintercept = 0, alpha = 0.6, linetype = "dashed") +
  #geom_smooth(aes(x = temp_sed_3day,
  #                y = exp(log_cch4_flux_mgm2d1) - 1,
  #                fill=veg_class), method = glm, alpha=0.3, show.legend = F)+
  scale_fill_manual(values = c("Open"="grey80",
                               "Vegetated" = "grey80") )+
  new_scale_fill() + 
  geom_point(aes(fill=veg_class), colour="black", shape=21, size = 2.5)+
  scale_fill_manual(values = c("Open"="dodgerblue2",
                               "Vegetated"="olivedrab3"))+
  scale_color_manual(values = c("Open"="grey70",
                                "Vegetated"="grey70"))+
  xlab("Sediment Temperature (°C)") +
  ylab(expression(CH[4]~Flux~"("*mg~C~m^{-2}~day^{-1}*")")) +
  labs(fill = NULL) +
  theme_classic(base_size = 15) +
  theme(legend.position = "bottom", 
        axis.title.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
# diff_sed_sig_plot

ebull_sed_sig_plotv2 <- no_clean_data %>%
  ggplot(aes(x=temp_sed_3day,
             y=ebull_flux_ch4, 
             colour=veg_class, fill = veg_class)) + 
  # geom_errorbar(aes(ymin=ebull_flux_ch4_mean-ebull_flux_ch4_se, 
  #                   ymax=ebull_flux_ch4_mean+ebull_flux_ch4_se),
  #               width = 0, alpha = 0.6) +
  # geom_errorbarh(aes(xmin=temp_sed_3day_mean-temp_sed_3day_se, 
  #                    xmax=temp_sed_3day_mean+temp_sed_3day_se), 
  #                width = 0, alpha = 0.6) +
  # geom_jitter(data = depth_modelling,
  #            aes(x = depth, y = cch4_flux_mgm2d1, fill = veg_class),
  #            alpha = 0.3, shape = 21) +
  geom_hline(yintercept = 0, alpha = 0.6, linetype = "dashed") +
  #geom_smooth(aes(x = temp_sed_3day,
  #                y = exp(log_ebull_flux_ch4) - 1,
  #                fill=veg_class), method = glm, alpha=0.3, show.legend = F)+
  scale_fill_manual(values = c("Open"="grey80",
                               "Vegetated" = "grey80") )+
  new_scale_fill() + 
  geom_point(aes(fill=veg_class), colour="black", shape=21, size = 2.5)+
  scale_fill_manual(values = c("Open"="dodgerblue2",
                               "Vegetated"="olivedrab3"))+
  scale_color_manual(values = c("Open"="white",
                                "Vegetated"="white"))+
  xlab("Sediment Temperature (°C)") +
  ylab(expression(CH[4]~Flux~"("*mg~C~m^{-2}~day^{-1}*")")) +
  labs(fill = NULL) +
  theme_classic(base_size = 15) +
  theme(legend.position = "bottom", 
        axis.title.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
# ebull_sed_sig_plot

sedv2 <- pm_sed_sig_plotv2 | diff_sed_sig_plotv2 | ebull_sed_sig_plotv2

## biomass ----
pm_biov2 <- no_veg_ch4 %>%
  filter(biomass_g < 200) %>%
  ggplot(aes(x=biomass_g, 
             y=plant_cch4_flux_mgm2d1, 
             colour=veg_class, fill = veg_class)) + 
  # geom_errorbar(aes(ymin=plant_cch4_flux_mgm2d1_mean-plant_cch4_flux_mgm2d1_se, 
  #                   ymax=plant_cch4_flux_mgm2d1_mean+plant_cch4_flux_mgm2d1_se),
  #               width = 0, alpha = 0.6) +
  # geom_errorbarh(aes(xmin=biomass_g_mean-biomass_g_se, 
  #                    xmax=biomass_g_mean+biomass_g_se), 
  #                width = 0, alpha = 0.6) +
  # geom_jitter(data = depth_modelling,
  #            aes(x = depth, y = cch4_flux_mgm2d1, fill = veg_class),
  #            alpha = 0.3, shape = 21) +
  geom_hline(yintercept = 0, alpha = 0.6, linetype = "dashed") +
  #geom_smooth(aes(x = biomass_g,
  #                y = exp(log_plant_cch4_flux_mgm2d1) - 1,
  #                fill=veg_class), method = glm, alpha=0.3, show.legend = F)+
  scale_fill_manual(values = c("Open"="grey80",
                               "Vegetated" = "grey80") )+
  new_scale_fill() + 
  geom_point(aes(fill=veg_class), colour="black", shape=23, size = 2.5)+
  scale_fill_manual(values = c("Open"="dodgerblue2",
                               "Vegetated"="darkolivegreen"))+
  scale_color_manual(values = c("Open"="grey70",
                                "Vegetated"="grey70"))+
  xlab(expression(Biomass~"("*g*")")) +
  ylab(expression(CH[4]~Flux~"("*mg~C~m^{-2}~day^{-1}*")")) +
  labs(fill = NULL) +
  theme_classic(base_size = 15) +
  theme(legend.position = "None", 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
# pm_bio

## aw vol ----
pm_abwvolv2 <- no_veg_ch4 %>%
  filter(biomass_g < 200) %>%
  ggplot(aes(x=abovewater_vol_cm3, 
             y=plant_cch4_flux_mgm2d1, 
             colour=veg_class, fill = veg_class)) + 
  # geom_errorbar(aes(ymin=plant_cch4_flux_mgm2d1_mean-plant_cch4_flux_mgm2d1_se,
  #                   ymax=plant_cch4_flux_mgm2d1_mean+plant_cch4_flux_mgm2d1_se),
  #               width = 0, alpha = 0.6) +
  # geom_errorbarh(aes(xmin=biomass_g_mean-biomass_g_se,
  #                    xmax=biomass_g_mean+biomass_g_se),
  #                width = 0, alpha = 0.6) +
  # geom_jitter(data = depth_modelling,
  #            aes(x = depth, y = cch4_flux_mgm2d1, fill = veg_class),
  #            alpha = 0.3, shape = 21) +
  geom_hline(yintercept = 0, alpha = 0.6, linetype = "dashed") +
  #geom_smooth(aes(x = abovewater_vol_cm3,
  #                y = exp(log_plant_cch4_flux_mgm2d1) - 1,
  #                fill=veg_class), method = glm, alpha=0.3, show.legend = F)+
  scale_fill_manual(values = c("Open"="grey80",
                               "Vegetated" = "grey70") )+
  new_scale_fill() + 
  geom_point(aes(fill=veg_class), colour="black", shape=23, size = 2.5)+
  scale_fill_manual(values = c("Open"="dodgerblue2",
                               "Vegetated"="darkolivegreen"))+
  scale_color_manual(values = c("Open"="grey70",
                                "Vegetated"="darkolivegreen"))+
  xlab(expression(Above~Water~Plant~Volume~"("*cm^{3}*")")) +
  ylab(expression(CH[4]~Flux~"("*mg~C~m^{-2}~day^{-1}*")")) +
  labs(fill = NULL) +
  theme_classic(base_size = 15) +
  theme(legend.position = "None", 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
# pm_abwvol

## GCC ----
pm_gcc_plotv2 <- no_veg_ch4 %>%
  filter(biomass_g < 200) %>%
  ggplot(aes(x=GCC,
             y=plant_cch4_flux_mgm2d1, 
             colour=veg_class, fill = veg_class)) + 
  # geom_errorbar(aes(ymin=plant_cch4_flux_mgm2d1_mean-plant_cch4_flux_mgm2d1_se, 
  #                   ymax=plant_cch4_flux_mgm2d1_mean+plant_cch4_flux_mgm2d1_se),
  #               width = 0, alpha = 0.6) +
  # geom_errorbarh(aes(xmin=GCC_mean-GCC_se, 
  #                    xmax=GCC_mean+GCC_se), 
  #                width = 0, alpha = 0.6) +
  # geom_jitter(data = depth_modelling,
  #            aes(x = depth, y = cch4_flux_mgm2d1, fill = veg_class),
  #            alpha = 0.3, shape = 21) +
  geom_hline(yintercept = 0, alpha = 0.6, linetype = "dashed") +
  #geom_smooth(aes(x = GCC,
  #                y = exp(log_plant_cch4_flux_mgm2d1) - 1,
  #                fill=veg_class), method = lm, alpha=0.3, show.legend = F)+
  scale_fill_manual(values = c("Open"="grey80",
                               "Vegetated" = "grey80") )+
  new_scale_fill() + 
  geom_point(aes(fill=veg_class), colour="black", shape=23, size = 2.5)+
  scale_fill_manual(values = c("Open"="dodgerblue2",
                               "Vegetated"="darkolivegreen"))+
  scale_color_manual(values = c("Open"="grey70",
                                "Vegetated"="grey70"))+
  xlab("Greenness Chromatic Coordinate") +
  ylab(expression(CH[4]~Flux~"("*mg~C~m^{-2}~day^{-1}*")")) +
  labs(fill = NULL) +
  theme_classic(base_size = 15) +
  theme(legend.position = "None", 
        axis.title.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
# pm_gcc_plot

## perc cover ----
pm_pc_plotv2 <- no_veg_ch4 %>%
  filter(biomass_g < 200) %>%
  ggplot(aes(x=perc_cover, 
             y=plant_cch4_flux_mgm2d1, 
             colour=veg_class, fill = veg_class)) + 
  # geom_errorbar(aes(ymin=plant_cch4_flux_mgm2d1_mean-plant_cch4_flux_mgm2d1_se, 
  #                   ymax=plant_cch4_flux_mgm2d1_mean+plant_cch4_flux_mgm2d1_se),
  #               width = 0, alpha = 0.6) +
  # geom_errorbarh(aes(xmin=perc_cover_mean-perc_cover_se, 
  #                    xmax=perc_cover_mean+perc_cover_se), 
  #                width = 0, alpha = 0.6) +
  # geom_jitter(data = depth_modelling,
  #            aes(x = depth, y = cch4_flux_mgm2d1, fill = veg_class),
  #            alpha = 0.3, shape = 21) +
  geom_hline(yintercept = 0, alpha = 0.6, linetype = "dashed") +
  #geom_smooth(aes(x = perc_cover,
  #                y = exp(log_plant_cch4_flux_mgm2d1) - 1,
  #                fill=veg_class), method = lm, alpha=0.3, show.legend = F)+
  scale_fill_manual(values = c("Open"="grey80",
                               "Vegetated" = "grey80") )+
  new_scale_fill() + 
  geom_point(aes(fill=veg_class), colour="black", shape=23, size = 2.5)+
  scale_fill_manual(values = c("Open"="dodgerblue2",
                               "Vegetated"="darkolivegreen"),
                    labels = c("Open" = "Open",
                               "Vegetated" = "Plant-Mediated"))+
  scale_color_manual(values = c("Open"="grey70",
                                "Vegetated"="grey70")) +
  xlab("Percent Cover (%)") +
  ylab(expression(CH[4]~Flux~"("*mg~C~m^{-2}~day^{-1}*")")) +
  labs(fill = NULL) +
  theme_classic(base_size = 15) +
  theme(legend.position = "None", 
        axis.title.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
# pm_pc_plot


## put final plot together ----
fig5v2_glmms <- depthv2 / sedv2 / (pm_biov2 | pm_pc_plotv2 | pm_gcc_plotv2)
fig5v2_glmms + plot_annotation(tag_levels = 'a',
                               tag_prefix = '(',
                               tag_suffix = ')') +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# FIGURE 6: PRODUCTION ----
## daily log production BAR plot ----
dailych4prod <- DailyNoOutIncProd %>%
  filter(treatment != "Control") %>%
  ggplot(aes(x = time_step,
             y = avg_production_dry,
             fill = factor(treatment),
             group = treatment)) +
  geom_bar(stat="identity",
           position = position_dodge(width = 0.9),
           alpha = 0.85) +
  geom_errorbar(aes(ymin = avg_production_dry - se_dry,
                    ymax = avg_production_dry + se_dry),
                width = 0.2,
                position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("Water + Benthic" = "olivedrab3",
                               "Water + Sediment" = "indianred3",
                               "Water + Sediment + Benthic" = "goldenrod"),
                    name = "Treatment") +
  xlab("Incubation Day") +
  ylab(expression("Potential CH"[4]*" Production ("*nmol~CH[4]~gdw^{-1}~d^{-1}*")")) +
  scale_y_continuous(transform = "pseudo_log",
                     limits = c(0, 300)) +
  facet_wrap(~ t_c, nrow = 2) +
  theme_minimal(base_size = 15) +
  theme(legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
# dailych4prod

## Q10 ----
AvgQ10Prod <- Q10IncProd %>%
  group_by(treatment) %>%
  summarise(
    mean_q10 = mean(q10_dry, na.rm = TRUE),
    sd_q10   = sd(q10_dry, na.rm = TRUE),
    se_q10 = sd_q10 / sqrt(sum(!is.na(q10_dry)))
  )
  
Q10_plot <- AvgQ10Prod %>%
  filter(treatment != "Control") %>%
  ggplot( aes(x = treatment, y = mean_q10,
              group = treatment)) +
  geom_point( aes(colour = treatment),
              size = 3, 
              position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = mean_q10 - se_q10,
                    ymax = mean_q10 + se_q10,
                    colour = treatment),
                width = 0.3, 
                position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = c("Water + Benthic" = "olivedrab3",
                                 "Water + Sediment" = "indianred3",
                                 "Water + Sediment + Benthic" = "goldenrod"),
                      name = "Treatment",
                      labels = c("Sediment Only",
                                 "Sediment & Benthic",
                                 "Benthic Only")) +
  ylab(expression("Average Q"[10]*"")) +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
# Q10_plot

## put final plot together ----
fig6_incubation <- dailych4prod | Q10_plot + 
  plot_layout(widths = c(3,1))
fig6_incubation + 
  plot_annotation(tag_levels = 'a',
                  tag_prefix = '(',
                  tag_suffix = ')') +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# FIGURE 7: OM + C:N ----
## CUMULATIVE + OM ----
cumulative_om <- CleanIncNutrient %>%
  ggplot( aes(x = avg_sed_om, y = avg_cumulative_prod,
              group = treatment_short)) + 
  geom_errorbar(aes(ymin=avg_cumulative_prod-se_cumulative_prod, 
                    ymax=avg_cumulative_prod+se_cumulative_prod)) +
  geom_errorbarh(aes(xmin=avg_sed_om-se_sed_om, 
                     xmax=avg_sed_om+se_sed_om), width = 0) +
  geom_point( aes(fill = treatment_short, shape = t_c), 
              size = 5, colour = "black") +
  # geom_smooth(aes(group=t_c), method = lm, alpha=0.3, se = F)+
  scale_shape_manual(name = "Temperature",
                     values = c("Average Field Temperature (20°C)" = 21,
                                "Potential Warming Temperature (30°C)" = 23),
                     labels = c("20°C",
                                "30°C")) +
  scale_fill_manual(name = "Treatments",
                    values = c("WB" = "olivedrab3",
                               "WS" = "indianred3",
                               "WSB" = "goldenrod"),
                    labels = c("Benthic Only", 
                               "Sediment Only", 
                               "Sediment & Benthic")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  xlab("Organic Matter Content (%)") +
  ylab(expression("Cumulative Potential CH"[4]*" Production ("*nmol~CH[4]~gdw^{-1}*")")) +
  theme_classic(base_size = 15) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
# cumulative_om

cumulative_cn <- CleanIncNutrient %>%
  ggplot( aes(x = avg_sed_cn, y = avg_cumulative_prod,
              group = treatment_short)) + 
  geom_errorbar(aes(ymin=avg_cumulative_prod-se_cumulative_prod, 
                    ymax=avg_cumulative_prod+se_cumulative_prod)) +
  geom_errorbarh(aes(xmin=avg_sed_cn-se_sed_cn, 
                     xmax=avg_sed_cn+se_sed_cn), width = 0) +
  geom_point( aes(fill = treatment_short, shape = t_c), 
              size = 5, colour = "black") +
  # geom_smooth(aes(group=t_c), method = lm, alpha=0.3, se = F)+
  scale_shape_manual(name = "Temperature",
                     values = c("Average Field Temperature (20°C)" = 21,
                                "Potential Warming Temperature (30°C)" = 23),
                     labels = c("20°C",
                                "30°C")) +
  scale_fill_manual(name = "Treatments",
                    values = c("WB" = "olivedrab3",
                               "WS" = "indianred3",
                               "WSB" = "goldenrod"),
                    labels = c("Benthic Only", 
                               "Sediment Only", 
                               "Sediment & Benthic")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  xlab("C:N") +
  ylab(expression("Cumulative Potential CH"[4]*" Production ("*nmol~CH[4]~gdw^{-1}*")")) +
  theme_classic(base_size = 15) +
  theme(legend.position = "right",
        axis.title.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
# cumulative_cn

## put final plot together ----
fig7_omcn <- cumulative_om | cumulative_cn 
fig7_omcn + plot_annotation(tag_levels = 'a',
                            tag_prefix = '(',
                            tag_suffix = ')')

# FIGURE 8: STUDY COMPS ----
open_observations <- observations %>%
  filter(veg_class == "Open")

open_data <- multistudy %>%
  filter(veg_class == "Open") %>%
  ggplot(aes(x = dataset, y = cch4_flux_mgm2d1)) +
  geom_jitter(aes(colour = ecosystem), alpha = 0.5,
              width = 0.35) +
  geom_boxplot(aes(fill = ecosystem), alpha = 0.7,
               outlier.shape = NA) +
  geom_text(data = open_observations, aes(x = dataset, y = Inf, label = n),
            inherit.aes = FALSE,
            vjust = 3, hjust = -0.5, size = 4) +
  scale_colour_manual(name = "Ecosystem",
                      values = c("Temperate" = 'steelblue2',
                                 "Arctic-boreal" = 'steelblue4',
                                 "Northern Boreal" = 'darkgoldenrod2')) +
  scale_fill_manual(name = "Ecosystem",
                    values = c("Temperate" = 'steelblue2',
                               "Arctic-boreal" = 'steelblue4',
                               "Northern Boreal" = 'darkgoldenrod2')) +
  xlab("Aquatic Carbon Dataset") +
  ylab(expression(CH[4]~Flux~"("*mg~C~m^{-2}~day^{-1}*")")) +
  # ylim(c(0, 1000)) +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "right") +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
# open_data

veg_observations <- observations %>%
  filter(veg_class == "Vegetated" &
           dataset != "BAWLD-CH4 (2021)")

veg_data <- multistudy %>%
  filter(veg_class == "Vegetated" & 
         dataset != "BAWLD-CH4 (2021)") %>%
  ggplot(aes(x = dataset, y = cch4_flux_mgm2d1)) +
  geom_jitter(aes(colour = ecosystem), alpha = 0.5,
              width = 0.35) +
  geom_boxplot(aes(fill = ecosystem), alpha = 0.7,
               outlier.shape = NA) +
  geom_text(data = veg_observations, aes(x = dataset, y = Inf, label = n),
            inherit.aes = FALSE,
            vjust = 3, hjust = -0.5, size = 4) +
  scale_colour_manual(name = "Ecosystem",
                      values = c("Temperate" = 'steelblue2',
                                 "Arctic-boreal" = 'steelblue4',
                                 "Northern Boreal" = 'darkgoldenrod2')) +
  scale_fill_manual(name = "Ecosystem",
                    values = c("Temperate" = 'steelblue2',
                               "Arctic-boreal" = 'steelblue4',
                               "Northern Boreal" = 'darkgoldenrod2')) +
  xlab("Aquatic Carbon Dataset") +
  ylab(expression(CH[4]~Flux~"("*mg~C~m^{-2}~day^{-1}*")")) +
  ylim(c(0, 1000)) +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title.y = element_blank(),
        legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
# veg_data

## put final plot together ----
fig8_studycomps <- open_data | veg_data
fig8_studycomps + plot_annotation(tag_levels = 'a',
                                  tag_prefix = '(',
                                  tag_suffix = ')') +
  plot_layout(guides = "collect")

# APPENDIX ----
## MODEL OUTPUTS ----
### bot water temp ----
bot_water <- model_plotting("temp_bottom_daily_", "Bottom Water Temperature (°C)")
# bot_water

### surf water temp ----
surf_water <- model_plotting("temp_surface_daily_", "Surface Water Temperature (°C)")
# surf_water

### air temp ----
air <- model_plotting("temp_air_daily_", "Air Temperature (°C)")
# air

### DO ----
do <- model_plotting("DO_", "Dissolved Oxygen (%)")
# do

### TDS ----
tds <- model_plotting("tds_mgL_", expression(Total~Dissolved~Solids~"("*mg~L^{-1}*")"))
# tds

### cond ----
cond <- model_plotting("cond_mScm_", expression(Conductivity~"("*mS~cm^{-1}*")"))
# cond

### pH ----
pH <- model_plotting("ph_", "pH")
# pH

### wind speed ----
wind_speed <- model_plotting("HOURLY_WIND_SPEED_", expression(Hourly~Wind~Speed~"("*m~s^{-2}*")"))
# wind_speed

## ALL LM PLOTS TOGETHER ----
extraglmms <- do / pH / wind_speed / air / surf_water / bot_water
extraglmms + plot_annotation(tag_levels = "a",
                             tag_prefix = "(",
                             tag_suffix = ")")

## Stack ebull, ch4 content, bubble vol. ----
ebull_stack <- ggplot(bubdata, 
                      aes(x = as.factor(sample_date), y = ebull_flux_ch4, 
                          fill = veg_class)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Open" = 'dodgerblue3',
                               "Vegetated" = 'olivedrab3')) +
  ylab(expression(Ebullitive~CH[4]~Flux~"("*mg~C~m^{-2}~day^{-1}*")")) +
  #  ggtitle(expression("Wheeler Lake Total"~CH[4]~Flux)) +
  theme_classic(base_size = 15) +
  theme(legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
# ebull_stack

percch4_stack <- ggplot(bubdata, 
                      aes(x = as.factor(sample_date), y = perc_ch4, 
                          fill = veg_class)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Open" = 'dodgerblue3',
                               "Vegetated" = 'olivedrab3')) +
  ylab(expression(Percent~CH[4]~"(%)")) +
  #  ggtitle(expression("Wheeler Lake Total"~CH[4]~Flux)) +
  theme_classic(base_size = 15) +
  theme(legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
# percch4_stack

bub_stack <- ggplot(bubdata, 
                        aes(x = as.factor(sample_date), y = bub_flux_lm2d1, 
                            fill = veg_class)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Open" = 'dodgerblue3',
                               "Vegetated" = 'olivedrab3')) +
  ylab(expression(Bubble~Flux~"("*L~m^{-2}~day^{-1}*")")) +
  xlab("Sample Date") + 
  theme_classic(base_size = 15) +
  theme(legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        axis.text.x = element_text(angle = 30, hjust = 1))
# bub_stack

stack <- ebull_stack / percch4_stack / bub_stack
stack +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")") +
  plot_layout(guides = "collect") &
  theme(legend.title = element_blank())

## flux time series ----
DiffTime <-no_clean_data %>%
  ggplot(aes(x = as.Date(sample_date), y = cch4_flux_mgm2d1,
             fill = veg_class)) +
  geom_smooth(aes(fill = veg_class, colour = veg_class), 
              method = "lm",
              se = TRUE, alpha = 0.2) +  # Match SE color to line color
  scale_fill_manual(values = c("Open"="grey80",
                               "Vegetated" = "grey80")) +
  new_scale_fill() + 
  geom_point(aes(fill=veg_class), colour="black", shape=21, size = 2.5)+
  scale_fill_manual(values = c("Open"="dodgerblue2",
                               "Vegetated"="olivedrab3")) +
  scale_color_manual(values = c("Open"="grey70",
                                "Vegetated"="grey70")) +
  coord_cartesian(ylim = c(0,300)) +
  theme_classic(base_size = 15) +
  ylab(expression(Chamber-Based~Flux~"("*mg~m^{-2}~day^{-1}*")")) +
  theme(
    legend.position = "None",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )
# DiffTime

EbullTime <- bubdata %>%
  ggplot(aes(x = as.Date(sample_date), y = ebull_flux_ch4,
             fill = veg_class)) +
  geom_smooth(aes(fill = veg_class, colour = veg_class), 
              method = "lm",
              se = TRUE, alpha = 0.2) +  # Match SE color to line color
  scale_fill_manual(values = c("Open"="grey80",
                               "Vegetated" = "grey80") )+
  new_scale_fill() + 
  geom_point(aes(fill=veg_class), colour="black", shape=21, size = 2.5)+
  scale_fill_manual(values = c("Open"="dodgerblue2",
                               "Vegetated"="olivedrab3"))+
  scale_color_manual(values = c("Open"="grey70",
                                "Vegetated"="grey70"))+
  coord_cartesian(ylim = c(0,300)) +
  theme_classic(base_size = 15) +
  xlab("Sample Date") +
  ylab(expression(Ebullitive~Flux~"("*mg~m^{-2}~day^{-1}*")")) +
  theme(
    legend.position = "None",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )
# EbullTime

DiffTime / EbullTime +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")")

## flux time series ----
GCCTime <-no_clean_data %>%
  ggplot(aes(x = as.Date(sample_date), y = GCC,
             fill = veg_class)) +
  geom_smooth(fill = "darkolivegreen", colour = "darkolivegreen", 
              method = "lm",
              se = TRUE, alpha = 0.2) +
  new_scale_fill() + 
  geom_point(fill="darkolivegreen", colour="black", shape=23, size = 2.5)+
  theme_classic(base_size = 15) +
  ylab("Green Chromatic Coordinate") +
  xlab("Date") +
  theme(
    legend.position = "None",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
GCCTime

## temperature ----
hobo_datasets <- CleanHOBOs(hobo_folderpath = '~/Desktop/masters/data/msc-aquatic-ch4/Wheeler_HOBO_data/hobo_csvs/',
                            hobo_rangemetafile = '~/Desktop/masters/data/msc-aquatic-ch4/Wheeler_HOBO_data/hobo_whe_2025.csv')

# separate out the hobo_data output into their respective dataframes
raw_hobo <- hobo_datasets$raw_hobo
avgtemps_byplacement <- hobo_datasets$byplacement

# plot a simple time series showing the raw temperatures by HOBO
tempr_plot <- raw_hobo %>%
  filter(!placement %in% c("AIR", "BOT", "SUR", "SUR_VEG") &
           temp_c < 29) %>%
  ggplot() +
  geom_line(aes(datetime, temp_c, colour = placement),
            alpha = 0.7) +
  xlab("Date") +
  ylab("Temperature (°C)") +
  labs(color = "Hobo Name") +
  scale_colour_manual(name = "HOBO Location",
                      labels = c("AIR" = "Air", "BOT" = "Bottom Water",
                                 "SED_OPEN" = "Sediment (Open)", "SED_VEG" = "Sediment (Vegetated)",
                                 "SUR" = "Surface Water (Open)", "SUR_VEG" = "Surface Water (Vegetated)"),
                      values = c("AIR" = "grey70", "BOT" = "dodgerblue4", 
                                 "SED_OPEN" = "sienna2", "SED_VEG" = "sienna",
                                 "SUR" = "steelblue3", "SUR_VEG" = "steelblue2")) +
  theme_classic(base_size = 15)
tempr_plot

# plot a simple time series showing the daily average temperature by placement
temp1_plot <- avgtemps_byplacement %>%
  ggplot() +
  geom_line(aes(date, temp_c_mean, colour = placement)) +
  scale_colour_manual(name = "HOBO Location",
                      labels = c("AIR" = "Air", "BOT" = "Bottom Water",
                                 "SED_OPEN" = "Sediment (Open)", "SED_VEG" = "Sediment (Vegetated)",
                                 "SUR" = "Surface Water (Open)", "SUR_VEG" = "Surface Water (Vegetated)"),
                      values = c("AIR" = "grey70", "BOT" = "dodgerblue4", 
                                 "SED_OPEN" = "sienna2", "SED_VEG" = "sienna", 
                                 "SUR" = "steelblue3", "SUR_VEG" = "steelblue2")) +
  scale_x_date(date_breaks = "1 week", date_labels = "%b %d") +
  xlab("Date") +
  ylab("Temperature (°C)") +
  labs(color = "Placement") +
  theme_classic(base_size = 15) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())
# temp1_plot

# plot a simple time series showing the three day average temperature by placement
temp3_plot <- avgtemps_byplacement %>%
  ggplot() +
  geom_line(aes(date, threeday_avg, colour = placement)) +
  scale_colour_manual(name = "HOBO Location",
                      labels = c("AIR" = "Air", "BOT" = "Bottom Water",
                                 "SED_OPEN" = "Sediment (Open)", "SED_VEG" = "Sediment (Vegetated)",
                                 "SUR" = "Surface Water (Open)", "SUR_VEG" = "Surface Water (Vegetated)"),
                      values = c("AIR" = "grey70", "BOT" = "dodgerblue4", 
                                 "SED_OPEN" = "sienna2", "SED_VEG" = "sienna", 
                                 "SUR" = "steelblue3", "SUR_VEG" = "steelblue2")) +
  scale_x_date(date_breaks = "1 week", date_labels = "%b %d") +
  xlab("Date") +
  ylab("Temperature (°C)") +
  labs(color = "Placement") +
  theme_classic(base_size = 15) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())
# temp3_plot

# plot a simple time series showing the seven day average temperature by placement
temp7_plot<- avgtemps_byplacement %>%
  ggplot() +
  geom_line(aes(date, sevenday_avg, colour = placement)) +
  scale_colour_manual(name = "HOBO Location",
                      labels = c("AIR" = "Air", "BOT" = "Bottom Water",
                                 "SED_OPEN" = "Sediment (Open)", "SED_VEG" = "Sediment (Vegetated)",
                                 "SUR" = "Surface Water (Open)", "SUR_VEG" = "Surface Water (Vegetated)"),
                      values = c("AIR" = "grey70", "BOT" = "dodgerblue4", 
                                 "SED_OPEN" = "sienna2", "SED_VEG" = "sienna", 
                                 "SUR" = "steelblue3", "SUR_VEG" = "steelblue2")) +
  scale_x_date(date_breaks = "1 week", date_labels = "%b %d") +
  xlab("Date") +
  ylab("Temperature (°C)") +
  labs(color = "Placement") +
  theme_classic(base_size = 15) 
# temp7_plot

# put temps all together
tempr_plot / temp1_plot / temp3_plot / temp7_plot +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")") 


## ebullitive flux histos ----
ehisto <- bubdata %>%
  ggplot() +
  geom_histogram(aes(x = ebull_flux_ch4), binwidth = 5,
                 fill = "orangered3", colour = "orangered4") + # see a very right skewed distribution
  xlab(expression(Ebullitive~CH[4]~Flux~"("*mg~C~m^{-2}~day^{-1}*")")) +
  ylab("Count") +
  theme_classic(base_size = 15)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))

logehisto <- bubdata %>%
  ggplot() +
  geom_histogram(aes(x = log(ebull_flux_ch4 + 1)), 
                 binwidth = 0.2,
                 fill = "orangered3", colour = "orangered4") + # see a very right skewed distribution
  xlab(expression(Log-Transformed~Ebullitive~CH[4]~Flux~"("*mg~C~m^{-2}~day^{-1}*")")) +
  ylab("Count") +
  theme_classic(base_size = 15)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))

histograms <- ehisto | logehisto
histograms + 
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")")
