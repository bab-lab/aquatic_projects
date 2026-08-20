# MSc Project: A little sedge goes a long way: emergent vegetation as a key driver of littoral methane flux

# Functions for running repetitive model setups -> lm; mlm; glmm, anovas

# 03/16/2026 (Last Edit: 08/17/2026)
# Author: Kelsey McGuire
# kmcgu@mail.ubc.ca; kmcguire.9@outlook.com

# load in regression functions ----
## simple linear regression function ----
runlm <- function(dataset, response, predictor) {
  lmlist <- list(model = list(), summary = list()) # list of all models
  sig_lms <- list(model = list(), summary = list()) # list of significant models
  response_var <- deparse(substitute(response)) # convert response variable into usable format
  
  for (i in 1:length(predictor)) {
    pred <- predictor[i] # get the specific predictor from the list of predictors
    
    formula <- as.formula(paste(response_var, "~", pred))
    lm_model <- lm(formula, data = dataset) #fit the LM model
    lm_sum <- summary(lm_model) # get the summary 
    
    lmlist$model[[pred]] <- lm_model # append the model
    lmlist$summary[[pred]] <- lm_sum # append the summary 
    
    coef <- lm_sum$coefficients
    print(coef)
    
    pval <- coef[2,4]
    
    if (!is.na(pval) && pval < 0.05) { # append only significant models in this list 
      sig_lms$model[[pred]] <- lm_model # append the model
      sig_lms$summary[[pred]] <- lm_sum # append the summary 
    }
  }
  
  models <- c(all_models = lmlist,
              sig_models = sig_lms)
  
  return(models)
}

## mixed linear regression function ----
runmlm <- function(dataset, response, predictor, random_effect) {
  lmlist <- list(model = list(), summary = list()) 
  sig_lms <- list(model = list(), summary = list())
  
  response_var <- deparse(substitute(response))
  
  for (pred in predictor) {
    # Build formula
    formula <- as.formula(
      paste0(response_var, " ~ ", pred, " + (1|", random_effect, ")")
    )
    
    # Filter out groups with <=1 obs
    df_filtered <- dataset %>%
      group_by(across(all_of(random_effect))) %>%
      filter(n() > 1) %>%
      ungroup()
    
    # Skip if less than 2 levels remain
    if (n_distinct(df_filtered[[random_effect]]) < 2) {
      message("Skipping predictor '", pred, "' because grouping factor has <2 levels after filtering")
      next
    }
    
    # Fit model
    lmm_model <- lmer(formula, data = df_filtered)
    lmm_sum <- summary(lmm_model)
    
    # Store all models
    lmlist$model[[pred]] <- lmm_model
    lmlist$summary[[pred]] <- lmm_sum
    
    # Extract p-value of the predictor
    pval <- lmm_sum$coefficients[2, "Pr(>|t|)"]
    
    # Store significant models
    if (!is.na(pval) && pval < 0.05) {
      sig_lms$model[[pred]] <- lmm_model
      sig_lms$summary[[pred]] <- lmm_sum
    }
  }
  
  list(all_models = lmlist, sig_models = sig_lms)
}

## generalized linear mixed model function ----
runglmm <- function(dataset, response, predictor, random_effect) {
  glmmlist <- list(model = list(), summary = list()) 
  
  response_var <- deparse(substitute(response))
  
  for (pred in predictor) {
    # Build formula
    formula <- as.formula(
      paste0(response_var, " ~ ", pred, " + zone + (1|", random_effect, ")")
    )
    
    # Fit model
    glmm_model <- glmmTMB(formula,
                          data = dataset,
                          family = Gamma(link = "log"))
    glmm_sum <- summary(glmm_model)
    
    # Store all models
    glmmlist$model[[pred]] <- glmm_model
    glmmlist$summary[[pred]] <- glmm_sum
  }
  
  glmmlist
}

## generalized linear mixed model function ----
runzaglmm <- function(dataset, response, predictor, random_effect) {
  glmmlist <- list(model = list(), summary = list()) 
  
  response_var <- deparse(substitute(response))
  
  for (pred in predictor) {
    # Build formula
    formula <- as.formula(
      paste0(response_var, " ~ ", pred, " + (1|", random_effect, ")")
    )
    
    # Fit model
    glmm_model <- glmmTMB(formula,
                          data = dataset,
                          family = Gamma(link = "log"))
    glmm_sum <- summary(glmm_model)
    
    # Store all models
    glmmlist$model[[pred]] <- glmm_model
    glmmlist$summary[[pred]] <- glmm_sum
  }
  
  glmmlist
}

## zero inflated generalized linear mixed model function ----
runziglm <- function(dataset, response, predictor, random_effect) {
  ziglmlist <- list(model = list(), summary = list()) 
  
  response_var <- deparse(substitute(response))
  
  for (pred in predictor) {
    # Build formula
    formula <- as.formula(
      paste0(response_var, " ~ ", pred, " + zone + (1|", random_effect, ")")
    )
    
    # Fit model
    ziglmm_model <- glmmTMB(formula,
                          data = dataset,
                          family = ziGamma(link = "log"), 
                          ziformula = ~ 1)
    ziglmm_sum <- summary(ziglmm_model)
    
    # Store all models
    ziglmlist$model[[pred]] <- ziglmm_model
    ziglmlist$summary[[pred]] <- ziglmm_sum
  }
  
  ziglmlist
}

## zero inflated generalized linear mixed model function ----
runzigzalm <- function(dataset, response, predictor, random_effect) {
  ziglmlist <- list(model = list(), summary = list()) 
  
  response_var <- deparse(substitute(response))
  
  for (pred in predictor) {
    # Build formula
    formula <- as.formula(
      paste0(response_var, " ~ ", pred, " + (1|", random_effect, ")")
    )
    
    # Fit model
    ziglmm_model <- glmmTMB(formula,
                            data = dataset,
                            family = ziGamma(link = "log"), 
                            ziformula = ~ 1)
    ziglmm_sum <- summary(ziglmm_model)
    
    # Store all models
    ziglmlist$model[[pred]] <- ziglmm_model
    ziglmlist$summary[[pred]] <- ziglmm_sum
  }
  
  ziglmlist
}

## model fitness plot function ----
modelfit <- function(model, xvar = NULL, main_label = "Model") {
  par(mfrow = c(2, 2), mai = c(0.6, 0.6, 0.6, 0.6), cex = 0.8)
  
  # Residual Plot
  plot(model$fitted.values,
       model$residuals,
       main = paste(main_label, "- Residual Plot"),
       xlab = "Fitted Values (ŷ)",
       ylab = "Residuals")
  abline(h = 0, col = "red")
  
  # Fitted Line Plot (requires x-variable)
  if (!is.null(xvar)) {
    plot(xvar,
         model$fitted.values,
         main = paste(main_label, "- Fitted Line Plot"),
         xlab = "x",
         ylab = "Fitted Values (ŷ)")
    
    abline(model, col = "red", lwd = 2)
  } else {
    plot(model$fitted.values,
         main = paste(main_label, "- Fitted Values"),
         ylab = "Fitted Values (ŷ)",
         xlab = "")
  }
  
  # Normal Q–Q Plot
  qqnorm(model$residuals,
         main = paste(main_label, "- Normal Q-Q Plot"))
  qqline(model$residuals, col = "red")
  
  # Histogram of Residuals
  hist(model$residuals,
       breaks = 8,
       density = 10,
       col = "green",
       border = "black",
       main = paste(main_label, "- Residual Distribution"))
  
  # Reset plotting parameters
  par(mfrow = c(1, 1), mai = c(1, 1, 1, 1), cex = 1)
}

## anova ----
runaov <- function(dataset, response, predictor) {
  aovlist <- list(model = list(), summary = list()) # list of all models
  sig_aov <- list(model = list(), summary = list()) # list of significant models
  response_var <- deparse(substitute(response)) # convert response variable into usable format
  
  for (i in 1:length(predictor)) {
    pred <- predictor[i] # get the specific predictor from the list of predictors
    
    formula <- as.formula(paste(response_var, "~", pred))
    aov_model <- aov(formula, data = dataset) # fit the AOV model
    aov_sum <- summary(aov_model) # get the summary 
    
    aovlist$model[[pred]] <- aov_model # append the model
    aovlist$summary[[pred]] <- aov_sum # append the summary 
    
    coef <- aov_sum$coefficients
    print(coef)
    
    pval <- coef[2,4]
    
    if (!is.na(pval) && pval < 0.05) { # append only significant models in this list 
      sig_aov$model[[pred]] <- aov_model # append the model
      sig_aov$summary[[pred]] <- aov_sum # append the summary 
    }
  }
  
  models <- c(all_models = aovlist,
              sig_models = sig_aov)
  
  return(models)
}

# model plotting for glms ----
model_plotting <- function(x_var, xlabel) {
  x_mean <- paste0(x_var, "mean")
  x_se <- paste0(x_var, "se")
  
  pm_plot <- 
    ggplot(no_veg_za, aes(x=.data[[x_mean]], 
                          y=plant_cch4_flux_mgm2d1_mean, 
                          colour=veg_class, fill = veg_class)) + 
    geom_errorbar(aes(ymin=plant_cch4_flux_mgm2d1_mean-plant_cch4_flux_mgm2d1_se, 
                      ymax=plant_cch4_flux_mgm2d1_mean+plant_cch4_flux_mgm2d1_se),
                  width = 0, alpha = 0.6) +
    geom_errorbarh(aes(xmin=.data[[x_mean]]-.data[[x_se]], 
                       xmax=.data[[x_mean]]+.data[[x_se]]), 
                   width = 0, alpha = 0.6) +
    # geom_jitter(data = depth_modelling,
    #            aes(x = depth, y = cch4_flux_mgm2d1, fill = veg_class),
    #            alpha = 0.3, shape = 21) +
    geom_hline(yintercept = 0, alpha = 0.6, linetype = "dashed") +
    geom_smooth(aes(x = .data[[x_mean]],
                    y = (exp(log_plant_cch4_flux_mgm2d1_mean) - 1),
                    fill=veg_class), method = glm, alpha=0.3)+
    scale_fill_manual(values = c("Open"="grey90",
                                 "Vegetated" = "grey90") )+
    new_scale_fill() + 
    geom_point(aes(fill=veg_class), colour="black", shape=23, size = 2.5)+
    scale_fill_manual(values = c("Open"="dodgerblue",
                                 "Vegetated"="darkolivegreen"))+
    scale_color_manual(values = c("Open"="grey80",
                                  "Vegetated"="grey80"))+
    xlab(xlabel) +
    ylab(expression(CH[4]~Flux~"("*mg~C~m^{-2}~day^{-1}*")")) +
    labs(fill = NULL) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none", 
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
  pm_plot
  
  diff_plot <- 
    ggplot(no_zone_avg, aes(x=.data[[x_mean]], 
                            y=cch4_flux_mgm2d1_mean, 
                            colour=veg_class, fill = veg_class)) + 
    geom_errorbar(aes(ymin=cch4_flux_mgm2d1_mean-cch4_flux_mgm2d1_se, 
                      ymax=cch4_flux_mgm2d1_mean+cch4_flux_mgm2d1_se),
                  alpha = 0.6) +
    geom_errorbarh(aes(xmin=.data[[x_mean]]-.data[[x_se]], 
                       xmax=.data[[x_mean]]+.data[[x_se]]), width = 0,
                   alpha = 0.6) +
    # geom_jitter(data = depth_modelling,
    #            aes(x = depth, y = cch4_flux_mgm2d1, fill = veg_class),
    #            alpha = 0.3, shape = 21) +
    geom_hline(yintercept = 0, alpha = 0.6, linetype = "dashed") +
    geom_smooth(aes(x = .data[[x_mean]],
                    y = (exp(log_cch4_flux_mgm2d1_mean) - 1),
                    fill=veg_class), method = glm, alpha=0.3)+
    scale_fill_manual(values = c("Open"="grey90",
                                 "Vegetated" = "grey90") )+
    new_scale_fill() + 
    geom_point(aes(fill=veg_class), colour="black", shape=21, size = 2.5)+
    scale_fill_manual(values = c("Open"="dodgerblue2",
                                 "Vegetated"="olivedrab3"))+
    scale_color_manual(values = c("Open"="grey80",
                                  "Vegetated"="grey80"))+
    xlab(xlabel) +
    ylab(expression(CH[4]~Flux~"("*mg~C~m^{-2}~day^{-1}*")")) +
    labs(fill = NULL) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none", 
          axis.title.y = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
  diff_plot
  
  ebull_plot <- 
    ggplot(no_zone_avg, aes(x=.data[[x_mean]], 
                            y=ebull_flux_ch4_mean, 
                            colour=veg_class, fill = veg_class)) + 
    geom_errorbar(aes(ymin=ebull_flux_ch4_mean-ebull_flux_ch4_se, 
                      ymax=ebull_flux_ch4_mean+ebull_flux_ch4_se),
                  width = 0, alpha = 0.6) +
    geom_errorbarh(aes(xmin=.data[[x_mean]]-.data[[x_se]], 
                       xmax=.data[[x_mean]]+.data[[x_se]]), 
                   width = 0, alpha = 0.6) +
    # geom_jitter(data = depth_modelling,
    #            aes(x = depth, y = cch4_flux_mgm2d1, fill = veg_class),
    #            alpha = 0.3, shape = 21) +
    geom_hline(yintercept = 0, alpha = 0.6, linetype = "dashed") +
    geom_smooth(aes(x = .data[[x_mean]],
                    y = (exp(log_ebull_flux_ch4_mean) - 1),
                    fill=veg_class), method = glm, alpha=0.3)+
    scale_fill_manual(values = c("Open"="grey90",
                                 "Vegetated" = "grey90") )+
    new_scale_fill() + 
    geom_point(aes(fill=veg_class), colour="black", shape=21, size = 2.5)+
    scale_fill_manual(values = c("Open"="dodgerblue2",
                                 "Vegetated"="olivedrab3"))+
    scale_color_manual(values = c("Open"="grey80",
                                  "Vegetated"="grey80"))+
    xlab(xlabel) +
    ylab(expression(CH[4]~Flux~"("*mg~C~m^{-2}~day^{-1}*")")) +
    labs(fill = NULL) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none", 
          axis.title.y = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
  ebull_plot
  
  ### plots together 
  combo_plot <- pm_plot | diff_plot | ebull_plot
  return(combo_plot)
}
