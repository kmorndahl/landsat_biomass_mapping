######################################################################################################
######################################################################################################

# CODE DESCRIPTION

# Produces plots to compare final model performance across data type and transform combinations as well as plant functional types
# Uses external validation data

# NOTE: output directory structure not hosted at github

library(dplyr)
library(tidyr)
library(yardstick)
library(ggplot2)
library(ggpmisc)
library(Metrics)

######################################################################################################
######################################################################################################

# 0. SET OUTPUT DIRECTORY AND GET DATA

output_results = FALSE

dir = '*/UAV_to_LS/results/FINAL_MODELS_nrmse_nmbe_corr_MC/externalValidation'

validation_FIELD_log = read.csv('external_biomass_data_sample_FIELD_log_nrmse_nmbe_corr_avg_MC_all.csv')
validation_FIELD_sqrt = read.csv('external_biomass_data_sample_FIELD_sqrt_nrmse_nmbe_corr_avg_MC_all.csv')
validation_UAV_log = read.csv('external_biomass_data_sample_UAV_log_nrmse_nmbe_corr_avg_MC_all.csv')
validation_UAV_sqrt = read.csv('external_biomass_data_sample_UAV_sqrt_nrmse_nmbe_corr_avg_MC_all.csv')

######################################################################################################
######################################################################################################

# 1. SET UP ------------------------------------------------------

# https://fishandwhistle.net/post/2018/modifying-facet-scales-in-ggplot2/
FacetEqualWrap <- ggproto(
  "FacetEqualWrap", FacetWrap,
  
  train_scales = function(self, x_scales, y_scales, layout, data, params) {
    
    # doesn't make sense if there is not an x *and* y scale
    if (is.null(x_scales) || is.null(x_scales)) {
      stop("X and Y scales required for facet_equal_wrap")
    }
    
    # regular training of scales
    ggproto_parent(FacetWrap, self)$train_scales(x_scales, y_scales, layout, data, params)
    
    # switched training of scales (x and y and y on x)
    for (layer_data in data) {
      match_id <- match(layer_data$PANEL, layout$PANEL)
      
      x_vars <- intersect(x_scales[[1]]$aesthetics, names(layer_data))
      y_vars <- intersect(y_scales[[1]]$aesthetics, names(layer_data))
      
      SCALE_X <- layout$SCALE_X[match_id]
      ggplot2:::scale_apply(layer_data, y_vars, "train", SCALE_X, x_scales)
      
      SCALE_Y <- layout$SCALE_Y[match_id]
      ggplot2:::scale_apply(layer_data, x_vars, "train", SCALE_Y, y_scales)
    }
    
  }
)

facet_wrap_equal <- function(...) {
  # take advantage of the sanitizing that happens in facet_wrap
  facet_super <- facet_wrap(...)
  
  ggproto(NULL, FacetEqualWrap,
          shrink = facet_super$shrink,
          params = facet_super$params
  )
}

# https://github.com/dgrtwo/drlib/blob/master/R/reorder_within.R

reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}

scale_x_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_x_discrete(labels = function(x) gsub(reg, "", x), ...)
}

scale_y_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_y_discrete(labels = function(x) gsub(reg, "", x), ...)
}

######################################################################################################
######################################################################################################

# 2. FIELD vs. UAV VALIDATION ------------------------------------------------------

# 2.1 CHOOSE LOG/SQRT ------------------------------------------------------

# Choice for each PFT based on cross validation, external validation, and visual inspection of maps in Google Earth Engine
# Deciduous shrubs: log
# Evergreen shrubs: sqrt
# Forbs: sqrt
# Graminoids: sqrt
# Lichens: sqrt

# 2.1.1 FIELD ------------------------------------------------------

decid_field =  dplyr::select(validation_FIELD_log, c("site", "region", "shrub_gm2_avg", "deciduousshrub_gm2_avg", "evergreenshrub_gm2_avg", "forb_gm2_avg", "graminoid_gm2_avg", "graminoid_live_gm2_avg", "lichen_gm2_avg", "biomass_decid_p50"))
eg_field = dplyr::select(validation_FIELD_sqrt, c("biomass_eg_p50"))
forb_field = dplyr::select(validation_FIELD_sqrt, c("biomass_forb_p50"))
gram_field = dplyr::select(validation_FIELD_sqrt, c("biomass_gram_p50"))
lichen_field = dplyr::select(validation_FIELD_sqrt, c("biomass_lichen_p50"))

validation_field = cbind(decid_field, eg_field, forb_field, gram_field, lichen_field)
validation_field$data_type = 'Field'

# 2.1.2 UAV ------------------------------------------------------

decid_uav =  dplyr::select(validation_UAV_log, c("site", "region", "shrub_gm2_avg", "deciduousshrub_gm2_avg", "evergreenshrub_gm2_avg", "forb_gm2_avg", "graminoid_gm2_avg", "graminoid_live_gm2_avg", "lichen_gm2_avg", "biomass_decid_p50"))
eg_uav = dplyr::select(validation_UAV_sqrt, c("biomass_eg_p50"))
forb_uav = dplyr::select(validation_UAV_sqrt, c("biomass_forb_p50"))
gram_uav = dplyr::select(validation_UAV_sqrt, c("biomass_gram_p50"))
lichen_uav = dplyr::select(validation_UAV_sqrt, c("biomass_lichen_p50"))

validation_uav = cbind(decid_uav, eg_uav, forb_uav, gram_uav, lichen_uav)
validation_uav$data_type = 'UAV'

# 2.2 COMBINE ------------------------------------------------------

validation_field_uav = rbind(validation_field, validation_uav)

# 2.3 TIDY ------------------------------------------------------

data = validation_field_uav

# Create live graminoid biomass only field
data$biomass_gram_live_p50 = data$biomass_gram_p50

# Create combined shrub biomass field
data$biomass_shrub_p50 = as.numeric(data$biomass_decid_p50 + data$biomass_eg_p50)

# Format shrub field
data$shrub_gm2_avg = as.numeric(gsub(',', '', data$shrub_gm2_avg))

# Divide into predicted and observed
data_predicted = dplyr::select(data, c("site", "region", "data_type", "biomass_shrub_p50", "biomass_eg_p50", "biomass_forb_p50", "biomass_gram_p50", "biomass_decid_p50", "biomass_lichen_p50", "biomass_gram_live_p50"))
data_observed = dplyr::select(data, c("site", "region", "data_type", "shrub_gm2_avg", "deciduousshrub_gm2_avg", "evergreenshrub_gm2_avg", "forb_gm2_avg", "graminoid_gm2_avg", "graminoid_live_gm2_avg", "lichen_gm2_avg"))

# Pivot longer - predicted
data_predicted = data_predicted %>%
  tidyr::pivot_longer(cols = starts_with("biomass_"), names_to = "pft", values_to = "predicted")

# Pivot longer - observed
data_observed = data_observed %>%
  tidyr::pivot_longer(!c(site, region, data_type), names_to = "pft", values_to = "observed")

# Rename - predicted
data_predicted$pft = gsub('biomass_gram_live_p50', 'Graminoids Live', data_predicted$pft)
data_predicted$pft = gsub('biomass_eg_p50', 'Evergreen Shrubs', data_predicted$pft)
data_predicted$pft = gsub('biomass_forb_p50', 'Forbs', data_predicted$pft)
data_predicted$pft = gsub('biomass_gram_p50', 'Graminoids', data_predicted$pft)
data_predicted$pft = gsub('biomass_decid_p50', 'Deciduous Shrubs', data_predicted$pft)
data_predicted$pft = gsub('biomass_lichen_p50', 'Lichens', data_predicted$pft)
data_predicted$pft = gsub('biomass_shrub_p50', 'Shrubs', data_predicted$pft)

# Rename - observed
data_observed$pft = gsub('graminoid_live_gm2_avg', 'Graminoids Live', data_observed$pft)
data_observed$pft = gsub('evergreenshrub_gm2_avg', 'Evergreen Shrubs', data_observed$pft)
data_observed$pft = gsub('forb_gm2_avg', 'Forbs', data_observed$pft)
data_observed$pft = gsub('graminoid_gm2_avg', 'Graminoids', data_observed$pft)
data_observed$pft = gsub('deciduousshrub_gm2_avg', 'Deciduous Shrubs', data_observed$pft)
data_observed$pft = gsub('lichen_gm2_avg', 'Lichens', data_observed$pft)
data_observed$pft = gsub('shrub_gm2_avg', 'Shrubs', data_observed$pft)

# Join predicted and observed
data_all_tidy = data_predicted %>% dplyr::left_join(data_observed, by = c("site", "region", "data_type", "pft"))

# Convert categorical fields to factor and tidy
data_all_tidy$pft = as.factor(data_all_tidy$pft)
data_all_tidy$site = as.factor(data_all_tidy$site)
data_all_tidy$region = as.factor(data_all_tidy$region)
data_all_tidy$region = gsub('AK_Seward', 'Seward Peninsula', data_all_tidy$region)
data_all_tidy$region = gsub('_', ' ', data_all_tidy$region)
data_all_tidy$region = gsub('AK', 'Alaska', data_all_tidy$region)

# Remove NAs
data_all_tidy = na.omit(data_all_tidy)

# Calculate grouped rmse and relative rmse
data_all_tidy = data.frame(data_all_tidy %>% 
                             dplyr::group_by(pft, data_type) %>%
                             dplyr::mutate(rmse = Metrics::rmse(observed, predicted), nrmse = mean(rmse)/mean(observed)))

# Calculate grouped mean bias error (mbe), normalized mean bias error, mean normalized bias error
data_all_tidy = data.frame(data_all_tidy %>% 
                             dplyr::group_by(pft, data_type) %>%
                             dplyr::mutate(mbe = mean(predicted - observed), nmbe = mean(mbe)/mean(observed)))

# Calculate correlation between actual and predicted
data_all_tidy = data.frame(data_all_tidy %>% 
                             dplyr::group_by(pft, data_type) %>%
                             dplyr::mutate(corr = cor(predicted, observed)))

# Get average of relative RMSE, mean normalized bias error, and correlation
# We need to convert correlation scores for use in average
# For average, we are trying to minimize the "error" score, so correlations need to be converted so that higher numbers are worse
# To do this: 1 - cor 
data_all_tidy$rmse_nmbe_corr_avg = (abs(data_all_tidy$nrmse) + abs(data_all_tidy$nmbe) + (1 - data_all_tidy$corr))/3

# Remove live+dead graminoids, rename live graminoids to graminoids
data_all_tidy = data_all_tidy[data_all_tidy$pft != 'Graminoids',]
data_all_tidy$pft = gsub('Graminoids Live', 'Gramioids', data_all_tidy$pft)

# Convert pft to factor
data_all_tidy$pft = as.factor(data_all_tidy$pft)

# 2.4 PLOT ------------------------------------------------------

# Set plotting parameters
colors = c('green4', 'cyan4', 'chartreuse3', 'darkseagreen3', 'yellow3', 'darkslategray')
textformula = y ~ x
title_text_size = 34
theme_text_size = title_text_size - 2
geom_text_size = theme_text_size / 4

plt = 
  ggplot(data_all_tidy, aes(x = observed,  y = predicted, col = pft, fill = pft))+
  geom_point(aes(shape = region), size = 8, alpha = 0.75, col = 'black')+
  geom_abline(slope = 1, intercept = 0, lty = 2)+
  geom_smooth(method='lm', formula= y~x, se = FALSE, col = 'black', size = 0.5)+
  facet_wrap_equal(pft ~ data_type, scales = 'free', nrow = 6, ncol = 2)+
  geom_label(data = data_all_tidy, aes(label=paste('RMSE = ', round(rmse, 2))), x = -Inf, y = Inf, hjust = -0.04, vjust = 1.4, color = 'grey20',  fill = 'white', size = geom_text_size-1, label.size = NA, alpha = 0.6, label.padding = unit(0.01, "lines"))+
  geom_label(data = data_all_tidy, aes(label=paste('nRMSE = ', round(nrmse, 2))), x = -Inf, y = Inf, hjust = -0.04, vjust = 2.7, color = 'grey20',  fill = 'white', size = geom_text_size-1, label.size = NA, alpha = 0.6, label.padding = unit(0.01, "lines"))+
  geom_label(data = data_all_tidy, aes(label=paste('MBE = ', round(mbe, 2))), x = -Inf, y = Inf, hjust = -0.04, vjust = 4.0, color = 'grey20',  fill = 'white', size = geom_text_size-1, label.size = NA, alpha = 0.6, label.padding = unit(0.01, "lines"))+
  geom_label(data = data_all_tidy, aes(label=paste('nMBE = ', round(nmbe, 2))), x = -Inf, y = Inf, hjust = -0.04, vjust = 5.3, color = 'grey20',  fill = 'white', size = geom_text_size-1, label.size = NA, alpha = 0.6, label.padding = unit(0.01, "lines"))+
  labs(y = expression(paste('Predicted Aboveground Biomass (g  ', m^-2, ')')), x = expression(paste('Observed Aboveground Biomass (g ', m^-2, ')')), shape = '')+
  theme_minimal()+
  theme(text = element_text(size=theme_text_size), 
        plot.margin = margin(10, 10, 10, 10),
        legend.position="bottom",
        strip.background = element_blank(),
        strip.text.x = element_blank())+
  guides(shape = guide_legend(override.aes = list(fill = 'black')), size = 30)+
  guides(color = 'none')+
  guides(fill = 'none')+
  scale_shape_manual(values=c(22, 24, 21))+
  scale_color_manual(values = colors)+
  scale_fill_manual(values = colors)+
  ggpmisc::stat_poly_eq(formula = textformula, geom = 'label', aes(label = paste(..rr.label.., sep = "~~~")), parse = TRUE, label.y = -Inf, label.x = Inf, hjust = 1.1, vjust = -0.5, color = 'grey20', fill = 'white', size = geom_text_size-1, label.size = NA, alpha = 0.6, label.padding = unit(0.01, "lines"))

plt

if(output_results){

  outName = paste0('externalValidation_field_vs_uav_bigger.png')
  
  ggsave(
    paste0(dir, outName),
    plt,
    width = 43.18,
    height = 55.88,
    units = 'cm',
    bg = 'white',
    dpi = 600
  )

}

######################################################################################################
######################################################################################################

# 3. OVERALL VALIDATION ------------------------------------------------------

# 3.1 CHOOSE LOG/SQRT ------------------------------------------------------

# Choice for each PFT based on cross validation, external validation, and visual inspection of maps in Google Earth Engine
# Deciduous shrubs: UAV log
# Evergreen shrubs: UAV sqrt
# Forbs: field sqrt
# Graminoids: UAV sqrt
# Lichens: field sqrt

decid =  dplyr::select(validation_UAV_log, c("site", "region", "shrub_gm2_avg", "deciduousshrub_gm2_avg", "evergreenshrub_gm2_avg", "forb_gm2_avg", "graminoid_gm2_avg", "graminoid_live_gm2_avg", "lichen_gm2_avg", "biomass_decid_p50"))
eg = dplyr::select(validation_UAV_sqrt, c("biomass_eg_p50"))
forb = dplyr::select(validation_FIELD_sqrt, c("biomass_forb_p50"))
gram = dplyr::select(validation_UAV_sqrt, c("biomass_gram_p50"))
lichen = dplyr::select(validation_FIELD_sqrt, c("biomass_lichen_p50"))

validation = cbind(decid, eg, forb, gram, lichen)

# 2.3 TIDY ------------------------------------------------------

data = validation

# Create live graminoid biomass only field
data$biomass_gram_live_p50 = data$biomass_gram_p50

# Create combined shrub biomass field
data$biomass_shrub_p50 = as.numeric(data$biomass_decid_p50 + data$biomass_eg_p50)

# Format shrub field
data$shrub_gm2_avg = as.numeric(gsub(',', '', data$shrub_gm2_avg))

# Divide into predicted and observed
data_predicted = dplyr::select(data, c("site", "region", "biomass_shrub_p50", "biomass_eg_p50", "biomass_forb_p50", "biomass_gram_p50", "biomass_decid_p50", "biomass_lichen_p50", "biomass_gram_live_p50"))
data_observed = dplyr::select(data, c("site", "region", "shrub_gm2_avg", "deciduousshrub_gm2_avg", "evergreenshrub_gm2_avg", "forb_gm2_avg", "graminoid_gm2_avg", "graminoid_live_gm2_avg", "lichen_gm2_avg"))

# Pivot longer - predicted
data_predicted = data_predicted %>%
  tidyr::pivot_longer(cols = starts_with("biomass_"), names_to = "pft", values_to = "predicted")

# Pivot longer - observed
data_observed = data_observed %>%
  tidyr::pivot_longer(!c(site, region), names_to = "pft", values_to = "observed")

# Rename - predicted
data_predicted$pft = gsub('biomass_gram_live_p50', 'Graminoids Live', data_predicted$pft)
data_predicted$pft = gsub('biomass_eg_p50', 'Evergreen Shrubs', data_predicted$pft)
data_predicted$pft = gsub('biomass_forb_p50', 'Forbs', data_predicted$pft)
data_predicted$pft = gsub('biomass_gram_p50', 'Graminoids', data_predicted$pft)
data_predicted$pft = gsub('biomass_decid_p50', 'Deciduous Shrubs', data_predicted$pft)
data_predicted$pft = gsub('biomass_lichen_p50', 'Lichens', data_predicted$pft)
data_predicted$pft = gsub('biomass_shrub_p50', 'Shrubs', data_predicted$pft)

# Rename - observed
data_observed$pft = gsub('graminoid_live_gm2_avg', 'Graminoids Live', data_observed$pft)
data_observed$pft = gsub('evergreenshrub_gm2_avg', 'Evergreen Shrubs', data_observed$pft)
data_observed$pft = gsub('forb_gm2_avg', 'Forbs', data_observed$pft)
data_observed$pft = gsub('graminoid_gm2_avg', 'Graminoids', data_observed$pft)
data_observed$pft = gsub('deciduousshrub_gm2_avg', 'Deciduous Shrubs', data_observed$pft)
data_observed$pft = gsub('lichen_gm2_avg', 'Lichens', data_observed$pft)
data_observed$pft = gsub('shrub_gm2_avg', 'Shrubs', data_observed$pft)

# Join predicted and observed
data_all_tidy = data_predicted %>% dplyr::left_join(data_observed, by = c("site", "region", "pft"))

# Convert categorical fields to factor and tidy
data_all_tidy$pft = as.factor(data_all_tidy$pft)
data_all_tidy$site = as.factor(data_all_tidy$site)
data_all_tidy$region = as.factor(data_all_tidy$region)
data_all_tidy$region = gsub('AK_Seward', 'Seward Peninsula', data_all_tidy$region)
data_all_tidy$region = gsub('_', ' ', data_all_tidy$region)
data_all_tidy$region = gsub('AK', 'Alaska', data_all_tidy$region)

# Remove NAs
data_all_tidy = na.omit(data_all_tidy)

# Calculate grouped rmse and relative rmse
data_all_tidy = data.frame(data_all_tidy %>% 
                             dplyr::group_by(pft) %>%
                             dplyr::mutate(rmse = Metrics::rmse(observed, predicted), nrmse = mean(rmse)/mean(observed)))

# Calculate grouped mean bias error (mbe), normalized mean bias error, mean normalized bias error
data_all_tidy = data.frame(data_all_tidy %>% 
                             dplyr::group_by(pft) %>%
                             dplyr::mutate(mbe = mean(predicted - observed), nmbe = mean(mbe)/mean(observed)))

# Calculate correlation between actual and predicted
data_all_tidy = data.frame(data_all_tidy %>% 
                             dplyr::group_by(pft) %>%
                             dplyr::mutate(corr = cor(predicted, observed)))

# Get average of relative RMSE, mean normalized bias error, and correlation
# We need to convert correlation scores for use in average
# For average, we are trying to minimize the "error" score, so correlations need to be converted so that higher numbers are worse
# To do this: 1 - cor 
data_all_tidy$rmse_nmbe_corr_avg = (abs(data_all_tidy$nrmse) + abs(data_all_tidy$nmbe) + (1 - data_all_tidy$corr))/3

# Remove live+dead graminoids, rename live graminoids to graminoids
data_all_tidy = data_all_tidy[data_all_tidy$pft != 'Graminoids',]
data_all_tidy$pft = gsub('Graminoids Live', 'Gramioids', data_all_tidy$pft)

# Convert pft to factor
data_all_tidy$pft = as.factor(data_all_tidy$pft)

# 2.4 PLOT ------------------------------------------------------

# Set plotting parameters
colors = c('green4', 'cyan4', 'chartreuse3', 'darkseagreen3', 'yellow3', 'darkslategray')
textformula = y ~ x
title_text_size = 34
theme_text_size = title_text_size - 2
geom_text_size = theme_text_size / 4

# Plot
plt = 
  ggplot(data_all_tidy, aes(x = observed,  y = predicted, col = pft, fill = pft))+
  geom_point(aes(shape = region), size = 6, alpha = 0.75, col = 'black')+
  geom_abline(slope = 1, intercept = 0, lty = 2)+
  geom_smooth(method='lm', formula= y~x, se = FALSE, col = 'black', size = 0.5)+
  facet_wrap_equal(~ pft, scales = 'free')+
  geom_label(data = data_all_tidy, aes(label=paste('RMSE = ', round(rmse, 2))), x = -Inf, y = Inf, hjust = -0.04, vjust = 1.4, color = 'grey20',  fill = 'white', size = geom_text_size-1, label.size = NA, alpha = 0.6, label.padding = unit(0.01, "lines"))+
  geom_label(data = data_all_tidy, aes(label=paste('nRMSE = ', round(nrmse, 2))), x = -Inf, y = Inf, hjust = -0.04, vjust = 2.7, color = 'grey20',  fill = 'white', size = geom_text_size-1, label.size = NA, alpha = 0.6, label.padding = unit(0.01, "lines"))+
  geom_label(data = data_all_tidy, aes(label=paste('MBE = ', round(mbe, 2))), x = -Inf, y = Inf, hjust = -0.04, vjust = 4.0, color = 'grey20',  fill = 'white', size = geom_text_size-1, label.size = NA, alpha = 0.6, label.padding = unit(0.01, "lines"))+
  geom_label(data = data_all_tidy, aes(label=paste('nMBE = ', round(nmbe, 2))), x = -Inf, y = Inf, hjust = -0.04, vjust = 5.3, color = 'grey20',  fill = 'white', size = geom_text_size-1, label.size = NA, alpha = 0.6, label.padding = unit(0.01, "lines"))+
  labs(y = expression(paste('Predicted Aboveground Biomass (g  ', m^-2, ')')), x = expression(paste('Observed Aboveground Biomass (g ', m^-2, ')')), shape = '')+
  theme_minimal()+
  theme(text = element_text(size=theme_text_size), 
        title = element_text(size=title_text_size), 
        plot.margin = margin(10, 10, 10, 10),
        legend.position="top")+
  guides(shape = guide_legend(override.aes = list(fill = 'black')), size = 30)+
  guides(color = 'none')+
  guides(fill = 'none')+
  scale_shape_manual(values=c(22, 24, 21))+
  scale_color_manual(values = colors)+
  scale_fill_manual(values = colors)+
  ggpmisc::stat_poly_eq(formula = textformula, geom = 'label', aes(label = paste(..rr.label.., sep = "~~~")), parse = TRUE, label.y = -Inf, label.x = Inf, hjust = 1.1, vjust = -0.5, color = 'grey20', fill = 'white', size = geom_text_size-1, label.size = NA, alpha = 0.6, label.padding = unit(0.01, "lines"))

plt

if(output_results){
  
  outName = paste0('externalValidation.png')
  
  ggsave(
    paste0(dir, outName),
    plt,
    width = 40,
    height = 30,
    units = 'cm',
    bg = 'white',
    dpi = 600
  )

}
