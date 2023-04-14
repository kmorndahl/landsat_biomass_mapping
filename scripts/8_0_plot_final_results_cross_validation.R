######################################################################################################
######################################################################################################

# CODE DESCRIPTION

# Produces plots to compare final model performance across data type and transform combinations as well as plant functional types
# Uses cross-validation data

# NOTE: output directory structure not hosted at github

library(dplyr)
library(tidyr)
library(yardstick)
library(ggplot2)
library(ggpmisc)

######################################################################################################
######################################################################################################

# 0. SET OUTPUT DIRECTORY AND GET DATA

output_results = FALSE

validation_FIELD_log = read.csv('*/UAV_to_LS/results/field/log/glmmLasso_BESTpredictions_nrmse_nmbe_corr_avg.csv')
validation_FIELD_sqrt = read.csv('*/UAV_to_LS/results/field/sqrt/glmmLasso_BESTpredictions_nrmse_nmbe_corr_avg.csv')
validation_UAV_log = read.csv('*/UAV_to_LS/results/UAV/log/glmmLasso_BESTpredictions_nrmse_nmbe_corr_avg.csv')
validation_UAV_sqrt = read.csv('*/UAV_to_LS/results/UAV/sqrt/glmmLasso_BESTpredictions_nrmse_nmbe_corr_avg.csv')

dir = '*/UAV_to_LS/results/FINAL_MODELS_nrmse_nmbe_corr_MC/crossValidation'

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

decid_field =  validation_FIELD_log[validation_FIELD_log$pft == 'DECIDUOUS_SHRUBS',]
eg_field = validation_FIELD_sqrt[validation_FIELD_sqrt$pft == 'EVERGREEN_SHRUBS',]
forb_field = validation_FIELD_sqrt[validation_FIELD_sqrt$pft == 'FORBS',]
gram_field = validation_FIELD_sqrt[validation_FIELD_sqrt$pft == 'GRAMINOIDS',]
lichen_field = validation_FIELD_sqrt[validation_FIELD_sqrt$pft == 'LICHENS',]

validation_field = rbind(decid_field, eg_field, forb_field, gram_field, lichen_field)
validation_field$data_type = 'Field'

# 2.1.2 UAV ------------------------------------------------------

decid_uav =  validation_UAV_log[validation_UAV_log$pft == 'DECIDUOUS_SHRUBS',]
eg_uav = validation_UAV_sqrt[validation_UAV_sqrt$pft == 'EVERGREEN_SHRUBS',]
forb_uav = validation_UAV_sqrt[validation_UAV_sqrt$pft == 'FORBS',]
gram_uav = validation_UAV_sqrt[validation_UAV_sqrt$pft == 'GRAMINOIDS',]
lichen_uav = validation_UAV_sqrt[validation_UAV_sqrt$pft == 'LICHENS',]

validation_uav = rbind(decid_uav, eg_uav, forb_uav, gram_uav, lichen_uav)
validation_uav$data_type = 'UAV'

# 2.2 COMBINE ------------------------------------------------------

validation_field_uav = rbind(validation_field, validation_uav)

# 2.3 TIDY ------------------------------------------------------

# Convert to factor
validation_field_uav$pft = as.factor(validation_field_uav$pft)

# Convert to g/m2
validation_field_uav$prediction_corrected = validation_field_uav$prediction_corrected/900
validation_field_uav$actual = validation_field_uav$actual/900

# Tidy names
validation_field_uav$pft = gsub('DECIDUOUS_SHRUBS', 'Deciduous Shrubs', validation_field_uav$pft)
validation_field_uav$pft = gsub('EVERGREEN_SHRUBS', 'Evergreen Shrubs', validation_field_uav$pft)
validation_field_uav$pft = gsub('FORBS', 'Forbs', validation_field_uav$pft)
validation_field_uav$pft = gsub('GRAMINOIDS', 'Graminoids', validation_field_uav$pft)
validation_field_uav$pft = gsub('LICHENS', 'Lichens', validation_field_uav$pft)

# 2.4 CALCULATE RMSE ------------------------------------------------------

RMSE.corrected = validation_field_uav %>%
  dplyr::group_by(pft, data_type) %>%
  yardstick::rmse(actual, prediction_corrected)

RMSE.corrected = dplyr::spread(RMSE.corrected, .metric, .estimate)
RMSE.corrected = subset(RMSE.corrected, select = -c(.estimator))
RMSE.corrected$rmse_corrected_label = paste0('RMSE = ', round(RMSE.corrected$rmse, 2),'g')
RMSE.corrected = dplyr::rename(RMSE.corrected, rmse_corrected = rmse)

validation_field_uav = dplyr::left_join(validation_field_uav, RMSE.corrected, by = c('pft', 'data_type'))

# 2.5 CALCULATE RELATIVE RMSE ------------------------------------------------------
# We standardized RMSE by mean AGB from field measurements (i.e. RSE/mean(AGB)) 

NRMSE.corrected = validation_field_uav %>%
  dplyr::group_by(pft, data_type) %>%
  summarise(nrmse_corrected = mean(rmse_corrected)/mean(actual))

NRMSE.corrected$nrmse_corrected_label = paste0('nRMSE = ', round(NRMSE.corrected$nrmse_corrected, 2))

validation_field_uav = dplyr::left_join(validation_field_uav, NRMSE.corrected, by = c('pft', 'data_type'))

# 2.6 CALCULATE MBE ------------------------------------------------------

MBE.corrected = validation_field_uav %>%
  dplyr::group_by(pft, data_type) %>%
  summarise(mbe_corrected = mean(prediction_corrected - actual))

MBE.corrected$mbe_corrected_label = paste0('MBE = ', round(MBE.corrected$mbe_corrected, 2), 'g')

validation_field_uav = dplyr::left_join(validation_field_uav, MBE.corrected, by = c('pft', 'data_type'))

# 2.5 CALCULATE RELATIVE MBE ------------------------------------------------------

NMBE.corrected = validation_field_uav %>%
  dplyr::group_by(pft, data_type) %>%
  summarise(nmbe_corrected = mean(mbe_corrected)/mean(actual))

NMBE.corrected$nmbe_corrected_label = paste0('nMBE = ', round(NMBE.corrected$nmbe_corrected, 2))

validation_field_uav = dplyr::left_join(validation_field_uav, NMBE.corrected, by = c('pft', 'data_type'))

# 4.4 SET TEXT SIZE ------------------------------------------------------

textformula = y ~ x
title_text_size = 34
theme_text_size = title_text_size - 2
geom_text_size = theme_text_size / 4

# 4.5 SET COLORS ------------------------------------------------------

colors = c('green4', 'cyan4', 'chartreuse3', 'darkseagreen3', 'yellow3')

# 4.6 PLOT ------------------------------------------------------

# Model
mod = lm(actual~prediction_corrected*pft, data = validation_field_uav)
validation_field_uav = cbind(validation_field_uav, predict(mod, interval='conf'))

# Plot
plt = 
  ggplot(validation_field_uav[validation_field_uav$pft != 'TOTAL',], aes(y = actual,  x = prediction_corrected, col = pft, fill = pft))+
  geom_point(shape = 21, col = 'black', size = 3, alpha = 0.5)+
  geom_abline(slope = 1, intercept = 0, lty = 2)+
  geom_smooth(method='lm', formula= y~x, se = FALSE, col = 'black', size = 0.5)+
  facet_wrap_equal(pft ~ data_type, scales = 'free', nrow = 5, ncol = 2)+
  geom_label(data = RMSE.corrected[RMSE.corrected$pft != 'TOTAL',], aes(label=rmse_corrected_label), x = -Inf, y = Inf, hjust = -0.04, vjust = 1.4, color = 'grey20',  fill = 'white', size = geom_text_size-1, label.size = NA, alpha = 0.6, label.padding = unit(0.01, "lines"))+
  geom_label(data = NRMSE.corrected[NRMSE.corrected$pft != 'TOTAL',], aes(label=nrmse_corrected_label), x = -Inf, y = Inf, hjust = -0.04, vjust = 2.7, color = 'grey20',  fill = 'white', size = geom_text_size-1, label.size = NA, alpha = 0.6, label.padding = unit(0.01, "lines"))+
  geom_label(data = MBE.corrected[MBE.corrected$pft != 'TOTAL',], aes(label=mbe_corrected_label), x = -Inf, y = Inf, hjust = -0.04, vjust = 4.0, color = 'grey20',  fill = 'white', size = geom_text_size-1, label.size = NA, alpha = 0.6, label.padding = unit(0.01, "lines"))+
  geom_label(data = NMBE.corrected[NMBE.corrected$pft != 'TOTAL',], aes(label=nmbe_corrected_label), x = -Inf, y = Inf, hjust = -0.04, vjust = 5.3, color = 'grey20',  fill = 'white', size = geom_text_size-1, label.size = NA, alpha = 0.6, label.padding = unit(0.01, "lines"))+
  labs(x = expression(paste('Predicted Aboveground Biomass (g ', m^-2, ')')), y = expression(paste('Observed Aboveground Biomass (g  ', m^-2, ')')), fill = 'Plant Functional Type')+
  theme_minimal()+
  theme(text = element_text(size=theme_text_size), 
        plot.margin = margin(10, 10, 10, 10),
        legend.position="bottom",
        strip.background = element_blank(),
        strip.text.x = element_blank())+
  guides(color = 'none')+
  guides(fill = 'none')+
  scale_color_manual(values = colors)+
  scale_fill_manual(values = colors)+
  ggpmisc::stat_poly_eq(formula = textformula, geom = 'label', aes(label = paste(..rr.label.., sep = "~~~")), parse = TRUE, label.y = -Inf, label.x = Inf, hjust = 1, vjust = 0, color = 'grey20', fill = 'white', size = geom_text_size-1, label.size = NA, alpha = 0.6, label.padding = unit(0.01, "lines"))

plt

if(output_results){

  outName = 'CVresults_field_vs_uav.png'
  
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

decid =  validation_UAV_log[validation_UAV_log$pft == 'DECIDUOUS_SHRUBS',]
eg = validation_UAV_sqrt[validation_UAV_sqrt$pft == 'EVERGREEN_SHRUBS',]
forb = validation_FIELD_sqrt[validation_FIELD_sqrt$pft == 'FORBS',]
gram = validation_UAV_sqrt[validation_UAV_sqrt$pft == 'GRAMINOIDS',]
lichen = validation_FIELD_sqrt[validation_FIELD_sqrt$pft == 'LICHENS',]

validation = rbind(decid, eg, forb, gram, lichen)

# 3.2 TIDY ------------------------------------------------------

# Convert to factor
validation$pft = as.factor(validation$pft)

# Convert to g/m2
validation$prediction_corrected = validation$prediction_corrected/900
validation$actual = validation$actual/900

# Tidy names
validation$pft = gsub('DECIDUOUS_SHRUBS', 'Deciduous Shrubs', validation$pft)
validation$pft = gsub('EVERGREEN_SHRUBS', 'Evergreen Shrubs', validation$pft)
validation$pft = gsub('FORBS', 'Forbs', validation$pft)
validation$pft = gsub('GRAMINOIDS', 'Graminoids', validation$pft)
validation$pft = gsub('LICHENS', 'Lichens', validation$pft)

# 3.3 CALCULATE RMSE ------------------------------------------------------

RMSE.corrected = validation %>%
  dplyr::group_by(pft) %>%
  yardstick::rmse(actual, prediction_corrected)

RMSE.corrected = dplyr::spread(RMSE.corrected, .metric, .estimate)
RMSE.corrected = subset(RMSE.corrected, select = -c(.estimator))
RMSE.corrected$rmse_corrected_label = paste0('RMSE = ', round(RMSE.corrected$rmse, 2),'g')
RMSE.corrected = dplyr::rename(RMSE.corrected, rmse_corrected = rmse)

validation = dplyr::left_join(validation, RMSE.corrected, by = c('pft'))

# 3.4 CALCULATE RELATIVE RMSE ------------------------------------------------------
# We standardized RSE by mean AGB from field measurements (i.e. RSE/mean(AGB)) 
# This statistic, which we refer to as RSE(%), controlled for differences in biomass magnitude between studies

NRMSE.corrected = validation %>%
  dplyr::group_by(pft) %>%
  summarise(nrmse_corrected = mean(rmse_corrected)/mean(actual))

NRMSE.corrected$nrmse_corrected_label = paste0('nRMSE = ', round(NRMSE.corrected$nrmse_corrected, 2))

validation = dplyr::left_join(validation, NRMSE.corrected, by = c('pft'))

# 3.5 CALCULATE MBE ------------------------------------------------------

MBE.corrected = validation %>%
  dplyr::group_by(pft) %>%
  summarise(mbe_corrected = mean(prediction_corrected - actual))

MBE.corrected$mbe_corrected_label = paste0('MBE = ', round(MBE.corrected$mbe_corrected, 2), 'g')

validation = dplyr::left_join(validation, MBE.corrected, by = c('pft'))

# 3.6 CALCULATE RELATIVE MBE ------------------------------------------------------

NMBE.corrected = validation %>%
  dplyr::group_by(pft) %>%
  summarise(nmbe_corrected = mean(mbe_corrected)/mean(actual))

NMBE.corrected$nmbe_corrected_label = paste0('nMBE = ', round(NMBE.corrected$nmbe_corrected, 2))

validation = dplyr::left_join(validation, NMBE.corrected, by = c('pft'))

# 3.7 SET TEXT SIZE ------------------------------------------------------

textformula = y ~ x
title_text_size = 34
theme_text_size = title_text_size - 2
geom_text_size = theme_text_size / 4

# 3.8 SET COLORS ------------------------------------------------------

colors = c('green4', 'cyan4', 'chartreuse3', 'darkseagreen3', 'yellow3')

# 3.9 PLOT ------------------------------------------------------

# Model
mod = lm(actual~prediction_corrected*pft, data = validation)
validation = cbind(validation, predict(mod, interval='conf'))

# Plot
plt = 
  ggplot(validation[validation$pft != 'TOTAL',], aes(y = actual,  x = prediction_corrected, col = pft, fill = pft))+
  geom_point(shape = 21, col = 'black', size = 3, alpha = 0.5)+
  geom_abline(slope = 1, intercept = 0, lty = 2)+
  geom_smooth(method='lm', formula= y~x, se = FALSE, col = 'black', size = 0.5)+
  facet_wrap_equal(~ pft, scales = 'free')+
  geom_label(data = RMSE.corrected[RMSE.corrected$pft != 'TOTAL',], aes(label=rmse_corrected_label), x = -Inf, y = Inf, hjust = -0.04, vjust = 1.4, color = 'grey20',  fill = 'white', size = geom_text_size-1, label.size = NA, alpha = 0.6, label.padding = unit(0.01, "lines"))+
  geom_label(data = NRMSE.corrected[NRMSE.corrected$pft != 'TOTAL',], aes(label=nrmse_corrected_label), x = -Inf, y = Inf, hjust = -0.04, vjust = 2.7, color = 'grey20',  fill = 'white', size = geom_text_size-1, label.size = NA, alpha = 0.6, label.padding = unit(0.01, "lines"))+
  geom_label(data = MBE.corrected[MBE.corrected$pft != 'TOTAL',], aes(label=mbe_corrected_label), x = -Inf, y = Inf, hjust = -0.04, vjust = 4.0, color = 'grey20',  fill = 'white', size = geom_text_size-1, label.size = NA, alpha = 0.6, label.padding = unit(0.01, "lines"))+
  geom_label(data = NMBE.corrected[NMBE.corrected$pft != 'TOTAL',], aes(label=nmbe_corrected_label), x = -Inf, y = Inf, hjust = -0.04, vjust = 5.3, color = 'grey20',  fill = 'white', size = geom_text_size-1, label.size = NA, alpha = 0.6, label.padding = unit(0.01, "lines"))+
  labs(x = expression(paste('Predicted Aboveground Biomass (g ', m^-2, ')')), y = expression(paste('Observed Aboveground Biomass (g  ', m^-2, ')')), fill = 'Plant Functional Type')+
  theme_minimal()+
  theme(text = element_text(size=theme_text_size), 
        title = element_text(size=title_text_size), 
        plot.margin = margin(10, 25, 10, 10))+
  guides(color = 'none')+
  guides(fill = 'none')+
  scale_color_manual(values = colors)+
  scale_fill_manual(values = colors)+
  ggpmisc::stat_poly_eq(formula = textformula, geom = 'label', aes(label = paste(..rr.label.., sep = "~~~")), parse = TRUE, label.y = -Inf, label.x = Inf, hjust = 1, vjust = 0, color = 'grey20', fill = 'white', size = geom_text_size-1, label.size = NA, alpha = 0.6, label.padding = unit(0.01, "lines"))

plt


if(output_results){
    
  outName = 'CVresults.png'
  
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

