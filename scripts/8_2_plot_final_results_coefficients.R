######################################################################################################
######################################################################################################

# CODE DESCRIPTION

# Produces plots to compare standardized coefficient values across plant functional types
# Uses external validation data

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

# Set output directory
dir = '*/UAV_to_LS/results/FINAL_MODELS_nrmse_nmbe_corr_MC/coef'

# Import coefficient model -----
coef_FIELD_log = read.csv('*/UAV_to_LS/results/field/FINAL_MODELS_nrmse_nmbe_corr_avg_log/glmmLasso_coefficients_finalModels_field_log_nrmse_nmbe_corr_avg_MC_all.csv')
coef_FIELD_sqrt = read.csv('*/UAV_to_LS/results/field/FINAL_MODELS_nrmse_nmbe_corr_avg_sqrt/glmmLasso_coefficients_finalModels_field_sqrt_nrmse_nmbe_corr_avg_MC_all.csv')
coef_UAV_log = read.csv('*/UAV_to_LS/results/UAV/FINAL_MODELS_nrmse_nmbe_corr_avg_log/glmmLasso_coefficients_finalModels_UAV_log_nrmse_nmbe_corr_avg_MC_all.csv')
coef_UAV_sqrt = read.csv('*/UAV_to_LS/results/UAV/FINAL_MODELS_nrmse_nmbe_corr_avg_sqrt/glmmLasso_coefficients_finalModels_UAV_sqrt_nrmse_nmbe_corr_avg_MC_all.csv')

# Import best model parameters (minimum loss function) -----
mins_FIELD_log = read.csv('*/UAV_to_LS/results/field/log/minimums_rmse_nmbe_corr_avg_log.csv')
mins_FIELD_sqrt = read.csv('*/UAV_to_LS/results/field/sqrt/minimums_rmse_nmbe_corr_avg_sqrt.csv')
mins_UAV_log = read.csv('*/UAV_to_LS/results/UAV/log/minimums_rmse_nmbe_corr_avg_log.csv')
mins_UAV_sqrt = read.csv('*/UAV_to_LS/results/UAV/sqrt/minimums_rmse_nmbe_corr_avg_sqrt.csv')

######################################################################################################
######################################################################################################

# 1. SET UP ------------------------------------------------------

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

# 2. CHOOSE FIELD/UAV LOG/SQRT ------------------------------------------------------

# Choice for each PFT based on cross validation, external validation, and visual inspection of maps in Google Earth Engine
# Deciduous shrubs: UAV log
# Evergreen shrubs: UAV sqrt
# Forbs: field sqrt
# Graminoids: UAV sqrt
# Lichens: field sqrt

# For each PFT, select the best data type/response transformation -- coefficients -----
decid_coef = coef_UAV_log[coef_UAV_log$pft=='DECIDUOUS_SHRUBS',]
eg_coef = coef_UAV_sqrt[coef_UAV_sqrt$pft=='EVERGREEN_SHRUBS',]
forb_coef = coef_FIELD_sqrt[coef_FIELD_sqrt$pft=='FORBS',]
gram_coef = coef_UAV_sqrt[coef_UAV_sqrt$pft=='GRAMINOIDS',]
lichen_coef = coef_FIELD_sqrt[coef_FIELD_sqrt$pft=='LICHENS',]

# For each PFT, select the best data type/response transformation -- model parameters -----
decid_mins = mins_UAV_log[mins_UAV_log$pft=='DECIDUOUS_SHRUBS',]
eg_mins = mins_UAV_sqrt[mins_UAV_sqrt$pft=='EVERGREEN_SHRUBS',]
forb_mins = mins_FIELD_sqrt[mins_FIELD_sqrt$pft=='FORBS',]
gram_mins = mins_UAV_sqrt[mins_UAV_sqrt$pft=='GRAMINOIDS',]
lichen_mins = mins_FIELD_sqrt[mins_FIELD_sqrt$pft=='LICHENS',]

# For each PFT, from the coefficient data, select the best data type/response transformation
decid = decid_coef[(decid_coef$cut_num == decid_mins$cut_num & decid_coef$LC == decid_mins$LC),]
eg = eg_coef[(eg_coef$cut_num == eg_mins$cut_num & eg_coef$LC == eg_mins$LC),]
forb = forb_coef[(forb_coef$cut_num == forb_mins$cut_num & forb_coef$LC == forb_mins$LC),]
gram = gram_coef[(gram_coef$cut_num == gram_mins$cut_num & gram_coef$LC == gram_mins$LC),]
lichen = lichen_coef[(lichen_coef$cut_num == lichen_mins$cut_num & lichen_coef$LC == lichen_mins$LC),]

# 3. CALCULATE COEFFICIENTS - AVERAGE ------------------------------------------------------

decid = decid %>%
  group_by(predictor) %>%
  summarise(coef_std_avg = mean(abs(coef_std), na.rm=TRUE),
            sign = mean(coef_std, na.rm=TRUE))

decid$pft = 'Deciduous Shrubs'

eg = eg %>%
  group_by(predictor) %>%
  summarise(coef_std_avg = mean(abs(coef_std), na.rm=TRUE),
            sign = mean(coef_std, na.rm=TRUE))

eg$pft = 'Evergreen Shrubs'

forb = forb %>%
  group_by(predictor) %>%
  summarise(coef_std_avg = mean(abs(coef_std), na.rm=TRUE),
            sign = mean(coef_std, na.rm=TRUE))

forb$pft = 'Forbs'

gram = gram %>%
  group_by(predictor) %>%
  summarise(coef_std_avg = mean(abs(coef_std), na.rm=TRUE),
            sign = mean(coef_std, na.rm=TRUE))

gram$pft = 'Graminoids'

lichen = lichen %>%
  group_by(predictor) %>%
  summarise(coef_std_avg = mean(abs(coef_std), na.rm=TRUE),
            sign = mean(coef_std, na.rm=TRUE))

lichen$pft = 'Lichens'

# Aggregate

coef_avg_all = rbind(decid, eg, forb, gram, lichen)

# 4. TIDY ------------------------------------------------------

# Remove NAs
coef_avg_all = na.omit(coef_avg_all)

# Convert to lower case
coef_avg_all$predictor = tolower(coef_avg_all$predictor)

# Remove intercept
coef_avg_all = coef_avg_all[coef_avg_all$predictor != '(intercept)',]

# Format land cover names
coef_avg_all$predictor = dplyr::recode(coef_avg_all$predictor, 
                             'as.factor(lc)1' = 'LC: evergreen forest', 
                             'as.factor(lc)2' = 'LC: deciduous forest', 
                             'as.factor(lc)3' = 'LC: mixed forest', 
                             'as.factor(lc)4' = 'LC: woodland', 
                             'as.factor(lc)5' = 'LC: low shrub', 
                             'as.factor(lc)6' = 'LC: tall shrub', 
                             'as.factor(lc)7' = 'LC: open shrubs', 
                             'as.factor(lc)8' = 'LC: herbaceous', 
                             'as.factor(lc)9' = 'LC: tussock tundra', 
                             'as.factor(lc)10' = 'LC: sparsely vegetated', 
                             'as.factor(lc)11' = 'LC: fen', 
                             'as.factor(lc)12' = 'LC: bog', 
                             'as.factor(lc)13' = 'LC: shallows/littoral', 
                             'as.factor(lc)14' = 'LC: barren',
                             'as.factor(lc)15' = 'LC: water')

# Format landsat names
coef_avg_all$predictor = gsub("landsatccdc_", "", coef_avg_all$predictor)
coef_avg_all$predictor = gsub("_", " ", coef_avg_all$predictor)
coef_avg_all$pft = as.factor(gsub("_", " ", coef_avg_all$pft))
coef_avg_all$predictor = gsub("changeratep", "p", coef_avg_all$predictor)
coef_avg_all$predictor = gsub("seasonalp", "p", coef_avg_all$predictor)

# Format cover names
coef_avg_all$predictor = gsub("predicted cover", "cover", coef_avg_all$predictor)

# Format seasonal predictor names
coef_avg_all$predictor = gsub("p0", " to ", coef_avg_all$predictor)
coef_avg_all$predictor = gsub("  to", "", coef_avg_all$predictor)
coef_avg_all$predictor = gsub("05", "early spring", coef_avg_all$predictor)
coef_avg_all$predictor = gsub("20", "spring", coef_avg_all$predictor)
coef_avg_all$predictor = gsub("35", "spring/summer", coef_avg_all$predictor)
coef_avg_all$predictor = gsub("50", "summer", coef_avg_all$predictor)
coef_avg_all$predictor = gsub("65", "summer/autumn", coef_avg_all$predictor)
coef_avg_all$predictor = gsub("80", "autumn", coef_avg_all$predictor)
coef_avg_all$predictor = gsub("95", "late autumn", coef_avg_all$predictor)

# Assign signs 
coef_avg_all$sign = with(coef_avg_all, ifelse(sign > 0, '(+)', '(-)'))

# Format signs
coef_avg_all$predictor_sign = paste0(coef_avg_all$predictor, ' ', coef_avg_all$sign)

# 5. PLOT ------------------------------------------------------

textformula = y ~ x
title_text_size = 26
theme_text_size = title_text_size - 2
geom_text_size = theme_text_size / 4

# PFT lookup
PFTcolors = data.frame(PFT = c("BRYOPHYTES", "DECIDUOUS SHRUBS", "EVERGREEN SHRUBS", "FORBS", "GRAMINOIDS", "LICHENS", "NON VEGETATED", "SHADOW", "TREES", "TOTAL"), color = c('darkorange', 'green4', 'cyan4', 'chartreuse3', 'darkseagreen3', 'yellow3', 'burlywood4', 'dimgray', 'mediumaquamarine', 'grey60'))

colors = c('green4', 'cyan4', 'chartreuse3', 'darkseagreen3', 'yellow3') 

plt = 
  ggplot(coef_avg_all, aes(x = coef_std_avg,  y = reorder_within(predictor_sign, coef_std_avg, pft), col = pft, fill = pft))+
  geom_point(shape = 21, col = 'black', size = 6)+
  facet_wrap(.~ pft, scales = 'free', ncol = 2)+
  labs(y = 'Predictor', x = '\nStandardized Coefficient')+
  theme_minimal()+
  theme(title = element_text(size=30),
        strip.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18))+
  guides(color = 'none')+
  guides(fill = 'none')+
  scale_color_manual(values = colors)+
  scale_y_reordered()+
  scale_fill_manual(values = colors)

plt

if(output_results){

  outName = 'finalModels_MC_avg_std_coefficients_vertical.png'
  
  ggsave(
    paste0(dir, outName),
    plt,
    width = 50,
    height = 60,
    units = 'cm',
    dpi = 600,
    bg = 'white'
  )

}
