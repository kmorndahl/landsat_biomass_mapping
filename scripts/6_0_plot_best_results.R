######################################################################################################
######################################################################################################

# CODE DESCRIPTION

# For each combination of data type and transform, gets the prediction and coefficient data for the best model across cut numbers and:
#   - Plots predicted vs. observed across cross-validation folds and reports RMSE
#   - Plots coefficient values

# NOTE: output directory structure not hosted at github

library(ggplot2)
library(dplyr)
library(yardstick)
library(tidyr)
library(qdapTools)
library(ggpmisc)

library(tidyverse)

######################################################################################################
######################################################################################################

# 1. SET UP ------------------------------------------------------

# 1.0 SET PARAMETERS ------------------------------------------------------

output_results = FALSE

# Set data type
data_type = 'field' # Choose 'field' or 'UAV'

# Set transformation
transform = 'sqrt' # Choose 'sqrt' or 'log'

# 1.1 GET FILES ------------------------------------------------------

# Set directory
dir = paste0('*/UAV_to_LS/results/', data_type, '/', transform, '/')
print(paste0('The output directory is: ', dir))
cat("\n")

# Get file names
coef_paths = list.files(dir, full.names = TRUE, pattern = '^glmmLasso_BESTcoefficients_.*\\.csv$')
pred_paths = list.files(dir, full.names = TRUE, pattern = '^glmmLasso_BESTpredictions_.*\\.csv$')

# 1.2 FUNCTIONS ------------------------------------------------------

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

# 2. PLOT PREDICTIONS ------------------------------------------------------

for(pred_path in pred_paths){
  
  df = read.csv(pred_path, header = TRUE)
  
  min_loss_criteria = gsub('^.*BESTpredictions_\\s*|\\s*_LC.*$', '', pred_path)
  
  # 2.1 TIDY ------------------------------------------------------
  
  df$pft = as.factor(df$pft)
  
  # 2.2 CALCULATE RMSE ------------------------------------------------------
  
  RMSE.corrected = df %>%
    dplyr::group_by(pft, LC_combo, loss_function_combo) %>%
    yardstick::rmse(actual, prediction_corrected)
  
  RMSE.corrected = tidyr::spread(RMSE.corrected, .metric, .estimate)
  RMSE.corrected = subset(RMSE.corrected, select = -c(.estimator))
  RMSE.corrected$rmse_corrected_label = paste0('RMSE = ', round(RMSE.corrected$rmse, 2),'g')
  RMSE.corrected = dplyr::rename(RMSE.corrected, rmse_corrected = rmse)
  
  df = dplyr::left_join(df, RMSE.corrected, by = c('pft', 'LC_combo', 'loss_function_combo'))
  
  # 2.3 CALCULATE RELATIVE RMSE ------------------------------------------------------
  # We standardized RSE by mean AGB from field measurements (i.e. RSE/mean(AGB)) 
  # This statistic, which we refer to as RSE(%), controlled for differences in biomass magnitude between studies
  
  NRMSE.corrected = df %>%
    dplyr::group_by(pft, LC_combo, loss_function_combo) %>%
    dplyr::summarise(nrmse_corrected = mean(rmse_corrected)/mean(actual))
  
  NRMSE.corrected$nrmse_corrected_label = paste0('Relative RMSE = ', round(NRMSE.corrected$nrmse_corrected, 2))
  
  df = dplyr::left_join(df, NRMSE.corrected, by = c('pft', 'LC_combo', 'loss_function_combo'))
  
  # 2.4 CALCULATE MBE ------------------------------------------------------
  
  MBE.corrected = df %>%
    dplyr::group_by(pft, LC_combo, loss_function_combo) %>%
    dplyr::summarise(mbe_corrected = mean(prediction_corrected - actual))
  
  MBE.corrected$mbe_corrected_label = paste0('MBE = ', round(MBE.corrected$mbe_corrected, 2), 'g')
  
  df = dplyr::left_join(df, MBE.corrected, by = c('pft', 'LC_combo', 'loss_function_combo'))
  
  # 2.5 CALCULATE NMBE ------------------------------------------------------

  NMBE.corrected = df %>%
    dplyr::group_by(pft, LC_combo, loss_function_combo) %>%
    dplyr::summarise(nmbe_corrected = mean(mbe_corrected)/mean(actual))

  NMBE.corrected$nmbe_corrected_label = paste0('NMBE = ', round(NMBE.corrected$nmbe_corrected, 2))

  df = dplyr::left_join(df, NMBE.corrected, by = c('pft', 'LC_combo', 'loss_function_combo'))

  # 2.6 SET TEXT SIZE ------------------------------------------------------
  
  textformula = y ~ x
  title_text_size = 26
  theme_text_size = title_text_size - 2
  geom_text_size = theme_text_size / 4
  
  # 2.7 SET UP LOOK UP TABLES ------------------------------------------------------
  
  # PFT lookup
  PFTcolors = data.frame(PFT = c("BRYOPHYTES", "DECIDUOUS_SHRUBS", "EVERGREEN_SHRUBS", "FORBS", "GRAMINOIDS", "LICHENS", "NON VEGETATED", "SHADOW", "TREES", "TOTAL"), color = c('darkorange', 'green4', 'cyan4', 'chartreuse3', 'darkseagreen3', 'yellow3', 'burlywood4', 'dimgray', 'mediumaquamarine', 'grey60'))
  
  # Set PFT colors
  PFTs = as.character(levels(df$pft))
  colors = PFTs %l% PFTcolors

  # 2.8 PLOT RESULTS ------------------------------------------------------
  
  mod = lm(actual~prediction_corrected*pft, data = df)
  df = cbind(df, predict(mod, interval='conf'))
  
  plt = 
    ggplot(df, aes(y = actual,  x = prediction_corrected, col = pft, fill = pft))+
    geom_point(shape = 21, col = 'black', size = 3, alpha = 0.5)+
    geom_abline(slope = 1, intercept = 0, lty = 2)+
    geom_line(aes(y = fit), col = 'black', size = 0.5)+
    facet_wrap_equal(~ pft, scales = 'free')+
    geom_label(data = NRMSE.corrected, aes(label=nrmse_corrected_label), x = -Inf, y = Inf, hjust = -0.04, vjust = 1, color = 'grey20',  fill = 'white', size = geom_text_size-1, label.size = NA, alpha = 0.6, label.padding = unit(0.01, "lines"))+
    geom_label(data = RMSE.corrected, aes(label=rmse_corrected_label), x = -Inf, y = Inf, hjust = -0.05, vjust = 2.3, color = 'grey20',  fill = 'white', size = geom_text_size-1, label.size = NA, alpha = 0.6, label.padding = unit(0.01, "lines"))+
    geom_label(data = MBE.corrected, aes(label=mbe_corrected_label), x = -Inf, y = Inf, hjust = -0.06, vjust = 3.6, color = 'grey20',  fill = 'white', size = geom_text_size-1, label.size = NA, alpha = 0.6, label.padding = unit(0.01, "lines"))+
    geom_label(data = NMBE.corrected, aes(label=nmbe_corrected_label), x = -Inf, y = Inf, hjust = -0.04, vjust = 4.9, color = 'grey20',  fill = 'white', size = geom_text_size-1, label.size = NA, alpha = 0.6, label.padding = unit(0.01, "lines"))+
    labs(x = expression(paste('Predicted Above-ground Biomass (g 900', m^-2, ')')), y = expression(paste('Observed Above-ground Biomass (g 900', m^-2, ')')), fill = 'Plant Functional Type')+
    theme_minimal()+
    theme(text = element_text(size=theme_text_size), 
          title = element_text(size=title_text_size), 
          plot.margin = margin(10, 10, 10, 10))+
    guides(color = 'none')+
    guides(fill = 'none')+
    scale_color_manual(values = colors)+
    scale_fill_manual(values = colors)+
    ggpmisc::stat_poly_eq(formula = textformula, geom = 'label', aes(label = paste(..rr.label.., sep = "~~~")), parse = TRUE, label.y = -Inf, label.x = Inf, hjust = 1, vjust = -1.2, color = 'grey20', fill = 'white', size = geom_text_size-1, label.size = NA, alpha = 0.6, label.padding = unit(0.01, "lines"))+
    ggpmisc::stat_poly_eq(formula = textformula, geom = 'label', aes(label = paste(..eq.label.., sep = "~~~")), parse = TRUE, label.y = -Inf, label.x = Inf, hjust = 1, vjust = 0, color = 'grey20', fill = 'white', size = geom_text_size-1, label.size = NA, alpha = 0.6, label.padding = unit(0.01, "lines"))
  
  plt
  
  if(output_results){
  
    outName = paste0('glmmLasso_results_', min_loss_criteria, '.png')
    
    ggsave(
      paste0(dir, outName),
      plt,
      width = 40,
      height = 30,
      units = 'cm'
    )
  
  }
  
}

######################################################################################################
######################################################################################################

# 2. PLOT COEFFICIENTS ------------------------------------------------------

for(coef_path in coef_paths){

  df = read.csv(coef_path, header = TRUE)
  
  min_loss_criteria = gsub('^.*BESTcoefficients_\\s*|\\s*_LC.*$', '', coef_path)
  
  # 2.1 TIDY ------------------------------------------------------
  
  df$pft = as.factor(df$pft)
  
  # 2.2 ABSOLUTE VALUE OF COEFFICIENTS FOR COMPARISON ------------------------------------------------------
  
  df$coef_std_abs = abs(df$coef_std)
  
  # 2.3 SUMMARIZE COEFFICIENTS ------------------------------------------------------
  
  df_summary = df %>% 
    dplyr::group_by(pft, predictor) %>%
    dplyr::summarise(coef_std_abs = mean(coef_std_abs))
  
  # 2.4 PLOT ------------------------------------------------------
  
  plt = 
    ggplot(df_summary[df_summary$coef_std_abs > 0,], aes(x = coef_std_abs,  y = reorder_within(predictor, coef_std_abs, pft), col = pft, fill = pft))+
    geom_point(shape = 21, col = 'black', size = 3, alpha = 0.5)+
    facet_wrap(.~ pft, scales = 'free_y')+
    labs(y = 'Predictor', x = 'Standardized Coefficient')+
    theme_minimal()+
    theme(title = element_text(size=20),
          strip.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12))+
    guides(color = 'none')+
    guides(fill = 'none')+
    scale_color_manual(values = colors)+
    scale_y_reordered()+
    scale_fill_manual(values = colors)
  
  plt
  
  if(output_results){
  
    outName = paste0('glmmLasso_coefficients_', min_loss_criteria, '.png')
    
    ggsave(
      paste0(dir, outName),
      plt,
      width = 60,
      height = 30,
      units = 'cm'
    )
  
  }
  
}