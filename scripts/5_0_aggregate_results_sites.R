library(ggplot2)
library(dplyr)
library(yardstick)
library(tidyr)
library(qdapTools)
library(ggpmisc)

######################################################################################################
######################################################################################################

# 1. SET UP ------------------------------------------------------

# 1.0 SET PARAMETERS ------------------------------------------------------

missing_cv = 'stop' # Choose what to do if less that 44 files are present i.e. missing file for one or more outer CV folds: 'warn': just warn and continue running, or 'stop': halt code execution
aggregated_skip = FALSE # Choose what to do if an aggregated file (aggregated across all outer CV folds) ending in '_ALL' is already present: TRUE: skip the cut number for which the aggregate file is already present, FALSE: don't skip, instead overwrite existing file

# Set data type
data_type = 'field' # Choose 'field' or 'UAV'

# Set transformation
transform = 'sqrt' # Choose 'sqrt' or 'log'

# 1.1 GET FOLDERS ------------------------------------------------------

dir = paste0('/scratch/kmo265/UAV_to_LS/results/', data_type, '/', transform, '/')

folders = list.files(dir, full.names = TRUE, recursive = FALSE, pattern = '*cut*')

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

######################################################################################################
######################################################################################################

# 2. AGGREGATE FILES ------------------------------------------------------

for(folder in folders){ ### START FOLDER LOOP
  
  print(paste0('===== The folder is: ', folder, ' ====='))
  cat("\n")
  
  # Get base name for files
  file_pieces = strsplit(folder, '/')[[1]]
  file_name = file_pieces[length(file_pieces)]

  # 2.1 GATHER FILES ------------------------------------------------------
  
  # Get files
  coef_paths = list.files(folder, '^glmmLasso_coefficients_.*\\.csv$')
  pred_paths = list.files(folder, '^glmmLasso_predictions_.*\\.csv$')
  
  # Get aggregated files ending with "_ALL.csv," if present
  aggregated_coef = coef_paths[(grepl("_ALL.csv$", coef_paths))]
  aggregated_pred = pred_paths[(grepl("_ALL.csv$", pred_paths))]
  
  # If indicated, skip cut number if aggregated file already present
  if(aggregated_skip & length(aggregated_coef) > 0 & length(aggregated_pred) > 0){
    print('Aggregated file already present for this folder, moving on to next folder...')
    cat("\n")
    next
  }
  
  # Remove aggregated files from file list, if present and if aggregated_present is not indicated
  if(length(aggregated_coef) > 0){coef_paths = coef_paths[coef_paths != aggregated_coef]}
  if(length(aggregated_pred) > 0){pred_paths = pred_paths[pred_paths != aggregated_pred]}
  
  # Check for 44 files
  if(length(coef_paths) < 44 | length(pred_paths) < 44){
    
    if(missing_cv == 'warn'){
      # Warn instead of stop if model fit issues at lower cut numbers
      print('========== WARNING: Not all outer CV files present for this cut number ==========')
    }else if(missing_cv  == 'stop'){
      stop('Not all files present, check folder and try again')
    }else{
      stop('Missing CV parameter not recognized, choose warn or stop')
    }
      
  }
  
  setwd(folder)
  
  # 2.2 AGGREGATE FILES ------------------------------------------------------
  
  # Loop through coefficient files
  coefTotal = data.frame()
  for(path in coef_paths){
    
    cvDF = try(read.csv(path, header = T))
    
    if(class(cvDF) == 'data.frame'){
    
      coefTotal = rbind(coefTotal, cvDF)
    
    }
    
  }
  
  cat("\n")
  
  # Loop through predictor files
  predTotal = data.frame()
  for(path in pred_paths){
    
    cvDF = try(read.csv(path, header = T))
    
    if(class(cvDF) == 'data.frame'){
      
      predTotal = rbind(predTotal, cvDF)
      
    }
    
  }
  
  # 2.3 SAVE AGGREGATED FILES ------------------------------------------------------
  
  coef_name = paste0('glmmLasso_coefficients_', file_name, '_ALL.csv')
  pred_name = paste0('glmmLasso_predictions_', file_name, '_ALL.csv')
  
  write.csv(coefTotal, coef_name, row.names = FALSE)
  write.csv(predTotal, pred_name, row.names = FALSE)

  print('Aggregated files saved')
  cat("\n")
  
  ######################################################################################################
  ######################################################################################################
  
  # 3. PREPARE DATA FOR GRAPHING ------------------------------------------------------
 
  # 3.1 TIDY ------------------------------------------------------
  
  coefTotal$pft = as.factor(coefTotal$pft)
  predTotal$pft = as.factor(predTotal$pft)
  
  # 3.2 CALCULATE RMSE ------------------------------------------------------

  RMSE.corrected = predTotal %>%
    dplyr::group_by(pft, loss_function) %>%
    yardstick::rmse(actual, prediction_corrected)
  
  RMSE.corrected = tidyr::spread(RMSE.corrected, .metric, .estimate)
  RMSE.corrected = subset(RMSE.corrected, select = -c(.estimator))
  RMSE.corrected$rmse_corrected_label = paste0('RMSE = ', round(RMSE.corrected$rmse, 2),'g')
  RMSE.corrected = dplyr::rename(RMSE.corrected, rmse_corrected = rmse)
  
  predTotal = dplyr::left_join(predTotal, RMSE.corrected, by = c('pft', 'loss_function'))
  
  print('RMSE calculated:')
  print(RMSE.corrected)
  cat("\n")
  
  # 3.2 CALCULATE RELATIVE RMSE ------------------------------------------------------
  # We standardized RSE by mean AGB from field measurements (i.e. RSE/mean(AGB)) 
  # This statistic, which we refer to as RSE(%), controlled for differences in biomass magnitude between studies

  NRMSE.corrected = predTotal %>%
    dplyr::group_by(pft, loss_function) %>%
    dplyr::summarise(rmse_corrected_relative = mean(rmse_corrected)/mean(actual))
  
  NRMSE.corrected$rmse_corrected_relative_label = paste0('Relative RMSE = ', round(NRMSE.corrected$rmse_corrected_relative, 2))
  
  predTotal = dplyr::left_join(predTotal, NRMSE.corrected, by = c('pft', 'loss_function'))
 
  print('Relative RMSE calculated:')
  print(NRMSE.corrected)
  cat("\n")
   
  # 3.3 SET TEXT SIZE ------------------------------------------------------
  
  textformula = y ~ x
  title_text_size = 26
  theme_text_size = title_text_size - 2
  geom_text_size = theme_text_size / 4
  
  # 3.4 SET UP LOOK UP TABLES ------------------------------------------------------
  
  # PFT lookup
  PFTcolors = data.frame(PFT = c("BRYOPHYTES", "DECIDUOUS_SHRUBS", "EVERGREEN_SHRUBS", "FORBS", "GRAMINOIDS", "LICHENS", "NON VEGETATED", "SHADOW", "TREES", "TOTAL"), color = c('darkorange', 'green4', 'cyan4', 'chartreuse3', 'darkseagreen3', 'yellow3', 'burlywood4', 'dimgray', 'mediumaquamarine', 'grey60'))
  
  # Set PFT colors
  PFTs = as.character(levels(predTotal$pft))
  colors = PFTs %l% PFTcolors

  # Manually set PFT colors
  colors = c('green4', 'cyan4', 'chartreuse3', 'darkseagreen3', 'yellow3', 'grey60')
  
  ######################################################################################################
  ######################################################################################################
  
  # 4. PLOT RMSE ------------------------------------------------------
  
  df = predTotal
  
  RMSE_corrected = 
    ggplot(df[df$loss_function == 'rmse_corrected',], aes(x = actual,  y = prediction_corrected, fill = pft))+
    geom_point(shape = 21, col = 'black', size = 3, alpha = 0.5)+
    geom_abline(slope = 1, intercept = 0, lty = 2)+
    facet_wrap_equal(~ pft, scales = 'free')+
    geom_label(data = NRMSE.corrected[NRMSE.corrected$loss_function == 'rmse_corrected',], aes(label=rmse_corrected_relative_label), x = -Inf, y = Inf, hjust = -0.04, vjust = 1.8, color = 'grey20',  fill = 'white', size = geom_text_size, label.size = NA, alpha = 0.6, label.padding = unit(0.01, "lines"))+
    geom_label(data = RMSE.corrected[RMSE.corrected$loss_function == 'rmse_corrected',], aes(label=rmse_corrected_label), x = -Inf, y = Inf, hjust = -0.06, vjust = 3.1, color = 'grey20',  fill = 'white', size = geom_text_size, label.size = NA, alpha = 0.6, label.padding = unit(0.01, "lines"))+
    labs(y = expression(paste('Predicted Above-ground Biomass (g 900', m^-2, ')')), x = expression(paste('Observed Above-ground Biomass (g 900', m^-2, ')')), fill = 'Plant Functional Type')+
    theme_minimal()+
    theme(text = element_text(size=theme_text_size), 
          title = element_text(size=title_text_size), 
          plot.margin = margin(10, 10, 10, 10))+
    guides(color = 'none')+
    guides(fill = 'none')+
    scale_fill_manual(values = colors)+
    ggpmisc::stat_poly_eq(formula = textformula, geom = 'label', aes(label = paste(..rr.label.., sep = "~~~")), parse = TRUE, label.y = -Inf, label.x = Inf, hjust = 1.1, vjust = -1.3, color = 'grey20', fill = 'white', size = geom_text_size, label.size = NA, alpha = 0.6, label.padding = unit(0.01, "lines"))
  
  RMSE_corrected
  
  outName = paste0('glmmLasso_preds_', file_name, '.png')
  
  ggsave(
    outName,
    RMSE_corrected,
    width = 40,
    height = 30,
    units = 'cm'
  )
  
  cat("\n")
  
}
