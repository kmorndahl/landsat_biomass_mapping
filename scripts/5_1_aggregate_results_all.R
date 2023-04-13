library(dplyr)
library(Metrics)
library(qdapTools)
library(ggplot2)

######################################################################################################
######################################################################################################

# 1. SET UP ------------------------------------------------------

# 1.0 SET PARAMETERS ------------------------------------------------------

# Set data type
data_type = 'field' # Choose 'field' or 'UAV'

# Set transformation
transform = 'sqrt' # Choose 'sqrt' or 'log'

# 1.1 GET FOLDERS ------------------------------------------------------

dir = paste0('/scratch/kmo265/UAV_to_LS/results/', data_type, '/', transform, '/')
print(paste0('The output directory is: ', dir))
setwd(dir)
cat("\n")

folders = list.files(dir, full.names = TRUE, recursive = FALSE, pattern = '*cut*')

# 1.1 GET DATA ------------------------------------------------------

if(data_type == 'field'){
  predictor_file = 'biomass_field_data_ls.csv'
  predictor_list_file = 'biomass_field_data_selected_features.csv'
  cut_number_list_file = 'biomass_field_data_cut_numbers.csv'
}else if(data_type == 'UAV'){
  predictor_file = 'biomass_UAV_data_ls.csv'
  predictor_list_file = 'biomass_UAV_data_selected_features.csv'
  cut_number_list_file = 'biomass_UAV_data_cut_numbers.csv'
}else{
  stop('Data type not recognized, please enter a valid data type')
}

# Get feature correlation data
corr_path = file.path('/scratch/kmo265/1_UAV_to_LS_final/data/', cut_number_list_file)

corr = read.csv(corr_path)

######################################################################################################
######################################################################################################

# 2. AGGREGATE FILES ------------------------------------------------------

coefTotal = data.frame()
predTotal = data.frame()

for(folder in folders){
  
  # Get aggregated coefficient and predictor tables ('_ALL')
  coef_path = list.files(folder, full.names = TRUE, pattern = '^glmmLasso_coefficients_.*ALL.*\\.csv$')
  pred_path = list.files(folder, full.names = TRUE, pattern = '^glmmLasso_predictions_.*ALL.*\\.csv$')
  
  print(coef_path)
  print(pred_path)
  cat("\n")
  
  # Read in csvs
  coef = read.csv(coef_path, header = T)
  pred = read.csv(pred_path, header = T)

  # Indicate LC or no LC
  if(grepl('noLC', folder)){ # No LC
    coef$LC = 'NO'
    pred$LC = 'NO'
  } else{
    coef$LC = 'YES'
    pred$LC = 'YES'
  }
  
  # Add to aggregated dataframe
  coefTotal = rbind(coefTotal, coef)
  predTotal = rbind(predTotal, pred)
  
}

write.csv(coefTotal, 'glmmLasso_coefficients_ALL.csv', row.names = FALSE)
write.csv(predTotal, 'glmmLasso_predictions_ALL.csv', row.names = FALSE)

######################################################################################################
######################################################################################################

# 3. TIDY ------------------------------------------------------

pfts = c("DECIDUOUS_SHRUBS", "EVERGREEN_SHRUBS", "FORBS", "GRAMINOIDS", "LICHENS", "TOTAL")

coefTotalTidy = coefTotal
predTotalTidy = predTotal

# Assign factors
coefTotalTidy$pft = as.factor(coefTotalTidy$pft)
predTotalTidy$pft = as.factor(predTotalTidy$pft)
coefTotalTidy$LC = as.factor(coefTotalTidy$LC)
predTotalTidy$LC = as.factor(predTotalTidy$LC)

# Join correlation data
corr = subset(corr, select=-c(X))
corr$response = gsub('_corr', '', corr$response)
predTotalTidy = dplyr::left_join(predTotalTidy, corr, by = c('pft' = 'response', 'cut_num'))

# Account for cut number = Inf i.e. only cover predictors, and LC if applicable
predTotalTidy$num_vars[predTotalTidy$cut_num == 'Inf'] = 0

# Remove lower cut numbers, for some PFTs they are subject to stalling out, non-convergence etc.
predTotalTidy = predTotalTidy[predTotalTidy$cut_num >= 4,]

# For each PFT, cut number, LC status, count the total number of predictions made
predCounts = predTotalTidy %>% dplyr::count(pft, cut_num, LC)

# Add the prediction counts to the full dataframe
predTotalTidy = dplyr::left_join(predTotalTidy, predCounts, by = c('pft', 'cut_num', 'LC'))

# Remove instances where the number of predictions is less than the total number of predictions expected for each PFT
for(pft in pfts){
  
  max_n = max(predTotalTidy[predTotalTidy$pft == pft,]$n)
  predTotalTidy = predTotalTidy[!(predTotalTidy$pft == pft & predTotalTidy$n < max_n),]
  
}

######################################################################################################
######################################################################################################

# 4. CALCULATE EVALUATION METRICS ------------------------------------------------------

# 4.1 INDIVIDUAL EVALUATION METRICS ------------------------------------------------------

# Calculate grouped rmse and relative rmse
predTotalTidy = data.frame(predTotalTidy %>% 
  dplyr::group_by(pft, cut_num, loss_function, LC) %>%
  dplyr::mutate(rmse_corrected = Metrics::rmse(actual, prediction_corrected), nrmse_corrected = mean(rmse_corrected)/mean(actual)))

# Calculate grouped mean bias error (mbe), normalized mean bias error
predTotalTidy = data.frame(predTotalTidy %>% 
  dplyr::group_by(pft, cut_num, loss_function, LC) %>%
  dplyr::mutate(mbe_corrected = mean(prediction_corrected - actual), nmbe_corrected = mean(mbe_corrected)/mean(actual)))

# Calculate correlation between actual and predicted
predTotalTidy = data.frame(predTotalTidy %>% 
  dplyr::group_by(pft, cut_num, loss_function, LC) %>%
  dplyr::mutate(corr_corrected = cor(prediction_corrected, actual)))

# 4.2 SUMMARIZE  ------------------------------------------------------

# Group and summarise by PFT
# No actual averaging done here, just gets dataframe down to the right size (gets rid of test site column)
# rmse/mbe values should be the same across pft/cut_num/loss_function/LC combos
predTotalSummaryPFT = data.frame(predTotalTidy %>% 
                                   dplyr::group_by(pft, cut_num, loss_function, LC, num_vars) %>% 
                                   dplyr::summarise(rmse_corrected = mean(rmse_corrected), nrmse_corrected = mean(nrmse_corrected), mbe_corrected = mean(mbe_corrected), nmbe_corrected = mean(nmbe_corrected), corr_corrected = mean(corr_corrected)))

# 4.3 COMBINED EVALUATION METRICS  ------------------------------------------------------

# Get average of relative RMSE, mean normalized bias error, and correlation
# We need to convert correlation scores for use in average
# For average, we are trying to minimize the "error" score, so correlations need to be converted so that higher numbers are worse
# To do this: 1 - cor
predTotalSummaryPFT$nrmse_nmbe_corr_avg = (abs(predTotalSummaryPFT$nrmse_corrected) + abs(predTotalSummaryPFT$nmbe_corrected) + (1 - predTotalSummaryPFT$corr_corrected))/3

######################################################################################################
######################################################################################################

# 5. IDENTIFY BEST MODELS ------------------------------------------------------
# Minimize nrmse, nmbe and corr

minimums_nrmse_nmbe_corr_avg = data.frame()

# Choose the best model considering all predictor combinations, and with corrected RMSE as the loss function for choosing lambda in the nested cross validation
combos = list(c('all', 'rmse_corrected'))

for(combo in combos){
  
  # Set up
  LC = combo[1]
  loss_function = combo[2]
  
  data_temp = predTotalSummaryPFT
  
  # Subset data
  if(LC != 'all'){
    data_temp = data_temp[data_temp$LC == LC,]
    print('Data subset by LC')
    cat('\n')
  }

  if(loss_function != 'all'){
    data_temp = data_temp[data_temp$loss_function == loss_function,]
    print('Data subset by loss function')
    cat('\n')
  }

  # Get minimum
  min_nrmse_nmbe_corr_avg_temp = data_temp %>%
    dplyr::group_by(pft) %>%
    dplyr::slice(which.min(nrmse_nmbe_corr_avg))
  
  min_nrmse_nmbe_corr_avg_temp$cut_num_label = paste0('Cut: ', min_nrmse_nmbe_corr_avg_temp$cut_num)
  min_nrmse_nmbe_corr_avg_temp$num_vars_label = paste0('Vars: ', min_nrmse_nmbe_corr_avg_temp$num_vars)
  
  min_nrmse_nmbe_corr_avg_temp$LC_combo = LC
  min_nrmse_nmbe_corr_avg_temp$loss_function_combo = loss_function
  
  minimums_nrmse_nmbe_corr_avg = dplyr::bind_rows(minimums_nrmse_nmbe_corr_avg, min_nrmse_nmbe_corr_avg_temp)
  
}

# Save chosen parameters
write.csv(minimums_nrmse_nmbe_corr_avg, paste0('minimums_nrmse_nmbe_corr_avg_', transform, '.csv'), row.names = FALSE)

######################################################################################################
######################################################################################################

# 6. PLOT BEST MODELS ------------------------------------------------------

title_text_size = 26
theme_text_size = title_text_size - 2
geom_text_size = theme_text_size / 4

responseLookup = data.frame(response = c("DECIDUOUS_SHRUBS", "EVERGREEN_SHRUBS", "FORBS", "GRAMINOIDS", "LICHENS", "TOTAL"), color = c('green4', 'cyan4', 'chartreuse3', 'darkseagreen3', 'yellow3', 'grey'))

response_vars = as.character(levels(as.factor(predTotalSummaryPFT$pft)))
colors = response_vars %l% responseLookup

for(combo in combos){
  
  # Set up
  LC = combo[1]
  loss_function = combo[2]
  
  data_plot = predTotalSummaryPFT
  
  # Subset data
  if(LC != 'all'){
    data_plot = data_plot[data_plot$LC == LC,]
    print('Data subset by LC')
  }
  
  if(loss_function != 'all'){
    data_plot = data_plot[data_plot$loss_function == loss_function,]
    print('Data subset by loss function')
  }
  
  # Get minimums
  min_nrmse_nmbe_corr_avg_plot = minimums_nrmse_nmbe_corr_avg[minimums_nrmse_nmbe_corr_avg$LC_combo == LC & minimums_nrmse_nmbe_corr_avg$loss_function_combo == loss_function,]
  
  # Plot
  plt = ggplot(data_plot, aes(y = nrmse_nmbe_corr_avg, x = num_vars, col = LC, shape = loss_function))+
    geom_point(size = 3)+
    geom_line(size = 1)+
    facet_wrap(~pft, scales = 'free')+
    geom_point(data = min_nrmse_nmbe_corr_avg_plot, color = 'black', size = 3)+
    geom_text(data = min_nrmse_nmbe_corr_avg_plot, aes(label = num_vars_label), hjust = -0.1, vjust = 1, color = 'black', size = 7)+
    geom_text(data = min_nrmse_nmbe_corr_avg_plot, aes(label = cut_num_label), hjust = -0.1, vjust = 0, color = 'black', size = 7)+
    theme_minimal()+
    scale_x_continuous(breaks = seq(0, 60, by = 5))+
    theme(text = element_text(size=theme_text_size),
          title = element_text(size=title_text_size))+
    labs(x = 'Number of variables', y = 'Relative RMSE/NMBE/Correlation Average', col = 'LC', shape = 'Loss Function')
  
  outName = 'glmmLasso_nrmse_nmbe_corr_avg_comparison.png'

  print(outName)
  cat('\n')
  
  ggsave(
    outName,
    plt,
    width = 40,
    height = 30,
    units = 'cm'
  )  
  
}

######################################################################################################
######################################################################################################

# 7. SAVE BEST MODELS ------------------------------------------------------

for(combo in combos){
  
  # Set up
  LC = combo[1]
  loss_function = combo[2]
  
  # Initialize output data frame
  bestPreds_nrmse_nmbe_corr_avg = data.frame()
  bestCoefs_nrmse_nmbe_corr_avg = data.frame()
  
  # Get minimums
  min_nrmse_nmbe_corr_avg_save = minimums_nrmse_nmbe_corr_avg[minimums_nrmse_nmbe_corr_avg$LC_combo == LC & minimums_nrmse_nmbe_corr_avg$loss_function_combo == loss_function,]
  
  # Get best predictions for each PFT
  for(pft in pfts){
    
    # Get PFT specific minimums
    min_nrmse_nmbe_corr_avg_pft = min_nrmse_nmbe_corr_avg_save[min_nrmse_nmbe_corr_avg_save$pft == pft,]
    
    # Save best predictions/coefficients
    preds_nrmse_nmbe_corr_avg_pft = predTotal[predTotal$pft == pft & predTotal$cut_num == min_nrmse_nmbe_corr_avg_pft$cut_num & predTotal$LC == min_nrmse_nmbe_corr_avg_pft$LC & predTotal$loss_function == min_nrmse_nmbe_corr_avg_pft$loss_function,]
    bestPreds_nrmse_nmbe_corr_avg = dplyr::bind_rows(bestPreds_nrmse_nmbe_corr_avg, preds_nrmse_nmbe_corr_avg_pft)
    
    coefs_nrmse_nmbe_corr_avg_pft = coefTotal[coefTotal$pft == pft & coefTotal$cut_num == min_nrmse_nmbe_corr_avg_pft$cut_num & coefTotal$LC == min_nrmse_nmbe_corr_avg_pft$LC & coefTotal$loss_function == min_nrmse_nmbe_corr_avg_pft$loss_function,]
    bestCoefs_nrmse_nmbe_corr_avg = dplyr::bind_rows(bestCoefs_nrmse_nmbe_corr_avg, coefs_nrmse_nmbe_corr_avg_pft)
      
  }
  
  # Save combo information
  bestPreds_nrmse_nmbe_corr_avg$LC_combo = LC
  bestPreds_nrmse_nmbe_corr_avg$loss_function_combo = loss_function
  bestCoefs_nrmse_nmbe_corr_avg$LC_combo = LC
  bestCoefs_nrmse_nmbe_corr_avg$loss_function_combo = loss_function

  # Save best predictions/coefficients
  outPred_nrmse_nmbe_corr_avg = 'glmmLasso_BESTpredictions_nrmse_nmbe_corr_avg.csv'
  outCoef_nrmse_nmbe_corr_avg = 'glmmLasso_BESTcoefficients_nrmse_nmbe_corr_avg.csv'
  
  print(outPred_nrmse_nmbe_corr_avg)
  print(outCoef_nrmse_nmbe_corr_avg)
  cat('\n')
  
  write.csv(bestPreds_nrmse_nmbe_corr_avg, outPred_nrmse_nmbe_corr_avg, row.names = FALSE)
  write.csv(bestCoefs_nrmse_nmbe_corr_avg, outCoef_nrmse_nmbe_corr_avg, row.names = FALSE)  
  
}
