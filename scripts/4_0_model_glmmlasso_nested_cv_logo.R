######################################################################################################
######################################################################################################

# CODE DESCRIPTION

# This script performs the glmmLasso modeling, including leave-one-group-out nested cross-validation

# NOTE: output directory structure not hosted at github

######################################################################################################
######################################################################################################

# 1. SET UP ------------------------------

# 1.1 LIBRARIES ------------------------------

library(dplyr)
library(caret)
library(glmmLasso)
library(ggplot2)
source('scripts/UAV_to_LS_fxns.R')

params = commandArgs(trailingOnly = TRUE)

# 1.2 SET PARAMETERS ------------------------------

output_results = FALSE

# Set data type
data_type = 'field' # Choose 'field' or 'UAV'

# Set transformation
transform = 'sqrt' # Choose 'sqrt' or 'log'

# Set land cover
incl_LC = TRUE # Choose TRUE or FALSE 

# Set loss function(s) to use for choosing lambda
loss_functions = c('rmse_corrected') # Choose one or more of: 'BIC', 'rmse_uncorrected', 'rmse_corrected', 'mae_uncorrected', 'mae_corrected'

# Set number of inner CV folds
inner_k = 10

# Set response variable
response_var = 'calval_biomassG'

# Set grouping (random effects) variable
group_var = 'calval_siteCode'

# Get outer CV fold number from bash script
fold_num = as.integer(params[1])

# Set to not use scientific notation
options(scipen = 10)

# Set seed
set.seed(1908)

# 1.3 LOAD DATA ------------------------------

# Get file names
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

# Set training data directory
trainingDataDir = file.path('data/')
print(paste0('The training data directory is: ', trainingDataDir))
cat("\n")

# Load predictor data
predictor_data_path = file.path(trainingDataDir, predictor_file)
predictor_data = read.csv(predictor_data_path, header = T)

# Load predictor list
predictor_list_all_path = file.path(trainingDataDir, predictor_list_file)
predictor_list_all = read.csv(predictor_list_all_path, header = T)

# Load cut number list
cut_number_list_path = file.path(trainingDataDir, cut_number_list_file)
cut_number_list = read.csv(cut_number_list_path, header = T)

######################################################################################################
######################################################################################################

# 2. TIDY PREDICTOR DATA ------------------------------

# Remove zeros, convert to factors
predictor_data = predictor_data[predictor_data[,response_var] > 0,]
predictor_data = predictor_data %>% dplyr::mutate_if(is.character, as.factor)
predictor_data$LC = factor(predictor_data$LC)

# Reorder LC factor levels to choose appropriate reference group: "Sparsely Vegetated"
predictor_data$LC = relevel(predictor_data$LC, '10')

print('Data tidyed')
cat("\n")

######################################################################################################
######################################################################################################

# 3. GET CUT NUMBERS FROM CORRELATION DATA ------------------------------

# Get TOTAL only (can choose any PFT, just need to choose one so data is not duplicated)
cut_numbers = cut_number_list[cut_number_list$response == 'TOTAL_corr',]

# Discard any observations that have the same 'num_vars' as another observation
# When there are matching 'num_vars,' keeps the first observation
# Result is a list of all cut numbers that produce unique numbers of variables
# These are the cut numbers to use in modeling
cut_numbers = dplyr::distinct(cut_numbers, num_vars, .keep_all = TRUE)

# Get chosen cut numbers as a list
cut_numbers = cut_numbers$cut_num

# OPTIONAL: Remove low cut numbers to improve model fitting
cut_numbers = cut_numbers[4:length(cut_numbers)]

cat('The cut numbers are: ', cut_numbers)
cat("\n")

######################################################################################################
######################################################################################################

# LOOP THROUGH CUT NUMBERS ------------------------------

for(cut_number in cut_numbers){ ### START CUT NUMBER LOOP
  
  print(paste0('========== The cut number is: ', cut_number, ' =========='))
  cat("\n")
  
  #################################################################################################################################
  # SET OUTPUT DIRECTORY ######################################################################################################
  #################################################################################################################################
  
  if(output_results){
  
    # Create path name
    if(!incl_LC){LCtext = '_noLC'}else{LCtext = ''}
    outPath = paste0('*/UAV_to_LS/results/', data_type, '/', transform, '/cut', cut_number, LCtext)
    
    # Create directory if it does not already exist
    dir.create(outPath)
    
    # Set as current working directory
    setwd(outPath)
    
    print(paste0('The current output directory is: ', outPath))
    cat("\n")
  
  }
  
  #################################################################################################################################
  #################################################################################################################################
  #################################################################################################################################
  
  # Subset to desired cut number
  predictor_list = predictor_list_all
  if(cut_number == 'Inf'){
    predictor_list = data.frame()
  } else{
    cut_number = as.numeric(cut_number)
    predictor_list$cut_num = as.numeric(predictor_list$cut_num)
    predictor_list = predictor_list[predictor_list$cut_num == cut_number,]
  }
  
  print('Correlation data subset to cut number:')
  print(dim(predictor_list))
  cat("\n")
  
  ######################################################################################################
  ######################################################################################################
  
  # LOOP THROUGH PFTS ------------------------------
  
  predictions = data.frame()
  coefficients = data.frame()
  
  pfts = c('DECIDUOUS_SHRUBS', 'EVERGREEN_SHRUBS', 'FORBS', 'GRAMINOIDS', 'LICHENS', 'TOTAL')

  for(pft in pfts){ ### START PFT LOOP
    
    print(paste0('========== The plant functional type is: ', pft, ' =========='))
    cat("\n")
    
    ######################################################################################################
    ######################################################################################################
    
    # 4. TIDY AND SUBSET ------------------------------
    
    # 4.1 SUBSET BY PFT ------------------------------
    
    # Predictor data
    df = predictor_data[predictor_data$calval_pft == pft,]
    
    # Correlation
    selectedFeatures = predictor_list[predictor_list$response == paste0(pft, '_corr'),]
    
    print('Data subset by plant functional type:')
    print(dim(df))
    print(head(selectedFeatures))
    cat("\n")
    
    # 4.2 GET RESPONSE VARIABLE ------------------------------
    
    # Get response variable
    y = dplyr::select(df, all_of(response_var))
    
    print('The response variable is:')
    print(head(y))
    cat("\n")
    
    # 4.3 GET SELECTED PREDICTOR VARIABLES ------------------------------
    
    # Get predictor names
    selectedFeatures = unique(as.character(selectedFeatures$predictor))
    
    # Include LC if specified
    if(incl_LC){
      print('Land cover categorical predictor included')
      cat("\n")
      selectedFeatures = append(selectedFeatures, 'LC')
    }
    
    cat('The predictor names are: ', selectedFeatures, "\n")
    cat("\n")
    
    X = dplyr::select(df, all_of(selectedFeatures))
    
    print('Response and predictor data prepared:')
    print(dim(y))
    print(dim(X))
    cat("\n")
    
    # 4.4 GET GROUPING VARIABLE ------------------------------
    
    # Get grouping variable
    groups = dplyr::select(df, all_of(group_var))
    
    print('The grouping variable is:')
    print(head(groups))
    cat("\n")
    
    # 4.5 REMOVE ZERO AND  ZERO VARIANCE PREDICTORS ------------------------------
    
    nzv = caret::nearZeroVar(X, uniqueCut = 1) # Find near zero variance predictors, uniqueCut default = 10

    if(length(nzv) > 0){

      print('Predictors with zero or near zero variance removed from dataset:')
      print(names(X)[nzv])
      cat("\n")

      X = X[, -nzv] # Remove near zero variance predictors from dataset

      }

    X_names = names(X)
    
    # 4.6 AGGREGATE DATA ------------------------------
    
    data = cbind(y, X, groups)
    
    # 4.7 SET UP OUTER CROSS VALIDATION ------------------------------
    
    # Outer CV is for testing accuracy of model
    
    # caret::groupKFold and logo functions produce a list of vectors
    # Length of list is k
    # Each vector within the list denotes which row indexes are included in the TRAINING SET for that fold
    
    # Leave One Group Out
    # logo() function is stable, produces same groupings each time, provided same input
    foldsLOGO_outer = logo(data[,group_var]) # Leave One Group Out
    
    print(paste0('---------- OUTER fold number ', fold_num, ' of ', length(foldsLOGO_outer),' ----------'))
    cat("\n")
    
    # If current PFT has fewer folds than current fold number, skip to next PFT gracefully
    # Some PFTs have fewer folds because when zeros are removed some sites are removed
    if(fold_num > length(foldsLOGO_outer)){
      print('Current fold number greater than number of folds for this plant functional type. Folds for this plant functional type completed. Skipping to next plant functional type.')
      cat("\n")
      next
    }
    
    # Fold denotes observations that should be in the training set
    inTraining = foldsLOGO_outer[[fold_num]]
    
    # Partition data
    train = data[inTraining,]
    test = data[-inTraining,]
    
    # Get held out test site
    out_site = unique(test[,group_var])
    
    print(paste0('The test site is: ', out_site))
    cat("\n")
    
    print('Data partitioned into training and testing for OUTER CV')
    print('Training:')
    print(dim(train))
    print('Testing:')
    print(dim(test))
    cat("\n")
    
    # 4.8 PREPARE DATA FOR MODELING ------------------------------
    
    # Convert to factors
    train[,group_var] = factor(train[,group_var])
    test[,group_var] = factor(test[,group_var])
    
    if(incl_LC){
      train$LC = factor(train$LC)
      test$LC = factor(test$LC)
    }
    
    # 5. MODELING ------------------------------
    
    # 5.1 SET UP INNER CV ------------------------------
    
    # Inner CV is for choosing best lambda value
    
    # caret::groupKFold and logo functions produce a list of vectors
    # Length of list is k
    # Each vector within the list denotes which row indexes are included in the TRAINING SET for that fold
    
    # Grouped K-fold
    set.seed(1908)
    folds10_inner = caret::groupKFold(train[,group_var], k = inner_k) # Grouped 10 fold
    
    # Create foldID vector from k fold object
    foldid = create_foldid(folds10_inner)
    foldid = foldid$foldid
    
    print(paste0(inner_k, ' inner folds created'))
    cat("\n")
    
    # Initialize list to store model output for each fold
    modListInnerCV = list()
    
    # 5.2 LOOP THROUGH INNER CV FOLDS ------------------------------
    
    # 5.2.1 FIRST, loop through and get lambda max for each fold ------------------------------
    
    lambda_max_list = c()
    
    for(fold_name in names(folds10_inner)){
      
      fold = folds10_inner[[fold_name]]
      
      # Get training and testing data for INNER fold
      train_train = train[fold,]
      train_test = train[-fold,]
      
      # Convert to factors
      train_train[,group_var] = factor(train_train[,group_var])
      train_test[,group_var] = factor(train_test[,group_var])
      
      if(incl_LC){
        train_train$LC = factor(train_train$LC)
        train_test$LC = factor(train_test$LC)
      }
      
      # Get predictor names for modeling
      X_train_train = dplyr::select(train_train, -c(calval_biomassG, calval_siteCode))
      X_train_train_names = names(X_train_train)
      
      # Check near zero variance again for inner fold to avoid fit issues
      # For these inner folds, only remove variable if it has zero variance (because it will cause fit issues)
      zeroVarCut = 1/nrow(X_train_train) * 100
      
      nzv = caret::nearZeroVar(X_train_train, uniqueCut = zeroVarCut) # Find near zero variance predictors, uniqueCut default = 10
      
      if(length(nzv) > 0){
        
        print('Predictors with zero or near zero variance removed from inner fold predictor names list:')
        print(names(X_train_train)[nzv])
        cat("\n")
        
        X_train_train_names = X_train_train_names[-nzv] # Remove near zero variance predictors from predictor names list
        
      } 
      
      # Arrange fixed effects formula
      fixed_eff = create_fixed_formula(names(y), X_train_train_names , 'no', transform)
      
      # Calculate lambda max
      lambda_max = calculate_lambda_max(data = train_train, fix = fixed_eff, response_var = response_var, transform = transform, scale = TRUE)
      
      # Add to list of lambda max
      lambda_max_list = c(lambda_max_list, lambda_max)
      
    }
    
    # 5.2.2 Get overall lambda max ------------------------------
    
    lambda_max = max(lambda_max_list)
    
    cat('The lambda max values for each fold are: ', lambda_max_list, "\n")
    print(paste0('The overall lambda max is: ', lambda_max))
    cat("\n")
    
    # 5.2.3 Using overall lambda max, loop through each fold and fit CV models ------------------------------
    
    for(fold_name in names(folds10_inner)){
      
      print(paste0('========== The inner CV fold number is: ', fold_name, ' =========='))
      cat("\n")
      
      fold = folds10_inner[[fold_name]]
      
      # Get training and testing data for INNER fold
      train_train = train[fold,]
      train_test = train[-fold,]
      
      # Report training and testing sets
      print('Data partitioned into training and testing for INNER CV')
      print('Training:')
      print(dim(train_train))
      print('Testing:')
      print(dim(train_test))
      cat("\n")    
  
      # Convert to factors
      train_train[,group_var] = factor(train_train[,group_var])
      train_test[,group_var] = factor(train_test[,group_var])
      
      if(incl_LC){
        train_train$LC = factor(train_train$LC)
        train_test$LC = factor(train_test$LC)
      }
      
      # Get predictor names for modeling
      X_train_train = dplyr::select(train_train, -c(calval_biomassG, calval_siteCode))
      X_train_train_names = names(X_train_train)
      
      # Check near zero variance again for inner fold to avoid fit issues
      # For these inner folds, only remove variable if it has zero variance (because it will cause fit issues)
      zeroVarCut = 1/nrow(X_train_train) * 100
      
      nzv = caret::nearZeroVar(X_train_train, uniqueCut = zeroVarCut) # Find near zero variance predictors, uniqueCut default = 10
      
      if(length(nzv) > 0){
        
        print('Predictors with zero or near zero variance removed from inner fold predictor names list:')
        print(names(X_train_train)[nzv])
        cat("\n")
        
        X_train_train_names = X_train_train_names[-nzv] # Remove near zero variance predictors from predictor names list
        
      }
      
      # Arrange fixed effects formula
      fixed_eff = create_fixed_formula(names(y), X_train_train_names , 'no', transform)
      print(paste0('The fixed effects formula is: ', paste0(deparse(fixed_eff), collapse = '')))
      cat("\n")
      
      # Define random effects formula
      random_eff = list(calval_siteCode = ~1)
      print(paste0('The random effects formula is: ', names(random_eff), random_eff))
      cat("\n")
      
      # Define sequence of lambdas to try
      lambda_seq = c(exp(seq(log(lambda_max), log(1), length.out = 99)), 0) # log scale sequence including zero 0
      cat('The sequence of lambdas to try is: ', lambda_seq, "\n")
      cat("\n")
      
      # Fit models
      # Output is a list containing:
      # $mods: a list of model fits for each value of lambda
      # $BIC: a vector of BIC values for each model fit, for each value of lambda
      fit_gauss = cv.glmmLasso.lmmen(dat = train_train,
                                     form.fixed = fixed_eff,
                                     form.rnd = random_eff,
                                     lambda = lambda_seq,
                                     switch.NR = FALSE,
                                     final.re = FALSE,
                                     family = gaussian(link="identity"))
      
      if(length(fit_gauss$mods) == 0){
        print('No successful model fits for this fold, skipping to next fold...')
        cat("\n")
        next
      }
      
      print(paste0('cv.glmmLasso models fit, total number of models fit successfully: ', length(fit_gauss$mods)))
      cat("\n")
      print('lambda max model:')
      print(fit_gauss$mods[[1]])
      cat("\n")

      # Prepare data for predictions
      # If there is only one LC level present in the test dataset, it will cause issues in prediction i.e. "contrasts can be applied only to factors with 2 or more levels" error
      if(incl_LC & length(levels(train_test$LC)) == 1){
        
        current_level = levels(train_test$LC) # Get the current level
        
        if(current_level == "10"){
          new_levels = c(current_level, "9") # Add a dummy level
        }else{
          new_levels = c(current_level, "10") # Add the reference level
        }
        
        levels(train_test$LC) = new_levels # Assign this as the new levels of the test dataset -- it will still only have one unique value, but now it will have two levels: the original level and a "ghost" level that is the reference level
        train_test$LC = relevel(train_test$LC, '10') # Reorder LC factor levels to choose appropriate reference group: "Sparsely Vegetated"
        
      }

      if(incl_LC){

        # Check factor order
        cat('The training LC factors are: ', levels(train_train$LC))
        cat("\n")
        cat('The testing LC factors are: ', levels(train_test$LC))
        cat("\n")
        
      }
      
      # Make predictions
      predict_gauss = cv.glmmLasso.lmmen.predict(fit_gauss, train_test, names(y), transform)
      print('Predictions made from cv.glmmLasso models')
      cat("\n")
      
      # 5.3 SAVE CURRENT INNER CV FOLD MODEL AND PREDICTIONS ------------------------------
      
      # Save results from this fold
      modListInnerCV = c(modListInnerCV, list(predict_gauss)) # Store model
      
      print('Model results saved, inner CV fold completed, moving on to next fold...')
      cat("\n")
      
    }
    
    print('All inner CV folds completed, aggregating results from inner CV...')
    cat("\n")
    
    # 6. INNER CV MODEL ASSESSMENT ------------------------------
    
    # 6.1 GATHER INNER CV DATA ------------------------------
    
    y_hat_final = data.frame()
    metrics_final = data.frame()
    
    for(i in 1:length(modListInnerCV)){
      
      # Get results
      results_fold = modListInnerCV[[i]]
      
      # Get lambda
      lambda_fold = results_fold$lambda
      
      # Get predictions
      y_hat_fold = dplyr::bind_rows(results_fold$y_hat, .id = "model_number") # Get y_hat
      y_hat_fold$model_number = as.numeric(y_hat_fold$model_number) # Convert model numbers to numeric
      y_hat_fold$lambda = lambda_fold[y_hat_fold$model_number] # Look up lambdas
      
      # Get loss metrics
      metrics_fold = data.frame(results_fold[c(2:3, 5:length(results_fold))]) # Select metric items from predict_gauss by position and convert to data frame
      metrics_fold_long = data.frame(metrics_fold %>% tidyr::pivot_longer(cols = !lambda, names_to = 'metric', values_to = 'value')) # Convert to long format
      metric_fold_names = metrics_fold_long$metric[1:7] # Get metric names in corrected order
      metrics_fold_long$metric = factor(metrics_fold_long$metric, levels = metric_fold_names)
      
      # Get fold number
      cv = i
      
      # Add to aggregated database
      y_hat_final = dplyr::bind_rows(y_hat_final, data.frame(cv = cv, y_hat_fold))
      metrics_final = dplyr::bind_rows(metrics_final, data.frame(cv = cv, metrics_fold_long))
      
    }
    
    print('Results gathered from all inner CV folds')
    cat("\n")
    
    # 6.2 SUMMARIZE INNER CV DATA ------------------------------
    
    # For each lambda, each loss metric take average across all inner CV folds
    metrics_summary = metrics_final %>% 
      dplyr::group_by(lambda, metric) %>% 
      dplyr::summarise(value = mean(value))
    metrics_summary$value = round(metrics_summary$value, 0)
    
    # Get minimum for each metric
    metrics_best_lambda = data.frame(metrics_summary %>% 
                                       dplyr::group_by(metric) %>% 
                                       dplyr::slice_min(order_by = value, with_ties = TRUE))     # If there is a tie, reports all tied lambda values
    metrics_best_lambda = data.frame(metrics_best_lambda %>% 
                                       dplyr::group_by(metric) %>% 
                                       dplyr::slice_max(order_by = lambda, with_ties = FALSE)) # From list of ties, selects highest lambda
    
    # Graph lambda tracking
    # Values for each lambda are averaged across inner CV folds
    lambda_tracking = ggplot(metrics_summary, aes(x = log(lambda), y = value))+
      facet_wrap(~metric, scales = 'free')+
      geom_point()+
      geom_line()+
      theme_minimal()+
      labs(x = 'log(Lambda)', y = 'Loss value') 
    
    if(output_results){
      
      outName = paste0('lambdaCV_cutNum', cut_number, '_', pft, '_test', out_site, '.png')
      
      ggsave(
        outName,
        lambda_tracking,
        width = 40,
        height = 30,
        units = 'cm'
      )
    
    }
    
    print('Results summarized from all inner CV folds')
    cat("\n")
    
    # 6.3 LOOP THROUGH LOSS FUNCTIONS ------------------------------
    # Get best lambda for each loss function and fit a final model using that lambda
    
    for(loss_function in loss_functions){ ### START LOSS FUNCTION LOOP
      
      print(paste0('========== Fitting final model! The loss function is: ', loss_function, ' =========='))
      cat("\n")
      
      # 6.4 FIT FINAL MODEL WITH BEST LAMBDA ------------------------------
      
      # Get lambda value using user provided loss function
      best_lambda = metrics_best_lambda[metrics_best_lambda$metric == loss_function,]$lambda
      
      # Set up fixed effects with all available predictors
      fixed_eff = create_fixed_formula(names(y), X_names , 'no', transform)
      
      # Check model formulas
      print(paste0('The fixed effects formula is: ', paste0(deparse(fixed_eff), collapse = '')))
      cat("\n")
      print(paste0('The random effects formula is: ', names(random_eff), random_eff))
      cat("\n")
      
      # Fit glmmLasso model on full training set using chosen lambda
      fit_final = try(glmmLasso::glmmLasso(data = train,
                                           fix = fixed_eff,
                                           rnd = random_eff,
                                           lambda = best_lambda,
                                           switch.NR = FALSE,
                                           final.re = FALSE,
                                           family = gaussian(link="identity")))
      
      if(class(fit_final) != "glmmLasso"){
        print('Best lambda model fitting failed, moving on to next loss function...')
        cat("\n")
        next
      }
      
      print('Best model fit on full outer CV training set')
      cat("\n")
      
      # 6.5 MAKE PREDICTIONS ------------------------------
      # Make predictions on the outer CV test set using the final model
      
      # Prepare data for predictions
      # If there is only one LC level present in the test dataset, it will cause issues in prediction i.e. "contrasts can be applied only to factors with 2 or more levels" error
      if(incl_LC & length(levels(test$LC)) == 1){
        
        current_level = levels(test$LC) # Get the current level
        
        if(current_level == "10"){
          new_levels = c(current_level, "9") # Add a dummy level
        }else{
          new_levels = c(current_level, "10") # Add the reference level
        }
        
        levels(test$LC) = new_levels # Assign this as the new levels of the test dataset -- it will still only have one unique value, but now it will have two levels: the original level and a "ghost" level that is the reference level
        test$LC = relevel(test$LC, '10') # Reorder LC factor levels to choose appropriate reference group: "Sparsely Vegetated"
        
      }
      
      if(incl_LC){
        # Check factor order
        cat('The training LC factors are: ', levels(train$LC))
        cat("\n")
        cat('The testing LC factors are: ', levels(test$LC))
        cat("\n")
      }  
    
      # Make predictions
      y_hat_mod = predict(fit_final, test)
      
      # 6.6 APPLY BIAS CORRECTIONS TO PREDICTIONS ------------------------------
      
      # Convert predictions to original scale and apply correction, if necessary
      # y_hat_uncorrected = uncorrected back-transformed predictions
      # y_hat_corrected = bias correction applied
      if(transform == 'log'){
        
        # Get Duan's smearing factor for log transform
        # https://www.jstor.org/stable/pdf/2288126.pdf?casa_token=QpztMyR1fj4AAAAA:lITwdyne-5DAxOoqcJB1qJUDcgmX6TjfHsPVpCKCJQzDYbp84ZZceXBw0nD-_Eaz6A7grNwjGvHaImLIW_po46K-W2gDeAs2Pe5R2mDuNb59IyldxuM5
        # http://139.70.23.11/people/newman_mc/pubs/Newman1993.pdf
        smear_factor = sum(exp(fit_final$y  - fit_final$y_hat))/length(fit_final$y)
        
        y_hat_uncorrected = exp(y_hat_mod) 
        
        y_hat_corrected = y_hat_uncorrected * smear_factor
        
      }else if(transform == 'sqrt'){
        
        # Get Duan's smearing factor for square root transform
        # https://www.jstor.org/stable/pdf/2288126.pdf?casa_token=QpztMyR1fj4AAAAA:lITwdyne-5DAxOoqcJB1qJUDcgmX6TjfHsPVpCKCJQzDYbp84ZZceXBw0nD-_Eaz6A7grNwjGvHaImLIW_po46K-W2gDeAs2Pe5R2mDuNb59IyldxuM5
        # http://139.70.23.11/people/newman_mc/pubs/Newman1993.pdf
        smear_factor = sum((fit_final$y  - fit_final$y_hat)^2)/length(fit_final$y)
        
        y_hat_uncorrected = y_hat_mod^2
        
        y_hat_corrected = y_hat_uncorrected + smear_factor
        
      }else if(transform == 'none'){
        
        y_hat_uncorrected = y_hat_mod
        y_hat_corrected = y_hat_mod
        
      }else{
        
        stop('Transformation indicator not recognized')
        
      }

      # 6.7 SAVE FINAL PREDICTIONS ------------------------------
      
      predict_final = data.frame(y_actual = test[,names(y)], y_hat_mod = y_hat_mod, y_hat_uncorrected = y_hat_uncorrected, y_hat_corrected = y_hat_corrected, smear_factor = smear_factor)
      
      # Add predictions to data frame
      best.mod.pred = data.frame(cut_num = cut_number,  pft = pft, fold_num = fold_num, test_site = out_site, loss_function = loss_function, lambda = best_lambda, prediction_uncorrected = predict_final$y_hat_uncorrected, prediction_corrected = predict_final$y_hat_corrected, actual = predict_final$y_actual, smear_factor = predict_final$smear_factor)
      predictions = dplyr::bind_rows(predictions, best.mod.pred)
      
      print('Predictions aggregated and saved')
      cat("\n")
      
      # 6.8 SAVE COEFFICENTS ------------------------------
      
      # Get coefficients
      # This are on the unstandardized (original) scale of the predictors
      # glmmLasso standardizes for model fitting, then unstandardizes to present the coefficients
      coef_final = stack(coef(fit_final))
      names(coef_final) = c('coef_orig', 'predictor')
      
      # Get coefficients in a space that allows for variable importance comparisons, Agresti method
      # https://www.jstor.org/stable/pdf/2684719.pdf
      train_matrix = model.matrix(fixed_eff, train)
      sds = apply(train_matrix, 2, sd)
      coef_final$coef_std = coef_final$coef_orig * sds
      
      # Add coefficients to data frame
      best.mod.coef = data.frame(cut_num = cut_number,  pft = pft, fold_num = fold_num, test_site = out_site, loss_function = loss_function, lambda = best_lambda, predictor = coef_final$predictor, coef_orig = coef_final$coef_orig, coef_std = coef_final$coef_std)
      coefficients = dplyr::bind_rows(coefficients, best.mod.coef)
      
      print('Coefficients aggregated and saved')
      cat("\n")
      
    } ### END LOSS FUNCTION LOOP
    
    print('Results summarized from all loss functions')
    cat("\n")
    
  } ### END PFT LOOP
  
  if(output_results){
    
    write.csv(predictions, paste0('glmmLasso_predictions_cutNum', cut_number, '_cv', fold_num, '.csv'), row.names = FALSE)
    write.csv(coefficients, paste0('glmmLasso_coefficients_cutNum', cut_number, '_cv', fold_num, '.csv'), row.names = FALSE)
    
    print('Final predictions and coefficients data frames saved as .csv')
    cat("\n")
  
  }
  
} ### END CUT NUMBER LOOP
