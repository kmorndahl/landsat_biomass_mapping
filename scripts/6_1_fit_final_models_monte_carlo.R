######################################################################################################
######################################################################################################

# CODE DESCRIPTION

# Calculates Monte Carlo permutations for best fitting model for each data type and transform combination

# NOTE: output directory structure not hosted at github

######################################################################################################
######################################################################################################

# 0. SET OUTPUT DIRECTORY

output_results = FALSE

# Set data type
data_type = 'field' # Choose 'field' or 'UAV'

# Set transformation
transform = 'sqrt' # Choose 'sqrt' or 'log'

# Set directory
dir = paste0('*/UAV_to_LS/results/', data_type, '/FINAL_MODELS_nrmse_nmbe_corr_avg_', transform)
print(paste0('The output directory is: ', dir))
setwd(dir)
cat("\n")

######################################################################################################
######################################################################################################

# 1. SET UP ------------------------------

# 1.1 LIBRARIES ------------------------------

library(dplyr)
library(caret)
library(glmmLasso)
library(Metrics)
source('scripts/UAV_to_LS_fxns.R')

# 1.2 SET PARAMETERS ------------------------------

# Choose which loss function to use for final model selection
final_loss_function = 'rmse_nmbe_corr_avg'

# Set loss function to use for picking lambda
loss_function = 'rmse_corrected'

# Set number of inner CV folds
inner_k = 10

# Set number of monte carlo iterations
mc_iter = 100

# Set response variable
response_var = 'calval_biomassG'

# Set grouping (random effects) variable
group_var = 'calval_siteCode'

# Set to not use scientific notation
options(scipen = 10)

# Set PFTs
pfts = c('DECIDUOUS_SHRUBS', 'EVERGREEN_SHRUBS', 'FORBS', 'GRAMINOIDS', 'LICHENS', 'TOTAL')

# 1.3 READ IN DATA ------------------------------

dataDir = paste0('*/UAV_to_LS/results/', data_type, '/', transform, '/')

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

# Get predictors
predictors = read.csv(file.path('data/', predictor_file), header = T)

# Get feature correlation data
corr_all = read.csv(file.path('data/', predictor_list_file), header = T)

# Get best model parameters

if(transform == 'log'){
  mins = read.csv(paste0(dataDir, 'minimums_', final_loss_function, '_log.csv'))
  mins$transform = 'log'
}else if(transform == 'sqrt'){
  mins = read.csv(paste0(dataDir, 'minimums_', final_loss_function, '_sqrt.csv'))
  mins$transform = 'sqrt'
}else{
  stop('Transform not recognized, please enter a valid transforme')
}

# 1.4 GET OVERALL MINIMUMS i.e. FINAL PARAMETERS ------------------------------

# Get correct LC/loss function combo
# Should only be one combo, but select to be sure
mins = mins[mins$LC_combo == 'all' & mins$loss_function_combo == 'rmse_corrected',]

# Select the proper parameters
params = dplyr::select(mins, c(pft, cut_num, LC, transform))

######################################################################################################
######################################################################################################

# 2. TIDY DATA ------------------------------

# 2.1 TIDY PREDICTOR DATA ------------------------------

# Remove zeros, convert to factors
predictors = predictors[predictors[,response_var] > 0,]
predictors = predictors %>% dplyr::mutate_if(is.character, as.factor)
predictors$LC = factor(predictors$LC)

# Reorder LC factor levels to choose appropriate reference group: "Sparsely Vegetated"
predictors$LC = relevel(predictors$LC, '10')

# 2.2 TIDY CORRELATION DATA ------------------------------

# Remove X column
corr_all = subset(corr_all, select=-c(X))

print('Data tidyed')
cat("\n")

######################################################################################################
######################################################################################################

# MONTE CARLO ------------------------------

for(k in 1:mc_iter){ ### START MONTE CARLO LOOP

  ######################################################################################################
  ######################################################################################################
  
  # 3. MONTE CARLO - PERMUTATE RESPONSE VARIABLE ------------------------------
  # Step 1: permute the response variable (biomass values) according to the standard error of either the input quadrats (field data) or the modeled, aggregated values (UAV data)
  # This captures uncertainty in sample selection (field data) or general uncertainty in the model (UAV data, although definitely an underestimate)
  
  # 3.0 INITIALIZE COEFFICIENTS DATAFRAME ------------------------------
  
  coefficients = data.frame()
  
  print(paste0('========== The Monte Carlo iteration is: ', k, ' =========='))
  cat("\n")
  
  # 3.1 CALCULATE STANDARD ERROR ------------------------------

  print('Original biomass values:')
  print(head(predictors$calval_biomassG))
  cat("\n")
  
  if(data_type == 'UAV'){
  
    # Get confidence intervals
    ci = dplyr::select(predictors, c('calval_biomassG_CI'))
    ci = as.numeric(ci[,1])
    
    print('The confidence intervals are:')
    print(head(ci))
    cat("\n")
  
    # For UAV data:  
    # We have 95% confidence intervals
    # Convert to standard error by dividing by 3.92

    # Convert confidence intervals to standard error
    se = ci/3.92
  
    print('The standard errors are:')
    print(head(se))
    cat("\n")
  
  }else if(data_type == 'field'){
    
    # Standard errors from quadrats used to calculate weighted means
    se = dplyr::select(predictors, c('calval_biomassGse'))
    se = as.numeric(se[,1])
    
    print('The standard errors are:')
    print(head(se))
    cat("\n")
    
  }
  
  # 3.2 PERMUTATION ------------------------------

  # Permute response variable
  set.seed(k) # We want Monte Carlo random functions to change with each Monte Carlo loop iteration
  permutation = rnorm(n = nrow(predictors), mean = 0, sd = 1) # Generate n random draws from standard normal distribution
  predictors$calval_biomassG_permuted = predictors$calval_biomassG + (se * permutation) # Permute
  
  # Ensure there are no negative or zero values (for log and sqrt transforms), assign a very small number instead
  predictors$calval_biomassG_permuted[predictors$calval_biomassG_permuted <= 0] = 0.01
  
  print('Biomass values permuted, the new biomass values are:')
  print(head(predictors$calval_biomassG_permuted))
  cat("\n")
  
  # Change response variable to select permuted biomass
  response_var = 'calval_biomassG_permuted'

  ######################################################################################################
  ######################################################################################################
  
  # LOOP THROUGH PFTs ------------------------------
  
  pfts = c('DECIDUOUS_SHRUBS', 'EVERGREEN_SHRUBS', 'FORBS', 'GRAMINOIDS', 'LICHENS', 'TOTAL')
  
  for(pft in pfts){ ### START PFT LOOP
  
    print(paste0('========== The pft is: ', pft, ' =========='))
    cat("\n")
    
    ######################################################################################################
    ######################################################################################################
    
    # 4. TIDY AND SUBSET ------------------------------
    
    # 4.1 GET PFT PARAMTERS ------------------------------
    
    pft_params = params[params$pft == pft,]
    cut_number = pft_params$cut_num
    LC = pft_params$LC
    transform = pft_params$transform

    # 4.2 SUBSET CORRELATION DATA BY CUT NUMBER ------------------------------
    
    corr = corr_all
    cut_number = as.numeric(cut_number)
    corr$cut_num = as.numeric(corr$cut_num)
    corr = corr[corr$cut_num == cut_number,]
    
    print('Correlation data subset to cut number:')
    print(dim(corr))
    cat("\n")
    
    # 4.3 SUBSET PREDICTOR AND CORRELATION DATA BY PFT ------------------------------
    
    # Predictors
    df = predictors[predictors$calval_pft == pft,]
    
    # Correlation
    selectedFeatures = corr[corr$response == paste0(pft, '_corr'),]
    
    print('Data subset by plant functional type:')
    print(dim(df))
    print(head(selectedFeatures))
    cat("\n")
    
    # 4.4 GET RESPONSE VARIABLE AND CI ------------------------------
    
    # Get response variable
    y = dplyr::select(df, all_of(response_var))
    
    print('The response variable is:')
    print(head(y))
    cat("\n")
    
    # 4.5 GET SELECTED PREDICTOR VARIABLES ------------------------------
    
    # Get predictor names
    selectedFeatures = unique(as.character(selectedFeatures$predictor))
    
    # Include LC if specified
    if(LC == 'YES'){
      print('Land cover categorical predictor included')
      cat("\n")
      selectedFeatures = append(selectedFeatures, 'LC')
    }
    
    cat('The predictor names are: ', selectedFeatures, "\n")
    cat("\n")
    
    # Select features
    X = dplyr::select(df, all_of(selectedFeatures))
  
    # Convert to factor
    if(LC == 'YES'){
      X$LC = factor(X$LC)
    }
    
    print('Response and predictor data prepared:')
    print(dim(y))
    print(dim(X))
    cat("\n")
    
    # 4.6 GET GROUPING VARIABLE ------------------------------
    
    # Get grouping variable
    groups = dplyr::select(df, all_of(group_var))
    
    print('The grouping variable is:')
    print(head(groups))
    cat("\n")
    
    # 4.7 REMOVE ZERO AND NEAR ZERO VARIANCE PREDICTORS ------------------------------
  
    nzv = caret::nearZeroVar(X, uniqueCut = 1) # Find near zero variance predictors, uniqueCut default = 10
    
    if(length(nzv) > 0){
      
      print('Predictors with zero or near zero variance removed from dataset:')
      print(names(X)[nzv])
      cat("\n")
      
      X = X[, -nzv] # Remove near zero variance predictors from dataset
      
    }
    
    X_names = names(X)
    
    # 4.7 AGGREGATE DATA ------------------------------
    
    data = cbind(y, X, groups)
    
    # 4.9 PREPARE DATA FOR MODELING ------------------------------
    
    # Convert to factors
    data[,group_var] = factor(data[,group_var])
    
    ######################################################################################################
    ######################################################################################################
    
    # 5. MONTE CARLO - SAMPLE WITH REPLACEMENT ------------------------------
    # Step 2: bootstrap the current sample with replacement
    # This captures uncertainty in the final biomass model parameters
    # https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-017-0437-y#Sec12
    # https://static-content.springer.com/esm/art%3A10.1186%2Fs12874-017-0437-y/MediaObjects/12874_2017_437_MOESM1_ESM.pdf
    
    set.seed(k) # We want Monte Carlo random functions to change with each Monte Carlo loop iteration
    data = dplyr::sample_n(data, nrow(data), replace = TRUE) # Draw n random samples with replacement from data set
  
    ######################################################################################################
    ######################################################################################################
    
    # 6. MODELING ------------------------------
    
    # 6.1 SET UP INNER CV ------------------------------
    
    # Inner CV is for choosing best lambda value
    
    # groupKFold and logo functions produce a list of vectors
    # Length of list is k
    # Each vector within the list denotes which row indexes are included in the TRAINING SET for that fold
    
    # Grouped K-fold
    set.seed(1908) # We want K-fold groups to stay the same
    folds10_inner = caret::groupKFold(data[,group_var], k = inner_k) # Grouped 10 fold
    
    # Create foldID vector from k fold object
    foldid = create_foldid(folds10_inner)
    foldid = foldid$foldid
    
    print(paste0(inner_k, ' inner folds created'))
    cat("\n")
    
    # Initialize list to store model output for each fold
    modListInnerCV = list()
    
    # 6.2 LOOP THROUGH INNER CV FOLDS ------------------------------
    
    # 6.2.1 FIRST, loop through and get lambda max for each fold ------------------------------
    
    lambda_max_list = c()
    
    for(fold_name in names(folds10_inner)){
      
      fold = folds10_inner[[fold_name]]
      
      # Get training and testing data for INNER fold
      train = data[fold,]
      test = data[-fold,]
      
      # Convert to factors
      train[,group_var] = factor(train[,group_var])
      test[,group_var] = factor(test[,group_var])
      
      if(LC == 'YES'){
        train$LC = factor(train$LC)
        test$LC = factor(test$LC)
      }
      
      # Arrange fixed effects formula
      fixed_eff = create_fixed_formula(names(y), X_names , 'no', transform)
      
      # Calculate lambda max
      lambda_max = calculate_lambda_max(data = train, fix = fixed_eff, response_var = response_var, transform = transform, scale = TRUE)
      
      # Add to list of lambda max
      lambda_max_list = c(lambda_max_list, lambda_max)
      
    }
    
    # 6.2.2 Get overall lambda max ------------------------------
    
    lambda_max = max(lambda_max_list)
    
    cat('The lambda max values for each fold are: ', lambda_max_list, "\n")
    print(paste0('The overall lambda max is: ', lambda_max))
    cat("\n")
    
    # 6.2.3 Using overall lambda max, loop through each fold and fit CV models ------------------------------
    
    for(fold_name in names(folds10_inner)){
      
      print(paste0('========== The inner CV fold number is: ', fold_name, ' =========='))
      cat("\n")
      
      fold = folds10_inner[[fold_name]]
      
      # Get training and testing data for INNER fold
      train = data[fold,]
      test = data[-fold,]
      
      # Convert to factors
      train[,group_var] = factor(train[,group_var])
      test[,group_var] = factor(test[,group_var])
      
      if(LC == 'YES'){
        
        train$LC = factor(train$LC)
        test$LC = factor(test$LC)

      }
      
      print('Data partitioned into training and testing for INNER CV')
      print('Training:')
      print(dim(train))
      print('Testing:')
      print(dim(test))
      cat("\n")
      
      # Arrange fixed effects formula
      fixed_eff = create_fixed_formula(names(y), X_names , 'no', transform)
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
      fit_gauss = cv.glmmLasso.lmmen(dat = train,
                                     form.fixed = fixed_eff,
                                     form.rnd = random_eff,
                                     lambda = lambda_seq,
                                     family = gaussian(link="identity"))
      
      if(length(fit_gauss$mods) == 0){
        print('No successful model fits for this fold, skipping to next fold...')
        cat("\n")
        next
      }
      
      print('cv.glmmLasso models fit, example model provided:')
      print(fit_gauss$mods[[1]])
      cat("\n")
      
      # Prepare data for predictions
      # If there is only one LC level present in the test dataset, it will cause issues in prediction i.e. "contrasts can be applied only to factors with 2 or more levels" error
      if(LC == 'YES' & length(levels(test$LC)) == 1){
        
        current_level = levels(test$LC) # Get the current level
        
        if(current_level == "10"){
          new_levels = c(current_level, "9") # Add a dummy level
        }else{
          new_levels = c(current_level, "10") # Add the reference level
        }
        
        levels(test$LC) = new_levels # Assign these as the new levels of the test data set -- it will still only have one unique value, but now it will have two levels: the original level and a "ghost" level that is the reference level
        test$LC = relevel(test$LC, '10') # Reorder LC factor levels to choose appropriate reference group: "Sparsely Vegetated"
        
      }
      
      if(LC == 'YES'){
        
        # Check factor order
        cat('The training LC factors are: ', levels(train$LC))
        cat("\n")
        cat('The testing LC factors are: ', levels(test$LC))
        cat("\n")
        
      }
      
      # Make predictions on test set
      predict_gauss = cv.glmmLasso.lmmen.predict(fit_gauss, test, names(y), transform)
      print('Predictions made from cv.glmmLasso models')
      cat("\n")
      
      # 6.3 SAVE CURRENT INNER CV FOLD MODEL AND PREDICTIONS ------------------------------
      
      # Save results from this fold
      modListInnerCV = c(modListInnerCV, list(predict_gauss)) # Store model
      
      print('Model results saved, inner CV fold completed, moving on to next fold...')
      cat("\n")
      
    }
    
    print('All inner CV folds completed, aggregating results from inner CV...')
    cat("\n")
    
    ######################################################################################################
    ######################################################################################################
    
    # 7. INNER CV MODEL ASSESSMENT - FINAL MODEL FIT ------------------------------
    
    # 7.1 GATHER INNER CV DATA ------------------------------
    
    metrics_final = data.frame()
    
    for(i in 1:length(modListInnerCV)){
      
      # Get results
      results_fold = modListInnerCV[[i]]
      
      # Get lambda
      lambda_fold = results_fold$lambda
      
      # Get loss metrics
      metrics_fold = data.frame(results_fold[c(2:3, 5:length(results_fold))]) # Select metric items from predict_gauss by position and convert to data frame
      metrics_fold_long = data.frame(metrics_fold %>% tidyr::pivot_longer(cols = !lambda, names_to = 'metric', values_to = 'value')) # Convert to long format
      metric_fold_names = metrics_fold_long$metric[1:7] # Get metric names in corrected order
      metrics_fold_long$metric = factor(metrics_fold_long$metric, levels = metric_fold_names)
      
      # Get fold number
      cv = i
      
      # Add to aggregated database
      metrics_final = dplyr::bind_rows(metrics_final, data.frame(cv = cv, metrics_fold_long))
      
    }
    
    print('Results gathered from all inner CV folds')
    cat("\n")
    
    # 7.2 SUMMARIZE INNER CV DATA ------------------------------
    
    # For each lambda, each loss metric, take average across all inner CV folds
    metrics_summary = metrics_final %>% dplyr::group_by(lambda, metric) %>% dplyr::summarise(value = mean(value))
    metrics_summary$value = round(metrics_summary$value, 0)
    
    # Get minimum for each metric
    metrics_best_lambda = data.frame(metrics_summary %>% dplyr::group_by(metric) %>% dplyr::slice_min(order_by = value, with_ties = TRUE)) # If there is a tie, reports all tied lambda values
    metrics_best_lambda = data.frame(metrics_best_lambda %>% dplyr::group_by(metric) %>% dplyr::slice_max(order_by = lambda, with_ties = FALSE)) # From list of ties, selects highest lambda
    
    print('Results summarized from all inner CV folds')
    cat("\n")
    
    # 7.3 SELECT LOSS FUNCTIONS ------------------------------
    # Get best lambda for selected loss function and fit a final model using that lambda
    
    print(paste0('========== The loss function is: ', loss_function, ' =========='))
    cat("\n")
    
    # 7.4 FIT FINAL MODEL WITH BEST LAMBDA ------------------------------
    
    # Get lambda value using user provided loss function
    best_lambda = metrics_best_lambda[metrics_best_lambda$metric == loss_function,]$lambda
    
    print(paste0('The best lambda value is: ', best_lambda))
    cat("\n")
    
    # Check model formulas
    print(paste0('The fixed effects formula is: ', paste0(deparse(fixed_eff), collapse = '')))
    cat("\n")
    print(paste0('The random effects formula is: ', names(random_eff), random_eff))
    cat("\n")
    
    if(LC == 'YES'){
      
      # Check factor order
      cat('The training LC factors are: ', levels(data$LC))
      cat("\n")
      
    }
    
    # Fit glmmLasso model on full training set using chosen lambda
    fit_final = try(glmmLasso::glmmLasso(data = data,
                                         fix = fixed_eff,
                                         rnd = random_eff,
                                         lambda = best_lambda,
                                         family = gaussian(link="identity")))

    if(class(fit_final) != "glmmLasso"){
      print('Best lambda model fitting failed, moving on to next loss function...')
      cat("\n")
      next
    }
    
    print('Best model fit on full outer CV training set')
    cat("\n")
    
    # 7.5 GET BIAS CORRECTION FACTOR ------------------------------
    
    y_hat_mod = fit_final$y_hat
    
    # Duan's smearing approach
    # https://www.jstor.org/stable/pdf/2288126.pdf?casa_token=QpztMyR1fj4AAAAA:lITwdyne-5DAxOoqcJB1qJUDcgmX6TjfHsPVpCKCJQzDYbp84ZZceXBw0nD-_Eaz6A7grNwjGvHaImLIW_po46K-W2gDeAs2Pe5R2mDuNb59IyldxuM5
    # http://139.70.23.11/people/newman_mc/pubs/Newman1993.pdf
    
    if(transform == 'log'){
      
      # Get Duan's smearing factor for log transform
      
      smear_factor = sum(exp(fit_final$y  - fit_final$y_hat))/length(fit_final$y)
      
      y_hat_uncorrected = exp(y_hat_mod) 
      
      y_hat_corrected = y_hat_uncorrected * smear_factor
      
    }else if(transform == 'sqrt'){
      
      # Get Duan's smearing factor for square root transform
      
      smear_factor = sum((fit_final$y  - fit_final$y_hat)^2)/length(fit_final$y)
      
      y_hat_uncorrected = y_hat_mod^2
      
      y_hat_corrected = y_hat_uncorrected + smear_factor
  
    }else if(transform == 'none'){
      
      y_hat_uncorrected = y_hat_mod
      
      y_hat_corrected = y_hat_mod
      
      print('No transformation applied, smearing factor not applicable')
      
    }else{
      
      stop('Transformation indicator not recognized')
      
    }
    
    # 7.6 GET PREDICTION ERROR (RMSE) ------------------------------
    
    y_actual = data[,response_var]
    
    prediction_error = Metrics::rmse(y_actual, y_hat_corrected)

    # 7.7 SAVE COEFFICENTS, SMEARING FACTOR, PREDICTION ERROR ------------------------------
    
    # Get coefficients
    # This are on the unstandardized (original) scale of the predictors
    # glmmLasso standardizes for model fitting, then unstandardizes to present the coefficients
    coef_final = stack(coef(fit_final))
    names(coef_final) = c('coef_orig', 'predictor')
    
    # Get coefficients in a space that allows for variable importance comparisons Agresti method
    # https://www.jstor.org/stable/pdf/2684719.pdf
    data_matrix = model.matrix(fixed_eff, data)
    sds = apply(data_matrix, 2, sd)
    coef_final$coef_std = coef_final$coef_orig * sds
    
    # Add coefficients, smearing factor, and MC iteration to data frame
    best.mod.coef = data.frame(loss_function = loss_function, transform = transform, cut_num = cut_number,  LC = LC, pft = pft, lambda = best_lambda, predictor = coef_final$predictor, coef_orig = coef_final$coef_orig, coef_std = coef_final$coef_std, smear_factor = smear_factor, MCiter = k, prediction_error = prediction_error)
    coefficients = dplyr::bind_rows(coefficients, best.mod.coef)
  
    print('Coefficients aggregated and saved')
    cat("\n")
  
  } ### END PFT LOOP

  if(output_results){write.csv(coefficients, paste0(dir, 'glmmLasso_coefficients_finalModels_MCiter', k, '.csv'), row.names = FALSE)}

} ### END MONTE CARLO LOOP