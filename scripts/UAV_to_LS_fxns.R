##############################################################################
#                                                                            #
# K.M.ORNDAHL 2021                                                           #
#                                                                            #
# FUNCTIONS TO SUPPORT GLMM/GLMNET/TIDYMODELS MODELING WITH CROSS VALIDATION #
#                                                                            #
##############################################################################

library(Metrics)
library(R.utils)

# ========== CROSS VALIDATION ==========

# Issue with groupKFold: https://github.com/topepo/caret/issues/1150
# Produces leave-one-group-out splits with the same formatting as groupKFold
# Length of list is k (in this case k = the number of groups)
# Each vector within the list denotes which row indexes are included in the TRAINING SET for that fold
logo <- function(group_list) {
  result <- list()
  fold_i <- 0
  for (left_out_group in unique(group_list)) {
    fold_i <- fold_i + 1
    result[[fold_i]] <- which(group_list != left_out_group, arr.ind = TRUE)
  }
  return(result)
}

# Converts a k fold group object to a foldid vector
# k fold group object:
# - Returned from groupKFold or logo functions
# - Produces list of vectors of length k where k = number of folds
# - Each vector within the list denotes which row indexes are included in the TRAINING SET for that fold
# foldid vector:
# - Required for cv.glmnet foldid parameter
# - Should be a vector the length of the input data frame
# - Values in the vector denote which fold the observation belongs to
# - 'Belonging' to a fold means the observation is part of the TEST SET for that fold
# - Each observation only belongs to ONE fold's test set
create_foldid = function(kFoldList){
  
  # Get maximum value from k fold split object
  # This represents the last index value of the data frame from which the k fold split object was derived i.e. the number of rows/observations in the original data frame
  n = max(unlist(kFoldList))
  
  # Initialize data frame to store fold ids
  foldid = data.frame(foldid = rep(0, n))
  
  # Initalize vector with row numbers
  idx = seq(1, n)
  
  # Initialize fold number
  i = 1
  
  # Loop through folds
  for(fold in kFoldList){
    
    # Get train and test row ids
    train_idx = fold # Get train row ids
    test_idx = idx[-train_idx] # Get test row ids
    
    # Where row id in data frame matches test set row id, assign the current fold number
    foldid$foldid[test_idx] = i
    
    i = i + 1
    
  }
  
  return(foldid)
  
}

# ========== FEATURE SELECTION ==========

# Takes a response variable and a list of predictors and constructs a formula
# y_name: name of response variable as string
# X_names: names of predictor variables as list of strings
# interactions: flag for which interactions should be added -- choose between 'no', 'lc' or 'all'
# transform: whether or not to transform the response variable -- choose 'log', 'sqrt' or 'none'
create_fixed_formula = function(y_name, X_names, interactions, transform){
  
  if(transform == 'log'){
    fixed = paste0('log(', y_name, ') ~')
  }else if(transform == 'sqrt'){
    fixed = paste0('sqrt(', y_name, ') ~')
  }else if(transform == 'none'){
    fixed = paste0(y_name, ' ~')
  }else{
    stop('Transformation indicator not recognized')
  }
  
  for(name in X_names){ # Loop through all predictors and add additive effects, using as.factor notation for categorical variables (in our case just LC)
    fixed = paste0(fixed, ' + ', name)
  }
  
  if(interactions == 'no'){
    # Do nothing for now
  } else if(interactions == 'lc'){
    # Loop through all predictors except LC and for each predictor add its interaction term with LC
    for(name in X_names[X_names != "LC"]){fixed = paste0(fixed, ' + LC:', name)} 
  } else if(interactions == 'all'){
    # Loop through all predictors and add all 2-way interactions
    combos = combn(X_names, 2) # Get all unique pairwise combinations
    for(i in 1:ncol(combos)){
      pair = combos[,i]
      fixed = paste0(fixed, ' + ', pair[1], ':', pair[2])
    }
  }
  
  fixed = gsub(' ~ + ', ' ~ ', fixed, fixed = TRUE) # Format first predictor correctly, get rid of extra '+'
  fixed = gsub('+ LC', '+ as.factor(LC)', fixed, fixed = TRUE) # Format LC correctly when additive or beginning of interactive
  fixed = gsub(':LC', ':as.factor(LC)', fixed, fixed = TRUE) # Format LC correctly when end of interactive
  fixed = gsub('LC + ', 'as.factor(LC) + ', fixed, fixed = TRUE) # Format LC correctly when first predictor
  
  fixed = as.formula(fixed) # Convert to formula
  
  environment(fixed) = globalenv() # Change environment to global R environment
  
  return(fixed)
  
}

# Calculates the maximum lambda value i.e. the lowest lambda that causes all predictors to be zero
# Calculating lambda max per 'Regularization Paths for Generalized Linear Models via Coordinate Descent' as described here:
# https://stats.stackexchange.com/questions/166630/glmnet-compute-maximal-lambda-value
# https://jerryfriedman.su.domains/ftp/glmnet.pdf
# data: full dataframe with response, predictors and random effects grouping variable (if applicable)
# fix: fixed effects formula
# response_var: response variable name, as a string
# transform: whether or not to transform the response variable -- choose 'log', 'sqrt' or 'none'
# scale: whether or not the predictors are scaled -- choose TRUE or FALSE
calculate_lambda_max = function(data, fix, response_var, transform, scale){
  
  alpha = 1 # alpha = 1 is LASSO 
  
  # Set up response
  y = data[,response_var]
  if(transform == 'log'){
    y = log(y)
  }else if(transform == 'sqrt'){
    y = sqrt(y)
  }else if(transform == 'none'){
    y = y
  }else{
    stop('Transformation indicator not recognized')
  }
  
  # Set up predictors 
  X = model.matrix(fix, data)
  if(scale){X = scale(X)}
  X[is.nan(X)] = 1
  
  lambda_max = (1/length(y))*max(abs(t(X)%*%y))
  
  if(transform == 'log'){
    lambda_max = exp(lambda_max)
  }else if(transform == 'sqrt'){
    lambda_max = (lambda_max)^2
  }else if(transform == 'none'){
    lambda_max = lambda_max
  }else{
    stop('Transformation indicator not recognized')
  }
  
  return(lambda_max)
  
}

# ========== CV.GLMMLASSO - LMMEN ==========
# https://github.com/yonicd/lmmen/blob/master/R/cv.glmmlasso.R
# https://www.rdocumentation.org/packages/lmmen/versions/1.0/topics/cv.glmmLasso
# https://raw.githubusercontent.com/cran/glmmLasso/master/demo/glmmLasso-soccer.r

#' @title Cross Validation for glmmLasso package
#' @description Cross Validation for glmmLasso package as shown in example xxx
#' @param dat data.frame, containing y,X,Z and subject variables
#' @param form.fixed formaula, fixed param formula, Default: NULL
#' @param form.rnd list, named list containing random effect formula, Default: NULL
#' @param lambda numeric, vector containing lasso penalty levels, Default: seq(500, 0, by = -5)
#' @param family family, family function that defines the distribution link of the glmm, Default: gaussian(link = "identity")
#' @return list of a fitted glmmLasso object and the cv BIC path
#' @examples
#' \dontrun{cv.glmmLasso(initialize_example(seed=1))}
#' @seealso 
#'  \code{\link[glmmLasso]{glmmLasso}}
#'  @references 
#'  Variable selection for generalized linear mixed models by ell 1-penalized estimation. 
#'  Statistics and Computing, pages 1-18, 2014.
#' @rdname cv.glmmLasso
#' @export 
#' @importFrom glmmLasso glmmLasso
#' @importFrom stats gaussian as.formula
cv.glmmLasso.lmmen=function(dat,
                            form.fixed=NULL,
                            form.rnd=NULL,
                            lambda=seq(500,0,by=-5),
                            switch.NR = FALSE,
                            final.re = FALSE,
                            family=stats::gaussian(link = "identity"),
                            timelimit = 900)
{
  
  # Convert to data frame if necessary
  if(inherits(dat,'matrix')) dat <- as.data.frame(dat)
  
  # Specifies size of Delta.start
  # d.size = # of fixed effects predictors + 1 for intercept + number of levels of random effect
  d.size = ncol(model.matrix(form.fixed, dat)) + length(levels(dat[,names(form.rnd)]))
  
  # Add column with unique identifier for each observation
  dat<-data.frame(subject=as.factor(row.names(dat)),dat,check.names = FALSE,row.names = NULL)
  
  # Create formulas if none given
  if(is.null(form.fixed)) form.fixed<-sprintf('y~%s',paste(grep('^X',names(dat),value = TRUE),collapse = '+'))
  if(is.null(form.rnd)) form.rnd<-eval(parse(text=sprintf('form.rnd<-list(subject=~1+%s)',paste(grep('^Z',names(dat),value = TRUE),collapse = '+'))))
  
  # Initialize vector to store BIC values for each lambda
  BIC_vec<-c()
  
  # Initialize vector to store lambda values
  lambda_vec<-c()
  
  # specify starting values for the very first fit; pay attention that Delta.start has suitable length! 
  Delta.start.base<-Delta.start<-as.matrix(t(rep(0,d.size)))
  Q.start.base<-Q.start<-0.1
  
  # Initialize list to store models for each lambda
  modList_lambdaJ <- list()
  
  # Loop through each value of lambda
  for(j in 1:length(lambda))
  {
    
    # Try to fit model
    suppressMessages({
      suppressWarnings({
        fn <- try(
            glmmLasso::glmmLasso(fix = stats::as.formula(form.fixed),
                                         rnd = form.rnd,
                                         data = dat,
                                         lambda = lambda[j],
                                         family = family,
                                         switch.NR = switch.NR,
                                         final.re = final.re,
                                         control = list(start=Delta.start[j,], q.start=Q.start[j]))
        ) # End try
      }) # End suppressWarnings      
    }) # End suppressMessages
    
    if(class(fn)!="try-error")
    { 
      # If modeling is successful...
      BIC_vec<-c(BIC_vec, fn$bic) # Store BIC
      lambda_vec<-c(lambda_vec, lambda[j]) # Store lambda
      Delta.start<-rbind(Delta.start,fn$Deltamatrix[fn$conv.step,]) # Get estimates of fixed and random effects from current model fit to use in next model fit
      Q.start<-c(Q.start,fn$Q_long[[fn$conv.step+1]]) # Get random effects variance-covariance parameters from current model fit to use in next model fit
      modList_lambdaJ = c(modList_lambdaJ, list(fn)) # Store model
      
    }else{
      Delta.start<-rbind(Delta.start,Delta.start.base)
      Q.start<-c(Q.start,Q.start.base)
    }
  }
  
  list(mods = modList_lambdaJ,BIC = BIC_vec, lambda = lambda_vec)
}

#' @title Predictions and RMSE for Cross Validation object from cv.glmmLasso.lmmen
#' @description Calculates predictions and RMSE for each model output from cv.glmmLasso.lmmen, will be used to identify best lambda value
#' @param fit_object fit object from cv.glmmLasso.lmmen, a list containing: $mods: a list of model fits for each value of lambda AND $BIC: a vector of BIC values for each model fit, for each value of lambda 
#' @param test test data set upon which to make predictions
#' @param response_var name of response variable
cv.glmmLasso.lmmen.predict = function(fit_object, test, response_var, transform){
  
  # Get models from model fit object
  models = fit_object$mods
  
  # Initialize list to store predictions for each model
  y_hat = list()
  rmse_mod = c()
  rmse_uncorrected = c()
  rmse_corrected = c()
  mae_mod = c()
  mae_uncorrected = c()
  mae_corrected = c()
  
  for(mod in models){
    
    # PREDICTIONS -----    
    
    # Get observed values of y from test set
    y_actual = test[,response_var]
    
    # Original predictions on transformed scale
    y_hat_mod = predict(mod, test)
    
    # Convert predictions to original scale and apply correction, if necessary
    # y_hat_uncorrected = uncorrected back-transformed predictions
    # y_hat_corrected = bias correction applied
    if(transform == 'log'){
      
      # Get Duan's smearing factor for log transform
      # https://www.jstor.org/stable/pdf/2288126.pdf?casa_token=QpztMyR1fj4AAAAA:lITwdyne-5DAxOoqcJB1qJUDcgmX6TjfHsPVpCKCJQzDYbp84ZZceXBw0nD-_Eaz6A7grNwjGvHaImLIW_po46K-W2gDeAs2Pe5R2mDuNb59IyldxuM5
      # https://www.sheffield.ac.uk/polopoly_fs/1.105589!/file/A_Jones_slides.pdf
      # http://139.70.23.11/people/newman_mc/pubs/Newman1993.pdf
      smear_factor = sum(exp(mod$y  - mod$y_hat))/length(mod$y)
      
      y_hat_uncorrected = exp(y_hat_mod) 
      
      y_hat_corrected = y_hat_uncorrected * smear_factor
      
    }else if(transform == 'sqrt'){
      
      # Get Duan's smearing factor for square root transform
      # https://www.jstor.org/stable/pdf/2288126.pdf?casa_token=QpztMyR1fj4AAAAA:lITwdyne-5DAxOoqcJB1qJUDcgmX6TjfHsPVpCKCJQzDYbp84ZZceXBw0nD-_Eaz6A7grNwjGvHaImLIW_po46K-W2gDeAs2Pe5R2mDuNb59IyldxuM5
      # https://www.sheffield.ac.uk/polopoly_fs/1.105589!/file/A_Jones_slides.pdf
      # http://139.70.23.11/people/newman_mc/pubs/Newman1993.pdf
      smear_factor = sum((mod$y  - mod$y_hat)^2)/length(mod$y)
      
      y_hat_uncorrected = y_hat_mod^2

      y_hat_corrected = y_hat_uncorrected + smear_factor
      
    }else if(transform == 'none'){
      
      y_hat_uncorrected = y_hat_mod
      y_hat_corrected = y_hat_mod
      
    }else{
      
      stop('Transformation indicator not recognized')
      
    }
    
    # Gather
    y_hat_all = data.frame(y_actual = y_actual, y_hat_mod = y_hat_mod, y_hat_uncorrected = y_hat_uncorrected, y_hat_corrected = y_hat_corrected)

    # Add to list of predictions
    y_hat = c(y_hat, list(y_hat_all)) # Store model

    # RMSE -----

    # Calcluate and store RMSE
    rmse_mod = c(rmse_mod, Metrics::rmse(y_actual, y_hat_mod))
    rmse_uncorrected = c(rmse_uncorrected, Metrics::rmse(y_actual, y_hat_uncorrected))
    rmse_corrected = c(rmse_corrected, Metrics::rmse(y_actual, y_hat_corrected))
    mae_mod = c(mae_mod, Metrics::mae(y_actual, y_hat_mod))
    mae_uncorrected = c(mae_uncorrected, Metrics::mae(y_actual, y_hat_uncorrected))
    mae_corrected = c(mae_corrected, Metrics::mae(y_actual, y_hat_corrected))
    
  }

  # Add predictions to model fit object
  fit_object$y_hat = y_hat
  fit_object$rmse_mod = rmse_mod
  fit_object$rmse_uncorrected = rmse_uncorrected
  fit_object$rmse_corrected = rmse_corrected
  fit_object$mae_mod = mae_mod
  fit_object$mae_uncorrected = mae_uncorrected
  fit_object$mae_corrected = mae_corrected

  return(fit_object)

}

