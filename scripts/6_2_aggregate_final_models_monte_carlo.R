######################################################################################################
######################################################################################################

# 1. GET FOLDERS AND DATA ------------------------------------------------------

# Set data type
data_type = 'field' # Choose 'field' or 'UAV'

# Set transformation
transform = 'sqrt' # Choose 'sqrt' or 'log'

# Set directory
dir = paste0('/scratch/kmo265/UAV_to_LS/results/', data_type, '/FINAL_MODELS_nrmse_nmbe_corr_avg_', transform)
print(paste0('The output directory is: ', dir))
setwd(dir)
cat("\n")

# Set output file name
outName = paste0('glmmLasso_coefficients_finalModels_', data_type, '_', transform, '_', 'nrmse_nmbe_corr_avg_MC_all.csv')
  
# Get folders in output directory
files = list.files(full.names = TRUE, recursive = FALSE, pattern = '*glmmLasso_coefficients_finalModels_MCiter*')

######################################################################################################
######################################################################################################

# 2. AGGREGATE ------------------------------------------------------

coefTotal = data.frame()

for(file in files){
  
  # Read in csvs
  coef = read.csv(file, header = T)

  # Add to aggregated dataframe
  coefTotal = rbind(coefTotal, coef)

}

######################################################################################################
######################################################################################################

# 3. TIDY ------------------------------------------------------

coefTidy = tidyr::complete(data = coefTotal, pft, MCiter) # Turn implicit missing values into explicit missing values

coefTidy$predictor[is.na(coefTidy$predictor)] = '(Intercept)' # Missing predictor names are the intercepts

coefTidy$coef_orig[is.na(coefTidy$coef_orig)] = -999999999 # Convert NAs for use in Google Earth Engine

coefTidy$smear_factor[is.na(coefTidy$smear_factor)] = -999999999 # Convert NAs for use in Google Earth Engine

# 4. SAVE ------------------------------------------------------

write.csv(coefTidy, outName, row.names = FALSE)
