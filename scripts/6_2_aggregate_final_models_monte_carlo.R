######################################################################################################
######################################################################################################

# CODE DESCRIPTION

# Aggregates coefficient data across Monte Carlo iterations

# NOTE: output directory structure not hosted at github

######################################################################################################
######################################################################################################

# 1. GET FOLDERS AND DATA ------------------------------------------------------

output_results = FALSE

# Set data type
data_type = 'field' # Choose 'field' or 'UAV'

# Set transformation
transform = 'sqrt' # Choose 'sqrt' or 'log'

# Set directory
dir = paste0('*/UAV_to_LS/results/', data_type, '/FINAL_MODELS_nrmse_nmbe_corr_avg_', transform)
print(paste0('The output directory is: ', dir))
cat("\n")

# Set output file name
outName = paste0('glmmLasso_coefficients_finalModels_', data_type, '_', transform, '_', 'nrmse_nmbe_corr_avg_MC_all.csv')
  
# Get files in output directory
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

if(output_results){write.csv(coefTidy, paste0(dir, outName), row.names = FALSE)}
