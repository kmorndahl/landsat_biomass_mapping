import os
import sys
import numpy
import numpy
import pandas

import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from scipy.cluster import hierarchy
from scipy.stats import spearmanr

# 1. ========== SET UP  ==========

# 1.1 ---------- Directories ----------

###########################################################
# SET OUTPUT DIRECTORY
###########################################################

outDir = os.path.join('/scratch/kmo265/UAV_to_LS/data/feature_correlation/')

os.chdir(outDir) # Sets the working directory

print('The output directory is: ' + outDir)
print('\n')

###########################################################

# Create training directory path ----------

trainingDir = os.path.join('/scratch/kmo265/1_UAV_to_LS_final/data/')

print('The training directory is: ' + trainingDir)
print('\n')

# 1.2 ---------- Set parameters ----------

# Choose data type ----------
data_type = 'field' # Choose 'field' or 'UAV'

# Select variables to exclude from modeling ----------
# These are variables for which the standard deviation for the training data was < 50% of the standard deviation across the full study area
exclVars = ["waterOccurrence",
"predicted_cover_bTREE",
"predicted_cover_talshr",
"landsatCcdc_ndsi_changeRatep080p095",
"phenology_firstSnowDoy",
"landsatCcdc_ndsi_seasonalp050",
"predicted_cover_cTREE",
"landsatCcdc_ndsi_seasonalp065",
"landsatCcdc_ndsi_seasonalp035",
"landsatCcdc_ndsi_median",
"landsatCcdc_ndsi_mean",
"landsatCcdc_ndsi_seasonalp080",
"landsatCcdc_ndsi_seasonalp095",
"landsatCcdc_ndsi_seasonalp020",
"landsatCcdc_ndsi_changeRatep065p080",
"landsatCcdc_ndsi_changeRatep050p065",
"landsatCcdc_ndmi_changeRatep080p095",
"permafrost_zonation_index",
"landsatCcdc_ndsi_seasonalp005",
"landsatCcdc_nbr_changeRatep080p095",
"landsatCcdc_ndsi_amplitude",
"topo_roughness",
"landsatCcdc_redcc_seasonalp065",
"landsatCcdc_greencc_seasonalp080",
"topo_surfaceAreaRatio",
"landsatCcdc_ndsi_changeRatep035p050",
"landsatCcdc_ndwi_changeRatep080p095",
"landsatCcdc_redcc_seasonalp080",
"landsatCcdc_ndwi_seasonalp035",
"landsatCcdc_ndvi_seasonalp095",
"landsatCcdc_ndwi_seasonalp050",
"landsatCcdc_ndvi_seasonalp080",
"landsatCcdc_ndwi_seasonalp065",
"landsatCcdc_ndwi_mean",
"landsatCcdc_ndwi_median",
"landsatCcdc_redcc_mean",
"landsatCcdc_ndvi_median",
"landsatCcdc_ndwi_changeRatep050p065",
"landsatCcdc_ndwi_seasonalp080",
"landsatCcdc_bluecc_changeRatep065p080",
"landsatCcdc_ndvi_mean",
"landsatCcdc_ndwi_seasonalp095",
"landsatCcdc_ndwi_seasonalp020",
"landsatCcdc_ndwi_changeRatep065p080",
"landsatCcdc_evi_changeRatep080p095",
"landsatCcdc_ndvi_seasonalp065",
"landsatCcdc_blue_seasonalp095",
"landsatCcdc_blue_seasonalp080",
"landsatCcdc_redcc_median",
"landsatCcdc_ndvi_seasonalp035",
"landsatCcdc_ndvi_seasonalp050",
"landsatCcdc_greencc_mean",
"landsatCcdc_redcc_seasonalp095",
"landsatCcdc_blue_median",
"climTerra_swi",
"landsatCcdc_tcwgd_seasonalp020",
"landsatCcdc_greencc_seasonalp065",
"landsatCcdc_ndvi_seasonalp020",
"landsatCcdc_tcwgd_seasonalp035",
"landsatCcdc_bluecc_seasonalp065",
"landsatCcdc_blue_seasonalp065",
"landsatCcdc_red_seasonalp080",
"landsatCcdc_bluecc_amplitude",
"landsatCcdc_green_changeRatep050p065",
"landsatCcdc_redcc_seasonalp050",
"landsatCcdc_ndmi_changeRatep005p020",
"landsatCcdc_bluecc_changeRatep050p065",
"landsatCcdc_blue_changeRatep080p095",
"landsatCcdc_blue_mean",
"landsatCcdc_evi_changeRatep005p020",
"landsatCcdc_red_seasonalp065",
"landsatCcdc_green_seasonalp080",
"landsatCcdc_green_median",
"landsatCcdc_green_seasonalp065",
"landsatCcdc_tcwgd_seasonalp050",
"landsatCcdc_bluecc_seasonalp080",
"landsatCcdc_ndwi_seasonalp005",
"landsatCcdc_red_median",
"landsatCcdc_tcb_seasonalp035",
"landsatCcdc_nbr_changeRatep005p020",
"landsatCcdc_blue_seasonalp050",
"landsatCcdc_green_seasonalp095"]

# Set predictor file ----------

if data_type == 'field':
    predictor_file = 'biomass_field_data_ls.csv'
elif data_type == 'UAV':
    predictor_file = 'biomass_UAV_data_ls.csv'
else:
    print('Data type not recognized, please enter a valid data type')
    sys.exit()

# 2. ========== DATA PREPARATION ==========

# 2.1 ---------- Load training points ----------

# Open training data
dfPath = os.path.join(trainingDir, predictor_file)

print('The predictor data path is: ' + dfPath)
print('\n')

df = pandas.read_csv(dfPath)

if data_type == 'UAV':

    print('UAV data, filtering by pixel percent...')
    print('\n')

    # Filter by pixel percent
    df = df[df['calval_pData'] >= 0.75]

# Report final full dataframe
print('The dataframe is:')
print(df)
print('\n')

# Get column indexes (multiple range slicing by label is not supported ...)
loc = df.columns.get_loc

# Select relevant predictors only
# Exclude categorical variables or variables that are not useful: sample_year, waterOccurrence, LC
# Include climate, landsat, phenology, topographic, permafrost, PFT cover
# Select columns using locations, add one because numpy slicing is right exclusive

X_full = df.iloc[:, numpy.r_[loc('climTerra_mjt'):loc('climTerra_tapmm')+1, loc('landsatCcdc_blue_amplitude'):loc('phenology_snowFreeDays')+1, loc('predicted_cover_DECIDUOUS_SHRUBS'):loc('predicted_cover_talshr')+1, loc('topo_elevation_m'):loc('topo_twi_100')+1]] 
print('PFT and climate included')
X_full['calval_pft'] = df['calval_pft']

# Exclude variables with small standard deviation over training data as compared to full study area
X_full = X_full.drop(columns=exclVars, errors='ignore')
print('The final predictors are:')
print(X_full.columns.tolist())
print('\n')

# Grab TOTAL only to avoid duplicate points
X = X_full[X_full['calval_pft'] == 'TOTAL']

# Remove pft column to prepare for creating correlation matrix
X = X.drop(columns=['calval_pft'])

# Report final dataframe
print('The final predictors dataframe is:')
print(X)
print('\n')

# 2.2 ---------- Scale data ----------

# Initialize scaler
scaler = StandardScaler()

# Scale dataset
numeric_cols = X.columns[X.dtypes.apply(lambda c: numpy.issubdtype(c, numpy.number))] # Get indices of numeric columns
X.loc[:, numeric_cols] = scaler.fit_transform(X.loc[:, numeric_cols]) # Scale and replace numeric columns only

print('The scaled predictor variables are:')
print(X)
print('\n')

# 3. ========== CREATE CORRELATION MATRIX ==========

corr_matrix = spearmanr(X).correlation

print('correlation matrix created:')
print(corr_matrix)
print('\n')

# Save correlation matrix
corr_matrix = pandas.DataFrame(corr_matrix, index = X.columns, columns = X.columns)
corr_matrix.to_csv('corr_matrix.csv')

# 4. ========== PERFORM HIERARCHICAL AGGLOMERATIVE CLUSTERING ==========
# https://scikit-learn.org/stable/auto_examples/inspection/plot_permutation_importance_multicollinear.html
# https://joernhees.de/blog/2015/08/26/scipy-hierarchical-clustering-and-dendrogram-tutorial/
# https://predictivehacks.com/hierarchical-clustering-in-python/

# 4.1 ---------- Create linkage ----------

# Calculate linkage
corr_linkage = hierarchy.ward(corr_matrix)
print('Linkages calculated:')
print(corr_linkage)
print('\n')

# 4.2 ---------- Create dendrogram ----------

# Create dendrogram
plt.figure(figsize=(30, 15))
plt.yticks(numpy.arange(0, 71, 1))
dendro = hierarchy.dendrogram(corr_linkage, labels=corr_matrix.columns, orientation = 'top', leaf_rotation=90)
print('Dendrogram created')
print('\n')

# Add potential cut points
for c in range(1, 70):
    plt.axhline(y = c, linestyle='--', linewidth = 1, c = 'black')

# Save dendrogram
plt.savefig('dendrogram_corr_ward_cut.png', format = 'png', dpi = 200, bbox_inches = 'tight')
print('Dendrogram saved')
print('\n')

# 4.3 ---------- Cluster ----------

# Get cluster IDs
clusters = pandas.DataFrame(corr_matrix.columns, columns = ['predictor'])

# Select cut point start and stop values
cut_nums = list(range(0, 70))

print('The cut numbers to consider are: ' + str(cut_nums))

for cut in cut_nums:
    clusters['clusterID_cut' + str(cut)] = hierarchy.fcluster(corr_linkage, cut, criterion = 'distance') # For each cut number, assign the correct cluster ID number to each predictor variable

print('Cluster IDs:')
print(clusters)
print('\n')

# 5. ========== CALCULATE CORRELATION METRICS BETWEEN PREDICTORS AND RESPONSES ==========

# Loop through PFTs and calculate correlation between each PFTs biomass predictions and each predictor variable
PFTs = df['calval_pft'].unique()
print('The PFTs are:' + str(PFTs))
PFT_counter = 1
for PFT in PFTs:

    # Report PFT progress
    print('Current response variable: ' + PFT)
    print('Response ' + str(PFT_counter) + ' of ' + str(len(PFTs)))

    # Subset data
    df_subset_PFT = df[df['calval_pft'] == PFT] # Subset full dataframe to include only observations pertaining to the current PFT
    X_subset_PFT = X_full[X_full['calval_pft'] == PFT] # Subset predictors only dataframe to include only observations pertaining to the current PFT

    # Grab variables
    response = df_subset_PFT['calval_biomassG'] # Grab the response variable i.e. the biomass predictions
    predictors = X_subset_PFT.drop(columns=['calval_pft']) # Grab the predictors, dropping the PFT column

    # Create correlation matrix
    corr = spearmanr(predictors, response).correlation # Creates 2D matrix with length and width equal to the number of columns in 'X' + the number of columns in 'response'
    corr = numpy.nan_to_num(corr) # NaNs caused when predictor variable is all one value, replace with zero
    corr = corr[:, -1].tolist() # Grab the last column, which should represent the correlations between 'response' and each predictor in 'X' and convert to list
    corr.pop() # Remove the last value, is which is the response variable's correlation with itself (and should be 1)

    # Add correlations to 
    clusters['calval_' + PFT + '_corr'] = corr
    PFT_counter += 1

print('Response variable correlation scores -- wide format:')
print(clusters)
print('\n')

# Save cluster information -- wide format
clusters.to_csv('clusters_corr_ward_wide.csv')

# Convert to long format
response_vars = ['calval_' + s + '_corr' for s in PFTs]
print(response_vars)
response_stubs = ['corr' + s for s in PFTs]
print(response_stubs)
clusters = pandas.wide_to_long(df = clusters, stubnames = ['clusterID_cut'], i = ['predictor'] + response_vars, j = 'cut_num').reset_index()
print(clusters)
clusters.rename(columns={'clusterID_cut': 'clusterID'}, inplace=True)
print(clusters)
clusters = pandas.wide_to_long(df = clusters, stubnames = ['calval_'], i = ['predictor', 'cut_num', 'clusterID'], j = 'response', suffix = '(\d+|\w+)').reset_index()
print(clusters)
clusters.rename(columns={'calval_': 'spearman_correlation'}, inplace=True)
print(clusters)
clusters.replace('BiomassG_corr', 'BiomassG', regex=True, inplace=True)

print('Response variable correlation scores -- long format:')
print(clusters)
print('\n')

# Save cluster information -- long format
clusters.to_csv('clusters_corr_ward_long.csv')

# 6. ========== SELECT FEATURES ==========

# Initialize data frame
selected_features = pandas.DataFrame()

responsePFTs = clusters.response.unique()

for cut in cut_nums:
    for responsePFT in responsePFTs:

        # Subset and tidy data
        print('The current cut number is: ' + str(cut))
        feature_data = clusters[clusters['cut_num'] == cut] # Subset to include only the current cut number
        print('The current PFT response is: ' + str(responsePFT))
        feature_data = feature_data[feature_data['response'] == responsePFT] # Subset to include only the current cut number
        print('The starting dataframe is:')
        print(feature_data)

        iteration = 1
        while not feature_data.empty:

            print('Iteration: ' + str(iteration))

            # Get overall maximum correlation, add it to selected features dataframe
            max_corr = max(feature_data.spearman_correlation.tolist(), key = abs)
            max_corr_row = feature_data[feature_data.spearman_correlation == max_corr]
            selected_features = selected_features.append(max_corr_row)
            print('The current maximum correlation with the response variable is:')
            print(max_corr_row)

            # Get cluster ID from current maximum correlation, remove that cluster from the dataframe
            max_corr_clusterID = max_corr_row.clusterID.tolist()[0]
            feature_data = feature_data[feature_data['clusterID'] != max_corr_clusterID]
            print('The current maximum correlation feature belongs to cluster ID ' + str(max_corr_clusterID) + ', removing that cluster from the dataframe: ')
            print(feature_data)

            # Get correlations between current maximum correlation feature and all other features
            current_feature_name = max_corr_row.predictor.unique()
            all_feature_names = feature_data.predictor.unique()
            check_correlations = corr_matrix.loc[all_feature_names, current_feature_name].abs()

            # Get all correlations above 0.7 and remove
            check_correlations = check_correlations[check_correlations[current_feature_name[0]] > 0.7]
            features_to_remove = check_correlations.index.values
            feature_data = feature_data[~feature_data['predictor'].isin(features_to_remove)]
            print('The current maximum correlation feature has high correlation with these other features: ' + str(features_to_remove) + '. Removing these features from the dataframe:')
            print(feature_data)

            iteration+=1

print('Features selected:')
print(selected_features)
print('\n')

# Save selected features -- long format
selected_features.to_csv('selected_features_corr_ward_long.csv')

# 7. ========== ADD CUT POINT SELECTION DATA ==========

# Get response variables
response_vars = selected_features['response'].unique()

# Set up dataframe
corr_scores = []

for var in response_vars:

    # Filter to response variable
    print('The current response variable is: ' + str(var))
    response_data = selected_features[selected_features['response'] == var] 

    for cut in cut_nums:

        # Filter to cut number
        print('The current cut number is: ' + str(cut))
        corr_data = response_data[response_data['cut_num'] == cut]

        # Get selected features
        features = corr_data['predictor'].unique()
        print('The number of features is: ' + str(len(features)))
       
        # Grab selected features from predictor dataframe
        X_sub = X[features]

        # Create correlation matrix
        corr_matrix_sub = spearmanr(X_sub).correlation
        corr_matrix_sub = pandas.DataFrame(corr_matrix_sub, index = X_sub.columns, columns = X_sub.columns)
        corr_matrix_sub_nans = corr_matrix_sub.mask(numpy.equal(*numpy.indices(corr_matrix_sub.shape))) # Mask out the values on the diagonal
        corr_matrix_sub_flat = corr_matrix_sub_nans.to_numpy().flatten() # Flatten to one dimensional array for ease of summary statistics
        corr_matrix_sub_flat = corr_matrix_sub_flat[~numpy.isnan(corr_matrix_sub_flat)] # Remove nans
        corr_matrix_sub_flat = numpy.absolute(corr_matrix_sub_flat) # Take absolute value

        # Calculate summary statistics
        corr_total = len(~numpy.isnan(corr_matrix_sub_flat))
        corr_avg = numpy.nanmean(corr_matrix_sub_flat)
        corr_median = numpy.nanmedian(corr_matrix_sub_flat)
        corr_count_gt05 = len(corr_matrix_sub_flat[corr_matrix_sub_flat > 0.5])
        corr_count_gt07 = len(corr_matrix_sub_flat[corr_matrix_sub_flat > 0.7])
        corr_count_gt09 = len(corr_matrix_sub_flat[corr_matrix_sub_flat > 0.9])

        # Append to list
        row = {'response': var, 'cut_num': cut, 'num_vars': corr_matrix_sub.shape[0], 'corr_total': corr_total, 'corr_avg': corr_avg, 'corr_median': corr_median, 'corr_count_gt05': corr_count_gt05, 'corr_percent_gt05': corr_count_gt05/corr_total, 'corr_count_gt07': corr_count_gt07, 'corr_percent_gt07': corr_count_gt07/corr_total, 'corr_count_gt09': corr_count_gt09, 'corr_percent_gt09': corr_count_gt09/corr_total}
        corr_scores.append(row)
        print('The data to append is: ' + str(row))
        print('\n')
  
corr_scores = pandas.DataFrame(corr_scores)
print('The correlation scores are:')
print(corr_scores)
print('\n')

# Save average correlation scores
corr_scores.to_csv('selected_features_corr_ward_cut_point_selection_criteria.csv')