######################################################################################################
######################################################################################################

# SET OUTPUT DIRECTORY

dir = '//minim.hpc.nau.edu/scratch/kmo265/1_UAV_to_LS_final/data/'
setwd(dir)

outName = 'biomass_field_data_ls.csv'

######################################################################################################
######################################################################################################

library(tidyverse)
library(sf)

######################################################################################################
######################################################################################################

# 1. IMPORT AND TIDY DATA --------------------------------------------------------------------------------------------------------------------------------

# 1.1 READ IN DATA AND TIDY --------------------------------------------------------------------------------------------------------------------------------

biomass_field = read.csv('biomass_field_data_lat_long.csv')
biomass_UAV = read.csv('pixel_centroids_predictors.csv') # Not quite uav_cover_biomass_satellite_predictors_coverPredictedWide_1m75_30m75_UPDATE.csv, has a few more plotIds bc this is from when Matt sent data that is only filtered by 75% at the 30m stage, the other is filtered by 75% at both 1m and 30m stage

# Select necessary columns
biomass_UAV = biomass_UAV %>% 
  dplyr::select(calval_siteCode, calval_plotId, latitude, longitude)

# Add suffixes to primary and secondary dataframes for data management
names(biomass_field) = paste0(colnames(biomass_field), '.field')
names(biomass_UAV) = paste0(colnames(biomass_UAV), '.UAV')

# 1.2 CONVERT TO SPATIAL POINTS --------------------------------------------------------------------------------------------------------------------------------

# Convert to spatial points
biomass_field_sp = sf::st_as_sf(biomass_field, coords = c("longitude.field", "latitude.field"), crs = 4326)
biomass_UAV_sp = sf::st_as_sf(biomass_UAV, coords = c("longitude.UAV", "latitude.UAV"), crs = 4326)

######################################################################################################
######################################################################################################

# 2. SELECT LANDSAT PIXEL CENTROIDS --------------------------------------------------------------------------------------------------------------------------------
# We only want to use Landsat pixel centroids when the centroid is within 50m of the site center

# 2.1 SUMMARIZE FIELD AND UAV DATA -- ONE OBSERVATION PER QUADRAT/PIXEL CENTROID --------------------------------------------------------------------------------------------------------------------------------

# Ensure there is only one observation per Landsat pixel centroid for UAV data
if(length(unique(biomass_UAV_sp$calval_plotId.UAV)) == nrow(biomass_UAV_sp)){
  UAV_pts_sp = biomass_UAV_sp
}else{
  stop('Duplicate landsat pixel centroids in UAV data, check data and try again')
}

# Reduce field data set to one observation per quadrat
field_quads_sp = biomass_field_sp %>% 
  dplyr::group_by(site_code.field, quadrat_num.field) %>% 
  dplyr::summarise(pft.field = first(pft.field))
field_quads_sp = dplyr::select(field_quads_sp, c(-pft.field))

# 3.2 FILTER UAV POINTS -- ONLY POINTS WITHIN 50 m OF PLOT CENTER --------------------------------------------------------------------------------------------------------------------------------

# Get only 50m quadrats, these denote the 'center' of the site
field_center_sp = field_quads_sp[field_quads_sp$quadrat_num.field == '50m',]

# The 50m field data points are the primary data frame
# The UAV data points are the secondary data frame
# For each 50m field data point, join to it any UAV points that are within 50m (half the distance of the transect)
# This identifies all of the UAV pixels that are close enough to the plot center to include
pts_within_50m_sp = st_join(field_center_sp, UAV_pts_sp, join = st_is_within_distance, dist = 50, left = FALSE)

# Grab the selected UAV points and drop geometry
selected_pts_UAV_df = st_drop_geometry(dplyr::select(pts_within_50m_sp, c(calval_siteCode.UAV, calval_plotId.UAV)))

# 3.3 ASSOCIATE FIELD QUADRATS WITH EACH SELECTED UAV POINT --------------------------------------------------------------------------------------------------------------------------------

# Join all field quadrats to each selected UAV point
selected_pts_quads_df = dplyr::left_join(selected_pts_UAV_df, st_drop_geometry(field_quads_sp), by = c('calval_siteCode.UAV' = 'site_code.field'))

# Separate into quadrat and points data sets
# Each data set will contain one observation for each unique Landsat pixel centroid/quadrat combination
# Note that for each data set, some observations will appear to be repeated, this is okay and necessary for joining Landsat pixel centroid/quadrat combinations
# These should be the same length and can be used to calculate pairwise distances
selected_quads_rep_df = dplyr::select(selected_pts_quads_df, c(calval_siteCode.UAV, quadrat_num.field))
selected_pts_rep_df = dplyr::select(selected_pts_quads_df, c(calval_siteCode.UAV, calval_plotId.UAV))

# 3.4 CALCULATE PAIRWISE DISTANCES --------------------------------------------------------------------------------------------------------------------------------

# Get coordinates from spatial data points to allow geometry information to be joined to selected, repeated quadrats and points data sets from above
UAV_lat_long = cbind(data.frame(st_drop_geometry(UAV_pts_sp)), st_coordinates(UAV_pts_sp))
field_lat_long = cbind(data.frame(st_drop_geometry(field_quads_sp)), st_coordinates(field_quads_sp))

# Left join lat/long information to selected, repeated points and quadrats data frames and convert back to spatial points
biomass_selected_quads_rep_sp = sf::st_as_sf(dplyr::left_join(selected_quads_rep_df, field_lat_long, by = c('calval_siteCode.UAV'='site_code.field', 'quadrat_num.field'='quadrat_num.field')), coords = c('X', 'Y'), crs = 4326)
biomass_selected_pts_rep_sp = sf::st_as_sf(dplyr::left_join(selected_pts_rep_df, UAV_lat_long, by = c('calval_plotId.UAV'='calval_plotId.UAV')), coords = c('X', 'Y'), crs = 4326)

# Calculate pairwise distances
distances = st_distance(biomass_selected_quads_rep_sp, biomass_selected_pts_rep_sp, by_element = TRUE)

# Put together full distance data frame
selected_quads_pts_distances = cbind(selected_pts_rep_df, dplyr::select(selected_quads_rep_df, c(quadrat_num.field)), distances)

######################################################################################################
######################################################################################################

# 3. AGGREGATE FIELD DATA --------------------------------------------------------------------------------------------------------------------------------
# Use inverse distance weighted averages
# https://www.geo.fu-berlin.de/en/v/soga/Geodata-analysis/geostatistics/Inverse-Distance-Weighting/index.html

# 3.1 CALCULATE IDW WEIGHTS --------------------------------------------------------------------------------------------------------------------------------

selected_quads_pts_distances$weight_idw_1 = as.numeric(selected_quads_pts_distances$distances^-1)

# 3.2 ADD IN BIOMASS DATA --------------------------------------------------------------------------------------------------------------------------------

biomass_selected_quads_pts_distances = dplyr::left_join(selected_quads_pts_distances, dplyr::select(biomass_field, c(site_code.field, quadrat_num.field, pft.field, biomass_g_900m2.field)), by = c('calval_siteCode.UAV' = 'site_code.field', 'quadrat_num.field' = 'quadrat_num.field'))

# 3.3 CALCULATE IDW MEAN BIOMASS --------------------------------------------------------------------------------------------------------------------------------

biomass_selected_quads_pts_distances$biomass_weighted_idw_1 = biomass_selected_quads_pts_distances$biomass_g_900m2.field * biomass_selected_quads_pts_distances$weight_idw_1

biomass_selected_quads_pts_distances_means = biomass_selected_quads_pts_distances %>% 
  dplyr::group_by(calval_plotId.UAV, pft.field) %>% 
  dplyr::summarise(biomass_g_900m2_idw_1 = sum(biomass_weighted_idw_1)/sum(weight_idw_1))

# 3.4 CALCULATE IDW MEAN STANDARD ERROR --------------------------------------------------------------------------------------------------------------------------------
# http://seismo.berkeley.edu/~kirchner/Toolkits/Toolkit_12.pdf

se_df = dplyr::left_join(biomass_selected_quads_pts_distances, biomass_selected_quads_pts_distances_means, by = c('calval_plotId.UAV', 'pft.field'))

plotId_list = unique(se_df$calval_plotId.UAV)
pft_list = unique(se_df$pft.field)

idw_biomass_se = data.frame()

for(plotId in plotId_list){
  for(pft in pft_list){
    
    # Get subset dataframe
    df = se_df[se_df$calval_plotId.UAV == plotId & se_df$pft.field == pft,]
    
    # Make sure there is data
    if(nrow(df) == 0){next}

    # Calculate standard error
    se = sqrt( ( ( sum(df$weight_idw_1 * (df$biomass_g_900m2.field)^2) / sum(df$weight_idw_1) ) - (mean(df$biomass_g_900m2_idw_1))^2 ) * ( sum(df$weight_idw_1^2) /  ( (sum(df$weight_idw_1)^2) - sum(df$weight_idw_1^2) ) ) )
    
    # If there is only one quadrat, force to zero
    if(nrow(df) == 1){se = 0}

    # Add to dataframe
    df$idw_1_se = se
    
    # Add to aggregating dataframe
    idw_biomass_se = dplyr::bind_rows(idw_biomass_se, df)
    
  }
}

# Remove observations where only one quadrat was used to calculate pixel centroid average
# If there is only one quadrat and zero biomass recorded in that quadrat, the PFT is likely not present at the site -- these observations can stay
idw_biomass_se = idw_biomass_se[!(idw_biomass_se$biomass_g_900m2_idw_1 > 0 & idw_biomass_se$idw_1_se == 0),]

# Summarise by pixel centroid
biomass_selected_quads_pts_distances_means_se = idw_biomass_se %>% 
  dplyr::group_by(calval_plotId.UAV, pft.field) %>% 
  dplyr::summarise(biomass_g_900m2_idw_1 = mean(biomass_g_900m2_idw_1),
            idw_1_se = mean(idw_1_se))

# Rename
names(biomass_selected_quads_pts_distances_means_se) = c("calval_plotId", "calval_pft", "calval_biomassG", "calval_biomassGse")

# 3.5 JOIN PREDICTORS --------------------------------------------------------------------------------------------------------------------------------

# Join
biomass_selected_quads_pts_distances_final = dplyr::left_join(biomass_selected_quads_pts_distances_means_se, biomass_UAV, by = 'calval_plotId')

# Replace spaces with underscores
biomass_selected_quads_pts_distances_final$calval_pft = gsub(' ', '_', biomass_selected_quads_pts_distances_final$calval_pft)

# 3.6 SAVE --------------------------------------------------------------------------------------------------------------------------------

write.csv(biomass_selected_quads_pts_distances_final, outName, row.names = FALSE)
