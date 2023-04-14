######################################################################################################
######################################################################################################

# CODE DESCRIPTION

# This script tidies the field biomass harvest data
#   - Removes outlier quadrats
#   - Adds location data
#   - Converts biomass units

######################################################################################################
######################################################################################################

# SET OUTPUT DIRECTORY

output_results = FALSE

outPath = 'data/'

outName = 'biomass_field_data_lat_long.csv'

######################################################################################################
######################################################################################################

# 1. IMPORT AND TIDY DATA -------------------------------------------------------------------------------------------------------------------------------

# 1.1 READ IN DATA --------------------------------------------------------------------------------------------------------------------------------

library(tidyverse)

# Read in field data
biomassData = read.csv(paste0(outPath, 'biomass_field_UAV_data_wide.csv'))

# Read in quadrat data
quads = read.csv(paste0(outPath, 'UAV_quadrats.csv'))

# 1.2 TIDY --------------------------------------------------------------------------------------------------------------------------------

# Select necessary columns
biomass_final = dplyr::select(biomassData, c(site_code, quadrat_num, PFT, weight_observed, quad_label))

# 1.3 REMOVE OBSERVATIONS WITH POOR DATA QUALITY

# TWELVEMILE 50m -- multi-bag sample and suspect missing bag based on comparison to cover/biomass from other quadrats
biomass_final = biomass_final[!(biomass_final$quad_label == 'TWELVEMILE 50m' & biomass_final$PFT == 'LICHENS'),]
biomass_final = biomass_final[!(biomass_final$quad_label == 'TWELVEMILE 50m' & biomass_final$PFT == 'TOTAL'),]

######################################################################################################
######################################################################################################

# 3. ADD LOCATION DATA -------------------------------------------------------------------------------------------------------------------------------

# Set up site/quadrat labels
quads$quad_label = paste(quads$site_code, quads$quadrat_num)

# Select necessary columns
quads = quads %>% dplyr::select(quad_label, latitude, longitude)

# Join to main data
biomass_final = dplyr::left_join(biomass_final, quads, by = "quad_label")

######################################################################################################
######################################################################################################

# 4. ADD BIOMASS UNITS -------------------------------------------------------------------------------------------------------------------------------

names(biomass_final) = c("site_code", "quadrat_num", "pft", "biomass_g_025m2", "quad_label", "latitude", "longitude")

# Convert from grams per 0.25m2 to more useful units
biomass_final$biomass_g_m2 = biomass_final$biomass_g_025m2 * 4 # g/m2
biomass_final$biomass_g_900m2 = biomass_final$biomass_g_m2 * 900 # g/900m2
biomass_final$biomass_kg_m2 = biomass_final$biomass_g_m2 / 1000 # kg/m2

######################################################################################################
######################################################################################################

# 4. SAVE -------------------------------------------------------------------------------------------------------------------------------

if(output_results){write.csv(biomass_final, paste(outPath, outName), row.names = FALSE)}
