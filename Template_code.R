# Template code -------------------------------
# Name of original paper
# doi
# This template code has 3 sections: #1 Preparation, #2 Model, #3 Model Results, #4 Niche comparison
#
# 1 Preparation ---------------------
#
## 1.1 Data sources ==================
#
# Source of occurence data:
# If it comes from a service such as Gbif remember to include the reference code of your query.
#  example:
#  Gbif
#  Reference code: GBIF.org (09 December 2022) GBIF Occurrence Download https://doi.org/10.15468/dl.8eqnqv
#
# If the original data is not available, state how the data set was recreated
#
# Source of climate data:
# Source, resolution, climate variables, link


## 1.2 Package library =============
# Packages always used:
library(tidyverse)
library(terra)
library(INLA)
library(pROC)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(ade4)
library(ecospat)


# Packages used, but not loaded:
# library(raster)

# Packages used for data prep:
library(here)


## 1.3 Function files ========
source(paste0(here(), "/Joe_functions.R"))
source(paste0(here(), "/Isak_functions.R"))

## 1.4 Import occurrence ========
sp_occurence <- read_csv(file = paste0(here(), "/data/sp_data.csv"))

# Clean the data if necessary
# Should include numeric colums for Latitude, Longitude and occurrenceStatus (which should always be 1)
sp_occurence_clean <- sp_occurence |>
  mutate(Latitude = as.double(decimalLatitude), Longitude = as.double(decimalLongitude)) |>
  select(Latitude, Longitude, occurrenceStatus) |>
  replace(is.character(sp_occurence$occurrenceStatus) == TRUE, sp_occurence$occurrenceStatus <- 1)

## 1.5 Import climate data ========
# This example data comes from Worldclim
bio1 <- rast("data/wc2.1_5m_bio/wc2.1_5m_bio_1.tif")
bio2 <- rast("data/wc2.1_5m_bio/wc2.1_5m_bio_2.tif")
bio3 <- rast("data/wc2.1_5m_bio/wc2.1_5m_bio_3.tif")
bio4 <- rast("data/wc2.1_5m_bio/wc2.1_5m_bio_4.tif")
bio5 <- rast("data/wc2.1_5m_bio/wc2.1_5m_bio_5.tif")
bio6 <- rast("data/wc2.1_5m_bio/wc2.1_5m_bio_6.tif")
bio7 <- rast("data/wc2.1_5m_bio/wc2.1_5m_bio_7.tif")
bio8 <- rast("data/wc2.1_5m_bio/wc2.1_5m_bio_8.tif")
bio9 <- rast("data/wc2.1_5m_bio/wc2.1_5m_bio_9.tif")
bio10 <- rast("data/wc2.1_5m_bio/wc2.1_5m_bio_10.tif")
bio11 <- rast("data/wc2.1_5m_bio/wc2.1_5m_bio_11.tif")
bio12 <- rast("data/wc2.1_5m_bio/wc2.1_5m_bio_12.tif")
bio13 <- rast("data/wc2.1_5m_bio/wc2.1_5m_bio_13.tif")
bio14 <- rast("data/wc2.1_5m_bio/wc2.1_5m_bio_14.tif")
bio15 <- rast("data/wc2.1_5m_bio/wc2.1_5m_bio_15.tif")
bio16 <- rast("data/wc2.1_5m_bio/wc2.1_5m_bio_16.tif")
bio17 <- rast("data/wc2.1_5m_bio/wc2.1_5m_bio_17.tif")
bio18 <- rast("data/wc2.1_5m_bio/wc2.1_5m_bio_18.tif")
bio19 <- rast("data/wc2.1_5m_bio/wc2.1_5m_bio_19.tif")



# Add all climate variable into a single terra object
climfull <- terra::`add<-`(bio1, c(bio2, bio3, bio4, bio5, bio6, bio7, bio8, bio9,
                                   bio10, bio11, bio12, bio13, bio14, bio15, bio16, bio17, bio18, bio19))

# Make a list of names for each climate variable (make sure to have the correct order)
clim_names <- c("Annual Mean Temperature (C)",
                "Mean Diurnal Range (Mean of monthly (max temp - min temp, C))",
                "Isothermality (BIO2/BIO7) (×100)", 
                "Temperature Seasonality (standard deviation ×100)",
                "Max Temperature of Warmest Month (C)",
                "Min Temperature of Coldest Month (C)",
                "Temperature Annual Range (BIO5-BIO6)",
                "Mean Temperature of Wettest Quarter (C)",
                "Mean Temperature of Driest Quarter (C)",
                "Mean Temperature of Warmest Quarter (C)",
                "Mean Temperature of Coldest Quarter",
                "Annual Precipitation (mm)",
                "Precipitation of Wettest Month (mm)",
                "Precipitation of Driest Month (mm)",
                "Precipitation Seasonality (Coefficient of Variation)",
                "Precipitation of Wettest Quarter",
                "Precipitation of Driest Quarter",
                "Precipitation of Warmest Quarter (mm)",
                "Precipitation of Coldest Quarter")

### 1.5Alt: Loading data with the sdmpredictors predictors package =========
# sdmpredictors is a package that allows you to download climate data from several common sources in a raster format
library(sdmpredictors)

# list_datasets() allows you to see which datasets are available
# list_layers() allows you to see all available climate variables with various details. You should save this to an object to avoid the print limit

# load_layers() loads raster layers of climate variables, the names of layers are found in list_layers()
climfull <- load_layers(c("BO_bathymean", "BO2_curvelmax_bdmean", "BO2_curvelrange_bdmean",
                          "BO2_dissoxmin_bdmean", "BO2_dissoxrange_bdmean", "BO2_tempmean_bdmean",
                          "BO2_temprange_bdmean", "BO2_salinitymin_bdmean", "BO2_salinityrange_bdmean"))

# For our purposes the object should be transformed into a spatraster
climfull <- climfull |> 
  rast()


## 1.6 Define expanded prediction range ==============

# This defines the geographical range of the predictions beyond the initial model prediction (W, E, S, N)
model_range_extended <- c(1, 2, 1, 2)

# This will be used in both native and invasive sections
# The function is available in the Isak_functions.R file
# (aggregation factor is only intended for trial runs, reducing calculation time and required memory.
# Default makes no aggregation, ignore warning message: [aggregate] all values in argument 'fact' are 1, nothing to do)
climate_extended_prediction_range <- define_range_and_aggregation(climfull, model_range_extended,
  aggregation_factor = 1
) |>
  raster::stack() |>
  as("SpatialGridDataFrame")


## 1.7 Native population ==============

# This defines the geographical range of the native model (W, E, S, N)
model_range_native <- c(-86, -72, 19, 23.5)

# clean occurrence data
sp_native <- sp_occurence_clean |>
  filter(
    Longitude > model_range_native[1], Longitude < model_range_native[2],
    Latitude > model_range_native[3], Latitude < model_range_native[4]
  )

# Prepare data
climate_native <- define_range_and_aggregation(climfull, model_range_native)

# The function is available in the Isak_functions.R file
climate_grid_native <- make_climate_grid(climate_native)

# The function is available in the Isak_functions.R file
occurence_grid_native <- make_occurence_grid(sp_native, climate_native)


## 1.8 Invasive population ==============

# This defines the geographical range of the invasive model (W, E, S, N)
model_range_invasive <- c(-100, -79.5, 24, 33.5)

# clean occurrence data
sp_invasive <- sp_occurence_clean |>
  filter(
    Longitude > model_range_invasive[1], Longitude < model_range_invasive[2],
    Latitude > model_range_invasive[3], Latitude < model_range_invasive[4]
  )

# Prepare data
climate_invasive <- define_range_and_aggregation(climfull, model_range_invasive)

# The function is available in the Isak_functions.R file
climate_grid_invasive <- make_climate_grid(climate_invasive)

# The function is available in the Isak_functions.R file
occurence_grid_invasive <- make_occurence_grid(sp_invasive, climate_invasive)


# ## 1.9 Plot occurrences ==============
# 
# coastline <- ne_countries(continent = "north america", returnclass = "sf", scale = 50) |>
#   st_crop(c(xmin = -105, xmax = -70, ymin = 10, ymax = 38))
# 
# ggplot() +
#   geom_sf(data = coastline, fill = "white") +
#   geom_point(data = sp_native, aes(Longitude, Latitude, color = "native")) +
#   geom_point(data = sp_invasive, aes(Longitude, Latitude, color = "invasive")) +
#   scale_x_continuous(expand = c(0, 0), limits = c(-102, -70)) +
#   scale_y_continuous(expand = c(0, 0), limits = c(13, 35)) +
#   labs(title = "Species occurences (up to 2011)", x = NULL, y = NULL, color = "Population") +
#   theme(panel.background = element_rect(
#     fill = "lightblue",
#     colour = "lightblue"
#   ), text = element_text(size = 20))


# 2 Model --------------

## 2.1 Native ============

# This function requires that the functions: sdmINLASpatRandEff, updateProgressBar,
# inlaMeshToSpatialPolygons, createProgressUpdater
# All required functions are available in the Joe_functions.R file
# Requires the INLA package to be loaded
mod_native <- sdmINLASpatRandEff(occurence_grid_native, climate_grid_native,
  meshParameters = list(max.n.strict = -1),
  createGeoTIFF = FALSE, outFolder = getwd(),
  myName = "sp_Native"
)

## 2.1 Invasive ============
mod_invasive <- sdmINLASpatRandEff(occurence_grid_invasive, climate_grid_invasive,
  meshParameters = list(max.n.strict = -1),
  createGeoTIFF = FALSE,
  outFolder = getwd(),
  myName = "sp_Invasive"
)

# 3 Model Results -------------

## 3.1 Native Results =========
# Retrive model results
obj_native <- readRDS(file = "sp_Native.rds")

### 3.1.1 Mean occurrence probability map #########
native_occurence_estimate <- as.data.frame(obj_native$spatialPredictions[1])

ggplot(native_occurence_estimate) +
  geom_raster(aes(s1, s2, fill = meanEst)) +
  coord_sf() +
  scale_fill_viridis_c(
    limits = c(0, 1), option = "magma",
    begin = 0.12, end = 1, name = "Occurence probability"
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(model_range_native[1], model_range_native[2])) +
  scale_y_continuous(expand = c(0, 0), limits = c(model_range_native[3], model_range_native[4])) +
  labs(title = "Native occurrence prediction", x = NULL, y = NULL) +
  theme_bw()


### 3.1.2 Spatial effect map #########
native_spatial_effect <- as.data.frame(obj_native$spatialPredictions[7])

ggplot(native_spatial_effect, aes(s1, s2)) +
  geom_raster(aes(fill = meanSpatialPred)) +
  scale_fill_viridis_c(
    option = "magma",
    end = 0.9, name = "Spatial effect"
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(model_range_native[1], model_range_native[2])) +
  scale_y_continuous(expand = c(0, 0), limits = c(model_range_native[3], model_range_native[4])) +
  coord_sf() +
  labs(title = "Spatial effect in native range", x = NULL, y = NULL) +
  theme_bw()

### 3.1.3 Predicted occurrence probability in an expanded range #########
predict_native <- predictSD(obj_native, climate_extended_prediction_range, origClimateData = climate_grid_native)
native_predicted_estimate <- as.data.frame(predict_native[1])

ggplot(native_predicted_estimate, aes(s1, s2)) +
  geom_raster(aes(fill = OccurrenceProb)) +
  coord_sf() +
  scale_fill_viridis_c(
    limits = c(0, 1), option = "magma",
    begin = 0.12, end = 1, name = "Suitability"
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(model_range_extended[1], model_range_extended[2])) +
  scale_y_continuous(expand = c(0, 0), limits = c(model_range_extended[3], model_range_extended[4])) +
  labs(title = "Native niche prediction", x = NULL, y = NULL) +
  theme_bw()

### 3.1.4 Individual climate response curves ########
# obj_native$responsePredictions contains each climate variable with the occurrence prediction (with upper and lower confidence interval) along the available climate gradient in the specified range
preds_native <- obj_native$responsePredictions

# This function gets the parameters of each response curve
# The function is available in the Joe_functions.R file.
pars_native <- getPars(obj_native, climate_grid_native)


# WARNING! The name after "preds_native$" will change depending on the original climate data,
# make sure to check preds_native beforehand.
# This function should be repeated for each climate variable
map2(preds_native, clim_names, plot_response_curve)


# Calculate standardized individual response curves
# This functions standardises the area under the curve to better show the change across the gradient regardless of effect size. As such disregard the actual response values
integral_native_bio5 <- calculate_integral_response_curve(
  climate_grid_native, "wc2.1_30s_bio_5",
  preds_native$wc2.1_30s_bio_5$covarVal,
  pars_native
)

plot_integral_response_curve(integral_native_bio5, "Name of the climate variable")

### 3.1.5 AUC ########
observations_native <- occurence_grid_native$occurrence
predictions_native <- obj_native$spatialPredictions$meanEst
roc_object_native <- roc(observations_native, predictions_native)
sp_AUC_native <- auc(roc_object_native)



## 3.2 Invasive Results =========
# Retrive model results
obj_invasive <- readRDS(file = "sp_invasive.rds")

### 3.2.1 Mean occurrence probability map #########
invasive_occurence_estimate <- as.data.frame(obj_invasive$spatialPredictions[1])

ggplot(invasive_occurence_estimate) +
  geom_raster(aes(s1, s2, fill = meanEst)) +
  coord_sf() +
  scale_fill_viridis_c(
    limits = c(0, 1), option = "magma",
    begin = 0.12, end = 1, name = "Occurence probability"
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(model_range_invasive[1], model_range_invasive[2])) +
  scale_y_continuous(expand = c(0, 0), limits = c(model_range_invasive[3], model_range_invasive[4])) +
  labs(title = "invasive occurrence prediction", x = NULL, y = NULL) +
  theme_bw()


### 3.2.2 Spatial effect map #########
invasive_spatial_effect <- as.data.frame(obj_invasive$spatialPredictions[7])

ggplot(invasive_spatial_effect, aes(s1, s2)) +
  geom_raster(aes(fill = meanSpatialPred)) +
  scale_fill_viridis_c(
    option = "magma",
    end = 0.9, name = "Spatial effect"
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(model_range_invasive[1], model_range_invasive[2])) +
  scale_y_continuous(expand = c(0, 0), limits = c(model_range_invasive[3], model_range_invasive[4])) +
  coord_sf() +
  labs(title = "Spatial effect in invasive range", x = NULL, y = NULL) +
  theme_bw()

### 3.2.3 Predicted occurrence probability in an expanded range #########
predict_invasive <- predictSD(obj_invasive, climate_extended_prediction_range, origClimateData = climate_grid_invasive)
invasive_predicted_estimate <- as.data.frame(predict_invasive[1])

ggplot(invasive_predicted_estimate, aes(s1, s2)) +
  geom_raster(aes(fill = OccurrenceProb)) +
  coord_sf() +
  scale_fill_viridis_c(
    limits = c(0, 1), option = "magma",
    begin = 0.12, end = 1, name = "Suitability"
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(model_range_extended[1], model_range_extended[2])) +
  scale_y_continuous(expand = c(0, 0), limits = c(model_range_extended[3], model_range_extended[4])) +
  labs(title = "invasive niche prediction", x = NULL, y = NULL) +
  theme_bw()

### 3.2.4 Individual climate response curves ########
# obj_invasive$responsePredictions contains each climate variable with the occurrence prediction (with upper and lower confidence interval) along the available climate gradient in the specified range
preds_invasive <- obj_invasive$responsePredictions

# This function gets the parameters of each response curve
# The function is available in the Joe_functions.R file.
pars_invasive <- getPars(obj_invasive, climate_grid_invasive)


# WARNING! The name after "preds_invasive$" will change depending on the original climate data,
# make sure to check preds_invasive beforehand.
# This function should be repeated for each climate variable
map2(preds_invasive, clim_names, plot_response_curve)


# Plot standardized individual response curves
# This functions standardises the area under the curve to better show the change across the gradient regardless of effect size. As such disregard the actual response values
integral_invasive_bio5 <- calculate_integral_response_curve(
  climate_grid_invasive, "wc2.1_30s_bio_5",
  preds_invasive$wc2.1_30s_bio_5$covarVal,
  pars_invasive
)

plot_integral_response_curve(integral_invasive_bio5, "Name of the climate variable")

### 3.2.5 AUC ########
observations_invasive <- occurence_grid_invasive$occurrence
predictions_invasive <- obj_invasive$spatialPredictions$meanEst
roc_object_invasive <- roc(observations_invasive, predictions_invasive)
sp_AUC_invasive <- auc(roc_object_invasive)

# 4 Niche comparison ---------

## 4.1 Univariate ==========


compare_response_curves(
  preds_native$wc2.1_30s_bio_1,
  preds_invasive$wc2.1_30s_bio_1
)
map2(preds_native, preds_invasive, compare_response_curves)

# Standardised curves

compare_response_curves(integral_native_bio5, integral_invasive_bio5)




# ## 4.2 Multivariate ===========
# 
# climate_background <- as.data.frame(climate_extended_prediction_range) |>
#   select(-s1, -s2) |>
#   rename(
#     "Annual Mean Temperature" = wc2.1_30s_bio_1,
#     "Mean Diurnal Range" = wc2.1_30s_bio_2,
#     "Max Temperature of Warmest Month" = wc2.1_30s_bio_5,
#     "Mean Temperature of Wettest Quarter" = wc2.1_30s_bio_8,
#     "Mean Temperature of Driest Quarter" = wc2.1_30s_bio_9,
#     "Mean Temperature of Warmest Quarter" = wc2.1_30s_bio_10,
#     "Annual Precipitation" = wc2.1_30s_bio_12,
#     "Precipitation of Wettest Month" = wc2.1_30s_bio_13,
#     "Precipitation of Driest Month" = wc2.1_30s_bio_14,
#     "Precipitation of Warmest Quarter" = wc2.1_30s_bio_18
#   )
# 
# occurence_native_data_frame <- as.data.frame(occurence_grid_native)
# 
# # Data frame showing each climate variable at each point in the native area
# native_data_frame <- as.data.frame(climate_grid_native) |>
#   left_join(occurence_native_data_frame) |>
#   select(-s1, -s2)
# 
# # Add predicted occurrence probability from the native model to extended climate area data frame
# native_prediction_data_frame <- as.data.frame(climate_extended_prediction_range) |>
#   left_join(native_predicted_estimate) |>
#   select(-s1, -s2)
# 
# occurence_invasive_data_frame <- as.data.frame(occurence_grid_invasive)
# 
# # Data frame showing each climate variable at each point in the invasive area
# invasive_data_frame <- as.data.frame(climate_grid_invasive) |>
#   left_join(occurence_invasive_data_frame) |>
#   select(-s1, -s2)
# 
# # Add predicted occurrence probability from the native model to extended climate area data frame
# invasive_prediction_data_frame <- as.data.frame(climate_extended_prediction_range) |>
#   left_join(invasive_predicted_estimate) |>
#   select(-s1, -s2)
# 
# # Make a pca of the climate in the extended prediction area
# pca.env <- dudi.pca(climate_background, scannf = F, nf = 2)
# 
# # Plots the contribution of each
# ecospat.plot.contrib(contrib = pca.env$co, eigen = pca.env$eig)
# 
# ## Convert
# # PCA scores for the whole study area
# scores.globclim <- pca.env$li
# # PCA scores for the species native distribution
# scores.sp.nat <- suprow(pca.env, native_data_frame[which(native_data_frame[, 11] == 1), 1:10])$li
# # PCA scores for the species invasive distribution
# scores.sp.inv <- suprow(pca.env, invasive_data_frame[which(invasive_data_frame[, 11] == 1), 1:10])$li
# # PCA scores for the whole native study area
# scores.clim.nat <- suprow(pca.env, native_data_frame[, 1:10])$li
# # PCA scores for the whole invaded study area
# scores.clim.inv <- suprow(pca.env, invasive_data_frame[, 1:10])$li
# # PCA scores for the occurrence prediction based on the native model
# scores.pred.nat <- suprow(pca.env, native_prediction_data_frame[, 1:10])$li
# # PCA scores for the occurrence prediction based on the invasive model
# scores.pred.inv <- suprow(pca.env, invasive_prediction_data_frame[, 1:10])$li
# 
# # Assigns the occurrence probability from the native model to coordinates in the pca
# clim_space_native <- bind_cols(scores.pred.nat, native_predicted_estimate) |>
#   select(-s1, -s2)
# # Assigns the occurrence probability from the invasive model to coordinates in the pca
# clim_space_invasive <- bind_cols(scores.pred.inv, invasive_predicted_estimate) |>
#   select(-s1, -s2)
# 
# 
# # Gridding the native niche
# grid.clim.nat <- ecospat.grid.clim.dyn(
#   glob = scores.globclim,
#   glob1 = scores.clim.nat,
#   sp = scores.sp.nat, R = 100,
#   th.sp = 0
# )
# 
# # Gridding the invasive niche
# grid.clim.inv <- ecospat.grid.clim.dyn(
#   glob = scores.globclim,
#   glob1 = scores.clim.inv,
#   sp = scores.sp.inv, R = 100,
#   th.sp = 0
# )
# 
# 
# # Plots the available climate and occurrences in the native and invasive ranges in a pca
# ecospat.plot.niche.dyn(grid.clim.nat, grid.clim.inv,
#   quant = 0.25, interest = 1,
#   title = "Niche Overlap", name.axis1 = "PC1 = 40.85%",
#   name.axis2 = "PC2 = 25.64%"
# )
# ecospat.shift.centroids(scores.sp.nat, scores.sp.inv, scores.clim.nat, scores.clim.inv)
# 
# # Save the coordinates of occurrence probabilities in the pca at several thresholds
# # No threshold saves all points from the extended prediction area
# coords_100 <- clim_space_native |>
#   select(-OccurrenceProb) |>
#   as.matrix()
# 
# coords_native_75 <- clim_space_native |>
#   filter(OccurrenceProb > 0.25) |>
#   select(-OccurrenceProb) |>
#   as.matrix()
# 
# coords_native_50 <- clim_space_native |>
#   filter(OccurrenceProb > 0.5) |>
#   select(-OccurrenceProb) |>
#   as.matrix()
# 
# coords_native_25 <- clim_space_native |>
#   filter(OccurrenceProb > 0.75) |>
#   select(-OccurrenceProb) |>
#   as.matrix()
# 
# 
# coords_invasive_75 <- clim_space_invasive |>
#   filter(OccurrenceProb > 0.25) |>
#   select(-OccurrenceProb) |>
#   as.matrix()
# 
# coords_invasive_50 <- clim_space_invasive |>
#   filter(OccurrenceProb > 0.5) |>
#   select(-OccurrenceProb) |>
#   as.matrix()
# 
# coords_invasive_25 <- clim_space_invasive |>
#   filter(OccurrenceProb > 0.75) |>
#   select(-OccurrenceProb) |>
#   as.matrix()
# 
# # Make hulls
# hull_100 <- inla.nonconvex.hull(coords_100, convex = -0.02, concave = -0.15, resolution = 300)
# hull_native_75 <- inla.nonconvex.hull(coords_native_75, convex = -0.02, concave = -0.15, resolution = 300)
# hull_native_50 <- inla.nonconvex.hull(coords_native_50, convex = -0.02, concave = -0.15, resolution = 300)
# hull_native_25 <- inla.nonconvex.hull(coords_native_25, convex = -0.02, concave = -0.15, resolution = 300)
# 
# hull_invasive_75 <- inla.nonconvex.hull(coords_invasive_75, convex = -0.02, concave = -0.15, resolution = 300)
# hull_invasive_50 <- inla.nonconvex.hull(coords_invasive_50, convex = -0.02, concave = -0.15, resolution = 300)
# hull_invasive_25 <- inla.nonconvex.hull(coords_invasive_25, convex = -0.02, concave = -0.15, resolution = 300)
# 
# ecospat.plot.niche(grid.clim.nat)
# lines(hull_100$loc, col = "red")
# lines(hull_native_75$loc, col = "red")
# lines(hull_native_50$loc, col = "red")
# lines(hull_native_25$loc, col = "red")
# # text() to ad values to the lines
# # polygon() for fill
# 
# 
# ecospat.plot.niche(grid.clim.inv)
# lines(hull_invasive_75$loc, col = "blue")
# 
# plot(hull_100$loc,
#   type = "lines",
#   xlab = "PC1 = 40.85%", ylab = "PC2 = 25.64%", main = "Predicted niches in climate space"
# )
# polygon(hull_invasive_75$loc, col = adjustcolor("blue", alpha.f = 0.1))
# polygon(hull_invasive_50$loc, col = adjustcolor("blue", alpha.f = 0.3))
# polygon(hull_invasive_25$loc, col = adjustcolor("blue", alpha.f = 0.5))
# polygon(hull_native_75$loc, col = adjustcolor("red", alpha.f = 0.1))
# polygon(hull_native_50$loc, col = adjustcolor("red", alpha.f = 0.3))
# polygon(hull_native_25$loc, col = adjustcolor("red", alpha.f = 0.5))
# points(scores.sp.inv, cex = 0.7, pch = 2, col = "red")
# points(scores.sp.nat, cex = 0.7, pch = 1, col = "blue")
# # look at legend()
# 
# 
# 
# ## 4.3 Overlap index =============
# A <- raster::raster(predict_invasive)
# B <- raster::raster(predict_native)
# 
# isocat::schoenersD(A / raster::cellStats(A, stat = "sum"), B / raster::cellStats(B, stat = "sum"))
# 
# a <- rast(predict_invasive)
# b <- rast(predict_native)
# 
# spatialEco::overlap(a, b)
