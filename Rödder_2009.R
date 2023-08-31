# Rödder 2009 -------------------------------
# https://dx.doi.org/10.1111/j.1466-8238.2009.00477.x
# This template code has 3 sections: #1 Preparation, #2 Model, #3 Model Results, #4 Niche comparison
#
# 1 Preparation ---------------------
#
## 1.1 Data sources ==================
#
# Source of occurence data:
#  Gbif
#  Reference code: GBIF.org (28 March 2023) GBIF Occurrence Download  https://doi.org/10.15468/dl.ek48jv
# 
# Add reference to Vertnet -------------
# http://www.vertnet.org/resources/datatoolscode.html#t-tab1

# Currently different datasets, needs comparison in discussion

# Georeferencing does not quite make sense (can be unambiguously assigned to a single grid cell)

#
# Source of climate data:
# WorldClim, 5 arc-minutes, all standard (19) WorldClim variables, https://www.worldclim.org/data/worldclim21.html
# via sdmpredict package


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

library(gridExtra)



# Packages used, but not loaded:
# library(raster)

# Packages used for data prep:
library(here)
library(readxl)


## 1.3 Function files ========
source(paste0(here(), "/Joe_functions.R"))
source(paste0(here(), "/Isak_functions.R"))
source(paste0(here(), "/resp_curves_update.R"))


## 1.4 Import occurrence ========
Hemidactylus_occurence <- read_xlsx(path = paste0(here(), "/data/Hemidactylus_turcicus_v3.xlsx"))

Vertnet <- read.csv(file = "data/VertNet_Reptilia_Oct2015/VertNet_Reptilia_Oct2015.csv")

Hemidactylus_occurence2 <- Vertnet |> 
  filter(scientificname == "Hemidactylus turcicus", year < 2010)

# Clean the data if necessary
# Should include numeric colums for Latitude, Longitude and occurrenceStatus (which should always be 1)
Hemidactylus_occurence_clean <- Hemidactylus_occurence |>
  mutate(Latitude = as.double(decimalLatitude), Longitude = as.double(decimalLongitude)) |>
  select(Latitude, Longitude, occurrenceStatus) |>
  replace(is.character(Hemidactylus_occurence$occurrenceStatus) == TRUE, Hemidactylus_occurence$occurrenceStatus <- 1)

Hemidactylus_occurence_clean2 <- Hemidactylus_occurence2 |>
  mutate(Latitude = as.double(decimallatitude), Longitude = as.double(decimallongitude)) |>
  select(Latitude, Longitude, occurrenceStatus) |>
  replace(is.character(Hemidactylus_occurence2$occurrenceStatus) == TRUE, Hemidactylus_occurence2$occurrenceStatus <- 1)

Hemidactylus_occurence_final <- full_join(Hemidactylus_occurence_clean, Hemidactylus_occurence_clean2)

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
bio12 <- rast("data/wc2.1_5m_bio/wc2.1_5m_bio_12.tif")
bio13 <- rast("data/wc2.1_5m_bio/wc2.1_5m_bio_13.tif")
bio14 <- rast("data/wc2.1_5m_bio/wc2.1_5m_bio_14.tif")
bio15 <- rast("data/wc2.1_5m_bio/wc2.1_5m_bio_15.tif")
bio18 <- rast("data/wc2.1_5m_bio/wc2.1_5m_bio_18.tif")
bio19 <- rast("data/wc2.1_5m_bio/wc2.1_5m_bio_19.tif")


# Add all climate variable into a single terra object
climfull <- c(bio1, bio2, bio3, bio4, bio5, bio6, bio7, bio8, bio9,
                                   bio12, bio13, bio14, bio15, bio18, bio19)


# Make a list of names for each climate variable (make sure to have the correct order)
clim_names <- c("Annual Mean Temperature (\u00B0C)",
                "Mean Diurnal Range (Mean of monthly (max temp - min temp, \u00B0C))",
                "Isothermality (BIO2/BIO7) (*100)", 
                "Temperature Seasonality (standard deviation *100)",
                "Max Temperature of Warmest Month (\u00B0C)",
                "Min Temperature of Coldest Month (\u00B0C)",
                "Temperature Annual Range (BIO5-BIO6)",
                "Mean Temperature of Wettest Quarter (\u00B0C)",
                "Mean Temperature of Driest Quarter (\u00B0C)",
                "Annual Precipitation (mm)",
                "Precipitation of Wettest Month (mm)",
                "Precipitation of Driest Month (mm)",
                "Precipitation Seasonality (Coefficient of Variation)",
                "Precipitation of Warmest Quarter (mm)",
                "Precipitation of Coldest Quarter")

clim_codes <- names(climfull)


## 1.6 Define expanded prediction range ==============

# This defines the geographical range of the predictions beyond the initial model prediction (W, E, S, N)
# model_range_extended <- c(-118.5, -76.5, 15, 37.5)
# 
# # This will be used in both native and invasive sections
# # The function is available in the Isak_functions.R file
# # (aggregation factor is only intended for trial runs, reducing calculation time and required memory.
# # Default makes no aggregation)
# climate_extended_prediction_range <- define_range_and_aggregation(climfull, model_range_extended,
#                                                                   aggregation_factor = 1
# ) |>
#   raster::stack() |>
#   as("SpatialGridDataFrame")

## 1.7 Native population ==============

# This defines the geographical range of the native model (W, E, S, N)
model_range_native <- c(-11, 39.3, 21, 46.7)

# clean occurrence data
Hemidactylus_native <- Hemidactylus_occurence_clean |>
  filter(
    Longitude > model_range_native[1], Longitude < model_range_native[2],
    Latitude > model_range_native[3], Latitude < model_range_native[4]
  )


# Prepare data
climate_native <- define_range_and_aggregation(climfull, model_range_native, aggregation_factor = 2)

# The function is available in the Isak_functions.R file
climate_grid_native <- make_climate_grid(climate_native)

# The function is available in the Isak_functions.R file
occurence_grid_native <- make_occurence_grid(Hemidactylus_native, climate_native)


## 1.8 Invasive population ==============

# This defines the geographical range of the invasive model (W, E, S, N)
model_range_invasive <- c(-120.2, -75.9, 15.3, 37.8)

# clean occurrence data
Hemidactylus_invasive <- Hemidactylus_occurence_clean |>
  filter(
    Longitude > model_range_invasive[1], Longitude < model_range_invasive[2],
    Latitude > model_range_invasive[3], Latitude < model_range_invasive[4]
  )

# Prepare data
climate_invasive <- define_range_and_aggregation(climfull, model_range_invasive, aggregation_factor = 2)

# The function is available in the Isak_functions.R file
climate_grid_invasive <- make_climate_grid(climate_invasive)

# The function is available in the Isak_functions.R file
occurence_grid_invasive <- make_occurence_grid(Hemidactylus_invasive, climate_invasive)


# ## 1.9 Plot occurrences ==============

native_coastline <- ne_countries(continent = c("europe", "africa", "asia"), returnclass = "sf", scale = 50) 


occ_nat <- ggplot() +
  geom_sf(data = native_coastline, fill = "white") +
  geom_point(data = Hemidactylus_native, aes(Longitude, Latitude), size = 1) +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(model_range_native[1],  model_range_native[2])) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(model_range_native[3], model_range_native[4])) +
  labs(title = "Native occurences", x = NULL, y = NULL)+
  theme(panel.background = element_rect(fill = "lightblue", colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))
occ_nat


invasive_coastline <- ne_countries(continent = "north america", returnclass = "sf", scale = 50)

occ_inv <- ggplot() +
  geom_sf(data = invasive_coastline, fill = "white") +
  geom_point(data = Hemidactylus_invasive, aes(Longitude, Latitude), size = 1) +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(model_range_invasive[1],  model_range_invasive[2])) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(model_range_invasive[3], model_range_invasive[4])) +
  labs(title = "Invasive occurences", x = NULL, y = NULL) +
  theme(panel.background = element_rect(fill = "lightblue", colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))
occ_inv


grid.arrange(occ_nat, occ_inv, nrow = 1)


# 2 Model --------------

## 2.1 Native ============

# This function requires that the functions: sdmINLASpatRandEff, updateProgressBar,
# inlaMeshToSpatialPolygons, createProgressUpdater
# All required functions are available in the Joe_functions.R file
# Requires the INLA package to be loaded
mod_native <- sdmINLASpatRandEff(occurence_grid_native, climate_grid_native,
                                 meshParameters = list(max.n.strict = -1),
                                 createGeoTIFF = FALSE, outFolder = getwd(),
                                 myName = "Hemidactylus_Native_agg_new"
)
#agg = 2

## 2.1 Invasive ============
mod_invasive <- sdmINLASpatRandEff(occurence_grid_invasive, climate_grid_invasive,
                                   meshParameters = list(max.n.strict = -1),
                                   createGeoTIFF = FALSE,
                                   outFolder = getwd(),
                                   myName = "Hemidactylus_Invasive_agg_new"
)
#agg = 2


# 3 Model Results -------------

## 3.1 Native Results =========
# Retrive model results
obj_native <- readRDS(file = "Hemidactylus_Native_agg_new.rds")

### 3.1.1 Mean occurrence probability map #########
# native_occurence_estimate <- as.data.frame(obj_native$spatialPredictions[1])
# 
# ggplot(native_occurence_estimate) +
#   geom_tile(aes(s1, s2, fill = meanEst)) +
#   coord_sf() +
#   scale_fill_viridis_c(
#     limits = c(0, 1), option = "magma",
#     begin = 0.12, end = 1, name = "Occurence probability"
#   ) +
#   scale_x_continuous(expand = c(0, 0), limits = c(model_range_native[1], model_range_native[2])) +
#   scale_y_continuous(expand = c(0, 0), limits = c(model_range_native[3], model_range_native[4])) +
#   labs(title = "Native occurrence prediction", x = NULL, y = NULL) +
#   theme_bw()

native_occurence_estimate <- as.data.frame(obj_native$spatialPredictions[1])
summary(native_occurence_estimate)


x_nat <- ggplot(native_occurence_estimate) +
  geom_tile(aes(s1, s2, fill = meanEst / max(native_occurence_estimate$meanEst))) +
  coord_sf() +
  scale_fill_viridis_c(
    # limits = c(0, .4),
    option = "magma",
    begin = 0.12, end = 1, name = "Relative \noccurence \nprobability"
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(model_range_native[1], model_range_native[2])) +
  scale_y_continuous(expand = c(0, 0), limits = c(model_range_native[3], model_range_native[4])) +
  labs(title = "a) Native distribution estimate", x = NULL, y = NULL) +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), panel.background = element_rect(fill = "lightblue"))
x_nat

native_suitability <- as.data.frame(obj_native$spatialPredictions[6]) |> 
  mutate(meanPred = invLogit(meanClimatePred))
summary(native_suitability)

y_nat <- ggplot(native_suitability) +
  geom_tile(aes(s1, s2, fill = meanPred / max(native_suitability$meanPred))) +
  coord_sf() +
  scale_fill_viridis_c(
    limits = c(0, .2),
    option = "magma",
    begin = 0.12, end = 1, name = "Relative \nsuitability"
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(model_range_native[1], model_range_native[2])) +
  scale_y_continuous(expand = c(0, 0), limits = c(model_range_native[3], model_range_native[4])) +
  labs(title = "b) Native habitat estimate", x = NULL, y = NULL) +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), panel.background = element_rect(fill = "lightblue"))
y_nat

native_uncertainty <- as.data.frame(obj_native$spatialPredictions[4])
summary(native_uncertainty)

z_nat <-
  ggplot(native_uncertainty) +
  geom_tile(aes(s1, s2, fill = uncertaintyEst)) +
  coord_sf() +
  scale_fill_viridis_c(
    limits = c(0, .7),
    option = "magma",
    begin = 0.12, end = 1, name = "95% Credible \ninterval width"
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(model_range_native[1], model_range_native[2])) +
  scale_y_continuous(expand = c(0, 0), limits = c(model_range_native[3], model_range_native[4])) +
  labs(title = "c) Native uncertainty", x = NULL, y = NULL) +
  theme(plot.margin=grid::unit(c(0,0,2,0), "mm"), panel.background = element_rect(fill = "lightblue"))
z_nat

library(gridExtra)
grid.arrange(x_nat, y_nat, z_nat, nrow = 3, padding = unit(0.2, "line"))

### 3.1.2 Spatial effect map #########
native_spatial_effect <- as.data.frame(obj_native$spatialPredictions[7])

ggplot(native_spatial_effect, aes(s1, s2)) +
  geom_tile(aes(fill = meanSpatialPred)) +
  scale_fill_viridis_c(
    option = "magma",
    end = 0.9, name = "Spatial effect"
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(model_range_native[1], model_range_native[2])) +
  scale_y_continuous(expand = c(0, 0), limits = c(model_range_native[3], model_range_native[4])) +
  coord_sf() +
  labs(title = "Spatial effect in native range", x = NULL, y = NULL) +
  theme_bw()

### 3.1.3 Predicted occurrence probability in an native and invasive ranges  #########
predict_native_native <- predictSD(obj_native, climate_grid_native, origClimateData = climate_grid_native)
native_native_predicted_estimate <- as.data.frame(predict_native_native[1])

ggplot(native_native_predicted_estimate, aes(s1, s2)) +
  geom_tile(aes(fill = OccurrenceProb)) +
  coord_sf() +
  scale_fill_viridis_c(
    limits = c(0, 1), option = "magma",
    begin = 0.12, end = 1, name = "Suitability"
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(model_range_native[1], model_range_native[2])) +
  scale_y_continuous(expand = c(0, 0), limits = c(model_range_native[3], model_range_native[4])) +
  labs(title = "Native niche prediction in the native range", x = NULL, y = NULL) +
  theme_bw()


predict_native_invasive <- predictSD(obj_native, climate_grid_invasive, origClimateData = climate_grid_native)
native_invasive_predicted_estimate <- as.data.frame(predict_native_invasive[1])

ggplot(native_invasive_predicted_estimate, aes(s1, s2)) +
  geom_tile(aes(fill = OccurrenceProb)) +
  coord_sf() +
  scale_fill_viridis_c(
    limits = c(0, 1), option = "magma",
    begin = 0.12, end = 1, name = "Suitability"
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(model_range_invasive[1], model_range_invasive[2])) +
  scale_y_continuous(expand = c(0, 0), limits = c(model_range_invasive[3], model_range_invasive[4])) +
  labs(title = "Native niche prediction in the invasive range", x = NULL, y = NULL) +
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


### 3.1.5 AUC ########
observations_native <- occurence_grid_native$occurrence
predictions_native <- obj_native$spatialPredictions$meanEst
roc_object_native <- roc(observations_native, predictions_native)
Hemidactylus_AUC_native <- auc(roc_object_native)
Hemidactylus_AUC_native


# 3.2 Invasive Results =========
# Retrive model results
obj_invasive <- readRDS(file = "Hemidactylus_Invasive_agg_new.rds")

### 3.2.1 Mean occurrence probability map #########
# invasive_occurence_estimate <- as.data.frame(obj_invasive$spatialPredictions[1])
# 
# ggplot(invasive_occurence_estimate) +
#   geom_tile(aes(s1, s2, fill = meanEst)) +
#   coord_sf() +
#   scale_fill_viridis_c(
#     limits = c(0, 1), option = "magma",
#     begin = 0.12, end = 1, name = "Occurence probability"
#   ) +
#   scale_x_continuous(expand = c(0, 0), limits = c(model_range_invasive[1], model_range_invasive[2])) +
#   scale_y_continuous(expand = c(0, 0), limits = c(model_range_invasive[3], model_range_invasive[4])) +
#   labs(title = "invasive occurrence prediction", x = NULL, y = NULL) +
#   theme_bw()

invasive_occurence_estimate <- as.data.frame(obj_invasive$spatialPredictions[1])
summary(invasive_occurence_estimate)

x_inv <-
  ggplot(invasive_occurence_estimate) +
  geom_tile(aes(s1, s2, fill = meanEst / max(invasive_occurence_estimate$meanEst))) +
  coord_sf() +
  scale_fill_viridis_c(
    # limits = c(0, .6),
    option = "magma",
    begin = 0.12, end = 1, name = "Relative \noccurence \nprobability"
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(model_range_invasive[1], model_range_invasive[2])) +
  scale_y_continuous(expand = c(0, 0), limits = c(model_range_invasive[3], model_range_invasive[4])) +
  labs(title = "a) Invasive distribution estimate", x = NULL, y = NULL) +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), panel.background = element_rect(fill = "lightblue"))
x_inv

invasive_suitability <- as.data.frame(obj_invasive$spatialPredictions[6]) |> 
  mutate(meanPred = invLogit(meanClimatePred))
summary(invasive_suitability)

y_inv <- ggplot(invasive_suitability) +
  geom_tile(aes(s1, s2, fill = meanPred / max(invasive_suitability$meanPred))) +
  coord_sf() +
  scale_fill_viridis_c(
    # limits = c(0, .3),
    option = "magma",
    begin = 0.12, end = 1, name = "Relative \nSuitability"
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(model_range_invasive[1], model_range_invasive[2])) +
  scale_y_continuous(expand = c(0, 0), limits = c(model_range_invasive[3], model_range_invasive[4])) +
  labs(title = "b) Invasive habitat estimate", x = NULL, y = NULL) +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), panel.background = element_rect(fill = "lightblue"))
y_inv

invasive_uncertainty <- as.data.frame(obj_invasive$spatialPredictions[4])
summary(invasive_uncertainty)

z_inv <-
  ggplot(invasive_uncertainty) +
  geom_tile(aes(s1, s2, fill = uncertaintyEst)) +
  coord_sf() +
  scale_fill_viridis_c(
    limits = c(0, .2),
    option = "magma",
    begin = 0.12, end = 1, name = "95% Credible \ninterval width"
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(model_range_invasive[1], model_range_invasive[2])) +
  scale_y_continuous(expand = c(0, 0), limits = c(model_range_invasive[3], model_range_invasive[4])) +
  labs(title = "Invasive uncertainty", x = NULL, y = NULL) +
  theme(plot.margin=grid::unit(c(0,0,2,0), "mm"), panel.background = element_rect(fill = "lightblue"))
z_inv

grid.arrange(x_inv, y_inv, z_inv, nrow = 3, padding = unit(0.2, "line"))

### 3.2.2 Spatial effect map #########
invasive_spatial_effect <- as.data.frame(obj_invasive$spatialPredictions[7])

ggplot(invasive_spatial_effect, aes(s1, s2)) +
  geom_tile(aes(fill = meanSpatialPred)) +
  scale_fill_viridis_c(
    option = "magma",
    end = 0.9, name = "Spatial effect"
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(model_range_invasive[1], model_range_invasive[2])) +
  scale_y_continuous(expand = c(0, 0), limits = c(model_range_invasive[3], model_range_invasive[4])) +
  coord_sf() +
  labs(title = "Spatial effect in invasive range", x = NULL, y = NULL) +
  theme_bw()

### 3.2.3 Predicted occurrence probability in an native and invasive ranges  #########
predict_invasive_native <- predictSD(obj_invasive, climate_grid_native,
                                     origClimateData = climate_grid_invasive)
invasive_native_predicted_estimate <- as.data.frame(predict_invasive_native[1])

ggplot(invasive_native_predicted_estimate, aes(s1, s2)) +
  geom_tile(aes(fill = OccurrenceProb)) +
  coord_sf() +
  scale_fill_viridis_c(
    limits = c(0, 1), option = "magma",
    begin = 0.12, end = 1, name = "Suitability"
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(model_range_native[1], model_range_native[2])) +
  scale_y_continuous(expand = c(0, 0), limits = c(model_range_native[3], model_range_native[4])) +
  labs(title = "Invasive niche prediction in the native range", x = NULL, y = NULL) +
  theme_bw()


predict_invasive_invasive <- predictSD(obj_invasive, climate_grid_invasive,
                                       origClimateData = climate_grid_invasive)
invasive_invasive_predicted_estimate <- as.data.frame(predict_invasive_invasive[1])

ggplot(invasive_invasive_predicted_estimate, aes(s1, s2)) +
  geom_tile(aes(fill = OccurrenceProb)) +
  coord_sf() +
  scale_fill_viridis_c(
    limits = c(0, 1), option = "magma",
    begin = 0.12, end = 1, name = "Suitability"
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(model_range_invasive[1], model_range_invasive[2])) +
  scale_y_continuous(expand = c(0, 0), limits = c(model_range_invasive[3], model_range_invasive[4])) +
  labs(title = "Invasive niche prediction in the invasive range", x = NULL, y = NULL) +
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

# Calculate standardized individual response curves
# This functions standardises the area under the curve to better show the change across the gradient regardless of effect size. As such disregard the actual response values

#NOTE!!! must define climate_data before running the function
climate_data <- climate_grid_invasive

standard_curves_invasive <- Map(calculate_integral_response_curve,
                              clim_codes,
                              preds_invasive, list(pars_invasive))

plot_integral_response_curve(integral_invasive_bio5, "Name of the climate variable")

### 3.2.5 AUC ########
observations_invasive <- occurence_grid_invasive$occurrence
predictions_invasive <- obj_invasive$spatialPredictions$meanEst
roc_object_invasive <- roc(observations_invasive, predictions_invasive)
Hemidactylus_AUC_invasive <- auc(roc_object_invasive)
Hemidactylus_AUC_invasive

# 4 Niche comparison ---------

## 4.1 Objects needed from other sections (section 1 must be run first) ==========
obj_native <- readRDS(file = "Hemidactylus_Native_agg_new.rds")
preds_native <- obj_native$responsePredictions
pars_native <- getPars(obj_native, climate_grid_native)

obj_invasive <- readRDS(file = "Hemidactylus_Invasive_agg_new.rds")
preds_invasive <- obj_invasive$responsePredictions
pars_invasive <- getPars(obj_invasive, climate_grid_invasive)


## 4.2 Niche overlap =============

### 4.2.1 Plot response curves =============

# #Raw response curves
# compare_response_curves(
#   preds_native$wc2.1_2.5m_bio_1,
#   preds_invasive$wc2.1_2.5m_bio_1,
#   "name"
# )
# 
# Map(compare_response_curves, preds_native, preds_invasive, clim_names)

### 4.2.2 Calculate standardised response curves =============

#Single curve
# integrated_clim01 <- calcOverlap(clim_codes[1], obj_native, climate_grid_native, obj_invasive, climate_grid_invasive)

#All curves
integrated_curves <- lapply(clim_codes, calcOverlap, obj_native, climate_grid_native, obj_invasive, climate_grid_invasive)
names(integrated_curves) <- clim_codes


### 4.2.3 Calculate niche overlap =============
# Overlap per climate variable
all_overlaps <- vapply(integrated_curves, function(x) x$overlap, FUN.VALUE = 1.0)
# Mean overlap
niche_overlap <- mean(all_overlaps)
niche_overlap
niche_sd <- sd(all_overlaps)
niche_sd

### 4.2.3 Occurrences across climate gradients =============

spatrast_climate_grid_native <- rast(climate_grid_native)
spatrast_occurence_grid_native <- rast(occurence_grid_native)

spatrast_climate_and_occurrence_native <- c(spatrast_climate_grid_native,spatrast_occurence_grid_native)

values_climate_and_occurrences_native <- terra::values(spatrast_climate_and_occurrence_native)
# str(j)
# head(j)
# table(j[ ,"occurrence"])

values_climate_at_occurrences_native <- 
  values_climate_and_occurrences_native[values_climate_and_occurrences_native[ ,"occurrence"] > 0.5, ] 

df_climate_at_occurrences_native <- as.data.frame(values_climate_at_occurrences_native)
summary(df_climate_at_occurrences_native)
climate_at_occurrences_native <- na.omit(df_climate_at_occurrences_native)
summary(climate_at_occurrences_native)


spatrast_climate_grid_invasive <- rast(climate_grid_invasive)
spatrast_occurence_grid_invasive <- rast(occurence_grid_invasive)

spatrast_climate_and_occurrence_invasive <- c(spatrast_climate_grid_invasive,spatrast_occurence_grid_invasive)

values_climate_and_occurrences_invasive <- terra::values(spatrast_climate_and_occurrence_invasive)
# str(j)
# head(j)
# table(j[ ,"occurrence"])

values_climate_at_occurrences_invasive <- 
  values_climate_and_occurrences_invasive[values_climate_and_occurrences_invasive[ ,"occurrence"] > 0.5, ] 

df_climate_at_occurrences_invasive <- as.data.frame(values_climate_at_occurrences_invasive)
summary(df_climate_at_occurrences_invasive)
climate_at_occurrences_invasive <- na.omit(df_climate_at_occurrences_invasive)
summary(climate_at_occurrences_invasive)




### 4.2.4 Plot standardised climate response curves (with occurrences) =============
library(vioplot)
# Plot overlap curves manually

# Uncomment the next 2 lines and the line with the function dev.off() to save plot as a png
# png(file= paste0(here(), "/", "test", ".png"),
#     width=730, height=460)

par(mfrow=c(2,1), mai = c(0, 1.3, 0.2, 1), cex = 1)

response_curve_native <- with(integrated_curves$wc2.1_2.5m_bio_17$native,
                              respCurve(xSeq, datSel, pars, clim_codes[10], integral))
response_curve_invasive <- with(integrated_curves$wc2.1_2.5m_bio_17$invasive,
                                respCurve(xSeq, datSel, pars, clim_codes[10], integral))

with(integrated_curves$wc2.1_2.5m_bio_17$invasive,
     plot(xSeq, response_curve_invasive,
          xlim = range(c(xSeq, integrated_curves$wc2.1_2.5m_bio_17$native$xSeq)), col = "red", 
          ylim = c(0, max(c(response_curve_native, response_curve_invasive))),
          type = "l", xlab = "",
          ylab = "",
          axes = FALSE))

with(integrated_curves$wc2.1_2.5m_bio_17$native, lines(xSeq, response_curve_native, col = "blue"))
abline(v = integrated_curves$wc2.1_2.5m_bio_17$cRange, lwd = 2, lty = 2)
legend("topleft", c("Native", "Invasive"), fill = c("blue","red"), bty = "o")
box()
axis(2, las = 1)
axis(2, at = mean(par("usr")[3:4]), labels = "Occurrence density", tick = FALSE, line = 2.5)


par(mai = c(1.3, 1.3, 0.1, 1))

vioplot(wc2.1_2.5m_bio_17~occurrence, data = climate_at_occurrences_native,
        side = "right", horizontal = TRUE, col = "blue",
        ylim = range(c(integrated_curves$wc2.1_2.5m_bio_17$invasive$xSeq,
                       integrated_curves$wc2.1_2.5m_bio_17$native$xSeq)),
        ylab = "", xlab = clim_names[10], xaxt = "n")

abline(h = 1, lwd = 2, lty = 2)

vioplot(wc2.1_2.5m_bio_17~occurrence, data = climate_at_occurrences_invasive,
        side = "left", horizontal = TRUE, col = "red", add = TRUE)

# dev.off()


# Single curve
niche_overlap_plot(integrated_curves$wc2.1_5m_bio_1, clim_codes[1], clim_names[1],
                   climate_at_occurrences_native, climate_at_occurrences_invasive, "test", all_overlaps[1])


# All curves
Map(niche_overlap_plot, integrated_curves, clim_codes, clim_names,
    list(climate_at_occurrences_native), list(climate_at_occurrences_invasive), "Rodder_2009_new", all_overlaps)


## 4.3  Climate/ Spatial effect contribution ==========
calcVariancePartitioning(obj_native)

calcVariancePartitioning(obj_invasive)

## 4.4 Shared climate coverage ===========
# shared climate range
clim_range_shared <- max(integrated_curves$BO2_curvelmax_bdmean$cRange) -
  min(integrated_curves$BO2_curvelmax_bdmean$cRange)

# total climate range
clim_range_total <- max(c(integrated_curves$BO2_curvelmax_bdmean$invasive$xSeq,
                          integrated_curves$BO2_curvelmax_bdmean$native$xSeq)) -
  min(c(integrated_curves$BO2_curvelmax_bdmean$invasive$xSeq,
        integrated_curves$BO2_curvelmax_bdmean$native$xSeq))

#Proportion
clim_shared_proportion <- clim_range_shared/clim_range_total
clim_shared_proportion


calculate_shared_climate_coverage(integrated_curves$BO2_curvelmax_bdmean)

all_climate_coverage <- Map(calculate_shared_climate_coverage, integrated_curves)
all_climate_coverage
mean(as.numeric(all_climate_coverage))
sd(as.numeric(all_climate_coverage))


## 4.5 Prevalence ===========
# Needs to run section 4.3

# Prevalence in native model
native_occurence_na_rm <- na.omit(native_occurences)
prevalence_native <- sum(native_occurence_na_rm == 1)/ sum(native_occurence_na_rm == 0)
prevalence_native

# Prevalence in invasive model
invasive_occurences <- replace(occurence_grid_invasive$occurrence, is.na(invasive_fundamental_estimate), NA)


invasive_occurence_na_rm <- na.omit(invasive_occurences)
prevalence_invasive <- sum(invasive_occurence_na_rm == 1)/ sum(invasive_occurence_na_rm == 0)
prevalence_invasive

## 4.6 Uncertainty ======
native_uncertainty <- as.data.frame(obj_native$spatialPredictions[4])
summary(native_uncertainty)
mean(native_uncertainty$uncertaintyEst)
sd(native_uncertainty$uncertaintyEst)

invasive_uncertainty <- as.data.frame(obj_invasive$spatialPredictions[4])
summary(invasive_uncertainty)
mean(invasive_uncertainty$uncertaintyEst)
sd(invasive_uncertainty$uncertaintyEst)
