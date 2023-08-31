# Battini et al. (2019) -------------------------------
# https://dx.doi.org/10.3354/meps12878
#
# 1 Preparation ---------------------
#
## 1.1 Data sources ==================
#
# Source of occurence data:
# The occurence data is available in the supplementary material for the paper
#
# Source of climate data:
# Bio-ORACLE, https://www.bio-oracle.org/index.php
# The data available in Bio-ORACLE are documented in two peer reviewed articles that you should cite:
#   
#   Tyberghein L, Verbruggen H, Pauly K, Troupin C, Mineur F, De Clerck O (2012) Bio-ORACLE: A global environmental dataset for marine species distribution modelling. Global Ecology and Biogeography, 21, 272–281.
# [access publication]   [supporting information]
# 
# Assis, J., Tyberghein, L., Bosh, S., Verbruggen, H., Serrão, E. A., & De Clerck, O. (2017). Bio-ORACLE v2.0: Extending marine data layers for bioclimatic modelling. Global Ecology and Biogeography.
# [access publication]   [supporting information]


## 1.2 Package library =============
# Packages always used:
library(tidyverse)
library(terra)
library(here)
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
library(parzer)
library(readxl)


library(gridExtra)



## 1.3 Function files ========
source(paste0(here(), "/Joe_functions.R"))
source(paste0(here(), "/Isak_functions.R"))
source(paste0(here(), "/resp_curves_update.R"))


## 1.4 Import occurrence ========
maculata_occurence <- read_xlsx(path = paste0(here(), "/data/Battini.xlsx"))

maculata_occurence_numeric = maculata_occurence |> 
  mutate(Lon = parzer::parse_lon(Longitude),
         Lat = parzer::parse_lat(Latitude))

# Clean the data if necessary
maculata_occurence_clean <- maculata_occurence_numeric |>
  mutate(Latitude = as.double(Lat), Longitude = as.double(Lon), occurrenceStatus = 1) |>
  select(Latitude, Longitude, occurrenceStatus)

# ## 1.5 Import climate data ========
# # This example data comes from Worldclim
# # C:\Users\Isak\OneDrive\Dokumenter\R_projects\Masters_niche_shift\data\Battini_clim
# bio1 <- rast(paste0(here(),"/data/Battini_clim/Present.Benthic.Mean.Depth.Current.Velocity.Max.tif.BOv2_1.tif"))
# bio2 <- rast(paste0(here(),"/data/Battini_clim/Present.Benthic.Mean.Depth.Current.Velocity.Range.tif.BOv2_1.tif"))
# bio3 <- rast(paste0(here(),"/data/Battini_clim/Present.Benthic.Mean.Depth.Dissolved.oxygen.Min.tif"))
# bio4 <- rast(paste0(here(),"/data/Battini_clim/Present.Benthic.Mean.Depth.Dissolved.oxygen.Range.tif"))
# bio5 <- rast(paste0(here(),"/data/Battini_clim/Present.Benthic.Mean.Depth.Salinity.Min.tif"))
# bio6 <- rast(paste0(here(),"/data/Battini_clim/Present.Benthic.Mean.Depth.Salinity.Range.tif"))
# bio7 <- rast(paste0(here(),"/data/Battini_clim/Present.Benthic.Mean.Depth.Temperature.Mean.tif"))
# bio8 <- rast(paste0(here(),"/data/Battini_clim/Present.Benthic.Mean.Depth.Temperature.Range.tif"))


# Load package
library(sdmpredictors)

# Explore datasets in the package
# list_datasets()

# Explore layers in a dataset
# x <-list_layers()

# Download specific layers to the current directory
clim <- load_layers(c("BO_bathymean", "BO2_curvelmax_bdmean", "BO2_curvelrange_bdmean", "BO2_dissoxmin_bdmean",
                      "BO2_dissoxrange_bdmean", "BO2_tempmean_bdmean", "BO2_temprange_bdmean",
                      "BO2_salinitymin_bdmean", "BO2_salinityrange_bdmean"))

clim <- clim |>
  rast()

# plot(clim$BO2_dissoxmin_bdmean)

# Check layer statistics
# layer_stats()

# Add all climate variable into a single terra object
# clim <- terra::`add<-`(bio1, c(bio2, bio5, bio8, bio9, bio10, bio12, bio13, bio14, bio18))

# Make a list of names for each climate variable (make sure to have the correct order)
clim_names <- c("Mean depth (m)", "Current velocity maximum (m/s)",
                "Current velocity range (m/s)", "Dissolved oxygen minimum (μmol/m^3)",
                "Dissolved oxygen range (μmol/m^3)", "Temperature mean (\u00B0C)", "Temperature range (\u00B0C)",
                "Salinity minimum (PSS)", "Salinity range (PSS)")

clim_codes <- names(clim)

## 1.6 Define expanded prediction range ==============
# 
# # This defines the geographical range of the predictions beyond the initial model prediction (W, E, S, N)
# model_range_extended <- c(-25.8, 34, 36, 70)
# 
# # This will be used in both invasive and invasive sections
# # The function is available in the Isak_functions.R file
# # (aggregation factor is only intended for trial runs, reducing calculation time and required memory.
# # Default makes no aggregation)
# climate_extended_prediction_range <- define_range_and_aggregation(clim, model_range_extended,
#                                                                   aggregation_factor = 1
# ) |>
#   raster::stack() |>
#   as("SpatialGridDataFrame")



## 1.7 native population ==============

# This defines the geographical range of the native model (W, E, S, N)
model_range_native <- c(139.7, 177.9, -47.9, -27.4)

# clean occurrence data
maculata_native <- maculata_occurence_clean |>
  filter(
    Longitude > model_range_native[1], Longitude < model_range_native[2],
    Latitude > model_range_native[3], Latitude < model_range_native[4]
  )

# Prepare data
climate_native <- define_range_and_aggregation(clim, model_range_native, aggregation_factor = 3)
#agg = 3
# X = 1

# The function is available in the Isak_functions.R file
climate_grid_native <- make_climate_grid(climate_native)

# The function is available in the Isak_functions.R file
occurence_grid_native <- make_occurence_grid(maculata_native, climate_native)

## 1.8 Invasive population ==============

# This defines the geographical range of the invasive model (W, E, S, N)
model_range_invasive <- c(-66.3, -53, -45, -34.1)

# clean occurrence data
maculata_invasive <- maculata_occurence_clean |>
  filter(
    Longitude > model_range_invasive[1], Longitude < model_range_invasive[2],
    Latitude > model_range_invasive[3], Latitude < model_range_invasive[4]
  )


# Prepare data
climate_invasive <- define_range_and_aggregation(clim, model_range_invasive)

# The function is available in the Isak_functions.R file
climate_grid_invasive <- make_climate_grid(climate_invasive)

# The function is available in the Isak_functions.R file
occurence_grid_invasive <- make_occurence_grid(maculata_invasive, climate_invasive)


## 1.9 Plot occurrences ==============

native_coastline <- ne_countries(continent = c("oceania"), returnclass = "sf", scale = 50) 

occ_nat <- ggplot() +
  geom_sf(data = native_coastline, fill = "white") +
  geom_point(data = maculata_native, aes(Longitude, Latitude), size = 1) +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(model_range_native[1],  model_range_native[2])) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(model_range_native[3], model_range_native[4])) +
  labs(title = "Native occurences", x = NULL, y = NULL)+
  theme(panel.background = element_rect(fill = "lightblue", colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))
occ_nat


invasive_coastline <- ne_countries(continent = "south america", returnclass = "sf", scale = 50)

occ_inv <- ggplot() +
  geom_sf(data = invasive_coastline, fill = "white") +
  geom_point(data = maculata_invasive, aes(Longitude, Latitude)) +
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

## 2.1 native ============

# This function requires that the functions: sdmINLASpatRandEff, updateProgressBar,
# inlaMeshToSpatialPolygons, createProgressUpdater
# All required functions are available in the Joe_functions.R file
# Requires the INLA package to be loaded
mod_native <- sdmINLASpatRandEff(occurence_grid_native, climate_grid_native,
                                 meshParameters = list(max.n.strict = -1),
                                 createGeoTIFF = FALSE, outFolder = getwd(),
                                 myName = "maculata_native_X"
)
#agg = 3
# X = 1
 
## 2.1 Invasive ============
mod_invasive <- sdmINLASpatRandEff(occurence_grid_invasive, climate_grid_invasive,
                                   meshParameters = list(max.n.strict = -1),
                                   createGeoTIFF = FALSE,
                                   outFolder = getwd(),
                                   myName = "maculata_Invasive_agg"
)
#agg = 1


# 3 Model Results -------------

## 3.1 Native Results =========
# Retrive model results
obj_native <- readRDS(file = "maculata_native_agg.rds")

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
#   labs(title = "native occurrence prediction", x = NULL, y = NULL) +
#   theme_bw()

native_occurence_estimate <- as.data.frame(obj_native$spatialPredictions[1])
summary(native_occurence_estimate)

x_nat <- ggplot(native_occurence_estimate) +
  geom_tile(aes(s1, s2, fill = meanEst / max(native_occurence_estimate$meanEst))) +
  coord_sf() +
  scale_fill_viridis_c(
    # limits = c(0, .7),
    option = "magma",
    begin = 0.12, end = 1, name = "Relative \noccurence \nprobability"
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(model_range_native[1], model_range_native[2])) +
  scale_y_continuous(expand = c(0, 0), limits = c(model_range_native[3], model_range_native[4])) +
  labs(title = "a) Native distribution estimate", x = NULL, y = NULL) +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))

x_nat

native_suitability <- as.data.frame(obj_native$spatialPredictions[6]) |> 
  mutate(meanPred = invLogit(meanClimatePred))
summary(native_suitability)

y_nat <- ggplot(native_suitability) +
  geom_tile(aes(s1, s2, fill = meanPred / max(native_suitability$meanPred))) +
  coord_sf() +
  scale_fill_viridis_c(
    limits = c(0, .4),
    option = "magma",
    begin = 0.12, end = 1, name = "Relative \n suitability"
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(model_range_native[1], model_range_native[2])) +
  scale_y_continuous(expand = c(0, 0), limits = c(model_range_native[3], model_range_native[4])) +
  labs(title = "b) Native habitat estimate", x = NULL, y = NULL) +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))
y_nat

native_uncertainty <- as.data.frame(obj_native$spatialPredictions[4])
summary(native_uncertainty)

z_nat <-
  ggplot(native_uncertainty) +
  geom_tile(aes(s1, s2, fill = uncertaintyEst)) +
  coord_sf() +
  scale_fill_viridis_c(
    # limits = c(0, .1),
    option = "magma",
    begin = 0.12, end = 1, name = "95% Credible \ninterval width"
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(model_range_native[1], model_range_native[2])) +
  scale_y_continuous(expand = c(0, 0), limits = c(model_range_native[3], model_range_native[4])) +
  labs(title = "c) Native uncertainty", x = NULL, y = NULL) +
  theme(plot.margin=grid::unit(c(0,0,2,0), "mm"), panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))
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

### 3.1.3 Predicted occurrence probability in an native and native ranges  #########
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
  labs(title = "native niche prediction in the native range", x = NULL, y = NULL) +
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
  labs(title = "native niche prediction in the invasive range", x = NULL, y = NULL) +
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
maculata_AUC_native <- auc(roc_object_native)
maculata_AUC_native


## 3.2 Invasive Results =========
# Retrive model results
obj_invasive <- readRDS(file = "maculata_invasive_agg.rds")

### 3.2.1 Mean occurrence probability map #########
# invasive_occurence_estimate <- as.data.frame(obj_invasive$spatialPredictions[1])
# 
# ggplot(invasive_occurence_estimate) +
#   geom_tile(aes(s1, s2, fill = meanEst)) +
#   coord_sf() +
#   scale_fill_viridis_c(
#     limits = c(0, 0.35), option = "magma",
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
    # limits = c(0, 0.9),
    option = "magma",
    begin = 0.12, end = 1, name = "Relative \noccurence \nprobability"
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(model_range_invasive[1], model_range_invasive[2])) +
  scale_y_continuous(expand = c(0, 0), limits = c(model_range_invasive[3], model_range_invasive[4])) +
  labs(title = "a) Invasive distribution estimate", x = NULL, y = NULL) +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))
# +
#   geom_point(data = maculata_invasive, aes(Longitude, Latitude), size = 0.7, color = "green")
x_inv


invasive_suitability <- as.data.frame(obj_invasive$spatialPredictions[6]) |> 
  mutate(meanPred = invLogit(meanClimatePred))
summary(invasive_suitability)

y_inv <- ggplot(invasive_suitability) +
  geom_tile(aes(s1, s2, fill = meanPred / max(invasive_suitability$meanPred))) +
  coord_sf() +
  scale_fill_viridis_c(
    # limits = c(0, .5),
    option = "magma",
    begin = 0.12, end = 1, name = "Relative \nSuitability"
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(model_range_invasive[1], model_range_invasive[2])) +
  scale_y_continuous(expand = c(0, 0), limits = c(model_range_invasive[3], model_range_invasive[4])) +
  labs(title = "b) Invasive habitat estimate", x = NULL, y = NULL) +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))
y_inv

invasive_uncertainty <- as.data.frame(obj_invasive$spatialPredictions[4])
summary(invasive_uncertainty)

z_inv <-
  ggplot(invasive_uncertainty) +
  geom_tile(aes(s1, s2, fill = uncertaintyEst)) +
  coord_sf() +
  scale_fill_viridis_c(
    # limits = c(0, .1),
    option = "magma",
    begin = 0.12, end = 1, name = "95% Credible \ninterval width"
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(model_range_invasive[1], model_range_invasive[2])) +
  scale_y_continuous(expand = c(0, 0), limits = c(model_range_invasive[3], model_range_invasive[4])) +
  labs(title = "Invasive uncertainty", x = NULL, y = NULL) +
  theme(plot.margin=grid::unit(c(0,0,2,0), "mm"), panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))
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



### 3.2.5 AUC ########
observations_invasive <- occurence_grid_invasive$occurrence
predictions_invasive <- obj_invasive$spatialPredictions$meanEst
#Error!!! roc_object_invasive <- roc(observations_invasive, predictions_invasive)
# Setting levels: control = 0, case = 1
# Error in roc.default(observations_invasive, predictions_invasive) : 
#   Response and predictor must be vectors of the same length.
roc_object_invasive <- roc(observations_invasive, predictions_invasive)
maculata_AUC_invasive <- auc(roc_object_invasive)

observations_invasive <- occurence_grid_invasive$occurrence
predictions_invasive <- obj_invasive$spatialPredictions$meanEst
roc_object_invasive <- roc(observations_invasive, predictions_invasive)
maculata_AUC_invasive <- auc(roc_object_invasive)
maculata_AUC_invasive

# 4 Niche comparison ---------

## 4.1 Objects needed from other sections (section 1 must be run first) ==========
obj_native <- readRDS(file = "maculata_native_agg.rds")
preds_native <- obj_native$responsePredictions
pars_native <- getPars(obj_native, climate_grid_native)

obj_invasive <- readRDS(file = "maculata_invasive_agg.rds")
preds_invasive <- obj_invasive$responsePredictions
pars_invasive <- getPars(obj_invasive, climate_grid_invasive)


## 4.2 Niche overlap =============

### 4.2.1 Plot response curves =============

#Raw response curves
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

response_curve_native <- with(integrated_curves$BO_bathymean$native,
                              respCurve(xSeq, datSel, pars, clim_codes[1], integral))
response_curve_invasive <- with(integrated_curves$BO_bathymean$invasive,
                                respCurve(xSeq, datSel, pars, clim_codes[1], integral))

with(integrated_curves$BO_bathymean$invasive,
     plot(xSeq, response_curve_invasive,
          xlim = range(c(xSeq, integrated_curves$BO_bathymean$native$xSeq)), col = "red", 
          ylim = c(0, max(c(response_curve_native, response_curve_invasive))),
          type = "l", xlab = "",
          ylab = "",
          axes = FALSE))

with(integrated_curves$BO_bathymean$native, lines(xSeq, response_curve_native, col = "blue"))
abline(v = integrated_curves$BO_bathymean$cRange, lwd = 2, lty = 2)
legend("topleft", paste0("Overlap = ", round(all_overlaps[2], digits = 3)), bty = "o")
# legend("topleft", c("Native", "Invasive"), fill = c("blue","red"), bty = "o")
box()
axis(2, las = 1)
axis(2, at = mean(par("usr")[3:4]), labels = "Occurrence density", tick = FALSE, line = 2.5)


par(mai = c(1.3, 1.3, 0.1, 1))

vioplot(BO_bathymean~occurrence, data = climate_at_occurrences_native,
        side = "right", horizontal = TRUE, col = "blue",
        ylim = range(c(integrated_curves$BO_bathymean$invasive$xSeq,
                       integrated_curves$BO_bathymean$native$xSeq)),
        ylab = "", xlab = clim_names[2], xaxt = "n")

abline(h = 1, lwd = 2, lty = 2)

vioplot(BO_bathymean~occurrence, data = climate_at_occurrences_invasive,
        side = "left", horizontal = TRUE, col = "red", add = TRUE)

# dev.off()


# Single curve
niche_overlap_plot(integrated_curves$BO_bathymean, clim_codes[1], clim_names[1], climate_at_occurrences_native, climate_at_occurrences_invasive, "test", all_overlaps[1])


# All curves
Map(niche_overlap_plot, integrated_curves, clim_codes, clim_names, list(climate_at_occurrences_native), list(climate_at_occurrences_invasive), "Battini_agg", all_overlaps)


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