# Bidinger 2012 -------------------------------
# https://dx.doi.org/10.1111/j.1439-0418.2010.01598.x
# This template code has 3 sections: #1 Preparation, #2 Model, #3 Model Results, #4 Niche comparison
#
# 1 Preparation ---------------------
#
## 1.1 Data sources ==================
#
# Source of occurence data:
# Appendix of the original paper
#
# Source of climate data:
# WorldClim, 2.5 arc-minutes, all standard (19) WorldClim variables, https://www.worldclim.org/data/worldclim21.html


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
library(vioplot)
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
Harmonia_occurence <- read_xlsx(path = paste0(here(), "/data/Bidinger_data_v2.xlsx"))

# Clean the data if necessary
# Should include numeric colums for Latitude, Longitude and occurrenceStatus (which should always be 1)
Harmonia_occurence_clean <- Harmonia_occurence |>
  mutate(Latitude = as.double(latitude), Longitude = as.double(longitude),
         occurrenceStatus = 1) |>
  select(Latitude, Longitude, occurrenceStatus, range)


## 1.5 Import climate data ========
# This example data comes from Worldclim
bio1 <- rast("data/wc2.1_2.5m_bio/wc2.1_2.5m_bio_1.tif")
bio3 <- rast("data/wc2.1_2.5m_bio/wc2.1_2.5m_bio_3.tif")
bio5 <- rast("data/wc2.1_2.5m_bio/wc2.1_2.5m_bio_5.tif")
bio6 <- rast("data/wc2.1_2.5m_bio/wc2.1_2.5m_bio_6.tif")
bio7 <- rast("data/wc2.1_2.5m_bio/wc2.1_2.5m_bio_7.tif")
bio8 <- rast("data/wc2.1_2.5m_bio/wc2.1_2.5m_bio_8.tif")
bio9 <- rast("data/wc2.1_2.5m_bio/wc2.1_2.5m_bio_9.tif")
bio15 <- rast("data/wc2.1_2.5m_bio/wc2.1_2.5m_bio_15.tif")
bio16 <- rast("data/wc2.1_2.5m_bio/wc2.1_2.5m_bio_16.tif")
bio17 <- rast("data/wc2.1_2.5m_bio/wc2.1_2.5m_bio_17.tif")



# Add all climate variable into a single terra object
climfull <- c(bio1, bio3, bio5, bio6, bio7, bio8, bio9,
                                   bio15, bio16, bio17)

clim_codes <- names(climfull)

# Make a list of names for each climate variable (make sure to have the correct order)
clim_names <- c("Annual Mean Temperature (\u00B0C)",
                "Isothermality (BIO2/BIO7) (×100)", 
                "Max Temperature of Warmest Month (\u00B0C)",
                "Min Temperature of Coldest Month (\u00B0C)",
                "Temperature Annual Range (BIO5-BIO6)",
                "Mean Temperature of Wettest Quarter (\u00B0C)",
                "Mean Temperature of Driest Quarter (\u00B0C)",
                "Precipitation Seasonality (Coefficient of Variation)",
                "Precipitation of Wettest Quarter (mm)",
                "Precipitation of Driest Quarter (mm)")


## 1.6 Define expanded prediction range ==============

# # This defines the geographical range of the predictions beyond the initial model prediction (W, E, S, N)
# model_range_extended <- c(1, 2, 1, 2)
# 
# # This will be used in both native and invasive sections
# # The function is available in the Isak_functions.R file
# # (aggregation factor is only intended for trial runs, reducing calculation time and required memory.
# # Default makes no aggregation, ignore warning message: [aggregate] all values in argument 'fact' are 1, nothing to do)
# climate_extended_prediction_range <- define_range_and_aggregation(climfull, model_range_extended,
#                                                                   aggregation_factor = 1
# ) |>
#   raster::stack() |>
#   as("SpatialGridDataFrame")


## 1.7 Native population ==============

# This defines the geographical range of the native model (W, E, S, N)
model_range_native <- c(86.7, 153.7, 21.9, 57.7)

# clean occurrence data
Harmonia_native <- Harmonia_occurence_clean |>
  filter(
    Longitude > model_range_native[1], Longitude < model_range_native[2],
    Latitude > model_range_native[3], Latitude < model_range_native[4]
  ) |> 
  select(-range)


# Prepare data
climate_native <- define_range_and_aggregation(climfull, model_range_native, aggregation_factor = 3)

# The function is available in the Isak_functions.R file
climate_grid_native <- make_climate_grid(climate_native)

# The function is available in the Isak_functions.R file
occurence_grid_native <- make_occurence_grid(Harmonia_native, climate_native)



## 1.8 Invasive population ==============

# This defines the geographical range of the invasive model (W, E, S, N)
model_range_invasive <- c(-7.1, 16.5, 42.4, 58.2)
#-7, 17
# clean occurrence data
Harmonia_invasive <- Harmonia_occurence_clean |>
  filter(
    Longitude > model_range_invasive[1], Longitude < model_range_invasive[2],
    Latitude > model_range_invasive[3], Latitude < model_range_invasive[4]
  ) |> 
  select(-range)


# Prepare data
climate_invasive <- define_range_and_aggregation(climfull, model_range_invasive, aggregation_factor = 2)

# The function is available in the Isak_functions.R file
climate_grid_invasive <- make_climate_grid(climate_invasive)

# The function is available in the Isak_functions.R file
occurence_grid_invasive <- make_occurence_grid(Harmonia_invasive, climate_invasive)


## 1.9 Plot occurrences ==============

native_coastline <- ne_countries(continent = c("asia", "europe"), returnclass = "sf", scale = 50) 

occ_nat <- ggplot() +
  geom_sf(data = native_coastline, fill = "white") +
  geom_point(data = Harmonia_native, aes(Longitude, Latitude)) +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(model_range_native[1],  model_range_native[2])) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(model_range_native[3], model_range_native[4])) +
  labs(title = "Native occurences", x = NULL, y = NULL)+
  theme(panel.background = element_rect(fill = "lightblue", colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
occ_nat
### NOTE: added black border, replicate in all datasets #######


invasive_coastline <- ne_countries(continent = "europe", returnclass = "sf", scale = 50)

occ_inv <- ggplot() +
  geom_sf(data = invasive_coastline, fill = "white") +
  geom_point(data = Harmonia_invasive, aes(Longitude, Latitude)) +
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
                                 myName = "Harmonia_Native_agg"
)
#agg = 3

## 2.1 Invasive ============
mod_invasive <- sdmINLASpatRandEff(occurence_grid_invasive, climate_grid_invasive,
                                   meshParameters = list(max.n.strict = -1),
                                   createGeoTIFF = FALSE,
                                   outFolder = getwd(),
                                   myName = "Harmonia_Invasive_agg"
)
#agg = 2


# 3 Model Results -------------

## 3.1 Native Results =========
# Retrive model results
obj_native <- readRDS(file = "Harmonia_Native_agg.rds")

### 3.1.1 Mean occurrence probability map #########
native_occurence_estimate <- as.data.frame(obj_native$spatialPredictions[1])
summary(native_occurence_estimate)

x_nat <- ggplot(native_occurence_estimate) +
  geom_tile(aes(s1, s2, fill = meanEst / max(native_occurence_estimate$meanEst))) +
  coord_sf() +
  scale_fill_viridis_c(
    # limits = c(0.1, .6),
    option = "magma",
    begin = 0.12, end = 1, name = "Relative \n occurence \n probability"
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
    # limits = c(0, .1),
    option = "magma",
    begin = 0.12, end = 1, name = "Relative \n suitability"
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
    limits = c(0, 0.05),
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



# ###### testing #################
# 
# library(biscale)
# x <- bi_class(stl_race_income, x = pctWhite, y = medInc, style = "quantile", dim = 3)
# y <- stl_race_income
# 
# z <- left_join(native_occurence_estimate, native_uncertainty)
# zz <- bi_class(z, x = meanEst, y = uncertaintyEst, style = "quantile", dim = 3)
# 
# map <- ggplot() +
#   geom_tile(data = zz, mapping = aes(s1, s2, fill = bi_class), color = "white", linewidth = 0.1, show.legend = FALSE) +
#   bi_scale_fill(pal = "GrPink", dim = 3) +
#   labs(
#     title = "occurrence with uncertainty",
#     subtitle = "Gray Pink (GrPink) Palette"
#   ) +
#   bi_theme()
# 
# legend <- bi_legend(pal = "GrPink",
#                     dim = 3,
#                     xlab = "meanEst",
#                     ylab = "uncertainty",
#                     size = 8)
# library(cowplot)
# # combine map with legend
# finalPlot <- ggdraw() +
#   draw_plot(map, 0, 0, 1, 1) +
#   draw_plot(legend, 0.2, .65, 0.2, 0.2)
# 
# ################


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
    limits = c(0, .001), option = "magma",
    begin = 0.12, end = 1, name = "Suitability"
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(model_range_native[1], model_range_native[2])) +
  scale_y_continuous(expand = c(0, 0), limits = c(model_range_native[3], model_range_native[4])) +
  labs(title = "Native niche prediction in the native range", x = NULL, y = NULL) +
  theme_bw()

#### testing ##############

mmm <- as.data.frame(obj_native$spatialPredictions[6]) |> 
  mutate(meanPred = invLogit(meanClimatePred))

ggplot(mmm) +
  geom_tile(aes(s1, s2, fill = meanPred)) +
  coord_sf() +
  scale_fill_viridis_c(
    limits = c(0, .001), option = "magma",
    begin = 0.12, end = 1, name = "Occurence probability"
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(model_range_native[1], model_range_native[2])) +
  scale_y_continuous(expand = c(0, 0), limits = c(model_range_native[3], model_range_native[4])) +
  labs(title = "Native occurrence prediction", x = NULL, y = NULL) +
  theme_bw()

max(native_native_predicted_estimate$OccurrenceProb)
max(mmm$meanPred)
summary(native_native_predicted_estimate$OccurrenceProb)
summary(mmm$meanPred)


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
################

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
Harmonia_AUC_native <- auc(roc_object_native)
summary(Harmonia_AUC_native)



## 3.2 Invasive Results =========
# Retrive model results
obj_invasive <- readRDS(file = "Harmonia_invasive_agg.rds")

### 3.2.1 Mean occurrence probability map #########
invasive_occurence_estimate <- as.data.frame(obj_invasive$spatialPredictions[1])

x_inv <-
  ggplot(invasive_occurence_estimate) +
  geom_tile(aes(s1, s2, fill = meanEst / max(invasive_occurence_estimate$meanEst))) +
  coord_sf() +
  scale_fill_viridis_c(
    # limits = c(0, 0.45),
    option = "magma",
    begin = 0.12, end = 1, name = "Relative \noccurence \nprobability"
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(model_range_invasive[1], model_range_invasive[2])) +
  scale_y_continuous(expand = c(0, 0), limits = c(model_range_invasive[3], model_range_invasive[4])) +
  labs(title = "a) Invasive realised niche", x = NULL, y = NULL) +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), panel.background = element_rect(fill = "lightblue"))
x_inv

invasive_suitability <- as.data.frame(obj_invasive$spatialPredictions[6]) |> 
  mutate(meanPred = invLogit(meanClimatePred))
summary(invasive_suitability)

y_inv <- ggplot(invasive_suitability) +
  geom_tile(aes(s1, s2, fill = meanPred / max(invasive_suitability$meanPred))) +
  coord_sf() +
  scale_fill_viridis_c(
    # limits = c(0, .45),
    option = "magma",
    begin = 0.12, end = 1, name = "Relative \nSuitability"
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(model_range_invasive[1], model_range_invasive[2])) +
  scale_y_continuous(expand = c(0, 0), limits = c(model_range_invasive[3], model_range_invasive[4])) +
  labs(title = "b) Invasive fundamental niche", x = NULL, y = NULL) +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), panel.background = element_rect(fill = "lightblue"))
y_inv

invasive_uncertainty <- as.data.frame(obj_invasive$spatialPredictions[4])
summary(invasive_uncertainty)

z_inv <-
  ggplot(invasive_uncertainty) +
  geom_tile(aes(s1, s2, fill = uncertaintyEst)) +
  coord_sf() +
  scale_fill_viridis_c(
    limits = c(0.0, .1), option = "magma",
    begin = 0.12, end = 1, name = "Confidence \nwidth"
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

# preds_ and pars_ will be used in the calculation of niche overlap later

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
roc_object_invasive <- roc(observations_invasive, predictions_invasive)
Harmonia_AUC_invasive <- auc(roc_object_invasive)
summary(Harmonia_AUC_invasive)

# 4 Niche comparison ---------

## 4.1 Objects needed from other sections (section 1 must be run first) ==========
obj_native <- readRDS(file = "Harmonia_Native_agg.rds")
preds_native <- obj_native$responsePredictions
pars_native <- getPars(obj_native, climate_grid_native)

obj_invasive <- readRDS(file = "Harmonia_invasive_agg.rds")
preds_invasive <- obj_invasive$responsePredictions
pars_invasive <- getPars(obj_invasive, climate_grid_invasive)


## 4.2 Niche overlap =============

### 4.2.1 Plot response curves =============
# 
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
niche_overlap_plot(integrated_curves$wc2.1_2.5m_bio_1, clim_codes[1], clim_names[1], climate_at_occurrences_native, climate_at_occurrences_invasive, "test", all_overlaps)


# All curves
Map(niche_overlap_plot, integrated_curves, clim_codes, clim_names, list(climate_at_occurrences_native), list(climate_at_occurrences_invasive), "Bidinger", all_overlaps)


## 4.3  Climate/ Spatial effect contribution ==========
calcVariancePartitioning(obj_native)

calcVariancePartitioning(obj_invasive)

## 4.4 Shared climate coverage ===========
# # shared climate range
# clim_range_shared <- max(integrated_curves$BO2_curvelmax_bdmean$cRange) -
#   min(integrated_curves$BO2_curvelmax_bdmean$cRange)
# 
# # total climate range
# clim_range_total <- max(c(integrated_curves$BO2_curvelmax_bdmean$invasive$xSeq,
#                           integrated_curves$BO2_curvelmax_bdmean$native$xSeq)) -
#   min(c(integrated_curves$BO2_curvelmax_bdmean$invasive$xSeq,
#         integrated_curves$BO2_curvelmax_bdmean$native$xSeq))
# 
# #Proportion
# clim_shared_proportion <- clim_range_shared/clim_range_total
# clim_shared_proportion
# 
# 
# calculate_shared_climate_coverage(integrated_curves$BO2_curvelmax_bdmean)

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