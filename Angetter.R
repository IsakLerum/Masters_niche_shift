# Angetter et al. (2011) -----------------
# https://doi.org/10.1111/j.1095-8312.2011.01780.x

# 1 Preparation ---------------------

# 1.1 Data sources ==================
# Source of occurence data:
# Gbif
# Reference code:
# GBIF.org (29 March 2023) GBIF Occurrence Download  https://doi.org/10.15468/dl.jh2rpq

# add vertnet reference

# Source of climate data:
# Worldclim, 30 arcsec, bio: 1, 2, 5, 8, 9, 10, 12, 13, 14, 18.
# link = https://www.worldclim.org/data/worldclim21.html

# original

# 1.2 Package library =============
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

library(hypervolume)

# library(tidylog)
# library(sp)
# # inla.upgrade() # for the stable version
# library(rgdal)
# library(vegan)
# library(ggfortify)
# library(mapproj)
# library(ade4)
# library(ecospat)
# # library(dismo)

# Packages used, but not loaded:
# library(raster)
# library(spatialEco)
# library(isocat)

# Packages used for data prep:
library(readxl)
library(here)

## 1.3 Function files ========
source(paste0(here(), "/Joe_functions.R"))
source(paste0(here(), "/Isak_functions.R"))
source(paste0(here(), "/resp_curves_update.R"))

#Note new occurrence data ----------
# range(-104.76345 9.42033,-68.8774 9.42033,-68.8774 43.83465,-104.76345 43.83465,-104.76345 9.42033)
# time (1911 - 2011)



## 1.4 Import occurrence ========
anolis_occurence <- read_xlsx(path = paste0(here(), "/data/anolis_sagrei_v2.xlsx"))


anolis_occurence_clean <- anolis_occurence |>
  mutate(Latitude = as.double(decimalLatitude), Longitude = as.double(decimalLongitude)) |>
  select(Latitude, Longitude, occurrenceStatus) |>
  replace(
    is.character(anolis_occurence$occurrenceStatus) == TRUE,
    anolis_occurence$occurrenceStatus <- 1
  )

Vertnet <- read.csv(file = "data/VertNet_Reptilia_Oct2015/VertNet_Reptilia_Oct2015.csv")

anolis_occurence2 <- Vertnet |> 
  filter(scientificname == "Anolis sagrei", year < 2011)

anolis_occurence_clean2 <- anolis_occurence2 |>
  mutate(Latitude = as.double(decimallatitude), Longitude = as.double(decimallongitude)) |>
  select(Latitude, Longitude, occurrenceStatus) |>
  replace(is.character(anolis_occurence2$occurrenceStatus) == TRUE, anolis_occurence2$occurrenceStatus <- 1)

anolis_occurence_final <- full_join(anolis_occurence_clean, anolis_occurence_clean2)

## 1.5 Import climate data ========
bio1 <- rast(paste0(here(), "/data/wc2.1_30s_bio/wc2.1_30s_bio_1.tif"))
bio2 <- rast(paste0(here(), "/data/wc2.1_30s_bio/wc2.1_30s_bio_2.tif"))
bio5 <- rast(paste0(here(), "/data/wc2.1_30s_bio/wc2.1_30s_bio_5.tif"))
bio8 <- rast(paste0(here(), "/data/wc2.1_30s_bio/wc2.1_30s_bio_8.tif"))
bio9 <- rast(paste0(here(), "/data/wc2.1_30s_bio/wc2.1_30s_bio_9.tif"))
bio10 <- rast(paste0(here(), "/data/wc2.1_30s_bio/wc2.1_30s_bio_10.tif"))
bio12 <- rast(paste0(here(), "/data/wc2.1_30s_bio/wc2.1_30s_bio_12.tif"))
bio13 <- rast(paste0(here(), "/data/wc2.1_30s_bio/wc2.1_30s_bio_13.tif"))
bio14 <- rast(paste0(here(), "/data/wc2.1_30s_bio/wc2.1_30s_bio_14.tif"))
bio18 <- rast(paste0(here(), "/data/wc2.1_30s_bio/wc2.1_30s_bio_18.tif"))


climfull <- c(bio1, bio2, bio5, bio8, bio9, bio10, bio12, bio13, bio14, bio18)


# Make a list of names for each climate variable (make sure to have the correct order)
clim_names <- c("Annual Mean Temperature (\u00B0C)",
                "Mean Diurnal Range (Mean of monthly (max temp - min temp, \u00B0C))",
                "Max Temperature of Warmest Month (\u00B0C)", "Mean Temperature of Wettest Quarter (\u00B0C)",
                "Mean Temperature of Driest Quarter (\u00B0C)", "Mean Temperature of Warmest Quarter (\u00B0C)",
                "Annual Precipitation (mm)", "Precipitation of Wettest Month (mm)",
                "Precipitation of Driest Month (mm)", "Precipitation of Warmest Quarter (mm)")

clim_codes <- names(climfull)

## 1.6 Define expanded prediction range ==============
# range(-104.76345 ,-68.8774, 9.42033, 43.83465)
model_range_extended <- c(-105, -67.5, 12, 40)
# -110, -67.5, 8, 44
# -105, -67.5, 12, 40
climate_extended_prediction_range <- define_range_and_aggregation(climfull, model_range_extended,
  aggregation_factor = 1
) |>
  raster::stack() |>
  as("SpatialGridDataFrame")


## 1.7 Native population ==============

model_range_native <- c(-86, -73.7, 19.4, 27.5)

# clean occurrence data
anolis_native <- anolis_occurence_final |>
  filter(Longitude > model_range_native[1], Longitude < model_range_native[2],
    Latitude > model_range_native[3], Latitude < model_range_native[4]
  )

native_exclusion_range1 <- c(-84, -79.5, 24, 27.5)
native_exclusion_range2 <- c(-81, -79, 19.4, 20)

anolis_native_exclusion1 <- anolis_occurence_final |>
  filter(Longitude > native_exclusion_range1[1], Longitude < native_exclusion_range1[2],
         Latitude > native_exclusion_range1[3], Latitude < native_exclusion_range1[4])

anolis_native_exclusion2 <- anolis_occurence_final |>
  filter(Longitude > native_exclusion_range2[1], Longitude < native_exclusion_range2[2],
         Latitude > native_exclusion_range2[3], Latitude < native_exclusion_range2[4])

anolis_native_merge <- anti_join(anolis_native, anolis_native_exclusion1)

anolis_native_final <- anti_join(anolis_native_merge, anolis_native_exclusion2)

# Prepare data
climate_native <- define_range_and_aggregation(climfull, model_range_native,
  aggregation_factor = 1
)

climate_native_exclusion1 <- define_range_and_aggregation(climfull, native_exclusion_range1)

climate_native_exclusion_null1 <- setValues(climate_native_exclusion1, values = NA)

climate_native_merge <- terra::merge(climate_native_exclusion_null1, climate_native)

climate_native_exclusion2 <- define_range_and_aggregation(climfull, native_exclusion_range2)

climate_native_exclusion_null2 <- setValues(climate_native_exclusion2, values = NA)

climate_native_final <- terra::merge(climate_native_exclusion_null2, climate_native_merge) |> 
  raster::aggregate(4)

# X 3
# base 4
climate_grid_native <- make_climate_grid(climate_native_final)

occurence_grid_native <- make_occurence_grid(anolis_native_final, climate_native_final)


## 1.8 Invasive population ==============

model_range_invasive <- c(-99.6, -75.6, 24.2, 37.9)

# clean occurrence data
anolis_invasive <- anolis_occurence_final |>
  filter(
    Longitude > model_range_invasive[1], Longitude < model_range_invasive[2],
    Latitude > model_range_invasive[3], Latitude < model_range_invasive[4]
  )

anolis_invasive_exclusion_range <- c(-79.5, -72, 24.2, 28)

anolis_native_exclusion <- anolis_occurence_final |>
  filter(Longitude > anolis_invasive_exclusion_range[1], Longitude < anolis_invasive_exclusion_range[2],
         Latitude > anolis_invasive_exclusion_range[3], Latitude < anolis_invasive_exclusion_range[4])


anolis_invasive_final <- anti_join(anolis_invasive, anolis_native_exclusion)


# Prepare data
climate_invasive <- define_range_and_aggregation(climfull, model_range_invasive,
  aggregation_factor = 1
)

climate_invasive_exclusion <- define_range_and_aggregation(climfull, anolis_invasive_exclusion_range)

climate_invasive_exclusion_null <- setValues(climate_invasive_exclusion, values = NA)

climate_invasive_new <- terra::merge(climate_invasive_exclusion_null, climate_invasive) |> 
  raster::aggregate(8)
# X 6
# 8

climate_grid_invasive <- make_climate_grid(climate_invasive_new)

occurence_grid_invasive <- make_occurence_grid(anolis_invasive, climate_invasive_new)


## 1.9 Plot occurrences ==============

coastline <- ne_countries(continent = "north america", returnclass = "sf", scale = 50)

ggplot() +
  geom_sf(data = coastline, fill = "white") +
  geom_point(data = anolis_native_final, aes(Longitude, Latitude, color = "Native")) +
  geom_point(data = anolis_invasive_final, aes(Longitude, Latitude, color = "Invasive")) +
  scale_x_continuous(expand = c(0, 0), limits = c(-102, -70)) +
  scale_y_continuous(expand = c(0, 0), limits = c(13, 35)) +
  labs(title = "Species occurences (2011)", x = NULL, y = NULL, color = "Population") +
  theme(panel.background = element_rect(fill = "lightblue", colour = "lightblue"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))

# 2 Model --------------

## 2.1 Native ============

mod_native <- sdmINLASpatRandEff(occurence_grid_native, climate_grid_native,
  meshParameters = list(max.n.strict = -1),
  createGeoTIFF = FALSE, outFolder = getwd(),
  myName = "Anolis_sagrei_Native_agg"
)
#agg_X = 3
#agg = 4

## 2.1 Invasive ============
mod_invasive <- sdmINLASpatRandEff(occurence_grid_invasive, climate_grid_invasive,
  meshParameters = list(max.n.strict = -1),
  createGeoTIFF = FALSE,
  outFolder = getwd(),
  myName = "Anolis_sagrei_Invasive_agg"
)
#agg_X = 6
#agg = 8

# 3 Model Results -------------

## 3.1 Native Results =========
# Retrive model results
obj_native <- readRDS(file = "Anolis_sagrei_Native_agg.rds")

### 3.1.1 Mean occurrence probability map #########
# native_occurence_estimate <- as.data.frame(obj_native$spatialPredictions[1])
# 
# ggplot(native_occurence_estimate) +
#   geom_tile(aes(s1, s2, fill = meanEst)) +
#   coord_sf() +
#   scale_fill_viridis_c(
#     limits = c(0, 0.35), option = "magma",
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
    limits = c(0, .6),
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
    limits = c(0, .6),
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
    limits = c(0., .4),
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

### 3.1.3 Predicted occurrence probability in an expanded range #########
predict_native <- predictSD(obj_native, climate_extended_prediction_range,
  origClimateData = climate_grid_native
)
native_predicted_estimate <- as.data.frame(predict_native[1])

ggplot(native_predicted_estimate, aes(s1, s2)) +
  geom_tile(aes(fill = OccurrenceProb)) +
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
# preds_ and pars_ will be used in the calculation of niche overlap later
preds_native <- obj_native$responsePredictions

pars_native <- getPars(obj_native, climate_grid_native)

# Plot individual response curves
map2(preds_native, clim_names, plot_response_curve)

### 3.1.5 AUC ########
observations_native <- occurence_grid_native$occurrence
predictions_native <- obj_native$spatialPredictions$meanEst
roc_object_native <- roc(observations_native, predictions_native)
anolis_AUC_native <- auc(roc_object_native)
anolis_AUC_native


## 3.2 Invasive Results =========
# Retrive model results
obj_invasive <- readRDS(file = "Anolis_sagrei_invasive_agg.rds")

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

x_inv <-
  ggplot(invasive_occurence_estimate) +
  geom_tile(aes(s1, s2, fill = meanEst / max(invasive_occurence_estimate$meanEst))) +
  coord_sf() +
  scale_fill_viridis_c(
    # limits = c(0, 0.3),
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
    # limits = c(0, .5),
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
    # limits = c(0.1, 0.4),
    option = "magma",
    begin = 0.12, end = 1, name = "95% Credible \ninterval width"
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(model_range_invasive[1], model_range_invasive[2])) +
  scale_y_continuous(expand = c(0, 0), limits = c(model_range_invasive[3], model_range_invasive[4])) +
  labs(title = "c) Invasive uncertainty", x = NULL, y = NULL) +
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

### 3.2.3 Predicted occurrence probability in an expanded range #########
predict_invasive <- predictSD(obj_invasive, climate_extended_prediction_range,
  origClimateData = climate_grid_invasive
)
invasive_predicted_estimate <- as.data.frame(predict_invasive[1])

ggplot(invasive_predicted_estimate, aes(s1, s2)) +
  geom_tile(aes(fill = OccurrenceProb)) +
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
# preds_ and pars_ will be used in the calculation of niche overlap later

preds_invasive <- obj_invasive$responsePredictions

pars_invasive <- getPars(obj_invasive, climate_grid_invasive)

# Plot individual response curves
map2(preds_invasive, clim_names, plot_response_curve)


### 3.2.5 AUC ########
observations_invasive <- occurence_grid_invasive$occurrence
predictions_invasive <- obj_invasive$spatialPredictions$meanEst
roc_object_invasive <- roc(observations_invasive, predictions_invasive)
anolis_AUC_invasive <- auc(roc_object_invasive)
anolis_AUC_invasive


# 4 Niche comparison ---------


## 4.1 Objects needed from other sections (section 1 must be run first) ==========
obj_native <- readRDS(file = "Anolis_sagrei_Native_agg.rds")
preds_native <- obj_native$responsePredictions
pars_native <- getPars(obj_native, climate_grid_native)

obj_invasive <- readRDS(file = "Anolis_sagrei_invasive_agg.rds")
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
niche_overlap_plot(integrated_curves$wc2.1_30s_bio_1, clim_codes[1], clim_names[1], climate_at_occurrences_native, climate_at_occurrences_invasive, "test", all_overlaps[1])


# All curves
Map(niche_overlap_plot, integrated_curves, clim_codes, clim_names, list(climate_at_occurrences_native), list(climate_at_occurrences_invasive), "Angetter", all_overlaps)




## 4.3  Climate/ Spatial effect contribution ==========
calcVariancePartitioning(obj_native)

calcVariancePartitioning(obj_invasive)

## 4.4 Shared climate coverage ===========
# shared climate range
clim_range_shared <- max(integrated_curves$wc2.1_30s_bio_1$cRange) -
  min(integrated_curves$wc2.1_30s_bio_1$cRange)

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