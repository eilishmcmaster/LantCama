# Example of building and visualizing an SDM
library(raster)
library(dismo)
library(sp)
library(sf)

# rarify to 5 km to prevent issues like overfitting and inaccurate inflation of model outcomes arising from spatial autocorrelation

adhikari_2024_bioclim <- c('Bio1','Bio2','Bio3','Bio12','Bio13','Bio14','Biome') #use these variables
# . Species distribution modeling of L. camara with RF was performed using the Biomod2 package, version 4.2–4 (Thuiller et al., 2023). 
# future scenarios SSP2–4.5 and SSP5–8.5 

occurrence_A_per <- m2[which(m2$cluster=="A"),c('lat','long')] %>% na.omit()

# Convert the data frame to an sf object with WGS84 (lat/lon) CRS
occurrence_sf <- st_as_sf(occurrence_A_per, coords = c("long", "lat"), crs = 4326)

# Transform to a projected CRS (e.g., UTM) for accurate distance calculations
occurrence_sf <- st_transform(occurrence_sf, crs = 3577)  # Replace 32633 with the UTM zone that covers your study area


# Round coordinates to the nearest 5000 meters (5 km)
occurrence_sf <- occurrence_sf %>%
  mutate(
    # Rounding longitude and latitude coordinates to the nearest 5000 meters (5 km)
    long_rounded = round(st_coordinates(.)[, 1] / 5000) * 5000,
    lat_rounded = round(st_coordinates(.)[, 2] / 5000) * 5000
  )

# Create a new sf object with the rounded coordinates
occurrence_thinned_sf <- occurrence_sf %>%
  distinct(long_rounded, lat_rounded, .keep_all = TRUE) %>%
  st_transform(crs = 4326)

# # If you want to convert back to a data frame for further use
# df_occurrence_thinned <- st_drop_geometry(occurrence_thinned_sf)
# 
# plot(occurrence_thinned_sf)

##
library(biomod2)
library(geodata)
library(raster)

# Download current climate data for Australia using WorldClim
bio_curr <- geodata::worldclim_country('Australia', var = 'bio', res = 5, path = 'worldclim/', download = FALSE)

# Extract only the desired bioclimatic variables (Bio1, Bio2, Bio3, Bio12, Bio13, Bio14, Biome)
selected_vars <- c('wc2.1_30s_bio_1', 'wc2.1_30s_bio_2', 'wc2.1_30s_bio_3', 'wc2.1_30s_bio_12', 'wc2.1_30s_bio_13', 'wc2.1_30s_bio_14')

# Subset the data to include only the chosen variables
bio_curr_selected <- subset(bio_curr, selected_vars)

# Sample occurrence data format (make sure to replace this with your actual data)
# Replace this with your actual data import step
occurrence_data <-occurrence_thinned_sf$geometry %>% st_coordinates() %>% as.data.frame()
colnames(occurrence_data) <- c('longitude','latitude')
coordinates(occurrence_data) <- ~longitude + latitude
proj4string(occurrence_data) <- CRS("+proj=longlat +datum=WGS84")


# Define the model options and training data
biomod_data <- biomod2::BIOMOD_FormatingData(
  resp.var = occurrence_data,
  expl.var = bio_curr_selected,
  resp.xy = coordinates(occurrence_data),
  resp.name = 'L_Camara',  # Name of the species
  PA.nb.rep = 1,
  PA.strategy = 'random'
)


# Run the SDM model
biomod_model <- biomod2::BIOMOD_Modeling(
  bm.format = biomod_data,
  models = c('RF'),  # Random Forest model
  CV.nb.rep = 2,
  CV.perc = 0.8,
  OPT.strategy = 'bigboss',
  metric.eval = c('TSS','ROC'),
  var.import = 2,
  seed.val = 42)


biomod_model

# Get evaluation scores & variables importance
get_evaluations(biomod_model)
get_variables_importance(biomod_model)

occurrence_predictions <- BIOMOD_Projection(
  bm.mod=biomod_model,
  nb.cp=4,
  new.env = bio_curr_selected,  # Current climate data
  proj.name = 'Current_Predictions'
)


predicted_rasters <- get_predictions(occurrence_predictions)

normalized_raster <- predicted_rasters[[3]] / 1000  # Assuming 1000 is the max value

# Define Australia's extent (approximately)
australia_extent <- ext(110, 155, -45, -10)  # Longitude/Latitude

# Crop and mask the raster to Australia's extent
australia_raster <- crop(normalized_raster, australia_extent)

# Plot using ggplot2
ggplot() +
  tidyterra::geom_spatraster(data = australia_raster) +
  scale_fill_gradientn(colors = rev(terrain.colors(100)), limits = c(0, 1), na.value = 'lightblue') +  # Use the inverted color scheme
  theme_few() +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  geom_point(data = occurrence_thinned_sf %>% st_coordinates(), 
             mapping = aes(x = X, y = Y), 
             col = 'red', 
             size = 0.5) +
  labs(title = "Predicted Occurrence of Cluster A (PER)", 
       fill = "Suitability")
