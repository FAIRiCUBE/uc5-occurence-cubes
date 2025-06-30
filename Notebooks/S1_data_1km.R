library(spatialEco)
library(terra)
library(sf)
library(tidyverse)
library(tictoc)
library(whitebox)


# r <- rast('C:/Users/chytryk/documents/GIS_db/DEM/DEM_90m/DEM90_Europe.tif')
# r3035 <- project(r, 'EPSG:3035')
# r3035_agg <- r3035 |> aggregate(16)
# r3035_agg |> writeRaster('C:/Users/chytryk/documents/GIS_db/DEM/DEM_90m/DEM90_Europe_res1002m_3035.tif')


r <- rast('C:/Users/chytryk/documents/GIS_db/DEM/DEM_90m/DEM90_Europe_res1002m_3035.tif')
plot(r)

corine <- rast("D:/GIS_db/Corine2019/clc2018_clc2018_v2018_20_raster100m/CLC2018_CLC2018_V2018_20.tif") 
corine1km <- aggregate(corine, 10, 'modal', na.rm = T)
r_to_cor <- project(r, corine1km)
r_to_cor |> writeRaster('C:/Users/chytryk/documents/GIS_db/DEM/DEM_90m/DEM90_Europe_res1000m_3035.tif')

countries <- read_sf("data/countries_mask.gpkg")
r <- rast('C:/Users/chytryk/documents/GIS_db/DEM/DEM_90m/DEM90_Europe_res1000m_3035.tif')

# rasters -----------------------------------------------------------------
path_name <- 'spatial_predictors_1km/'
attri <- '.tif'

#' whitebox TWI
dem_path <- 'C:/Users/chytryk/documents/GIS_db/DEM/DEM_90m/DEM90_Europe_res1000m_3035.tif'
flow_accum_fd8_path <- "metadata/fd8_flow_accumulation.tif"
slope_path <- "metadata/slope_radians.tif"
wbt_fd8_flow_accumulation(dem = dem_path, output = flow_accum_fd8_path, exponent = 1.0) 
wbt_slope(dem = dem_path, output = slope_path, units = "radians")
flow_accum_rast <- rast(flow_accum_fd8_path)
slope_rast <- rast(slope_path)
slope_rast[slope_rast == 0] <- 0.0001 # Or a more robust handling of flat areas
twi_rast <- log(flow_accum_rast / tan(slope_rast))
twi_rast |>
  setNames("twi") |>
  crop(countries) |>
  mask(countries) |>
  project(r) |>
  writeRaster(paste0(path_name, "twi", attri), overwrite = T)


layers <- list.files(r'(C:\Users\chytryk\Dropbox\PC (6)\Documents\GIS_db\FCrasters)', pattern = '.tif$', full.names = T)

chelsa_narrow <- rast(layers[-1]) |>
  crop(st_transform(countries, 4326)) |>
  project(r) |> 
  mask(countries) |>
  setNames(c("temperature", "precipitation", "tempseasonality", "precseasonality"))

for (i in names(chelsa_narrow)) {
  chelsa_narrow[[i]] |>
    writeRaster(paste0(path_name, i, attri), overwrite = T)
}

rast(layers[1]) |>
  crop(st_transform(countries, 4326)) |>
  project(r) |> 
  mask(st_transform(countries)) |>
  setNames('AI')|>
  writeRaster(paste0(path_name, 'AI', attri), overwrite = T)

hli(r) |>
  setNames("hli") |>
  crop(countries) |>
  mask(countries) |>
  project(r) |> 
  writeRaster(paste0(path_name, "hli", attri), overwrite = T)

tri(r) |>
  setNames("tri") |>
  crop(countries) |>
  mask(countries) |>
    project(r) |> 
  writeRaster(paste0(path_name, "tri", attri), overwrite = T)

terrain(r, v = "slope", unit = "degrees") |>
  setNames("slope") |>
  crop(countries) |>
  mask(countries) |>
    project(r) |> 
  writeRaster(paste0(path_name, "slope", attri), overwrite = T)

# tpi(r, scale = 5, "total") |>
#   setNames("tpi") |>
#   crop(countries) |>
#   mask(countries) |>
#     project(r) |> 
#   writeRaster(paste0(path_name, "tpi", attri), overwrite = T)

r |>
  setNames("elevation") |>
  crop(countries) |>
  mask(countries) |>
  project(r) |> 
  writeRaster(paste0(path_name, "elevation", attri), overwrite = T)

ifel(corine1km == 999, NA, corine1km) |> 
  crop(countries) |>
  mask(countries) |>
  setNames("corine") |>
  project(r) |> 
  writeRaster(paste0(path_name, "corine", attri), overwrite = T)


# -------------------------------------------------------------------------
# species data ------------------------------------------------------------

list.files("species_occurences/", full.names = T) |>
  map(read_csv) |>
  bind_rows() |>
  filter(precision < 500) |>
  st_as_sf(coords = c("x", "y"), crs = 4326) |>
  write_sf("data/species_occurences.gpkg")
