library(spatialEco)
library(terra)
library(sf)
library(tidyverse)
library(tictoc)
library(whitebox)

path_name <- 'spatial_predictors_100m/Europe/'
attri <- '.tif'

countries <- read_sf("data/countries_mask.gpkg")
# r <- rast('C:/Users/chytryk/documents/GIS_db/DEM/DEM_90m/DEM90_Europe.tif')
# corine <- rast("D:/GIS_db/Corine2019/clc2018_clc2018_v2018_20_raster100m/CLC2018_CLC2018_V2018_20.tif")
# tic()
# r_100 <- project(r, corine)
# toc()
# dir.create("spatial_predictors_100m/Europe")
# dir.create('spatial_predictors_100m/Alps')
#r_100 |> writeRaster('spatial_predictors_100m/dem.tif')

r_100 <- rast('spatial_predictors_100m/dem.tif')

tic()
rast("D:/GIS_db/Corine2019/clc2018_clc2018_v2018_20_raster100m/CLC2018_CLC2018_V2018_20.tif") |>
  crop(countries) |>
  mask(countries) |>
  project(r_100) |>
  setNames("corine") |>
  writeRaster(paste0(path_name, "corine", attri), overwrite = T)
toc()

 
#' whitebox TWI
tic()
dem_path <- 'spatial_predictors_100m/dem.tif'
flow_accum_fd8_path <- "metadata/fd8_flow_accumulation_100m.tif"
slope_path <- "metadata/slope_radians_100m.tif"
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
  project(r_100) |>
  writeRaster(paste0(path_name, "twi", attri), overwrite = T)
toc() 


layers <- list.files(r'(C:\Users\chytryk\Dropbox\PC (6)\Documents\GIS_db\FCrasters)', pattern = '.tif$', full.names = T)

chelsa_narrow <- rast(layers[-1]) |>
  crop(st_transform(countries, 4326)) |>
  project(r_100) |> 
  mask(countries) |>
  setNames(c("temperature", "precipitation", "tempseasonality", "precseasonality"))

for (i in names(chelsa_narrow)) {
  chelsa_narrow[[i]] |>
    writeRaster(paste0(path_name, i, attri), overwrite = T)
}

rast(layers[1]) |>
crop(st_transform(countries, 4326)) |>
  project(r_100) |> 
  mask(st_transform(countries)) |>
  setNames('AI')|>
  writeRaster(paste0(path_name, 'AI', attri), overwrite = T)

hli(r_100) |>
  setNames("hli") |>
  crop(countries) |>
  mask(countries) |>
  project(r_100) |> 
  writeRaster(paste0(path_name, "hli", attri), overwrite = T)

tri(r_100) |>
  setNames("tri") |>
  crop(countries) |>
  mask(countries) |>
  project(r_100) |>
  writeRaster(paste0(path_name, "tri", attri), overwrite = T)

terrain(r_100, v = "slope", unit = "degrees") |>
  setNames("slope") |>
  crop(countries) |>
  mask(countries) |>
  project(r_100) |> 
  writeRaster(paste0(path_name, "slope", attri), overwrite = T)

# tpi(r, scale = 5, "total") |>
#   setNames("tpi") |>
#   crop(countries) |>
#   mask(countries) |>
#   project(r) |>
#   writeRaster(paste0(path_name, "tpi", attri), overwrite = T)

r_100 |>
  setNames("elevation") |>
  crop(countries) |>
  mask(countries) |>
  project(r_100) |> 
  writeRaster(paste0(path_name, "elevation", attri), overwrite = T)

#' ========================================================================
#' cropping to Alps
#' ========================================================================
alps <- rast('data/alps.tiff')

for (i in list.files("spatial_predictors_100m/Europe/", full.names = T, pattern = ".tif$")) {
  rast(i) |>
    crop(alps) |>
    mask(alps) |>
    writeRaster(gsub('Europe', 'Alps', i), overwrite = T)
}

