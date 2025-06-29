#' Script S2: Create presence-absence data
#' 
#' this script reads individual species ocurrences as provided by the user stored in `species_ocrrences` folder
#' it also reads the spatial predictors which are rasters stored in `spatial_predictors` folder
#' next, in a species-specific loop, it creates a buffer around the species ocurrences 
#' to exclude this area from pseudo-absence selections
#' lastly it stores presence-and-pseudo-absence data in `species_presence_absence` folder for each species

#' user can choose the buffer size (in metres)
#' buffer_size <- 1e3

#' on potentially can use different sets of species or predictors if these are exchanged in given folders

#' script uses EPSG:3035 projection for all operations

#' ========================================================================
#' script starts here
#' ========================================================================

library(sf)
library(terra)
library(tidyverse)

if(!exists('buffer_size')){buffer_size <- 100}


#' ========================================================================
#' species data
#' ========================================================================
read_csv("data/S22df_datasource_grouped.csv") |>
  filter(precision < 500) |>
  st_as_sf(coords = c("x", "y"), crs = 4326) |>
  st_transform(3035) |>
  select(spp_name, data_source_group) |> 
  write_sf("data/species_occurences.gpkg")

#' ========================================================================
#' predictors
#' ========================================================================

#' 100 m use
rs <- rast(list.files('spatial_predictors_100m/Europe', full.names = T, pattern = '.tif$'))
#' 1 km use
#' rs <- rast(list.files('spatial_predictors_1km', full.names = T, pattern = '.tif$'))

species <- unique(read_sf('data/species_occurences.gpkg')$spp_name)
#' i <- 'Dryas octopetala'

for (i in species) {
print(i)  
  species_iter <- read_sf('data/species_occurences.gpkg') |> 
    filter(spp_name == i) |> 
    mutate(layer = 1) |> 
    select(layer) |>
    mutate(Y = round(st_coordinates(geom)[,2]/1000),
    X = round(st_coordinates(geom)[,1]/1000)) |> 
    group_by(X, Y) |> 
    slice(1) 
  
   valid_area <- as.polygons(!is.na(rs[['elevation']] |> 
                                     setNames('validarea')), dissolve = T) |> 
    st_as_sf() |> 
    filter(validarea == 1)
  
  spp_buf <- st_union(st_buffer(species_iter, buffer_size)) #' user-defined buffer size
  valid_raster <- rasterize(st_difference(valid_area, spp_buf), rs[[1]])
  spp_buf |> write_sf(paste0('test/', i, '_spp_buf.gpkg'), overwrite = T)
  valid_raster |> writeRaster(paste0('test/',i, '_valid_raster.tif'), overwrite = T)

  random_samples <- spatSample(valid_raster, size = nrow(species_iter)*10, as.points = T) |> 
    st_as_sf() |> 
    filter(layer == 1) |> 
    sample_n(size = nrow(species_iter)) |> 
    rename(geom = geometry) |> 
    st_transform(3035) |> 
    mutate(layer = 0)
  
  pas <- bind_rows(species_iter, random_samples) |>
    mutate(layer = replace_na(layer, 0)) |> 
    select(layer)
  
  bind_cols(pas, terra::extract(rs, pas)) |>
    drop_na() |>
    as_tibble() |>
    select(-geom, -ID) |>
    write_csv(paste0("species_presence_absence_100m/", i, ".csv"))
  
  print('over')
}
