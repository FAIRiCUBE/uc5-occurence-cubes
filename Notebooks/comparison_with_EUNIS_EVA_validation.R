library(terra)
library(sf)
library(tidyverse)


all_sp <- rast("outputs_100m/all_spp.tif")
eunis_narrow <- rast("data/EUNIS_Prob_S22_100m.tif") |>
    crop(all_sp) |>
    mask(all_sp)
eunis_narrow_cat <- eunis_narrow > 0
species_cat <-  (all_sp > 6)

code_both <- 1
code_only_species <- 2
code_only_eunis <- 3
code_nothing <- 4

combined_numeric <- ifel(eunis_narrow_cat & species_cat, code_both,
                         ifel(species_cat, code_only_species,
                              ifel(eunis_narrow_cat, code_only_eunis, code_nothing))
                       )

category_levels <- data.frame(
  ID = 1:4,
  category = c("both", "only species", "only eunis", "nothing")
)

levels(combined_numeric) <- category_levels

combined_numeric |>
    writeRaster('outputs_100m/combined_numeric.tif',
                overwrite = TRUE)

#' ========================================================================
#' 1 km data

all_sp_1km <- rast("outputs_1km/all_spp.tif")
eunis_1km <- rast("data/EUNIS_Prob_S22_100m.tif") |>
    aggregate(10) |>
    crop(all_sp_1km) |>
    mask(all_sp_1km)
    

eunis_cat_1km <- eunis_1km > 0
species_cat_1km <- (all_sp_1km > 3)

code_both <- 1
code_only_species <- 2
code_only_eunis <- 3
code_nothing <- 4

combined_numeric <- ifel(eunis_cat_1km & species_cat_1km, code_both,
                         ifel(species_cat_1km, code_only_species,
                              ifel(eunis_cat_1km, code_only_eunis, code_nothing))
                       )

category_levels <- data.frame(
  ID = 1:4,
  category = c("both", "only species", "only eunis", "nothing")
)

levels(combined_numeric) <- category_levels

combined_numeric |>
    writeRaster('outputs_1km/combined_numeric.tif',
                overwrite = TRUE)




#' ========================================================================

plot(combined_numeric)

eunis_narrow <- crop(eunis, all_sp) |>
    mask(all_sp)

both <- c(eunis_narrow |> setNames('eunis'), all_sp |> setNames('species')) |>
    as.data.frame(xy = T)


both |>
    sample_n(100000) |>
    ggplot(aes(factor(species), eunis)) +
    geom_boxplot(fill = "red", alpha = .2) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(
        y = "EUNIS probability",
        x = "Number of species for which suitable"
    )

ggsave(paste0("EUNIS-Species_comparison_Alps_100m", ".png"), width = 10, height = 8)

#' ========================================================================
#' 

countries <- read_sf("data/countries_mask.gpkg")

clipped <- read_delim("data/eva-export_2025-5-29.tsv") |>
    select(Latitude, Longitude) |>
    st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) |>
    st_transform(3035) |>
    mutate(
        X = (round(st_coordinates(geometry)[, 1]/1000)),
        Y = (round(st_coordinates(geometry)[, 2]/1000))
    ) |>     
    group_by(X, Y) |>
    slice(1) |>
    st_intersection(countries)


clipped |>
    mutate(terra::extract(combined_numeric, clipped)) |>
    drop_na() |>
    arrange(category) |>
    ggplot() +
    geom_sf(data = countries, fill = "white") +
    geom_sf(size = 1, aes(colour = category), fill = "grey80") +
    coord_sf(xlim = c(2964041, 5620692), ylim = c(1660330, 4269174)) +
    theme_bw() +
    theme(legend.position = c(0, 1), legend.justification = c(0, 1), legend.background = element_blank()) -> m

ggsave('outputs_1km/clipped_categories.png', m, width = 10, height = 8)

clipped |>
    mutate(terra::extract(combined_numeric, clipped)) |>
    drop_na() |>
    arrange(category) |>
    count(category)

bind_cols(terra::extract(combined_numeric, clipped)) |>
    drop_na() |>
    ggplot() +
    geom_sf()

    group_by(category) |> 
    count()
