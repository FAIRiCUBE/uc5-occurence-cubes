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


create_presence_absence <- function(occurrence_data = one_spp, buffer_size = 100, stratification = 1000) {

  require(terra)
  require(tidyverse)
  require(sf)


  if(!('data.frame' %in% class(occurrence_data))){
    occurrence_data <- read_csv(occurrence_data)
  }

  #' ========================================================================
  #' species data
  #' ========================================================================
  occurrence_data |>
    filter(precision < 500) |>
    st_as_sf(coords = c("x", "y"), crs = 4326) |>
    st_transform(3035) |>
    mutate(layer = 1) |>
    mutate(
      Y = round(st_coordinates(geometry)[, 2] / stratification),
      X = round(st_coordinates(geometry)[, 1] / stratification)
    ) |>
    group_by(X, Y) |>
    slice(1) |>
    select(layer) -> species

  #' ========================================================================
  #' predictors
  #' ========================================================================

  rs <- rast(list.files("spatial_predictors_1km/", full.names = T, pattern = ".tif$"))

  valid_area <- as.polygons(!is.na(rs[["elevation"]] |>
    setNames("validarea")), dissolve = T) |>
    st_as_sf() |>
    filter(validarea == 1)

  spp_buf <- st_union(st_buffer(species, buffer_size)) #' user-defined buffer size
  valid_raster <- rasterize(st_difference(valid_area, spp_buf), rs[[1]])

  # spp_buf |> write_sf(paste0('test/', i, '_spp_buf.gpkg'), overwrite = T)
  # valid_raster |> writeRaster(paste0('test/',i, '_valid_raster.tif'), overwrite = T)

  random_samples <- spatSample(valid_raster, size = nrow(species) * 10, as.points = T) |>
    st_as_sf() |>
    filter(layer == 1) |>
    sample_n(size = nrow(species)) |>
    st_transform(3035) |>
    mutate(layer = 0)

  pas <- bind_rows(species, random_samples) |>
    mutate(layer = replace_na(layer, 0)) |>
    select(layer) #|>
 #   mutate(spp_name = gsub("_", " ", str_remove(basename(occurrence_data), ".csv"))) |>
 #   relocate(spp_name)

  bind_cols(pas, terra::extract(rs, pas)) |>
    drop_na() |>
    as_tibble() |>
    select(-geometry, -ID) -> out
  
  return(out)
}

    

#' ========================================================================


#' Script S3: Fit species distribution models
#' 
#' current implementation uses GLM, GAM and RF models and then builds an ensemble model
#' ensemble is made as as a weighted mean of the three models (their predicted probabilities)
#' weighting is done using TSS of individual models
#' TSS values are calculated based on independent 10 fold cross-validation using occurrences as strata

library(irr)
library(terra)
library(furrr)
library(foreach)
library(tidyverse)
library(sf)
library(tidymodels)
library(rsample)
library(pROC)
library(mgcv)
library(randomForest)
library(tictoc)

fit_gam <- function(df) {
  gam_model <- gam(
    layer ~ s(AI) +
      s(elevation) +
      s(twi) +
      s(hli) +
      s(precipitation) +
      s(precseasonality) +
      s(slope) +
      s(temperature) +
      s(tempseasonality) +
      s(tri) +
      corine,
    data = df, family = binomial
  )
  return(gam_model)
}

fit_glm <- function(df) {
  glm(
    layer ~ poly(AI, 2) +
      poly(elevation, 2) +
      poly(twi, 2) +
      poly(hli, 2) +
      poly(precipitation, 2) +
      poly(precseasonality, 2) +
      poly(slope, 2) +
      poly(temperature, 2) +
      poly(tempseasonality, 2) +
      poly(tri, 2) +
      corine,
    data = df, 
    family = binomial
  )
}

fit_rf <- function(df) {
  rf_model <- randomForest(layer ~ .,
    data = df |> mutate(layer = factor(layer)), family = binomial()
  )
  return(rf_model)
}

get_final_tss <- function(df) {
  auc <- as.numeric(roc(df$observed, df$predicted, quiet = T)$auc)
  vals <- coords(roc(df$observed, df$predicted, quiet = T), "best")
  predicted_binary <- ifelse(df$predicted > vals$threshold[[1]], 1, 0)
  kappa_value <- kappa2(tibble(predicted = predicted_binary, observed = df$observed))$value
  accuracy <- mean(predicted_binary == df$observed)
  tss <- vals$sensitivity[1] + vals$specificity[1] - 1

  out <- tibble(
    threshold = vals$threshold[[1]],
    TSS = tss, 
    Kappa = kappa_value, 
    sensitivity = vals$sensitivity,
    specificity = vals$specificity,
    auc = auc,
    accuracy = accuracy
  )

  return(out)
}

get_tss <- function(m, df) {
  if (class(m)[2] == "randomForest") {
    prediction <- data.frame(predict(m, df, type = "prob"))$X1
  } else {
    prediction <- predict(m, df, type = "response")
  }

  auc <- as.numeric(roc(df$layer, prediction, quiet = T)$auc)
  vals <- coords(roc(df$layer, prediction, quiet = T), "best")[1, ]
  predicted_binary <- ifelse(prediction > vals$threshold[[1]], 1, 0)
  kappa_value <- kappa2(tibble(predicted = predicted_binary, observed = df$layer))$layer
  accuracy <- mean(predicted_binary == df$layer)
  tss <- vals$sensitivity[1] + vals$specificity[1] - 1


  stats <- tibble(
    predicted = prediction,
    observed = df$layer,
    TSS = tss,
    Kappa = kappa_value,
    sensitivity = vals$sensitivity,
    specificity = vals$specificity,
    threshold = vals$threshold,
    accuracy = accuracy,
    auc = auc
  )

  return(stats)
}


fit_and_predict_ensemble <- function(presence_absence, predict = F){

#' set of rasters to predict on
#' converting corine to forest/non-forest
rs <- rast(list.files("spatial_predictors_1km/", full.names = T, pattern = ".tif$"))
rs_df <- as.data.frame(rs, xy = T) |> drop_na() |> 
  mutate(corine = factor(corine %in% 311:313, levels = c(TRUE, FALSE), labels = c("forest", "non-forest")))

data_to_work_with <- presence_absence |>
  #select(-spp_name) |> 
    mutate(corine = factor(corine %in% 311:313, levels = c(TRUE, FALSE), labels = c("forest", "non-forest")))

  vfold_cv(data_to_work_with, strata = layer, v = 2) |>
    mutate(output = map(splits, \(mdf){
      list(
        get_tss(fit_glm(analysis(mdf)), assessment(mdf)) |>
          mutate(model = "glm"),
        get_tss(fit_gam(analysis(mdf)), assessment(mdf)) |>
          mutate(model = "gam"),
        get_tss(fit_rf(analysis(mdf)), assessment(mdf)) |>
          mutate(model = "rf")
      ) |>
        bind_rows()
    })) |>
    select(output) |>
    unnest() -> out

  out |>
    group_by(model) |>
    mutate(id = row_number()) |>
    ungroup() |>
    group_by(id) |>
    summarise(
      predicted = weighted.mean(predicted, w = TSS),
      observed = mean(observed)
    ) |>
    get_final_tss() |>
    select(-Kappa) |> 
    mutate(model = 'ensemble') -> final_TSS

  out |>
    summarise(across(TSS:auc, mean), .by = model) |>
    bind_rows(final_TSS) -> stats 

  if(predict == T){

  glm_prediction <- rs_df[c("x", "y")] |> bind_cols(predict(fit_glm(data_to_work_with), rs_df, type = "response"))
  glm_prob <- rast(glm_prediction, crs = "EPSG:3035") |> setNames(paste("GLM"))
  glm_class <- ifel(glm_prob > stats$threshold[stats$model == "glm"], 1, 0)
  #glm_class |> writeRaster(paste0("outputs_100m/", i, "_GLM.tif"), overwrite = T)

  m_gam <- fit_gam(data_to_work_with)
  prediction <- rs_df |>
    select(-x, -y) |>
    group_by(chunk = ntile(row_number(), 32)) |> #
    group_split() |>
    future_map(\(x) predict(m_gam, x, type = "response")) |>
    unlist() |>
    as_tibble()
  gam_prob <- bind_cols(rs_df[c("x", "y")], prediction) |>
    rast(crs = "EPSG:3035") |>
    setNames('GAM')
  gam_class <- ifel(gam_prob > stats$threshold[stats$model == "gam"], 1, 0)
  #gam_class |> writeRaster(paste0("outputs_100m/", i, "_GAM.tif"), overwrite = T)

  rf_prediction <- rs_df[c("x", "y")] |> bind_cols(predict(fit_rf(data_to_work_with), rs_df, type = "prob")[, 2])
  rf_prob <- rast(rf_prediction, crs = "EPSG:3035") |> setNames("RF")
  rf_class <- ifel(rf_prob > stats$threshold[stats$model == "rf"], 1, 0)
  #rf_class |> writeRaster(paste0("outputs_100m/", i, "_RF.tif"), overwrite = T)

  ensemble_prediction <- terra::weighted.mean(
    x = c(
      glm_prob,
      gam_prob,
      rf_prob
    ),
    w = c(
      stats$threshold[stats$model == "glm"],
      stats$threshold[stats$model == "gam"],
      stats$threshold[stats$model == "rf"]
    )
  )

  ensemble_prediction_class <- ifel(ensemble_prediction > stats$threshold[stats$model == "ensemble"], 1, 0)

  return_ls <- list(ensemble = setNames(ensemble_prediction_class, 'prediction'), 
  stats = stats)
  } else {
    return_ls <- list(stats = stats)
  }

}
