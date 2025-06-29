#' Script S3: Fit species distribution models
#' 
#' current implementation uses GLM, GAM and RF models and then builds an ensemble model
#' ensemble is made as as a weighted mean of the three models (their predicted probabilities)
#' weighting is done using TSS of individual models
#' TSS values are calculated based on independent 10 fold cross-validation using occurrences as strata

#' outputs of the script are model statistics (csv) and predictions (rasters), both stored in the `outputs` folder

library(irr)
library(terra)
library(furrr)
library(randomForest)
library(foreach)
library(tidyverse)
library(sf)
library(tidymodels)
library(rsample)
library(pROC)
library(mgcv)
library(randomForest)
# library(doParallel)
library(tictoc)
# library(car)


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


#' set of rasters to predict on
#' 100 m resolution
rs <- rast(list.files("spatial_predictors_100m/Alps", full.names = T, pattern = ".tif$"))

#' 1km resolution
#' rs <- rast(list.files("spatial_predictors_1km/", full.names = T, pattern = ".tif$"))

rs_df <- as.data.frame(rs, xy = T) |> drop_na()

#' converting corine to forest/non-forest
rs_df <- rs_df |>
  mutate(corine = factor(corine %in% 311:313, levels = c(TRUE, FALSE), labels = c("forest", "non-forest")))

#' ========================================================================
#' this is parallel processing implementation
#' if desired, uncomment the lines below
#' however, watch out for memory issues
#' ========================================================================
#' numCores <- detectCores()
# cl <- makeCluster(numCores)
# registerDoParallel(cl)

# foreach(i = 1:nrow(df), .packages = c("tidymodels", "tidyverse", "pROC", "mgcv", "randomForest")) %dopar% {
   #data_to_work_with <- df$data[[i]]

goover <- list.files('species_presence_absence_100m/', full.names = T)
#' y <- goover[[1]]

for (y in goover) {
  
  i <- gsub(".csv", "", basename(y))
  print(paste(Sys.time(), i))

  sp <- read_csv(y)

  data_to_work_with <- sp |>
    #' mutate(corine = factor(corine, levels = levels_corine[[1]]))
    mutate(corine = factor(corine %in% 311:313, levels = c(TRUE, FALSE), labels = c("forest", "non-forest")))

  #' mdf <- (vfold_cv(data_to_work_with, strata = layer)$splits[[2]])


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
    get_final_tss() -> final_TSS

  out |>
    summarise(TSS = mean(TSS), threshold = mean(threshold), .by = model) |>
    bind_rows(tibble(model = "ensemble", TSS = final_TSS$TSS, threshold = final_TSS$threshold)) -> stats

  stats |> write_csv(paste0("outputs_100m/", i, "_model_statistics.csv"))

  glm_prediction <- rs_df[c("x", "y")] |> bind_cols(predict(fit_glm(data_to_work_with), rs_df, type = "response"))
  glm_prob <- rast(glm_prediction, crs = "EPSG:3035") |> setNames(paste(i, "GLM"))
  glm_class <- ifel(glm_prob > stats$threshold[stats$model == "glm"], 1, 0)
  glm_class |> writeRaster(paste0("outputs_100m/", i, "_GLM.tif"), overwrite = T)

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
    setNames(i)
  gam_class <- ifel(gam_prob > stats$threshold[stats$model == "gam"], 1, 0)
  gam_class |> writeRaster(paste0("outputs_100m/", i, "_GAM.tif"), overwrite = T)

  rf_prediction <- rs_df[c("x", "y")] |> bind_cols(predict(fit_rf(data_to_work_with), rs_df, type = "prob")[, 2])
  rf_prob <- rast(rf_prediction, crs = "EPSG:3035") |> setNames(paste(i, "RF"))
  rf_class <- ifel(rf_prob > stats$threshold[stats$model == "rf"], 1, 0)
  rf_class |> writeRaster(paste0("outputs_100m/", i, "_RF.tif"), overwrite = T)

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

  plot(ensemble_prediction_class)

  ensemble_prediction_class |> writeRaster(paste0("outputs_100m/", i, "_ensemble.tif"), overwrite = T)

}
      


list.files("outputs_100m/", pattern = "ensemble.tif", full.names = T) |>
  map(rast)


all_spp <- app(rast(list.files('outputs_100m/', pattern = 'ensemble.tif', full.names = T)), 'sum')

all_spp |> writeRaster('outputs_100m/all_spp.tif', overwrite = T)

