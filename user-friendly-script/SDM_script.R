library(tidyverse)

#' loads functions 
source('scripts/Functions.R')

#' species dasta as obtained from github stored in occurrences folder
i <- read_csv('occurrences/herbs/Artemisia_vulgaris.csv')

#' create presence absence data using disc selection with buffer of 1 km size
human_record <- i |>
  filter(basisOfRecord == "HUMAN_OBSERVATION") |>
  create_presence_absence()

#' fit ensemble model and derive statistics using 10-fold cross validation
#' by default, raster of predicted suitabilities is not generated
model_output <- fit_and_predict_ensemble(human_record)

#' if you want to generate raster of predicted suitabilities, set predict = T
#' only results of the ensemble model are returned
model_output_prediction <- fit_and_predict_ensemble(human_record, predict = T)
model_output_prediction$stats #' model statistics
plot(model_output_prediction$ensemble) #' ensemble model prediction

#' store model statistics
model_output$stats |>
  write_csv('outputs/model_statistics.csv')

#' store raster data
model_output_prediction$ensemble |>
  writeRaster('outputs/ensemble_prediction.tif', overwrite = TRUE)
  