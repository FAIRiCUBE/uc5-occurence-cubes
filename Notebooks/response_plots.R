library(terra)
library(sf)
library(tidyverse)

rs <- rast(list.files("spatial_predictors_1km/", full.names = T, pattern = ".tif$"))
rs_df <- as.data.frame(rs, xy = T) |> as_tibble() |> drop_na()

species <- rast(list.files("outputs_1km", pattern = "ensemble", full.names = T))
names(species) <- gsub(".csv_ensemble.tif", "", basename(sources(species)))
species_df <- as.data.frame(species, xy = T) |> as_tibble()

for(i in names(species_df)[-c(1, 2)]) {
    
    tibble(response = species_df[[i]]) |>
        bind_cols(rs_df[-c(1, 2)]) |>
        group_by(response) |>
        sample_n(2000) |>
        select(-corine) |>
        pivot_longer(-response) |>
        ggplot(aes(value, response)) +
        facet_wrap(~name, scales = "free") +
        geom_smooth(colour = 'red', alpha = .1, fill = 'red') + 
        #geom_smooth(
        #    method = "glm", method.args = list(family = "binomial"),
        #    formula = y ~ poly(x, 2), # y and x correspond to your mapped variables
        #    se = FALSE, color = "blue"
        #) +
        theme_bw() +
        labs(
            x = "Value of the predictor", y = "Probability of occurrence",
            title = i
        ) +
        theme(
            panel.grid = element_blank(),
            strip.text = element_text(size = 12, face = "bold", hjust = 0),
            strip.background = element_blank(),
            plot.title = element_text(face = 'italic', size = 16)
        )
        ggsave(paste0("response_plots_1km/", i, ".png"), width = 10, height = 8)

}
