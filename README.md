# UC 5: Validation of Phytosociological Methods through Occurrence Cubes


# About this Use Case
Use case 5 aims to refine and validate habitat prediction methods by incorporating species occurrence data from scientific, museum collections and citizen science sources, with environmental variables, through machine learning models. By integrating plant species occurrences from the Global Biodiversity Information Facility (GBIF) with climatic and topographic data from Earth Observation sources, UC5 seeks to improve habitat mapping accuracy within Europe and highlight the usefulness of data sources as GBIF. Our approach, built over the EUNIS habitat study case 'S22', focuses on comparing predicted species distributions with existing habitat prediction maps, particularly from the European Nature Information System (EUNIS), and exploring the effectiveness of data cubes in improving habitat classification accuracy. The ultimate goal is to improve our understanding of habitat dynamics and promote a more data-driven approach to habitat prediction, with applications in biodiversity monitoring, conservation planning, and museum collections management.


# Research Questions
1. Which environmental factors influence the distribution and community formation of plant species in European habitats (Habitat case study S22)?
2. Can integrating species occurrence data from databases like GBIF, including scientific and citizen science sources, and environmental factors retrieved from EO sources improve habitat prediction accuracy, especially within the context of the EUNIS system?
3. How do the predictions of species distribution compare with the existing EUNIS habitat probability map, and where do discrepancies arise?
   

# Methodology
Use Case 5 integrates species occurrence data, environmental predictors, and machine learning (ML) models to improve habitat prediction accuracy across Europe. The EUNIS habitat type S22 serves as the study case on which the entire modelling framework is built.
To test our approach, we collected occurrence data from the Global Biodiversity Information Facility (GBIF) for the eight diagnostic species (including taxonomic synonyms) associated with habitat S22. These records were combined with a suite of topographic and climatic variables derived from Earth Observation sources—including elevation, slope, TRI, TWI, HLI, Aridity Index (AI), temperature, seasonal temperature, precipitation, and seasonal precipitation—sourced from platforms such as Copernicus, WorldClim, and CHELSA Bioclim.
All datasets were cleaned, harmonised, and integrated into spatially aligned data cubes at two spatial resolutions: 1 km across Europe and 100 m for the Alps region. To enable supervised learning, pseudo-absence data were generated using a 1 km disk buffer method and combined with GBIF presence data.

Species distribution models (SDMs) were developed using an ensemble ML pipeline composed of Generalised Linear Models (GLM), Generalised Additive Models (GAM), and Random Forest (RF). These models produce species occurrence probability maps, which are then combined into a single ensemble prediction. Each model's contribution is weighted by its True Skill Statistic (TSS), calculated through 10-fold cross-validation, stratified by presence–absence data.

Model predictions were validated using the EUNIS-ESy habitat distribution maps derived from the EVA database, which served as an independent reference for habitat S22 at 1 km resolution. 
Lastly, UC5 results were compared to the official EUNIS probability map for the habitat S22 to assess the alignment and potential improvements offered by the species-based modelling approach.
This two-tiered evaluation—against both observational (EVA) and modelled (EUNIS probability) references—provides insights into the accuracy and robustness of the predictions and highlights areas where the ensemble model may refine or complement existing habitat mapping efforts.


# UC5 Outputs
UC5 outputs include:
1. Aggregated species distributions built over individual models (GLM, GAM, RF) and weighted through TSS values, combined in an ensemble model.  
2. A comparative analysis of predicted species distributions against existing EUNIS habitat maps, with recommendations from the UC5 experience for improving habitat classification accuracy.
3. R scripts produced, which can be used for further research and decision-making in biodiversity conservation.



# Getting Started with Our Approach 
R scripts available in the UC5 GitHub directory 'Notebook': 
1. Download and Clean GBIF Occurrence Data (GETDATA_GBIF.R)
2. Integrate GBIF occurrence data with environmental predictors (elevation, slope, TPI, TWI, HLI, AI, temperature, seasonal temperature, precipitation, seasonal precipitation) from EO sources at 1 km (Europe) and 100 m (Alps) resolutions (S1_data_100m.R, S1_data_1km.R)
3. Generate Pseudo-absence data based on GBIF presence only occurrence data, using a disk buffer method (1 km radius) (create_pseudo-absence_based_on_occurrence_presence_data.R)
4. Predict Species Distributions through an Ensemble model built over individual models (GLM, GAM, RF) weighted with TSS scores (Ensemble_model_workflow.R)
5. Compare the UC5 predicted areas with the EUNIS probability map of the Habitat S22 (compare_with_EUNIS.R)


