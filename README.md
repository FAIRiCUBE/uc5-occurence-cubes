# uc5-occurence-cubes

UC 5: Validation of Phytosociological Methods through Occurrence Cubes

1. Objectives
Phytosociology is a branch of vegetation science that studies vegetation communities and classifies the species occurring in such communities in hierarchical vegetation units. A classic phytosociological method is based on the cover that the species have on ground and translates this pattern into units of species combinations forming the community.
However, this approach does not explain the reason why they occur on the spot and form such communities since the environmental conditions are only partially or not at all considered.
In this context, FAIRiCUBE provides the ideal framework for ground-truthing, as we can compare known environmental conditions with the areas where these plant communities occur.
In particular, the main objective of this use case is to validate the traditional methods applied in phytosociology to characterize and classify plant communities. This will be approached by linking distribution data of plant species based on records from human observations and collection samples from the Global Biodiversity Information Facility (GBIF, www.gbif.org) and an online collaboration platform for botanical collections (JACQ).
Moreover, this use case aims also to develop a new phytosociological approach to characterize and predict the presence of plant communities in unknown localities, based on satellite and occurrence data of corresponding known communities.

2.	Research questions
For this use case, the research questions addressed are as follows:
-	Which factors abiotic or biotic contribute to the distribution of taxa to form vegetation communities?
-	Do occurrences of taxa vary along environmental gradients and what are the driving forces behind the observed distribution patterns.
-	Is it possible to predict the presence of vegetation communities in unknown locations based on known occurrences of corresponding communities and contributing environmental factors?

3. Methods
First and foremost, a list of habitats will be chosen from the EUNIS classification (European Nature Information System) of Habitat types. The diagnostic taxa related to the habitat types will be obtained together with their occurrences from the Global Biodiversity Information Facility (GBIF) and the Virtual Herbaria System (JACQ). 
The rest of the taxa comprised in the vegetation units of the Habitats chosen will be obtained from vegetation units present in Mucina et al. (1993)¹. 
Furthermore, known distributions of plant communities will be obtained as raster data based on the vegetation units of Europe to produce Community Cubes. 
A second set of Data Cubes, based on occurrences of taxa, will be produced by combining biotic and abiotic data from Rasdaman services together with taxon occurrence data using the tool Wormpicker developed by UC3. The tool will retrieve EO point estimates from the Rasdaman interface based on point coordinates derived from taxa occurrences. 
Once we obtain Occurrence Cubes and Community Cubes, we will investigate the distribution patterns of the taxa where plant communities have been identified.
Lastly, we will use ML and AI approaches to identify relations between identified communities and EO data, determine locations with favourable environmental conditions and predict possible presences at these locations of plant communities corresponding to the ones investigated.

4. Expected results
The expected results of this use case include the creation of Data Cubes on the basis of habitat and vegetation types, taxon occurrences and earth observation data, the investigation of the distribution patterns of vegetation communities based on environmental factors and predictions of vegetation communities’ presence in unknown localities. 
The results of this use case will be in line with F.A.I.R. principles and available for use by a broader scientific audience. In particular, the data cubes produced will be disseminated and identified with Digital Object Identifiers (DOIs). Moreover, API specifications will be provided for data delivery from JACQ and GBIF together with Notebooks for collection data integrated into data cube processes and analysis. Taxon distribution patterns in relation to biotic and abiotic factors will be published. Lastly, guidance for legacy data preparation will be provided. 


