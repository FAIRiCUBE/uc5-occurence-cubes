library(pacman)
pacman::p_load(here, tidyverse, rgbif, lubridate, dplyr, CoordinateCleaner, sf)

#link to packages of rGBIF: https://docs.ropensci.org/rgbif/reference/index.html

setwd("P:/FAIRiCUBE/Workpackages_FiC/WP2_use_cases/Occurrence_Cubes/DATA/TAXA_Occurrence_Data/habitat_S26_Rho_hir")

# Setup access to gbif
GBIF_USER = 
GBIF_PWD = 
GBIF_EMAIL = 

#Get taxon key with name_backbone and usageKey. Lookup names in the GBIF backbone taxonomy
taxonKey <- name_backbone("Calamagrostis villosa")$usageKey

basionyms <- name_usage(key=4106690)$data # look at basionym column
synonyms <- name_usage(key=4106690,data="synonyms")$data # this will get all the synonyms 

all_taxon_keys <- c(taxonKey, basionyms$key, synonyms$key)


#Download requests
request_taxonKey <- occ_download(pred_in("taxonKey", c(4106690, 11256675,  8042708)), #all taxa codes
                                 pred("hasCoordinate", TRUE), 
                                 pred("hasGeospatialIssue", FALSE),
                                 pred("occurrenceStatus","PRESENT"), 
                                 pred_not(pred_in("basisOfRecord",c("FOSSIL_SPECIMEN","LIVING_SPECIMEN", "MACHINE_OBSERVATION"))), #there are different options available not included here. e.g., PRESERVED_SPECIMEN (stands for museum collection data)
                                 pred("CONTINENT", "EUROPE"),
                                 pred_gte("year", 1980),
                                 pred_lte("year", 2024),
                                 pred_gte("month", 1),
                                 pred_lte("month", 9),
                                 format = "DWCA",
                                 user = GBIF_USER,
                                 pwd = GBIF_PWD,
                                 email = GBIF_EMAIL)


#download data and clean from centroids
occ_download_wait(request_taxonKey) 
df_inv_species <- occ_download_get(request_taxonKey) %>%
  occ_download_import() 


# Convert the eventDate column to consistent yyyy-mm-dd format
if ("eventDate" %in% colnames(df_inv_species)) {
  # Parse the eventDate column to a consistent format
  df_inv_species$eventDate <- parse_date_time(df_inv_species$eventDate, 
                                              orders = c("Ymd", "Ym", "Y", "Ymd HMS", "Ymd HM", "mdY", "dmY"))
  # Convert to date only (yyyy-mm-dd)
  df_inv_species$eventDate <- as.Date(df_inv_species$eventDate)
}


colnames(df_inv_species)
unique_countryCode <- unique(df_inv_species$countryCode)
print(unique_countryCode)

Cal_vil <- df_inv_species %>%
  filter(coordinateUncertaintyInMeters <= 500) %>%
  select(1, 5, 7, 9, 11, 12, 13, 14, 15, 16, 17, 18, 40, 63, 85, 98, 99, 100, 149, 139) %>%
  filter(!countryCode %in% c("CH", "IT", "AT", "LI", "FR", "DE", "SI")) %>%
  rename(
    y = decimalLatitude, 
    x = decimalLongitude, 
    time = eventDate, 
    spp_name = scientificName, 
    precision = coordinateUncertaintyInMeters
  ) %>%
  mutate(spp_name = "Calamagrostis villosa")

write.csv(Cal_vil, "P:/FAIRiCUBE/Workpackages_FiC/WP2_use_cases/Occurrence_Cubes/DATA/TAXA_Occurrence_Data/habitat_S26_Rho_hir/Cal_vil_raw.csv", row.names = FALSE)

cat("\nMissing Values Per Column:\n")
print(colSums(is.na(Cal_vil)))
# Remove rows with missing latitude or longitude
Cal_vil <- Cal_vil[!is.na(Cal_vil$x) & !is.na(Cal_vil$y), ]

# Check Coordinate Reference System (CRS)
data_sf <- st_as_sf(Cal_vil, coords = c("x", "y"), crs = 4326)
cat("\nCoordinate Reference System:\n")
print(st_crs(data_sf))

Cal_vil_cleaned <- st_coordinates(data_sf) %>%
  as.data.frame() %>%
  cbind(st_drop_geometry(data_sf))

colnames(Cal_vil_cleaned)
Calamagrostis_villosa <- subset(Cal_vil_cleaned, select = c(1, 12, 13, 14, 15, 16, 17, 18, 19))

# Save the formatted dataset to a CSV file
write.csv(Calamagrostis_villosa, "P:/FAIRiCUBE/Workpackages_FiC/WP2_use_cases/Occurrence_Cubes/DATA/TAXA_Occurrence_Data/habitat_S26_Rho_hir/Calamagrostis_villosa.csv", row.names = FALSE)


#GBIF citation
print(gbif_citation(occ_download_meta(request_taxonKey))$download)


