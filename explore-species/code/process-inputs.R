#' ---
#' title: DOB Cal Coast <br> <strong> Explore Species - Process Inputs
#' ---
#' 
#' ## Load/process environmental predictor data
#' #### West Coast Regional Ocean Modeling System (ROMS)
#' #### Extract ROMS summary data for the US West Coast
#' *Temperature*
temp_summary_rasters <- purrr::map(list.files("explore-species/data/", pattern = "^temp_.*\\.grd$*", full.names = TRUE), raster)
names(temp_summary_rasters) <- gsub(".grd", "", list.files("explore-species/data/", pattern = "^temp_.*\\.grd$*", full.names = FALSE))
names(temp_summary_rasters) <- gsub("temp_", "", names(temp_summary_rasters))
#' *Salinity*
salt_summary_rasters <- purrr::map(list.files("explore-species/data/", pattern = "^salt_.*\\.grd$*", full.names = TRUE), raster)
names(salt_summary_rasters) <- gsub(".grd", "", list.files("explore-species/data/", pattern = "^salt_.*\\.grd$*", full.names = FALSE))
names(salt_summary_rasters) <- gsub("salt_", "", names(salt_summary_rasters)) 
#' *Bathymetry*
h_summary_rasters <- purrr::map(list.files("explore-species/data/", pattern = "^h_.*\\.grd$*", full.names = TRUE), raster)
names(h_summary_rasters) <- gsub(".grd", "", list.files("explore-species/data/", pattern = "^h_.*\\.grd$*", full.names = FALSE))
names(h_summary_rasters) <- gsub("h_", "", names(h_summary_rasters)) 
#' *Sea surface height*
zeta_summary_rasters <- purrr::map(list.files("explore-species/data/", pattern = "^zeta_.*\\.grd$*", full.names = TRUE), raster)
names(zeta_summary_rasters) <- gsub(".grd", "", list.files("explore-species/data/", pattern = "^zeta_.*\\.grd$*", full.names = FALSE))
names(zeta_summary_rasters) <- gsub("zeta_", "", names(zeta_summary_rasters)) 
#'
#' ### West Coast Environmental Sensitivity Index (ESI) Shore Types
#' #### Extract ESI data for the US West Coast
west_coast_ESI_raster <- readRDS("explore-species/data/west_coast_ESI_raster.rds")
west_coast_ESI_raster_fct <- as.factor(west_coast_ESI_raster)
#'
#' ### Calculate measures of sampling effort based on iNat observations across the study area
west_coast_effort_rasters <- purrr::map(2012:2019, function(yr){
  get_effort_rasters(selected_year = yr)
})
names(west_coast_effort_rasters) <- as.character(2012:2019)
#'
#' ### Data from MARINe surveys
#' #### Data from biodiversity surveys
survey_data_biodiversity <- readRDS("explore-species/data/survey_data_biodiversity.rds")
survey_data_biodiversity_presabs <- survey_data_biodiversity %>% dplyr::mutate(presence_absence = as.numeric(total > 0))
survey_data_biodiversity_presabs <- survey_data_biodiversity_presabs %>% 
  dplyr::group_by(year, latitude, longitude, species_code) %>%
  dplyr::summarise(presence_absence = max(presence_absence, na.rm = TRUE)) 
#' #### Data from long-term counts for mobile species
survey_data_lt_counts <- readRDS("explore-species/data/survey_data_lt_counts.rds")
survey_data_lt_counts_presabs <- survey_data_lt_counts %>% dplyr::mutate(presence_absence = as.numeric(total > 0))
survey_data_lt_counts_presabs <- survey_data_lt_counts_presabs %>% 
  dplyr::group_by(year, latitude, longitude, species_code) %>%
  dplyr::summarise(presence_absence = max(presence_absence, na.rm = TRUE))
#' #### Data from long-term photoplots
survey_data_lt_plots <- readRDS("explore-species/data/survey_data_lt_cov.rds")
survey_data_lt_presabs <- survey_data_lt_plots %>% dplyr::mutate(presence_absence = as.numeric(total > 0))
survey_data_lt_presabs <- survey_data_lt_presabs %>% 
  dplyr::group_by(year, latitude, longitude, species_code) %>%
  dplyr::summarise(presence_absence = max(presence_absence, na.rm = TRUE))
#' #### Combine MARINe data from long-term surveys
survey_data_lt <- rbind(survey_data_lt_counts, survey_data_lt_plots %>% dplyr::filter(species_code %in% setdiff(species_code, survey_data_lt_counts %>% distinct(species_code) %>% pull())))


