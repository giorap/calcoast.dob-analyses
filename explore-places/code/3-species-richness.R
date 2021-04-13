#' ---
#' title: DOB Cal Coast <br> <strong> Explore Places - Species Richness
#' ---
#' ### Calculate species richness using Chao1
#' #### All species
#' #### MPAs
richness_estimates_byMPA <- split(cal_coast_count_data_byMPA, cal_coast_count_data_byMPA$MPA_ID) %>%
  purrr::map(estimate_richness_from_counts) %>% 
  clean_richness_estimates(ID_fields = cal_coast_count_data_byMPA["MPA_ID"], suffix = "all")
#' #### Watersheds
richness_estimates_byWatershed <- split(cal_coast_count_data_byWatershed, cal_coast_count_data_byWatershed$watershed_ID) %>%
  purrr::map(estimate_richness_from_counts) %>% 
  clean_richness_estimates(ID_fields = cal_coast_count_data_byWatershed["watershed_ID"], suffix = "all")
#' #### Counties
richness_estimates_byCounty <- split(cal_coast_count_data_byCounty, cal_coast_count_data_byCounty$county_ID) %>%
  purrr::map(estimate_richness_from_counts) %>% 
  clean_richness_estimates(ID_fields = cal_coast_count_data_byCounty["county_ID"], suffix = "all")
#' #### Hexagons
richness_estimates_byHexagon <- split(cal_coast_count_data_byHexagon, cal_coast_count_data_byHexagon$hexagon_ID) %>%
  purrr::map(estimate_richness_from_counts) %>% 
  clean_richness_estimates(ID_fields = cal_coast_count_data_byHexagon["hexagon_ID"], suffix = "all")
#'
#' #### Individual taxonomic groups
taxonomic_groups <- c("Mollusca", "Arthropoda", "Echinodermata", "Cnidaria", "Asteroidea", "Nudibranchia", "Chromista")
#' #### MPAs
richness_estimates_byMPA_byTaxon <- purrr::map(taxonomic_groups, function(taxon) {
  split(cal_coast_count_data_byMPA, cal_coast_count_data_byMPA$MPA_ID) %>%
    purrr::map(estimate_richness_from_counts, focal_taxon = taxon) %>% 
    clean_richness_estimates(ID_fields = cal_coast_count_data_byMPA["MPA_ID"], suffix = taxon)
})
#' #### Watersheds
richness_estimates_byWatershed_byTaxon <- purrr::map(taxonomic_groups, function(taxon) {
  split(cal_coast_count_data_byWatershed, cal_coast_count_data_byWatershed$watershed_ID) %>%
    purrr::map(estimate_richness_from_counts, focal_taxon = taxon) %>% 
    clean_richness_estimates(ID_fields = cal_coast_count_data_byWatershed["watershed_ID"], suffix = taxon)
})
#' #### Counties
richness_estimates_byCounty_byTaxon <- purrr::map(taxonomic_groups, function(taxon) {
  richness_estimates_byCounty <- split(cal_coast_count_data_byCounty, cal_coast_count_data_byCounty$county_ID) %>%
    purrr::map(estimate_richness_from_counts, focal_taxon = taxon) %>% 
    clean_richness_estimates(ID_fields = cal_coast_count_data_byCounty["county_ID"], suffix = taxon)
})
#' #### Hexagons
richness_estimates_byHexagon_byTaxon <- purrr::map(taxonomic_groups, function(taxon) {
  split(cal_coast_count_data_byHexagon, cal_coast_count_data_byHexagon$hexagon_ID) %>%
    purrr::map(estimate_richness_from_counts, focal_taxon = taxon) %>% 
    clean_richness_estimates(ID_fields = cal_coast_count_data_byHexagon["hexagon_ID"], suffix = taxon)
})
#' #### Rename all lists
names(richness_estimates_byMPA_byTaxon) <- names(richness_estimates_byWatershed_byTaxon) <- names(richness_estimates_byCounty_byTaxon) <- names(richness_estimates_byHexagon_byTaxon) <- taxonomic_groups