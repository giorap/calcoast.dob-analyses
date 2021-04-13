#' ---
#' title: DOB Cal Coast <br> <strong> Explore Places - Process Outputs
#' ---
#' ### Generate outputs objects to be included within the Places Shiny app
#' 
#' ### Species richness
#' #### Store output within spatial boudary polygons
#' #### MPAs
MPA_polygons@data <- MPA_polygons@data %>% left_join(richness_estimates_byMPA, by = c("ID" = names(richness_estimates_byMPA)[grep("ID", names(richness_estimates_byMPA))])) %>%
  left_join(join_all(richness_estimates_byMPA_byTaxon), by = c("ID" = "MPA_ID"))
#' #### Watersheds
watershed_polygons@data <- watershed_polygons@data %>% left_join(richness_estimates_byWatershed, by = c("ID" = names(richness_estimates_byWatershed)[grep("ID", names(richness_estimates_byWatershed))])) %>%
  left_join(join_all(richness_estimates_byWatershed_byTaxon), by = c("ID" = "watershed_ID"))
#' #### Counties
county_polygons@data <- county_polygons@data %>% left_join(richness_estimates_byCounty, by = c("ID" = names(richness_estimates_byCounty)[grep("ID", names(richness_estimates_byCounty))])) %>%
  left_join(join_all(richness_estimates_byCounty_byTaxon), by = c("ID" = "county_ID"))
#' #### Hexagons
hexagon_polygons@data <- hexagon_polygons@data %>% left_join(richness_estimates_byHexagon, by = c("ID" = names(richness_estimates_byHexagon)[grep("ID", names(richness_estimates_byHexagon))])) %>%
  left_join(join_all(richness_estimates_byHexagon_byTaxon), by = c("ID" = "hexagon_ID"))
#'
#' ### Biodiversity uniqueness
#' #### Store output within spatial polygons 
#' #### MPAs
MPA_polygons@data <- MPA_polygons@data %>% left_join(uniqueness_estimates_byMPA, by = c("ID" = names(uniqueness_estimates_byMPA)[grep("ID", names(uniqueness_estimates_byMPA))])) %>%
  left_join(join_all(uniqueness_estimates_byMPA_byTaxon), by = c("ID" = "MPA_ID"))
#' #### Watersheds
watershed_polygons@data <- watershed_polygons@data %>% left_join(uniqueness_estimates_byWatershed, by = c("ID" = names(uniqueness_estimates_byWatershed)[grep("ID", names(uniqueness_estimates_byWatershed))])) %>%
  left_join(join_all(uniqueness_estimates_byWatershed_byTaxon), by = c("ID" = "watershed_ID"))
#' #### Counties
county_polygons@data <- county_polygons@data %>% left_join(uniqueness_estimates_byCounty, by = c("ID" = names(uniqueness_estimates_byCounty)[grep("ID", names(uniqueness_estimates_byCounty))])) %>%
  left_join(join_all(uniqueness_estimates_byCounty_byTaxon), by = c("ID" = "county_ID"))
#' #### Hexagons
hexagon_polygons@data <- hexagon_polygons@data %>% left_join(uniqueness_estimates_byHexagon, by = c("ID" = names(uniqueness_estimates_byHexagon)[grep("ID", names(uniqueness_estimates_byHexagon))])) %>%
  left_join(join_all(uniqueness_estimates_byHexagon_byTaxon), by = c("ID" = "hexagon_ID"))
#'
#' ### Species rarity
#' #### Store output within spatial polygons 
#' #### MPAs
MPA_polygons@data <- MPA_polygons@data %>% left_join(rarity_byMPA$WE, by = c("ID" = names(rarity_byMPA$WE)[grep("ID", names(rarity_byMPA$WE))]))
#' #### Watersheds
watershed_polygons@data <- watershed_polygons@data %>% left_join(rarity_byWatershed$WE, by = c("ID" = names(rarity_byWatershed$WE)[grep("ID", names(rarity_byWatershed$WE))]))
#' #### Counties
county_polygons@data <- county_polygons@data %>% left_join(rarity_byCounty$WE, by = c("ID" = names(rarity_byCounty$WE)[grep("ID", names(rarity_byCounty$WE))]))
#' #### Hexagons
hexagon_polygons@data <- hexagon_polygons@data %>% left_join(rarity_byHexagon$WE, by = c("ID" = names(rarity_byHexagon$WE)[grep("ID", names(rarity_byHexagon$WE))]))
#'
#' ### Biodiversity importance
#' #### Store output within spatial polygons 
#' #### MPAs
MPA_polygons@data <- MPA_polygons@data %>% left_join(biodiversity_importance_byMPA %>% dplyr::select(MPA_ID, rank_priority), by = c("ID" = "MPA_ID"))
#' #### Watersheds
watershed_polygons@data <- watershed_polygons@data %>% left_join(biodiversity_importance_byWatershed %>% dplyr::select(watershed_ID, rank_priority), by = c("ID" = "watershed_ID"))
#' #### Counties
county_polygons@data <- county_polygons@data %>% dplyr::mutate(rank_priority = NA)
#' #### Hexagons
hexagon_polygons@data <- hexagon_polygons@data %>% left_join(biodiversity_importance_byHexagon %>% dplyr::select(hexagon_ID, rank_priority), by = c("ID" = "hexagon_ID"))
#'
#' ### Percentage area protected
#' #### MPAs
MPA_polygons@data <- MPA_polygons@data %>% left_join(proportion_protected_byMPA, by = "ID")
#' #### Watersheds
watershed_polygons@data <- watershed_polygons@data %>% left_join(proportion_protected_byWatershed, by = "ID")
#' #### Counties
county_polygons@data <- county_polygons@data %>% left_join(proportion_protected_byCounty, by = "ID")
#' #### Hexagons
hexagon_polygons@data <- hexagon_polygons@data %>% left_join(proportion_protected_byHexagon, by = "ID")
#'
#' ### Output spatial polygon objects
saveRDS(MPA_polygons, "explore-places/output/MPA_polygons.rds")
saveRDS(watershed_polygons, "explore-places/output/watershed_polygons.rds")
saveRDS(county_polygons, "explore-places/output/county_polygons.rds")
saveRDS(hexagon_polygons, "explore-places/output/hexagon_polygons.rds")
#'
#' ### Temporal stability
#' #### MPAs
saveRDS(temporal_stability_byMPA, "explore-places/outputs/temporal_stability_byMPA.rds")
#' #### Watersheds
saveRDS(temporal_stability_byWatershed, "explore-places/outputs/temporal_stability_byWatershed.rds")
#' #### Counties
saveRDS(temporal_stability_byCounty, "explore-places/outputs/temporal_stability_byCounty.rds")
#' #### Hexagons
saveRDS(temporal_stability_byHexagon, "explore-places/outputs/temporal_stability_byHexagon.rds")

