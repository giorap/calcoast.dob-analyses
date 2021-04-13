#' ---
#' title: DOB Cal Coast <br> <strong> Explore Places - Species rarity
#' ---
#' ### Calculate rarity-weighted richness (weighted endemism) of each place
#' #### All species
#' #### MPAs
rarity_byMPA <- prepare_weighted_endemism_inputs(count_data = cal_coast_count_data_byMPA, polys = MPA_polygons) %>% 
  weighted.endemism(records = "site", deg.resolution = 0.01)
#' #### Watersheds
rarity_byWatershed <- prepare_weighted_endemism_inputs(count_data = cal_coast_count_data_byWatershed, polys = watershed_polygons) %>% 
  weighted.endemism(records = "site", deg.resolution = 0.25)
#' #### Counties
rarity_byCounty <- prepare_weighted_endemism_inputs(count_data = cal_coast_count_data_byCounty, polys = county_polygons) %>% 
  weighted.endemism(records = "site", deg.resolution = 0.01)
#' #### Hexagons
rarity_byHexagon <- prepare_weighted_endemism_inputs(count_data = cal_coast_count_data_byHexagon, polys = hexagon_polygons) %>% 
  weighted.endemism(records = "site", deg.resolution = 0.01)
