#' ---
#' title: DOB Cal Coast <br> <strong> Explore Places - Biodiversity Importance
#' ---
#' ### Load priority ranks from a run of the Core Area Zonation algorithm at the network scale
#' #### MPAs
biodiversity_importance_byMPA <- readRDS("explore-places/data/MPA_priorities.rds")
#' #### Watersheds
biodiversity_importance_byWatershed <- readRDS("explore-places/data/watershed_priorities.rds")
#' #### Hexagons
biodiversity_importance_byHexagon <- readRDS("explore-places/data/hexagon_priorities.rds")
#'
#' ### Calculate proportion area protected
#' #### Loop through all MPAs
proportion_protected_byMPA <- data.frame(ID = MPA_polygons@data$ID, matrix(NA, ncol = 8, nrow = nrow(MPA_polygons), dimnames = list(1:nrow(MPA_polygons), paste0("percentage_", c(as.character(unique(MPA_polygons$Type)), "protected")))))
proportion_protected_byMPA$percentage_protected <- 1
for (i in 1:nrow(MPA_polygons)){
  proportion_protected_byMPA[i, grep(paste0("^percentage_", as.character(unique(MPA_polygons$Type[i])), "$"), names(proportion_protected_byMPA))] <- 1
}
#' #### Loop through all Watersheds
proportion_protected_byWatershed <- data.frame(ID = watershed_polygons@data$ID, matrix(NA, ncol = 8, nrow = nrow(watershed_polygons), dimnames = list(1:nrow(watershed_polygons), paste0("percentage_", c(as.character(unique(MPA_polygons$Type)), "protected")))))
for (i in 1:nrow(watershed_polygons)){
  proportion_protected_byWatershed[i, -1] <- calculate_proportion_protected(watershed_polygons[i, ])
}
#' #### Loop through all Counties
proportion_protected_byCounty <- data.frame(ID = county_polygons@data$ID, matrix(NA, ncol = 8, nrow = nrow(county_polygons), dimnames = list(1:nrow(county_polygons), paste0("percentage_", c(as.character(unique(MPA_polygons$Type)), "protected")))))
for (i in 1:nrow(county_polygons)){
  proportion_protected_byCounty[i, -1] <- calculate_proportion_protected(county_polygons[i, ])
}
#' #### Loop through all Hexagons
proportion_protected_byHexagon <- data.frame(ID = hexagon_polygons@data$ID, matrix(NA, ncol = 8, nrow = nrow(hexagon_polygons), dimnames = list(1:nrow(hexagon_polygons), paste0("percentage_", c(as.character(unique(MPA_polygons$Type)), "protected")))))
for (i in 1:nrow(hexagon_polygons)){
  proportion_protected_byHexagon[i, -1] <- calculate_proportion_protected(hexagon_polygons[i, ])
}
