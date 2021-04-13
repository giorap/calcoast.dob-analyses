#' ---
#' title: DOB Cal Coast <br> <strong> Load and process spatial boundaries
#' ---
#' #### This code loads and processes spatial layers corresponding to the various spatial boundaries of interest
#' 
#' ## West coast
#' #### Load west coast polygon
west_coast_polygon <- readRDS("get-data/data/west_coast_polygon.rds")
#' #### Additionally, generate a ~10km raster from the polygon
west_coast_grid <- west_coast_polygon %>%
  spsample(type = "regular", cellsize = 0.0833) %>%
  points2grid() %>% 
  raster()
west_coast_grid[] <- 1
west_coast_grid <- west_coast_grid %>% mask(west_coast_polygon) %>% trim()
#' #### and a data.frame of grid cell centroids from the grid
west_coast_points <- rasterToPoints(west_coast_grid, spatial = FALSE, )[, 1:2] %>% as.data.frame()
west_coast_points$cellID <- cellFromXY(west_coast_grid, west_coast_points)
#' 
#' ## California State boundary
#' #### Load California polygon
california_polygon <- readRDS("get-data/data/california_polygon.rds")
#'
#' ## California's Marine Protected Area Network
#' #### Load MPA polygons
MPA_polygons <- readRDS("get-data/data/MPA_polygons.rds")
#'
#' ## California's North, Central, and South coastal regions
#' #### Load coastal region polygons
region_polygons <- readRDS("get-data/data/region_polygons.rds")
#'
#' ## California's coastal counties 
county_polygons <- readRDS("get-data/data/county_polygons.rds")
#'
#' ## California's coastal watersheds
watershed_polygons <- readRDS("get-data/data/watershed_polygons.rds")
#'
#' ## California's coastal hexagonal grid 
hexagon_polygons <- readRDS("get-data/data/hexagon_polygons.rds")