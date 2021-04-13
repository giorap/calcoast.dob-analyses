#' ---
#' title: DOB Cal Coast <br> <strong> Master 
#' ---
#' 
#' #### Master code file to run all individual code files to build the DOB Cal Coast
#' 
#' ## Get Data
#' ### Load and process spatial boundary layers
#' #### Duration (s): 0.498
system.time(source("get-data/code/process-spatial-boundaries.R"))
#' ### Load and process iNaturalist observations
#' #### Duration (s): 1917.348
system.time(source("get-data/code/process-iNat-observations.R"))
#' ### Cache .RData
cache_RData("get-data", "get-data/cache")
#'
#'
#' ## Explore Places
#' ### Load most recently cached data workspace
load(paste0("get-data/cache/", list.files("get-data/cache/")[which.max(as.numeric(regmatches(list.files("get-data/cache/"), gregexpr("[[:digit:]]+", list.files("get-data/cache/")))))]))
#' ### Load functions for place-based analyses
#' #### Duration (s): 0.041  
system.time(source("explore-places/code/functions.R"))
#' ### Process inputs for place-based analyses
#' #### Duration (s): 123.906
system.time(source("explore-places/code/process-inputs.R"))
#' ### Estimate species richness
#' #### Duration (s): 208.592
system.time(source("explore-places/code/species-richness.R"))
#' ### Estimate biodiversity uniqueness
#' #### Duration (s): 31.535
system.time(source("explore-places/code/spatial-beta-diversity.R"))
#' ### Estimate species rarity
#' #### Duration (s): 25.793
system.time(source("explore-places/code/species-rarity.R"))
#' ### Estimate network-wide biodiversity importance
#' #### Duration (s): 0.004
system.time(source("explore-places/code/biodiversity-importance.R"))
#' ### Estimate temporal stability
#' #### Duration (s): 1769.712
system.time(source("explore-places/code/temporal-stability.R"))
#' ### Process outputs to de displayed in places app
#' #### Duration: 8.5
system.time(source("explore-places/code/process-outputs.R"))
#' ### Cache .RData
cache_RData("explore-places", "explore-places/cache")
#' ### Clear workspace
rm(list = ls())
#'
#'
#' ## Explore Species
#' ### Load most recently cached data workspace
load(paste0("get-data/cache/", list.files("get-data/cache/")[which.max(as.numeric(regmatches(list.files("get-data/cache/"), gregexpr("[[:digit:]]+", list.files("get-data/cache/")))))]))
#' ### Load functions for species-based analyses
#' #### Duration (s): 0.011  
system.time(source("explore-species/code/functions.R"))
#' ### Process inputs for species-based analyses
#' #### Duration (s): 73.771
system.time(source("explore-species/code/process-inputs.R"))
#' ### Identify species associations in iNaturalist data
#' #### Duration (s): 3.195
system.time(source("explore-species/code/species-associations.R"))
#' ### Identify temporal trends from iNaturalist data
#' #### Duration (s): 4000
system.time(source("explore-species/code/temporal-trends.R"))
#' ### Run species distribution models for a subset of species
#' #### Duration (s): 30
system.time(source("explore-species/code/species-distribution-models.R"))
#' ### Cache .RData
cache_RData("explore-species", "explore-species/cache")


