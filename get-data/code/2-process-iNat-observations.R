#' ---
#' title: DOB Cal Coast <br> <strong> Load and process iNaturalist observations
#' ---
#' 
#' ## Load/process species presence-only observations from iNaturalist 
#' #### This code identifies all iNaturalist observations of the broader taxa of interest overlapping spatial boundaries of interest
#'
#' ### Load iNat observations downloaded from GBIF dump via a temporary file
#' #### Download inaturalist-dwca-with-taxa.zip folder
iNat_observations_files <- list.files("get-data/data", "iNat_observations")
if (length(iNat_observations_files) > 0){ # iNat_observations are already downloaded, load them 
  iNat_observations_file_dates <- as.numeric(regmatches(iNat_observations_files, gregexpr("[[:digit:]]+", iNat_observations_files)))
  ### Load most recently downloaded .csv with iNat observations
  iNat_observations <- read_csv(paste0("get-data/data/", iNat_observations_files[which.max(iNat_observations_file_dates)]))
  } else { # Load iNat observations downloaded from GBIF dump via a temporary file
    temp_file <- tempfile()
    download.file("http://inaturalist.org/observations/inaturalist-dwca-with-taxa.zip", temp_file) # Download inaturalist-dwca-with-taxa.zip folder
    iNat_observations <- unz(temp_file, "observations.csv") %>% read.csv(header = TRUE)
    unlink(temp_file) 
}
#' ### Filter out iNat observations without georeferenced coordinates and turn into spatial feature object
iNat_observations_withCoordinates <- iNat_observations %>% 
  dplyr::filter(complete.cases(decimalLongitude, decimalLatitude)) %>% 
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>% 
  dplyr::filter(informationWithheld == "")
#' ### Identify west coast observations
west_coast_observations <- iNat_observations_withCoordinates %>% 
  st_transform("+init=epsg:2163") %>%
  dplyr::filter(st_intersects(iNat_observations_withCoordinates %>% 
                                st_transform("+init=epsg:2163"), 
                              west_coast_polygon %>% 
                                st_as_sf() %>% 
                                st_transform("+init=epsg:2163"),
                              sparse = FALSE
                              )
                ) %>%
  st_transform("+init=epsg:4326") 
#' ### Filter out taxa that are not of interest. In this case, we exclude phyla Plantae (except the three species below), Fungi, and Protozoa, and classes Insecta, Arachnida, Amphibia, and Reptilia
west_coast_observations_filtered <- west_coast_observations %>%
  dplyr::filter(kingdom %in% c("Animalia", "Chromista"), 
                !(class %in% c("Insecta", "Arachnida", "Amphibia", "Reptilia", "Aves", "Mammalia"))  
  ) %>%
  rbind(west_coast_observations %>% dplyr::filter(scientificName %in% c("Phyllospadix scouleri", "Phyllospadix torreyi", "Zostra marina"))) 
#' #### Label "Phyllospadix scouleri", "Phyllospadix torreyi", "Zostra marina" as "Chromista" for the purposes of these analyses
west_coast_observations_filtered[west_coast_observations_filtered$scientificName %in% c("Phyllospadix scouleri", "Phyllospadix torreyi", "Zostra marina"), "kingdom"] <- "Chromista"
#' 
west_coast_observations_filtered <- west_coast_observations_filtered %>%
  cbind(st_coordinates(west_coast_observations_filtered))
#' #### Remove duplicate observations which are likely to be of the same individual(s)
west_coast_observations_filtered <- west_coast_observations_filtered %>% 
  distinct(scientificName, eventDate, X, Y, .keep_all = TRUE)
#' ### Identify sampling visits 
#' #### Sampling visit are unique observation events in which a single observer generated a species list within a discrete space and time
#' #### Here, we define sampling visits as all observations made by a single user on the same date within an area contained within a single grid cell
west_coast_observations_filtered <- west_coast_observations_filtered %>% 
  dplyr::mutate(year = as.numeric(substr(eventDate, 1, 4))) %>% 
  dplyr::mutate(cellID = cellFromXY(west_coast_grid, as(west_coast_observations_filtered, "Spatial"))) %>% #### Generate cellID field based on the cell each observation falls into
  dplyr::filter(complete.cases(cellID)) %>% 
  dplyr::mutate(eventID = paste(cellID, year, sep = "_")) %>% #### Generate eventID field: an event is a combination of site by year
  dplyr::mutate(visitID = paste(eventID, recordedBy, substr(eventDate, 6, 10), sep = "_")) #### Generate visitID field: a visit is a unique day/time by user combination
#" ### Include a field for quarter within a year
west_coast_observations_filtered <- west_coast_observations_filtered %>% 
  dplyr::mutate(yearmonth = as.Date(paste0(substr(west_coast_observations_filtered$eventDate, 1, 7), "-01"), format = "%Y-%m-%d"))
west_coast_observations_filtered <- west_coast_observations_filtered %>% 
  dplyr::mutate(quarter = ifelse(grepl("01-01|02-01|03-01", west_coast_observations_filtered$yearmonth), paste0(substr(west_coast_observations_filtered$yearmonth, 1, 4), "-01"), 
                                 ifelse(grepl("04-01|05-01|06-01", west_coast_observations_filtered$yearmonth), paste0(substr(west_coast_observations_filtered$yearmonth, 1, 4), "-04"),
                                        ifelse(grepl("07-01|08-01|09-01", west_coast_observations_filtered$yearmonth), paste0(substr(west_coast_observations_filtered$yearmonth, 1, 4), "-07"),
                                               ifelse(grepl("10-01|11-01|12-01", west_coast_observations_filtered$yearmonth), paste0(substr(west_coast_observations_filtered$yearmonth, 1, 4), "-10"),
                                                      NA
                                               )
                                        )
                                 )
  )
  )
#' Create a data.frame object
west_coast_observations_filtered_df <- west_coast_observations_filtered %>% 
  dplyr::mutate(longitude = west_coast_observations_filtered %>% st_coordinates() %>% as.data.frame() %>% pull(X), latitude = west_coast_observations_filtered %>% st_coordinates() %>% as.data.frame() %>% pull(Y)) %>%
  st_set_geometry(NULL)
#' ### Extract out taxonomy data 
west_coast_taxa <- west_coast_observations_filtered %>%
  dplyr::select(scientificName, genus, family, order, class, phylum, kingdom) %>%
  unique(.)
#' ### Isolate California observations
#' #### Using stateProvince field
cal_coast_observations1 <- west_coast_observations_filtered %>% 
  dplyr::filter(stateProvince == "California")
#' #### Using spatial boundary overlap
cal_coast_observations2 <- west_coast_observations_filtered %>% 
  st_transform("+init=epsg:2163") %>%
  dplyr::filter(st_intersects(west_coast_observations_filtered %>% 
                                st_transform("+init=epsg:2163"), 
                              region_polygons %>% 
                                st_as_sf() %>% 
                                st_transform("+init=epsg:2163") %>% 
                                st_combine(),
                              sparse = FALSE
  )
  ) %>%
  st_transform("+init=epsg:4326")
#' ### Bind the two results and keep only distinct rows
cal_coast_observations <- rbind(cal_coast_observations1, cal_coast_observations2) %>%
  dplyr::distinct(.keep_all = TRUE)
#' ### Identify which iNat observations fall within which polygon 
#' #### Attribute each observation to the corresponding spatial boundary polygon 
cal_coast_observations <- cal_coast_observations %>% 
  st_transform("+init=epsg:2163") %>%
  st_join(region_polygons %>% st_as_sf() %>% st_transform("+init=epsg:2163") %>% dplyr::select(ID) %>% dplyr::transmute(region_ID = ID), join = st_intersects) %>%
  st_join(MPA_polygons %>% st_as_sf() %>% st_transform("+init=epsg:2163") %>% dplyr::select(ID) %>% dplyr::transmute(MPA_ID = ID), join = st_intersects) %>%
  st_join(county_polygons %>% st_as_sf() %>% st_transform("+init=epsg:2163") %>% dplyr::select(ID) %>% dplyr::transmute(county_ID = ID), join = st_intersects) %>%
  st_join(watershed_polygons %>% st_as_sf() %>% st_transform("+init=epsg:2163") %>% dplyr::select(ID) %>% dplyr::transmute(watershed_ID = ID), join = st_intersects) %>%
  st_join(hexagon_polygons %>% st_as_sf() %>% st_transform("+init=epsg:2163") %>% dplyr::select(ID) %>% dplyr::transmute(hexagon_ID = ID), join = st_intersects) %>%
  st_transform("+init=epsg:4326")
#' #### Add a within_mpa versus outside_mpa
cal_coast_observations <- cal_coast_observations %>% 
  dplyr::mutate(within_MPA = factor(ifelse(!is.na(MPA_ID), "within_MPA", "outside_MPA"))) %>% 
  dplyr::mutate(within_MPA_byRegion = paste(within_MPA, region_ID, sep = "_"))
#' ### Create a data.frame object
cal_coast_observations_df <- cal_coast_observations %>% 
  dplyr::mutate(longitude = cal_coast_observations %>% st_coordinates() %>% as.data.frame() %>% pull(X), latitude = cal_coast_observations %>% st_coordinates() %>% as.data.frame() %>% pull(Y)) %>% 
  st_set_geometry(NULL)
#' ### Output useful objects 
# saveRDS(west_coast_observations, "get-data/data/west_coast_observations.rds")
# saveRDS(west_coast_observations_filtered, "get-data/data/west_coast_observations_filtered.rds")
# saveRDS(cal_coast_observations, "get-data/data/cal_coast_observations.rds")
# saveRDS(west_coast_taxa, "get-data/data/west_coast_taxa.rds")
#' ### Remove no longer needed objects
rm(iNat_observations, iNat_observations_withCoordinates, west_coast_observations, cal_coast_observations1, cal_coast_observations2)
