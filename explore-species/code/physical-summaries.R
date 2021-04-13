#' ---
#' title: DOB Cal Coast <br> <strong> Explore Species - Physical Summaries
#' ---
#' 
#' ## Generate statewide and regional summaries from ROMS data
#' ### Temperature
temp_means_byRegion <- purrr::map(list("North", "Central", "South", c("North", "Central", "South")), function(region){
  temp_list <- purrr::map(temp_summary_rasters, function(rast){
    raster::extract(rast, region_polygons[region_polygons$Region %in% region, ]) %>% unlist() %>% mean(na.rm = TRUE)
    })
  temp_df <- data.frame(temp_min = temp_list[grep("min", names(temp_list))] %>% unlist(),
                        temp_median = temp_list[grep("median", names(temp_list))] %>% unlist(),
                        temp_max = temp_list[grep("max", names(temp_list))] %>% unlist()
                        )
})
names(temp_means_byRegion) <- c("North", "Central", "South", "Statewide")
#' ### Salinity
salt_means_byRegion <- purrr::map(list("North", "Central", "South", c("North", "Central", "South")), function(region){
  salt_list <- purrr::map(salt_summary_rasters, function(rast){
    raster::extract(rast, region_polygons[region_polygons$Region %in% region, ]) %>% unlist() %>% mean(na.rm = TRUE)
  })
  salt_df <- data.frame(salt_min = salt_list[grep("min", names(salt_list))] %>% unlist(),
                        salt_median = salt_list[grep("median", names(salt_list))] %>% unlist(),
                        salt_max = salt_list[grep("max", names(salt_list))] %>% unlist()
  )
})
names(salt_means_byRegion) <- c("North", "Central", "South", "Statewide")
#' ### Sea level
zeta_means_byRegion <- purrr::map(list("North", "Central", "South", c("North", "Central", "South")), function(region){
  zeta_list <- purrr::map(zeta_summary_rasters, function(rast){
    raster::extract(rast, region_polygons[region_polygons$Region %in% region, ]) %>% unlist() %>% mean(na.rm = TRUE)
  })
  zeta_df <- data.frame(zeta_min = zeta_list[grep("min", names(zeta_list))] %>% unlist(),
                        zeta_median = zeta_list[grep("median", names(zeta_list))] %>% unlist(),
                        zeta_max = zeta_list[grep("max", names(zeta_list))] %>% unlist()
  )
})
names(zeta_means_byRegion) <- c("North", "Central", "South", "Statewide")

