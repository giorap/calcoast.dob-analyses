#' ---
#' title: DOB Cal Coast <br> <strong> Explore Places - Functions
#' ---
#'
#' ### Process-data
#' #### Extract observations over a given polygon
get_observations_over_polygon <- function(selected_poly){
  selected_poly <- selected_poly %>% st_as_sf() %>% st_transform("+init=epsg:2163")
  selected_observations <- cal_coast_iNat_observations_filtered_sf[which(st_intersects(cal_coast_iNat_observations_filtered_sf %>% st_transform("+init=epsg:2163"), selected_poly, sparse = FALSE)), ]
  return(selected_observations)
}
#' #### Generate place x species abundance data matrix
### Generate a place x species table indicating species observations counts in each visit
get_count_data <- function(observations, group_field = "MPA_ID", extra_grouping_var = NULL){
  
  if (!is.null(extra_grouping_var)){
    grouping_vars <- c(group_field, "scientificName", extra_grouping_var)
  } else {
    grouping_vars <- c(group_field, "scientificName")
  }
    count_data <- observations %>%
      dplyr::filter(taxonRank %in% c("species", "subspecies")) %>% 
      dplyr::group_by_at(grouping_vars) %>%
      dplyr::count(scientificName) %>%
      dplyr::mutate(species_name = gsub(" ", "_", scientificName)) %>%
      dplyr::ungroup() %>%
      dplyr::select(-scientificName) %>%
      tidyr::spread(key = species_name, value = n, fill = 0) 
  
    count_data <- count_data %>% dplyr::filter(complete.cases(.))
    count_data <- count_data %>% dplyr::select(-names(which(colSums(count_data[, -1]) == 0)))
    
  return(count_data)
  
}
#' ### Species Richness
#' #### Wrapper function for running ChaoSpecies to estimate diversity
estimate_richness_from_counts <- purrr::safely(function(count_data, focal_taxon = NULL){
  
  count_data <- count_data[setdiff(names(count_data), c(names(count_data)[grep("ID", names(count_data))], "period"))] 
  
    if (!is.null(focal_taxon)){
      selected_species <- west_coast_taxa %>% dplyr::filter_at(vars(kingdom, phylum, class, order, family, genus), any_vars(str_detect(., pattern = paste0("\\b", focal_taxon, "\\b")))) %>% pull("scientificName") %>% as.character() 
      occurring_species <- intersect(names(count_data), gsub(" ", "_", selected_species))
      count_data <- count_data[occurring_species] 
    } 
    
    ChaoSpecies_out <- count_data %>% 
      ChaoSpecies(datatype = "abundance")
    
    diversity_estimates <- data.frame(
      num_observations = as.numeric(as.character(ChaoSpecies_out$Basic_data_information[1, 2])),
      obs_species = as.numeric(as.character(ChaoSpecies_out$Basic_data_information[2, 2])),
      obs_species_proportion = as.numeric(as.character(ChaoSpecies_out$Basic_data_information[2, 2]))/ncol(count_data),
      sample_completeness = as.numeric(as.character(ChaoSpecies_out$Basic_data_information[3, 2])),
      chao1_estimate = as.numeric(as.character(ChaoSpecies_out$Species_table[3, 1])),
      chao1_error = as.numeric(as.character(ChaoSpecies_out$Species_table[3, 2])),
      chao1_lower = as.numeric(as.character(ChaoSpecies_out$Species_table[3, 3])),
      chao1_upper = as.numeric(as.character(ChaoSpecies_out$Species_table[3, 4]))
    )
    
    diversity_estimates <- diversity_estimates %>% 
      dplyr::mutate(chao1_uncertainty = (chao1_upper * 100)/chao1_estimate - (chao1_lower * 100)/chao1_estimate,
                    chao1_confidence = cut(chao1_uncertainty, breaks = c(0, 10, 30, 50, 100, max(chao1_uncertainty, na.rm = TRUE)), labels = c("very high", "high", "medium", "low", "very low")) 
      )
    
    return(diversity_estimates)
    
  })
#' #### Clean up diversity estimates object
clean_richness_estimates <- function(diversity_estimates, ID_fields, suffix = NULL){
  
  empty_result <- data.frame(
    num_observations = NA,
    obs_species = NA,
    obs_species_proportion = NA,
    sample_completeness = NA,
    chao1_estimate = NA,
    chao1_error = NA,
    chao1_lower = NA,
    chao1_upper = NA,
    chao1_uncertainty = NA,
    chao1_confidence = NA
  ) 
  
  diversity_estimates <- lapply(diversity_estimates, "[[", "result")
  
  for (i in which(unlist(lapply(diversity_estimates, is.null)))){
    diversity_estimates[[i]] <- empty_result
  }
  
  for (i in which(unlist(lapply(diversity_estimates, "[", "num_observations")) < 10)){
    diversity_estimates[[i]] <- empty_result
  }
  
  diversity_estimates <- do.call("rbind", diversity_estimates)
  
  if (!is.null(suffix)) names(diversity_estimates) <- paste(names(diversity_estimates), suffix, sep = "_")
  
  diversity_estimates <- data.frame(ID_fields, diversity_estimates)
  
  return(diversity_estimates)
  
}
#' ### Beta diversity
#' #### Spatial beta diversity
estimate_beta_diversity_from_counts <- function(count_data = cal_coast_count_data_byMPA, min_species = 10, focal_taxon = NULL, suffix = NULL){
  
  count_data_filtered <- count_data[setdiff(names(count_data), c(names(count_data)[grep("ID", names(count_data))], "period"))]
  
  if (!is.null(focal_taxon)){
    selected_species <- west_coast_taxa %>% dplyr::filter_at(vars(kingdom, phylum, class, order, family, genus), any_vars(str_detect(., pattern = paste0("\\b", focal_taxon, "\\b")))) %>% pull("scientificName") %>% as.character() 
    occurring_species <- intersect(names(count_data_filtered), gsub(" ", "_", selected_species))
    count_data_filtered <- count_data_filtered[occurring_species] 
  } 
  
  rows_to_keep <- apply(count_data_filtered, 1, function(x) sum(x > 0) >= min_species)
  count_data_filtered <- count_data_filtered[rows_to_keep, ]
  
  pairwise_dissimilarities <- vegdist(count_data_filtered, method = "horn", na.rm = TRUE) %>% as.matrix() %>% as.data.frame()
  
  uniqueness_estimates <- count_data[rows_to_keep, intersect(names(count_data), c(names(count_data)[grep("ID", names(count_data))], "period"))]
  uniqueness_estimates$uniqueness_median <- apply(pairwise_dissimilarities, 1, function(x) 1 - quantile(x[x > 0], .5, na.rm = TRUE))
  uniqueness_estimates$uniqueness_lower <- apply(pairwise_dissimilarities, 1, function(x) 1 - quantile(x[x > 0], .675, na.rm = TRUE))
  uniqueness_estimates$uniqueness_upper <- apply(pairwise_dissimilarities, 1, function(x) 1 - quantile(x[x > 0], .375, na.rm = TRUE))
  
  if (!is.null(suffix)) names(uniqueness_estimates)[-1] <- paste(names(uniqueness_estimates)[-1], suffix, sep = "_")
  
  return(uniqueness_estimates)
}
#' 
#' #### Convenience function to calculate confidence intervals
#' #### Function taken from package Rmisc 
CI <- function (x, ci = 0.95){
  a <- mean(x)
  s <- sd(x)
  n <- length(x)
  error <- qt(ci + (1 - ci)/2, df = n - 1) * s/sqrt(n)
  return(c(upper = a + error, mean = a, lower = a - error))
}
#' #### Temporal stability
estimate_temporal_stability <- function(observations = cal_coast_observations, boundary_ID = "watershed_ID", temporal_resolution = c("year", "quarter"), min_num_observations = 100, num_iterations = 100){

  temporal_resolution <- match.arg(temporal_resolution)
  
  st_geometry(observations) <- NULL
  
  observations <- observations %>%
    dplyr::mutate(time_period = cal_coast_observations %>% pull(temporal_resolution))
  
  observations_byBoundary <- observations %>% 
    split(observations[, boundary_ID]) 
  
  observations_byBoundary_count <- purrr::map(observations_byBoundary, function(obs){
    obs %>% 
      dplyr::group_by(time_period) %>% 
      dplyr::summarise(count = n()) %>% 
      dplyr::arrange(time_period) %>% 
      dplyr::filter(count >= min_num_observations) %>%
      dplyr::filter(c(diff(time_period), 1) == 1)
  })
  
  temporal_stability_byBoundary <- purrr::map2(observations_byBoundary, observations_byBoundary_count, function(obs, obs_count){
    
    stability <- purrr::map(1:num_iterations, function(iterations){
      count_mat <- obs %>% 
        dplyr::group_by(time_period) %>% 
        dplyr::filter(time_period %in% obs_count$time_period) %>%
        dplyr::sample_n(min(obs_count$count)) %>% 
        get_count_data(group_field = "time_period") %>%
        dplyr::arrange(time_period) 
      dist_mat <- count_mat %>%
        vegan::vegdist(method = "horn", na.rm = TRUE) %>% as.matrix() %>% as.data.frame()
      data.frame(time_period = count_mat$time_period[-1], dist = dist_mat[row(dist_mat) == (col(dist_mat) + 1)])    
    }) %>% join_all("time_period")
    
    names(stability)[-1] <- paste0(names(stability)[-1], 1:(ncol(stability)-1))
    
    stability
    
  })
  
  temporal_stability <- do.call("rbind", temporal_stability_byBoundary) 
  temporal_stability <- temporal_stability %>%
    dplyr::mutate(place = unlist(lapply(strsplit(row.names(temporal_stability), "\\."), "[[", 1))) %>%
    dplyr::select(place, time_period, contains("dist"))
  
  return(temporal_stability)
  
}
#' Function to clean temporal stability estimates
clean_temporal_stability <- function(stability){
  stability <- stability %>% dplyr::mutate_at(vars(contains("dist")), function(y) 1 - y)
  global_stability <- stability %>% dplyr::select(contains("dist")) %>% unlist() %>% CI()
  stability <- stability %>% 
    data.frame(apply(stability %>% dplyr::select(contains("dist")), 1, CI) %>% t()) %>%
    dplyr::select(place, time_period, upper, mean, lower)
  
  stability <- list(place_stability = stability, global_stability = global_stability)
  return(stability)
}

#' ### Weighted endemism
weighted.endemism <- function(we_data, records="single", species="SPECIES", longitude="LONGITUDE", latitude="LATITUDE", frame.raster, deg.resolution=c(0.25,0.25), extent.vector, type="weighted", plot.raster=TRUE, own.weights, weight.type="cell", geo.type="cell", geo.calc="max.dist", outlier_pct=95, verbose=TRUE, own.grid.matrix)
  #Description --
  #
  #Calculates (taxonomic / species) weighted endemism (species richness inversely weighted by species ranges) across gridded maps using single or site-based point records.
  #
  #Usage --
  #
  #For example:
  #
  #require(vegan)
  #
  #endemism_mydata <- weighted.endemism(mite, site.coords = mite.xy, records="site")
#
#endemism_mydata2 <- weighted.endemism(mite, site.coords = mite.xy, records="site", weight.type="geo", own.grid.matrix = endemism_mydata$grid.matrix, frame.raster=endemism_mydata$WE_raster)
#
#endemism_mydata3 <- weighted.endemism(mite, site.coords = mite.xy, records="site", own.weights = endemism_mydata2$weights)
#
#Arguments --
#
#species_records: a data.frame, either with:
#a) rows as individual species records, and columns that include fields for species name, longitude and latitude (see species, longitude, latitude below); or
#b) rows as sites and columns as species, in which case site.coords (below) must also be supplied
#
#records: are the species_records in single/long format (the default, "single") or in site-based/short format (records="site")
#
#site.coords: for site-based data (records="site"), a data.frame with the sites (/field plots) that match the column names of species_records and their longlat coordinates
#
#species: for records="single"; what colname in the supplied species_records contains species names?
#
#latitude: for records="single"; what colname in the supplied species_records contains latitude values?
#
#longitude: for records="single"; what colname in the supplied species_records contains longitude values?
#
#frame.raster: an existing raster object the user can optionally elect to supply as the frame for calculations and mapping
#
#deg.resolution/extent.vector: arguments specifying the map resolution and extent in degrees the user wishes the calculations and mapping to use. If no frame is specified, an arbitrary resolution and extent is supplied. If a frame.raster is specified, these arguments are ignored (function bases mapping on the supplied raster)
#
#type: either "weighted" (default; weighted endemism), or "corrected" (corrected weighted endemism - the 'per-species' weighted endemism as per Crisp et al. (2001). This is provided for convenience, but is not particularly recommended - insterad, use the outputs of weighted.endemism in endemism.null.test function to compare to null expectations of endemism given the species richness)
#
#plot.raster: whether or not to plot the output raster with endemism scores. Either way, the raster object is stored in the output.
#
#own.weights: an optional user-supplied numeric vector of species weights for calculating endemism. Values must have names that are a complete and exact match for the species names in species_records as each species must have a weight. This optional argument is intended mainly so that the more time consuming calculation of geographic weights can be done once and the result stored and used for subsequent re-runs of the endemism calculations (see example under Usage above).
#
#weight.type: default is "cell" (cell-based range weights), while "geo" will calculate geographic range weights. Weight.type "richness" sets the weights to 1, which is equivalent to calculating simple species richness (note setting type="corrected" in this case will give each cell a score of 1). Argument is ignored if own.weights is supplied.
#
#geo.type: default is "cell". This argument is only considered when weight.type is set to "geo" (above), in which case geo.type="cell" calculates geographic ranges based on map grid cell centroids. This can optionally be set to geo.type="point", in which case the geographic range weightings are calculated based on the point locations of each species - this is supplied for reference, but is not especially recommended as it is slower and doesnât provide much more information at the resolution and scale of most analyses. 
#
#geo.calc: default is "max.dist". This argument is only considered when when weight.type is set to "geo" and applies to both geo.type="cell" and geo.type="point". If geo.calc="max.dist" the default maximum geographic distance (âspanâ) between cell centroids/point locations is calculated for a species range. If geo.calc is set to "polygon", the function weights species ranges based on the area of a MCP (minimum convex polygon) that contains all points. Further arguments to this function can be included, such as changing the default outlier_pct=95 (removes outlying locations). Additionally, the geo.type="point" option is not recommended for the "polygon" method, as it is more likely to lead to errors where nearby point locations do not allow drawing of a spanning polygon. In this case, cell-centroid based calculations ensure that multiple records are spatially separated (in different cells) and that occurrences within a single cell are returned as the area of that cell.
#
#outlier_pct: for the calculation of range span or area via convex polygons (at least 5 records of the species), this argument can be used to remove outliers on a percentage basis via the mcp function in package adehabitat, i.e. 95 means 5% most outlying points are removed. Default is 95.
#
#own.grid.matrix: user can supply a binary matrix of species against grid cell numbers, rather than this being generated within the function. The purpose of this argument is that the step can be time consuming for large datasets, so the user can return the matrix that is returned from the function in subsequent runs with different setting to improve speed. If this is supplied, a frame.raster must also be supplied that has cell numbers which match the row.names of the own.grid.matrix
#
#Details --
#
#This implementation of weighted endemism allows alternative calculation of weights for species ranges as well as the option of user-supplied weights. Weights can be calculated based on the frequency of occurrence in grid cells, or alternatively by the geographical size of the species range, calculated in one (span) or two (area) dimensions.
#
#Value --
#
#Returns a list of length 4:
#
#$WE (/$CWE) : vector of weighted endemism scores
#
#$WE_Raster (/$CWE) : raster map with endemism scores
#
#$weights : a named numeric vector of weights used to calculate endemism (equivalent to range size in metres if weight.type="geo", range size in cells if weight.type="cell" (default), or the user supplied weights if own.weights was supplied)
#
#$grid.matrix : a binary data.frame of species against grid cell numbers used in the function which is returned so that it can be re-used in subsequent runs to save time
#
#Required packages --
#
#simba, geosphere, adehabitat, raster
#
#Authors --
#
#Greg R. Guerin & Lasse Ruokolainen
#
#References --
#
#Guerin, G.R., Ruokolainen, L. & Lowe, A.J. (2015) A georeferenced implementation of weighted endemism. Methods in Ecology and Evolution. DOI: 10.1111/2041-210X.12361
#
#License --
#
#GPL-3
#
{
  
  require(raster)
  require(simba)
  require(adehabitatHR)
  require(geosphere)
  
  species_records <- we_data$species.records %>% dplyr::select(-ends_with("_ID"))
  
  site.coords <- we_data$site.coords
  
  if(outlier_pct > 99 | outlier_pct < 1) {
    stop("Outlier_pct should be a percentage")
  }
  
  if(class(species_records) != "data.frame") {
    stop("Species data must be in a data.frame")
  }
  
  if(records == "site") {
    convert <- function(an.occurrence.matrix, site.coords) {
      dat <-  data.frame(SPECIES = "hold",LONGITUDE = 0,LATITUDE = 0)
      nam <-  names(an.occurrence.matrix)
      for(ii in 1:ncol(an.occurrence.matrix)){
        w <-  an.occurrence.matrix[,ii]>0
        dat <-  rbind(dat, setNames(data.frame(rep(nam[ii],sum(w)),site.coords[w,]), names(dat)))
      }
      return(dat[-1,])
    }
    species_records <- convert(species_records, site.coords)
  }
  
  if(records == "single") {
    species_records <- species_records[,c(species, longitude, latitude)]
    colnames(species_records) <- c("SPECIES", "LONGITUDE", "LATITUDE")
  }
  
  if(!("SPECIES" %in% colnames(species_records))) {stop("Cannot locate species data")}
  if(!("LATITUDE" %in% colnames(species_records))) {stop("Cannot locate latitude data")}
  if(!("LONGITUDE" %in% colnames(species_records))) {stop("Cannot locate longitude data")}
  if(any(is.na(species_records$LONGITUDE))) {
    species_records <- species_records[-which(is.na(species_records$LONGITUDE)),] 
  }
  if(any(is.na(species_records$LATITUDE))) {
    species_records <- species_records[-which(is.na(species_records$LATITUDE)),] 
  }
  
  coordinates(species_records) <- c("LONGITUDE", "LATITUDE")
  
  if(missing(frame.raster)) {
    if(!(missing(own.grid.matrix))) {stop("You must supply a frame.raster with cells that match own.grid.matrix")}
    frame.raster <- raster()
    if(missing(extent.vector)) {
      extent(frame.raster)@xmin <- floor(extent(species_records)@xmin)
      extent(frame.raster)@ymin <- floor(extent(species_records)@ymin)
      extent(frame.raster)@xmax <- ceiling(extent(species_records)@xmax)
      extent(frame.raster)@ymax <- ceiling(extent(species_records)@ymax)
    }
    if(!(missing(extent.vector))) {
      extent(frame.raster) <- extent.vector
    }
    res(frame.raster) <- deg.resolution
    cat("Generating frame raster at ", deg.resolution, " resolution and extent defined by: ", extent(frame.raster)@xmin, extent(frame.raster)@xmax, extent(frame.raster)@ymin, extent(frame.raster)@ymax,"\n")
  } 
  
  
  if(!(extent(species_records)@xmin >= extent(frame.raster)@xmin & extent(species_records)@xmax <= extent(frame.raster)@xmax & extent(species_records)@ymin >= extent(frame.raster)@ymin & extent(species_records)@ymax <= extent(frame.raster)@ymax)) {
    cat("Some point locations lie outside the frame raster -- trimming these records", "\n")
    species_record_COORDS <- as.data.frame(coordinates(species_records))
    species_records <- species_records[-which(species_record_COORDS$LONGITUDE < extent(frame.raster)@xmin | species_record_COORDS$LONGITUDE > extent(frame.raster)@xmax | species_record_COORDS$LATITUDE < extent(frame.raster)@ymin | species_record_COORDS$LATITUDE > extent(frame.raster)@ymax),]
  }
  
  
  if(!(missing(own.grid.matrix))) {
    cat("Reading user-defined gridded occurrence matrix", "\n")
    if(class(own.grid.matrix) != "data.frame") {stop("Supplied own.grid.matrix must be a data.frame")}
    cell_occur_matrix <- own.grid.matrix
  } #clse if(!(missing(own.grid...
  
  if(missing(own.grid.matrix)) {
    cat("Generating the gridded occurrence matrix", "\n")
    cell_numbers <- cellFromXY(frame.raster, species_records)
    cell_occur_matrix_prep <- data.frame(cell=cell_numbers, species=species_records$SPECIES, presence=rep(1, length(cell_numbers)))
    cell_occur_matrix_prep$species <- factor(cell_occur_matrix_prep$species)
    if(any(duplicated(cell_occur_matrix_prep))) {
      cell_occur_matrix_prep <- cell_occur_matrix_prep[-which(duplicated(cell_occur_matrix_prep)),]
    }
    if(any(is.na(cell_occur_matrix_prep$cell))) {
      cell_occur_matrix_prep <- cell_occur_matrix_prep[-which(is.na(cell_occur_matrix_prep$cell)),]
    }
    cell_occur_matrix <- mama(cell_occur_matrix_prep)
    cat("Occurrence matrix generated with dimensions: ", dim(cell_occur_matrix), "\n")
  } #cls if(missing(own.grid...
  
  
  
  
  if(!(missing("own.weights"))) {
    if(class(own.weights) != "numeric") {stop("Supplied species weights must be numeric")}
    if(!(all(names(own.weights) %in% colnames(cell_occur_matrix)) & all(colnames(cell_occur_matrix) %in% names(own.weights)))) {stop("Species names for supplied weights are not a complete match for species in supplied records")} #true if all in it both ways
    if(length(own.weights) != length(colnames(cell_occur_matrix))) {stop("Supplied weights vector is different length than number of species - must supply a weighting for each species")}
    
    cat("Calculating user supplied weights", "\n")
    inv_rang_cell_occur_mat <- cell_occur_matrix
    for (i in 1:ncol(inv_rang_cell_occur_mat)) {inv_rang_cell_occur_mat[,i] <- inv_rang_cell_occur_mat[,i]/own.weights[which(names(own.weights) == colnames(inv_rang_cell_occur_mat)[i])]}
    ranges <- own.weights
  } #cls if(!(missing...
  
  
  
  if(missing("own.weights")) {
    
    
    if(weight.type=="cell") {
      cat("Calculating cell-based range weights", "\n")
      inv_rang_cell_occur_mat <- apply(cell_occur_matrix, 2, function(x) {x/sum(x)})
      ranges <- apply(cell_occur_matrix, 2, function(x) {sum(x)}) #can we just do colSums?
    } #cls if(weight.type = cell...
    
    if(weight.type=="richness") {
      inv_rang_cell_occur_mat <- cell_occur_matrix
      ranges <- apply(cell_occur_matrix, 2, function(x) {x = 1})
    } #cls if(w.t = richness...
    
    if(weight.type=="geo") {	
      
      CalcDists <- function(longlats) { #modified from CalcDists.R, see https://gist.githubusercontent.com/sckott/931445/raw/9db1d432b2308a8861f6425f38aaabbce44eb994/CalcDists.R
        name <- list(rownames(longlats), rownames(longlats))
        n <- nrow(longlats)
        z <- matrix(0, n, n, dimnames = name)
        for (i in 1:n) {
          for (j in 1:n) z[i, j] <- distCosine(c(longlats[j, 1], longlats[j, 2]), c(longlats[i, 1], longlats[i, 2]))
        }
        z <- as.dist(z)
        return(z)
      } #cls CalcDists
      
      
      if(geo.type == "cell") {
        
        cell_centroids <- as.data.frame(coordinates(frame.raster))
        colnames(cell_centroids) <- c("LONGITUDE", "LATITUDE")
        cell_centroids$cells <- row.names(cell_centroids)
        
        species_ranges <- function(x) {
          v <- rep(0, ncol(x))
          for (i in 1:ncol(x)) {
            temp <- as.data.frame(x[,i])
            colnames(temp) <- "species_i"
            row.names(temp) <- row.names(x)
            temp$delete <- temp$species_i
            if(any(temp$species_i == 0)) {
              temp <- temp[-which(temp$species_i == 0),]
            } #cls if(any(temp...
            temp$cells <- row.names(temp)
            temp <- merge(temp, cell_centroids, by="cells")
            temp <- data.frame(LONGITUDE=temp$LONGITUDE, LATITUDE=temp$LATITUDE)
            
            if(geo.calc == "max.dist") {
              cell_dimensions <- (mean(values(area(frame.raster)))*1000000)^(0.5) #area in km^2 so converting to m^2 to match areaPolygon, then find square root to get back to 1-dimension
              if(nrow(temp) < 2) {
                v[i] <- cell_dimensions
                warning("Only one record, returning the span of single grid cell for ", colnames(x[i]))
              } else {
                if(nrow(temp) < 5) {
                  v[i] <- max(CalcDists(temp))
                  if(v[i] < cell_dimensions) {
                    v[i] <- cell_dimensions
                    warning("Range is less than spatial grain of frame.raster, returning span of single grid cell for ", colnames(x[i]))
                  } #cls if(v[i] < cell_dimensions)...
                } #cls if nrow(temp) < 5...
                if(nrow(temp) > 4) {
                  spp_i_range_polygon <- try(mcp(temp, id=rep(1, nrow(temp)), percent=outlier_pct))
                  if(class(spp_i_range_polygon)[1] == "try-error") {
                    v[i] <- max(CalcDists(temp))
                  } #cls if(class(spp...
                  if(!class(spp_i_range_polygon)[1] == "try-error") {
                    v[i] <- max(CalcDists(as.data.frame(spp_i_range_polygon[,2:3])))
                  }	 #cls if(!class(spp...
                  if(v[i] < cell_dimensions) {
                    v[i] <- cell_dimensions
                    warning("Range is less than spatial grain of frame.raster, returning span of single grid cell for ", colnames(x[i]))
                  } #cls if(v[i] <...
                }	#cls if(nrow(temp) > 4)...						
              } #cls else...
            } #cls if geo.calc = max.dist...
            
            if(geo.calc == "polygon") {
              cell_size <- mean(values(area(frame.raster)))*1000000 #area in km^2 so convert to m^2 to match areaPolygon
              if(nrow(temp) < 5) {
                polygon_area <- try(areaPolygon(temp))
                if(class(polygon_area) == "try-error") {
                  v[i] <- cell_size
                  warning("Less than minimum number of required points to compute polygon, returning approximate area of single grid cell for ", colnames(x[i]))
                } #cls if(class(polyon...
                if(class(polygon_area) == "numeric") {
                  if(polygon_area < cell_size) {
                    v[i] <- cell_size
                    warning("Polygon area is less than spatial grain of frame.raster, returning approximate area of a single grid cell for ", colnames(x[i]))
                  } #cls if(polygon_area <...
                  if(polygon_area >= cell_size) {
                    v[i] <- polygon_area
                  } #CLS if(class(polygon..
                } #cls if(class(polygon...
              } #cls if(nrow(temp) < 5)...
              if(nrow(temp) > 4) {
                spp_i_range_polygon <- try(mcp(temp, id=rep(1, nrow(temp)), percent=outlier_pct))
                if(class(spp_i_range_polygon)[1] == "try-error") {
                  v[i] <- cell_size
                  warning("Unable to compute polygon, returning approximate area of single grid cell for ", colnames(x[i]))
                } #cls if(class(spp...
                if(class(spp_i_range_polygon)[1] != "try-error") {
                  plot(spp_i_range_polygon, main=paste(polygon_area))
                  polygon_area <- try(areaPolygon(as.data.frame(spp_i_range_polygon[,2:3])))
                  if(class(polygon_area) != "numeric") {
                    v[i] <- cell_size
                    warning("Cannot compute polygon area, returning approximate area of single grid cell for ", colnames(x[i]))
                  } #cls if(class(polyon...	
                  if(class(polygon_area) == "numeric") {
                    if(polygon_area < cell_size) {
                      v[i] <- cell_size
                      warning("Polygon area is less than spatial grain of frame.raster, returning approximate area of single grid cell for ", colnames(x[i]))
                    } #cls if(polygon_area <...
                    if(polygon_area >= cell_size) {
                      v[i] <- polygon_area
                    } #cls if(polygon_area >=...
                  } #cls if(class(polygon...
                } #cls if(!class(spp...		
              } #cls if(nrow(temp) > 4)...
            } #cls if(geo.calc = polygon)
            
            if(verbose) cat("spp complete:", i, v[i], colnames(x[i]), "\n")
            
          } #cls for (i in...
          
          names(v) <- colnames(cell_occur_matrix)
          return(v)
        } #cls species_ranges function
        
        ranges <- species_ranges(cell_occur_matrix)	
        
      } #cls if(geo.type = cell...
      
      
      if(geo.type == "point") {
        
        spp_ranges <- function(x) {
          n <- 0
          v <- rep(0, length(colnames(cell_occur_matrix)))
          x$SPECIES <- sub(pattern = " ", replacement = ".", x = x$SPECIES, fixed=TRUE)
          x$SPECIES <- sub(pattern = "-", replacement = ".", x = x$SPECIES, fixed=TRUE)
          for (i in colnames(cell_occur_matrix)) {
            n <- n + 1
            temp <- x[which(x$SPECIES == i),]
            colnames(temp) <- colnames(x) #?necessary
            temp <- data.frame(LONGITUDE=temp$LONGITUDE, LATITUDE=temp$LATITUDE)
            
            if(geo.calc == "max.dist") {
              cell_dimensions <- (mean(values(area(frame.raster)))*1000000)^(0.5) #area in km^2 so convert to m^2 to match areaPolygon, then find square root to get back to 1-dimension
              if(nrow(temp) < 2) {
                v[n] <- cell_dimensions
                warning("Only one record, returning the size of a single grid cell (1-dimension) for ", i)
              } else {
                if(nrow(temp) < 5) {
                  v[n] <- max(CalcDists(temp))
                  if(v[n] < cell_dimensions) {
                    v[n] <- cell_dimensions
                    warning("Range span less than spatial grain of frame.raster, returning span of single grid cell for ", i)
                  } #cls if(v[n] < cell_dimensions...
                } #cls if(nrow(temp) < 5)...
                if(nrow(temp) > 4) {
                  spp_i_range_polygon <- try(mcp(temp, id=rep(1, nrow(temp)), percent=outlier_pct))
                  if(class(spp_i_range_polygon)[1] == "try-error") {
                    v[n] <- max(CalcDists(temp))
                  } #cls if(class...
                  if(!class(spp_i_range_polygon)[1] == "try-error") {
                    v[n] <- max(CalcDists(as.data.frame(spp_i_range_polygon[,2:3])))
                  } #cls if(!(class...
                  if(v[n] < cell_dimensions) {
                    v[n] <- cell_dimensions
                    warning("Range span less than spatial grain of frame.raster, returning span of single grid cell for ", i)
                  } #cls if(v[n] <...
                } #cls if(nrow(temp) > 4)...
              } #cls else...
            } #cls if(geo.calc == "max.dist")... [three open now...]
            
            if(geo.calc == "polygon") {
              cell_size <- mean(values(area(frame.raster)))*1000000 #returns size in m^2 (convert from km^2)
              if(nrow(temp) < 5) {
                polygon_area <- try(areaPolygon(temp))
                if(class(polygon_area) == "try-error") {	
                  v[n] <- cell_size
                  warning("Cannot compute polygon, returning approximate area of single grid cell for ", i)
                } #cls if(class(polygon_area)...
                if(class(polygon_area) == "numeric") {
                  if(polygon_area < cell_size) {
                    v[n] <- cell_size
                    warning("Polygon area less than spatial grain of frame.raster, returning approximate area of single grid cell for ", i)
                  } #cls if(polygon_area <...
                  if(polygon_area >= cell_size) {
                    v[n] <- polygon_area
                  } #cls if (polygon_area >=...
                } #cls if(class(polyon_area) == "numeric...
              } #cls if(nrow(temp) < 5)...
              if(nrow(temp) > 4) {
                spp_i_range_polygon <- try(mcp(temp, id=rep(1, nrow(temp)), percent=outlier_pct))
                if(class(spp_i_range_polygon)[1] == "try-error") {
                  v[n] <- cell_size
                  warning("Unable to compute a polygon, returning approximate area of single grid cell for ", i)
                } #cls if(class(spp...
                if(class(spp_i_range_polygon)[1] != "try-error") {
                  plot(spp_i_range_polygon, main=paste(polygon_area))
                  polygon_area <- try(areaPolygon(as.data.frame(spp_i_range_polygon[,2:3])))
                  if(class(polygon_area) == "try-error") {
                    v[n] <- cell_size
                    warning("Cannot compute polygon area, returning approximate area of single grid cell for ", i)
                  } #cls if(class(polygon_area)...
                  if(class(polygon_area) == "numeric") {
                    v[n] <- polygon_area
                    if(v[n] < cell_size) {
                      v[n] <- cell_size
                      warning("Polygon area is less than spatial grain of frame.raster, returning approximate area of single grid cell for ", i)
                    } #cls if(v[n] < cell_size)...
                  } #cls if(class(poygon_area == numeric...	
                } #cls if(!class(spp...
              } #cls if(nrow(temp) > 4)...
            } #cls if geo.calc = polygon
            
            if(verbose) cat(n, "species complete:", i, v[n], "\n")
            
          } #cls for (i in colnames...
          
          names(v) <- colnames(cell_occur_matrix)
          return(v)
        } #cls spp_ranges function
        
        ranges <- spp_ranges(species_records)
        
      } #cls if geo.type = point...
      
      
      
      cat("Calculating geographic range weights", "\n")
      inv_rang_cell_occur_mat <- cell_occur_matrix
      for(i in 1:ncol(inv_rang_cell_occur_mat)) {
        inv_rang_cell_occur_mat[,i] <- inv_rang_cell_occur_mat[,i]/ranges[which(names(ranges) == colnames(inv_rang_cell_occur_mat)[i])]
      } #cls for(i in 1:ncol...
      
      
    } #close if(weight.type = geo...
  } #close if(missing(own.weights... 		
  
  cat("Calculating weighted endemism", "\n")
  rawEndemism <- rowSums(inv_rang_cell_occur_mat) 
  
  
  if(type=="weighted") {
    WE_raster <- frame.raster
    WE_raster[] <- NA
    WE_raster[as.numeric(names(rawEndemism))] <- rawEndemism
    if(plot.raster==TRUE) {raster::plot(WE_raster, main="Weighted Endemism")}
    
    rawEndemism <- data.frame(we_data$species.records %>% dplyr::select(ends_with("_ID")),
                              species_rarity = rawEndemism[match(cellFromXY(WE_raster, site.coords), names(rawEndemism))]
                     )
    
    outputs <- list(WE = rawEndemism, WE_raster = WE_raster, weights = ranges, grid.matrix = cell_occur_matrix)
    return(outputs)
  } #cls if(type = w...
  
  if(type=="corrected"){
    richness4endemism <- rowSums(cell_occur_matrix)
    corrected.WE <- rawEndemism/richness4endemism
    corrected.WE_raster <- frame.raster 
    corrected.WE_raster[] <- NA
    corrected.WE_raster[as.numeric(names(corrected.WE))] <- corrected.WE
    if(plot.raster==TRUE) {
      raster::plot(corrected.WE_raster, main="Corrected Weighted Endemism")
    } #cls if(plot...

    outputs <- list(CWE = corrected.WE, corrected.WE_raster = corrected.WE_raster, weights = ranges, grid.matrix = cell_occur_matrix)
    return(outputs)
  } #cls if(type = c...
  
  
} #cls function

#' ### Prepare inputs for weighted endemism function
prepare_weighted_endemism_inputs <- function(count_data, polys){
  count_data <- count_data %>% as.data.frame()
  rows_to_keep <- apply(count_data %>% dplyr::select(-ends_with("_ID")), 1, function(x) sum(x > 0) >= 10)
  count_data <- count_data[rows_to_keep, ]
  IDs <- count_data %>% dplyr::select(ends_with("_ID")) %>% pull()
  polys <- polys[polys$ID %in% IDs, ]
  poly.xy <- do.call("rbind", lapply(slot(polys, "polygons"), function(i) slot(i, "labpt"))) %>% as.data.frame()
  names(poly.xy) <- c("x", "y")
  row.names(poly.xy) <- polys$ID
  count_data <- count_data[match(polys$ID, IDs), ]
  out <- list(species.records = count_data, site.coords = poly.xy)
  return(out)
}

#######################################
##### -- Conservation priority -- #####
#######################################
##### -- Zonation analyses -- #####
##### -- Rasters from counts -- #####
rasters_from_counts <- function(polys = MPA_polys, count_data = cal_coast_count_data_byMPA, min_species = 10, focal_taxon = NULL, out_dir = "output", suffix = "_preseabs_byMPA"){
  
  ID_field <- names(count_data)[grep("ID", names(count_data))]
  
  if (!is.null(focal_taxon)){
    selected_species <- cal_coast_taxa %>% dplyr::filter_at(vars(kingdom, phylum, class, order, family, genus), any_vars(str_detect(., pattern = paste0("\\b", focal_taxon, "\\b")))) %>% pull("scientificName") %>% as.character() 
    occurring_species <- intersect(names(count_data), gsub(" ", "_", selected_species))
    count_data <- count_data[c(ID_field, occurring_species)] 
  } 
  
  rows_to_keep <- apply(count_data %>% dplyr::select(-ID_field), 1, function(x) sum(x > 0) >= min_species)
  count_data <- count_data[rows_to_keep, ]
  
  poly_names <- intersect(polys@data[, "ID"], count_data[, ID_field] %>% pull())
  poly_centroids <- data.frame(poly_names, do.call("rbind", lapply(slot(polys[polys@data[, "ID"] %in% poly_names, ], "polygons"), function(i) slot(i, "labpt"))))
  names(poly_centroids) <- c(ID_field, "centroid_x", "centroid_y")
  poly_centroids <- poly_centroids %>% dplyr::mutate(cellID = cellFromXY(california_coast_grid, poly_centroids[c("centroid_x", "centroid_y")]))
  poly_counts <- left_join(poly_centroids, count_data, by = ID_field)
  species_names <- which(count_data %>% dplyr::select(-ID_field) %>% colSums() > 0) %>% names()
  
  for (i in 1:length(species_names)){
    species_grid <- california_coast_grid
    species_grid[] <- NA
    species_grid[poly_counts$cellID] <- poly_counts[, species_names[i]]
    species_grid[] <- as.numeric(species_grid[] > 0)
    writeRaster(species_grid, filename = paste0(getwd(), "/", out_dir, "/", species_names[i], suffix, sep = ""), format = "GTiff", overwrite = TRUE)
  }

  return(poly_centroids)
}

##### This code requires installation of the "zonator" R package by Atte Moilanen, Joona Lehtomaki et al. 
##### -- create_zproject_edit() -- #####
##### Modified version of zonator::create_zproject, to call alternative templates files to run either core area Zonation (CAZ) or additive benefit function (ABF)
create_zproject_edit <- function(name, dir, variants, dat_template_file = NULL,
                                 spp_template_file = NULL, spp_template_dir = NULL,
                                 overwrite = FALSE, debug = FALSE, ...) {
  if (!file.exists(dir)) {
    stop("Directory ", dir, " provided does not exist.")
  }
  
  # Create the new location
  project_dir <- file.path(dir, name)
  if (file.exists(project_dir)) {
    if (overwrite) {
      if (debug) message("Removing existing directory ", project_dir)
      unlink(project_dir, recursive = TRUE, force = TRUE)
    } else {
      stop("Project ", project_dir, " already exists and overwrite is off")
    }
  }
  
  if (debug) message("Creating a project directory ", project_dir)
  dir.create(project_dir)
  
  # Create an empty README file for the project
  if (debug) message("Creating an empty README file")
  file.create(file.path(project_dir, "README.md"), showWarnings = FALSE)
  
  # Create the variant subfolders with content
  for (variant in variants) {
    variant_dir <- file.path(project_dir, variant)
    if (debug) message("Creating a variant directory ", variant_dir)
    dir.create(variant_dir)
    
    # If no templates are provided, use the ones shipped with zonator. Change
    # the filenames to match the variant.
    if (is.null(dat_template_file)) {
      dat_template_file <- system.file("extdata", paste("template_ABF.dat", sep = ""),
                                       package = "zonator")
    }
    dat_to <- file.path(variant_dir, paste0(variant, ".dat"))
    
    # If no templates are provided, use the ones shipped with zonator. Change
    # the filenames to match the variant.
    if (is.null(spp_template_file) & is.null(spp_template_dir)) {
      spp_template_file <- system.file("extdata", "template.spp",
                                       package = "zonator")
    }
    
    # Define the target variant spp file path
    spp_to <- file.path(variant_dir, paste0(variant, ".spp"))
    
    # Copy the templates to the new variant folder
    if (debug) message("Copying template dat-file ", dat_template_file,
                       " to variant directory ", variant_dir)
    if (file.exists(dat_template_file)) {
      file.copy(from = dat_template_file, to = dat_to, overwrite = TRUE)
    } else {
      stop("dat-file template ", dat_template_file, " not found")
    }
    # Work out the details depending if using a template file or a
    # directory of input rasters.
    if (!is.null(spp_template_dir)) {
      # We may have multiple directories
      if (all(sapply(spp_template_dir, function(x) file.exists(x)))) {
        if (debug) {
          if (length(spp_template_dir) > 1) {
            dir_msg <- paste("Creating a spp file from rasters in directories ",
                             paste(spp_template_dir, collapse = ", "))
          } else{
            dir_msg <- paste("Creating a spp file from rasters in directory ",
                             spp_template_dir)
          }
          message(dir_msg)
        }
        create_spp(filename = spp_to, spp_file_dir = spp_template_dir, ...)
      } else {
        stop("Spp template dir ", spp_template_dir, " not found.")
      }
    } else if (!is.null(spp_template_file)) {
      if (file.exists(spp_template_file)) {
        if (debug) {
          message("Copying template spp-file  ", spp_template_file,
                  " to variant directory ", variant_dir)
        }
        file.copy(from = spp_template_file, to = spp_to, overwrite = TRUE)
      } else {
        stop("Input template spp-file ", spp_template_file, " not found!")
      }
    }
    
    # Create to output folder
    output_dir <- file.path(variant_dir, paste0(variant, "_out"))
    if (debug) message("Creating an output directory ", output_dir)
    dir.create(output_dir, recursive = TRUE)
    # Create a bat file, first read the template content
    bat_from <- system.file("extdata", "template.bat", package = "zonator")
    cmd_sequence <- scan(file = bat_from, "character", sep = " ",
                         quiet = TRUE)
    # Replace tokens with actual (relative) paths
    dat_relative <- gsub(paste0(project_dir, .Platform$file.sep), "", dat_to)
    spp_relative <- gsub(paste0(project_dir, .Platform$file.sep), "", spp_to)
    output_dir_relative <- gsub(paste0(project_dir, .Platform$file.sep), "",
                                output_dir)
    cmd_sequence <- gsub("INPUT_DAT", dat_relative, cmd_sequence)
    cmd_sequence <- gsub("INPUT_SPP", spp_relative, cmd_sequence)
    cmd_sequence <- gsub("OUTPUT", file.path(output_dir_relative,
                                             paste0(variant, ".txt")),
                         cmd_sequence)
    # Write bat-file
    bat_to <- file.path(project_dir, paste0(variant, ".bat"))
    if (debug) message("Writing bat file ", bat_to)
    cat(paste0(paste(cmd_sequence, collapse = " "), "\n"), file = bat_to)
  }
  
  return(invisible(NULL))
}

##### -- run_zonation() -- #####
##### Wrapper function for running multiple iterations (select number of runs using "runs" argument) of the Zonation software
##### using species (input = "species"), trait categories (input = "traits") or phylogenetic nodes (input = "phylo") as input rasters
##### and either the CAZ or ABF algorithms ("algorithm" argument).
##### The "period" argument specifies which period ("present", "DDnt" or "DDt") the input rasters should be selected from
##### This wrapper function uses the create_zproject() (the "_edit" version above) and run_bat() functions of the zonator R package
run_zonation <- function(id = "MPA",
                         runs = 10, 
                         working_dir = getwd(), 
                         zonation_dir = "output/zonation_runs",
                         ...){
  #### Set working directory
  setwd(working_dir)
  #### Create directory to store outputs
  dir.create(zonation_dir, showWarnings = FALSE)
  #### Define project settings
  zonation_project <- create_zproject_edit(name = id, dir = zonation_dir, 
                                           variants = paste(id, 1:runs, sep = ""),
                                           spp_template_dir = paste0(working_dir, "/data/input_rasters/", id),
                                           overwrite = FALSE,
                                           ...)
  #### Run Zonation
  for (run in 1:runs){
    temp_dir <- paste(working_dir, "/output/zonation_runs/", id, sep = "")
    setwd(temp_dir)
    run_bat(paste(id, run, ".bat", sep = ""), exe = paste(working_dir, "/data/zig4", sep = ""))
  }
}

##### -- extract_zonation_output() -- #####
##### Function to extract output rasters for spatial prioritizations generated from Zonation runs
extract_zonation_ranks <- function(zonation_dir, variant_name, runs = 100){
  zonation_rasters <- vector("list", runs)
  ### Extract Zonation rasters from n runs
  for (run in 1:runs){
    zonation_rasters[[run]] <- raster(paste(zonation_dir, "/", variant_name, run, "/", variant_name, run, "_out/", variant_name, run, ".rank.compressed.tif", sep = ""))
  }
  ### Extract ranks from
  zonation_ranks <- lapply(zonation_rasters, function(x){
    target_rank_df <- x %>% as.data.frame(xy = TRUE) %>% mutate(cellID = 1:length(x[])) %>% dplyr::filter(complete.cases(.))
    names(target_rank_df)[3] <- "rank"
    target_rank_df <- target_rank_df[order(target_rank_df$rank, decreasing = TRUE), ]
    target_rank_df
  })
  zonation_ranks <- data.frame(zonation_ranks[[1]], do.call("cbind", lapply(zonation_ranks[-1], "[[", 4)))
  names(zonation_ranks)[-c(1:2)] <- c("rank_priority", paste("rank", 1:(runs), sep = ""))
  ### Merge zonation ranks with the full reference coordinates: this will enable assigning the correct rows with the target matrix later on
  zonation_ranks <- list(zonation_ranks = zonation_ranks, reference_coordinates = as.data.frame(zonation_rasters[[1]], xy = TRUE)[, 1:2])
  return(zonation_ranks) 
}


#######################################
##### -- Conservation priority -- #####
#######################################
##### -- calculate_proportion_protected() -- #####
#### Calculate area of overlap between a polygon and MPAs of different types
calculate_proportion_protected <- function(focal_polygon, min_area = 100){
  
  focal_polygon <- focal_polygon %>% 
    st_as_sf() %>% 
    st_transform("+init=epsg:2163")
  coastal_regions <- region_polygons %>% 
    st_as_sf() %>% 
    st_transform("+init=epsg:2163") %>% 
    st_combine()
  MPA_polys <- MPA_polygons %>% 
    st_as_sf() %>% 
    st_transform("+init=epsg:2163")
    
  # Intersect polygon with coastal area
  focal_polygon <- st_intersection(focal_polygon, coastal_regions) 
  
  proportion_protected_area <- matrix(NA, ncol = length(unique(MPA_polygons$Type)), nrow = 1, dimnames = list(NA, as.character(unique(MPA_polys$Type)))) %>% as.data.frame()
  if (length(focal_polygon$geometry) != 0){
    if (as.numeric(st_area(focal_polygon)) >= min_area){
      for (k in 1:length(unique(MPA_polys$Type))){
        protected_polygon <- st_intersection(focal_polygon, MPA_polys[MPA_polys$Type == unique(MPA_polys$Type)[k], ])
        proportion_protected_area[, k] <- ifelse(!is.null(protected_polygon), st_area(protected_polygon)/st_area(focal_polygon), NA)
      }
      proportion_protected_area <- proportion_protected_area %>% dplyr::mutate(Overall = rowSums(., na.rm = TRUE))  
    }
  } 

  return(proportion_protected_area)
  
}
