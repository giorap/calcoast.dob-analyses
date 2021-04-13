#' ---
#' title: Species-based biodiversity indicators from iNaturalist observations <br> <strong> Functions
#' ---
#' 
#' ## Load/process environmental predictor data
#' 
#' ### fill_na
#' #### Function to fill na values in extended ROMS rasters
#' 
fill_na <- function(x, i = 13) {
  if( is.na(x)[i] ) {
    return(mean(x, na.rm=TRUE))
  } else {
    return(round(x[i], 0))
  }
}
#' ### fill_na_habitat
#' #### Function to fill na values in ESI habitat
#' 
fill_na_habitat <- function(x, i = 13) {
  if( is.na(x)[i] ) {
    return(raster::modal(x, na.rm=TRUE))
  } else {
    return(round(x[i], 0))
  }
}
#'
#' #### raster_from_nc
#' ### Function to generate extended rasters from loaded nc files
#'
raster_from_nc <- function(nc_file, var = "temp"){
  print(nc_file)
  f_nc <- nc_open(nc_file)
  f_df <- data.frame(
    x = ncvar_get(f_nc, "lon_rho") %>% matrix(ncol = 1),
    y = ncvar_get(f_nc, "lat_rho") %>% matrix(ncol = 1),
    var = ncvar_get(f_nc, var) %>% matrix(ncol = 1)
  )
  names(f_df)[3] <- var
  nc_close(f_nc)
  f_r <- rasterFromXYZ(f_df, crs = CRS("+init=epsg:4326"))
  f_r[which(is.na(f_r[]))] <- NA
  f_r <- focal(f_r, w = matrix(1, 5, 5), fun = fill_na, pad = TRUE, na.rm = FALSE)
  f_r <- f_r %>% mask(accessible_area_polygon)
  f_r
}
#' ### extract_roms
#' #### Function to extract ROMS nc files and generate extended rasters from them
#' #### the function is run in safe mode using purrr::safely() to run despite errors but make a record of them
#' 
extract_roms <- function(var = "temp", year = "2011", rm_files = c(1, 93, 264, 364, 951, 1289, 1817)){
  ##### Get ROMS raster files
  roms_files <- list.files("X://grapacciuolo/ROMS_download/ROMS_ncdf/", full.names = TRUE, pattern = paste0("West-Coast-ROMS_var=", var))
  roms_files <- roms_files[grep(paste0("time=", year), roms_files)]
  roms_rasters <- purrr::map(roms_files, purrr::safely(raster_from_nc), var = var)
  roms_rasters <- lapply(roms_rasters, "[[", "result")
  names(roms_rasters) <- roms_files
  roms_rasters <- roms_rasters[unlist(lapply(roms_rasters, function(x) !is.null(x)))]
  
  saveRDS(roms_rasters, paste0("C://Users/grapacciuolo/Dropbox/cal_academy-citizen_science/data/CeNCOOS-ROMS/rasters/", var, "_", year, "_rasters_west_coast.rds"))
  
}
#'
#' ### get_summary_rasters
get_summary_rasters <- function(var = "temp"){
  
  years <- c(as.list(2011:2019), list(2011:2012), list(2013:2014), list(2014:2015), list(2015:2016), list(2016:2017), list(2017:2018), list(2018:2019), list(2011:2019))  
  
  summary_rasters <- purrr::map(years, function(year){
    
    if (length(year) == 1){
      roms_rasters <- readRDS(paste0("C://Users/grapacciuolo/Dropbox/cal_academy-citizen_science/data/CeNCOOS-ROMS/rasters/", var, "_", year, "_rasters_west_coast.rds"))
    } else {
      roms_rasters <- purrr::map(year, function(multiple){
        readRDS(paste0("C://Users/grapacciuolo/Dropbox/cal_academy-citizen_science/data/CeNCOOS-ROMS/rasters/", var, "_", multiple, "_rasters_west_coast.rds"))
      }
      ) 
    }
    
    roms_raster_stack <- roms_rasters %>% unlist() %>% raster::stack()
    raster_median <- roms_raster_stack %>% raster::calc(function(x){quantile(x, probs = 0.50, na.rm = TRUE)})
    raster_min <- roms_raster_stack %>% raster::calc(function(x){quantile(x, probs = 0.10, na.rm = TRUE)})
    raster_max <- roms_raster_stack %>% raster::calc(function(x){quantile(x, probs = 0.90, na.rm = TRUE)})
    
    list(median = raster_median, min = raster_min, max = raster_max)
    
  })
  
  names(summary_rasters) <- map(years, paste0, collapse = "-") %>% unlist()
  
  return(summary_rasters)
  
}
#' ### get_observations_count_raster
#' #### Function to count iNat observations in each grid cell and generate raster
get_effort_rasters <- function(observations = west_coast_observations_filtered, background_raster = west_coast_grid, selected_year = NULL, logged = TRUE){
  
  observations <- observations
  
  if (!is.null(selected_year)) observations <- observations %>% dplyr::filter(year %in% selected_year)
  
  # Total number of observations
  num_observations <- observations %>% dplyr::count(cellID)
  num_observations <- data.frame(xyFromCell(background_raster, num_observations$cellID), num_observations)
  num_observations_raster <- rasterFromXYZ(num_observations[c("x", "y", "n")], crs = CRS("+init=epsg:4326"))
  num_observations_raster <- raster::extend(num_observations_raster, extent(background_raster))
  num_observations_raster[is.na(num_observations_raster[])] <- 0
  
  # Number of observers
  num_observers <- observations %>% dplyr::group_by(cellID) %>% dplyr::summarize(n = n_distinct(recordedBy))
  num_observers <- data.frame(xyFromCell(background_raster, num_observers$cellID), num_observers)
  num_observers_raster <- rasterFromXYZ(num_observers[c("x", "y", "n")], crs = CRS("+init=epsg:4326"))
  num_observers_raster <- raster::extend(num_observers_raster, extent(background_raster))
  num_observers_raster[is.na(num_observers_raster[])] <- 0
  
  # Number of visits
  num_visits <- observations %>% dplyr::group_by(cellID) %>% dplyr::summarize(n = n_distinct(visitID))
  num_visits <- data.frame(xyFromCell(background_raster, num_visits$cellID), num_visits)
  num_visits_raster <- rasterFromXYZ(num_visits[c("x", "y", "n")], crs = CRS("+init=epsg:4326"))
  num_visits_raster <- raster::extend(num_visits_raster, extent(background_raster))
  num_visits_raster[is.na(num_visits_raster[])] <- 0
  
  # Number of comprehensive visits
  high_effort_visits <- observations %>% group_by(visitID) %>% dplyr::summarize(count = n()) %>% dplyr::filter(count >= 10) %>% pull(visitID)
  high_effort_observations <- observations %>% dplyr::filter(visitID %in% high_effort_visits)
  num_high_effort_visits <- high_effort_observations %>% dplyr::group_by(cellID) %>% dplyr::summarize(n = n_distinct(visitID, .drop = FALSE))
  num_high_effort_visits <- data.frame(xyFromCell(background_raster, num_high_effort_visits$cellID), num_high_effort_visits)
  num_high_effort_visits_raster <- rasterFromXYZ(num_high_effort_visits[c("x", "y", "n")], crs = CRS("+init=epsg:4326"))
  num_high_effort_visits_raster <- raster::extend(num_high_effort_visits_raster, extent(background_raster))
  num_high_effort_visits_raster[is.na(num_high_effort_visits_raster[])] <- 0
  
  output_rasters <- list(num_observations = num_observations_raster, num_observers = num_observers_raster, num_visits = num_visits_raster, num_high_effort_visits = num_high_effort_visits_raster)
  
  if (isTRUE(logged)) output_rasters <- purrr::map(output_rasters, function(r) log(r + 1))
  
  return(output_rasters)
  
}
#'
#' ### Get association matrix
get_species_association_matrix <- function(observations = cal_coast_observations_df){
  
  common_species <- observations %>%
    dplyr::group_by(scientificName) %>%
    dplyr::count(scientificName) %>% 
    dplyr::filter(n >= 100) %>%
    pull(scientificName)
  
  count_data <- observations %>%
    dplyr::filter(scientificName %in% common_species) %>%
    dplyr::group_by(visitID, scientificName) %>%
    dplyr::count(scientificName) %>%
    dplyr::mutate(species_name = gsub(" ", "_", scientificName)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-scientificName) %>%
    tidyr::spread(key = species_name, value = n, fill = 0) 
  
  #### Generate a reduced set that includes only visits with at least two species and species with at least 50 observations across those visits
  count_data_reduced <- count_data[, -1]
  count_data_reduced <- count_data_reduced[which(rowSums(count_data_reduced) >= 2), ]
  count_data_reduced <- count_data_reduced[, which(colSums(count_data_reduced) >= 50)]
  count_data_reduced <- count_data_reduced[which(rowSums(count_data_reduced) >= 2), ]
  #' #### Convert count data to presence absence
  presabs <- count_data_reduced %>% 
    dplyr::mutate_all(~as.numeric(. > 0))
  #' #### Keep only visits with at least two species observed
  presabs <- presabs[which(rowSums(presabs) >= 2), ]
  #' #### Flip site x species matrix to obtain species x sites matrix
  #' #### This is necessary to estimate associations among species, not sites, within vegdist
  presabs_flipped <- t(presabs) %>% as.data.frame() 
  
  return(presabs_flipped)
  
}

#' ### multispeciesPP_edit
multispeciesPP_edit <- 
  function(sdm.formula, bias.formula, PA, PO, BG, species = names(PO), 
           species.PA = species, species.PO = species, quadrat.size = 1, 
           region.size = 1, start = NULL, inverse.hessian = FALSE, penalty.l2.sdm = 0.1, 
           penalty.l2.bias = 0.1, penalty.l2.intercept = 1e-04, weights = rep(1, n.species * nrow(x)), control = list()){
    control <- do.call("glm.control", control)
    species <- union(species.PO, species.PA)
    sdm.formula <- update(sdm.formula, ~. + 1)
    bias.formula <- update(bias.formula, ~. - 1)
    sdm.mf <- model.frame(sdm.formula, data = BG)
    bias.mf <- model.frame(bias.formula, data = BG)
    sdm.BG.model.matrix <- model.matrix(terms(sdm.mf), BG)
    sdm.means <- c(0, apply(sdm.BG.model.matrix[, -1, drop = FALSE], 
                            2, mean))
    sdm.BG.model.matrix <- sweep(sdm.BG.model.matrix, 2, sdm.means, 
                                 "-")
    sdm.sds <- c(1, apply(sdm.BG.model.matrix[, -1, drop = FALSE], 
                          2, sd))
    sdm.BG.model.matrix <- sweep(sdm.BG.model.matrix, 2, sdm.sds, 
                                 "/")
    sdm.standardize <- function(mat) sweep(sweep(mat, 2, sdm.means, 
                                                 "-"), 2, sdm.sds, "/")
    bias.BG.model.matrix <- model.matrix(terms(bias.mf), BG)
    bias.means <- apply(bias.BG.model.matrix, 2, mean)
    bias.BG.model.matrix <- sweep(bias.BG.model.matrix, 2, bias.means, 
                                  "-")
    bias.sds <- apply(bias.BG.model.matrix, 2, sd)
    bias.BG.model.matrix <- sweep(bias.BG.model.matrix, 2, bias.sds, 
                                  "/")
    bias.standardize <- function(mat) sweep(sweep(mat, 2, bias.means, 
                                                  "-"), 2, bias.sds, "/")
    BG.good.rows <- intersect(rownames(sdm.BG.model.matrix), 
                              rownames(bias.BG.model.matrix))
    sdm.PA.model.matrix <- sdm.standardize(model.matrix(terms(sdm.mf), 
                                                        PA))
    PA.good.rows <- rownames(sdm.PA.model.matrix)
    if (!is.null(species.PO)) {
      sdm.PO.model.matrices <- lapply(as.list(species.PO), 
                                      function(sp) sdm.standardize(model.matrix(terms(sdm.mf), 
                                                                                PO[[sp]])))
      names(sdm.PO.model.matrices) <- species.PO
      bias.PO.model.matrices <- lapply(as.list(species.PO), 
                                       function(sp) bias.standardize(model.matrix(terms(bias.mf), 
                                                                                  PO[[sp]])))
      names(bias.PO.model.matrices) <- species.PO
      PO.good.rows <- lapply(as.list(species.PO), function(sp) intersect(rownames(sdm.PO.model.matrices[[sp]]), 
                                                                         rownames(bias.PO.model.matrices[[sp]])))
      names(PO.good.rows) <- species.PO
    }
    n.species <- length(species)
    p.sdm <- ncol(sdm.BG.model.matrix) - 1
    p.bias <- ncol(bias.BG.model.matrix)
    sdm.margins.ab <- matrix(0, n.species, p.sdm + 1, dimnames = list(species, 
                                                                      colnames(sdm.BG.model.matrix)))
    sdm.margins.gamma <- matrix(0, n.species, 1, dimnames = list(species, 
                                                                 "isPO"))
    bias.margins <- matrix(0, 1, p.bias, dimnames = list(NULL, 
                                                         colnames(bias.BG.model.matrix)))
    for (sp in species.PO) {
      k <- match(sp, species)
      sdm.margins.ab[k, ] <- colSums(sdm.PO.model.matrices[[sp]][PO.good.rows[[sp]], 
                                                                 , drop = FALSE])
      sdm.margins.gamma[k, ] <- length(PO.good.rows[[sp]])
      bias.margins <- bias.margins + colSums(bias.PO.model.matrices[[sp]][PO.good.rows[[sp]], 
                                                                          , drop = FALSE])
    }
    abcd.from.all.coef <- function(all.coef) {
      sdm.coef <- matrix(all.coef[1:(n.species * (p.sdm + 2))], 
                         p.sdm + 2, n.species)
      alpha <- sdm.coef[1, ]
      beta <- t(sdm.coef[2:(p.sdm + 1), , drop = FALSE])
      gamma <- sdm.coef[p.sdm + 2, ]
      delta <- all.coef[-(1:(n.species * (p.sdm + 2)))]
      names(alpha) <- names(gamma) <- species
      colnames(beta) <- colnames(sdm.margins.ab)[-1]
      rownames(beta) <- species
      names(delta) <- colnames(bias.BG.model.matrix)
      return(list(alpha = alpha, beta = beta, gamma = gamma, 
                  delta = delta))
    }
    all.coef.from.abcd <- function(alpha, beta, gamma, delta) {
      c(rbind(alpha, beta, gamma), delta)
    }
    n.PA <- length(PA.good.rows)
    n.BG <- length(BG.good.rows)
    subsamp.PA.offset <- 0
    subsamp.BG.offset <- 0
    n.sites <- n.BG + n.PA
    x <- cbind(rbind(sdm.margins.ab, 0, sdm.PA.model.matrix[PA.good.rows, 
                                                            , drop = FALSE], sdm.BG.model.matrix[BG.good.rows, , 
                                                                                                 drop = FALSE]), c(sdm.margins.gamma, rep(0:1, c(1 + n.PA, 
                                                                                                                                                 n.BG))))
    x <- rbind(x, diag(sqrt(c(penalty.l2.intercept, rep(penalty.l2.sdm, 
                                                        p.sdm), penalty.l2.intercept))), matrix(0, p.bias, p.sdm + 
                                                                                                  2))
    z <- rbind(matrix(0, n.species, p.bias), bias.margins, matrix(0, 
                                                                  n.PA, p.bias), bias.BG.model.matrix[BG.good.rows, , drop = FALSE], 
               matrix(0, p.sdm + 2, p.bias), sqrt(penalty.l2.bias/n.species) * 
                 diag(p.bias))
    y <- rep(0, nrow(x) * n.species)
    offset <- rep(0, nrow(x) * n.species)
    for (k in 1:n.species) {
      yk <- rep(0, nrow(x))
      yk[1:n.species] <- 1 * (1:n.species == k)
      yk[1 + n.species] <- 1 * (1 == k)
      if (species[k] %in% species.PA) {
        yk[1 + n.species + (1:n.PA)] <- PA[PA.good.rows, 
                                           species[k]]
      }
      else {
        yk[1 + n.species + (1:n.PA)] <- NA
      }
      if (species[k] %in% species.PO) {
        yk[1 + n.species + n.PA + (1:n.BG)] <- 0
      }
      else {
        yk[1 + n.species + n.PA + (1:n.BG)] <- NA
      }
      yk[1 + n.species + n.sites + (1:(p.sdm + 2 + p.bias))] <- 0
      y[(k - 1) * nrow(x) + 1:nrow(x)] <- yk
      offk <- rep(0, nrow(x))
      offk[1 + n.species + (1:n.PA)] <- log(quadrat.size)
      offk[1 + n.species + n.PA + (1:n.BG)] <- log(region.size) - 
        log(n.BG)
      offset[(k - 1) * nrow(x) + 1:nrow(x)] <- offk
    }
    which.PA <- (2 + n.species):(1 + n.species + n.PA) + rep((0:(n.species - 
                                                                   1)) * nrow(x), each = n.PA)
    which.BG <- (2 + n.species + n.PA):(1 + n.species + n.PA + 
                                          n.BG) + rep((0:(n.species - 1)) * nrow(x), each = n.BG)
    if (is.null(start)) {
      start.alpha <- start.gamma <- rep(0, n.species)
      for (k in 1:n.species) {
        if ((species[k] %in% species.PA) && sum(!is.na(PA[PA.good.rows, 
                                                          species[k]]) > 0)) 
          start.alpha[k] <- log((1 + sum(PA[PA.good.rows, 
                                            species[k]], na.rm = TRUE))/n.PA/quadrat.size)
        if (species[k] %in% species.PO) 
          start.gamma[k] <- log1p(sdm.margins.gamma[k, 
                                                    ]) - start.alpha[k] - log(region.size)
      }
      start <- all.coef.from.abcd(start.alpha, matrix(0, p.sdm, 
                                                      n.species), start.gamma, rep(0, p.bias))
    }
    fit <- block.glm.fit(x, z, y, weights = weights, start = start, 
                         offset = offset, families = list(linear(), binomial(link = "cloglog"), 
                                                          poisson(), gaussian()), row.families = rep(rep(1:4, 
                                                                                                         c(1 + n.species, n.PA, n.BG, p.sdm + p.bias + 2)), 
                                                                                                     n.species), control = control)
    all.coef <- fit$coefficients
    eta <- fit$linear.predictors
    mu <- fit$fitted.values
    names(all.coef)[1:(n.species * (p.sdm + 2))] <- paste(rep(species, 
                                                              each = p.sdm + 2), c(colnames(sdm.BG.model.matrix)[1:(p.sdm + 
                                                                                                                      1)], "isPO"), sep = ":")
    names(all.coef)[-(1:(n.species * (p.sdm + 2)))] <- paste("isPO:", 
                                                             colnames(bias.BG.model.matrix), sep = "")
    std.errs <- fit$fit$std.errs
    names(std.errs) <- names(all.coef)
    species.coef <- matrix(all.coef[1:(n.species * (p.sdm + 2))], 
                           p.sdm + 2, n.species, dimnames = list(c(colnames(sdm.margins.ab), 
                                                                   "isPO"), species))
    bias.coef <- all.coef[-(1:(n.species * (p.sdm + 2)))]
    names(bias.coef) <- colnames(bias.BG.model.matrix)
    fit.PA <- linear.fit.PA <- matrix(NA, nrow(PA), length(species), 
                                      dimnames = list(dimnames(PA)[[1]], species))
    linear.fit.PA[PA.good.rows, ] <- eta[which.PA]
    fit.PA[PA.good.rows, ] <- mu[which.PA]
    fit.BG <- linear.fit.BG <- bias.fit.BG <- linear.bias.fit.BG <- matrix(NA, 
                                                                           nrow(BG), length(species), dimnames = list(dimnames(BG)[[1]], 
                                                                                                                      species))
    linear.fit.BG[BG.good.rows, ] <- matrix(eta[which.BG], ncol = n.species) + 
      log(n.BG) - log(region.size)
    fit.BG[BG.good.rows, ] <- matrix(mu[which.BG], ncol = n.species) * 
      n.BG/region.size
    linear.bias.fit.BG[BG.good.rows, ] <- c(bias.BG.model.matrix[BG.good.rows, 
                                                                 , drop = FALSE] %*% bias.coef)
    bias.fit.BG[BG.good.rows, ] <- exp(linear.bias.fit.BG[BG.good.rows, 
                                                          ])
    fitted.sdm.margins.gamma <- colSums(fit.BG[BG.good.rows, 
                                               , drop = FALSE]) * region.size/n.BG
    fitted.bias.margins <- colSums(t(fit.BG[BG.good.rows, species.PO, 
                                            drop = FALSE]) %*% bias.BG.model.matrix[BG.good.rows, 
                                                                                    , drop = FALSE] * region.size/n.BG)
    score.check.gamma <- fitted.sdm.margins.gamma - sdm.margins.gamma + 
      penalty.l2.intercept * species.coef[p.sdm + 2, ]
    score.check.gamma <- score.check.gamma[species %in% species.PO]
    score.check.bias <- fitted.bias.margins - bias.margins + 
      penalty.l2.bias * bias.coef
    if (length(score.check.gamma) > 0) 
      stopifnot(mean((score.check.gamma/fit$deviance)^2) < 
                  control$epsilon)
    stopifnot(mean((score.check.bias/fit$deviance)^2) < control$epsilon)
    sd.normalizer <- c(rep(c(sdm.sds, 1), n.species), bias.sds)
    unstandardized.coef <- all.coef/sd.normalizer
    gamma.adjust <- sum(unstandardized.coef[-(1:(n.species * 
                                                   (p.sdm + 2)))] * bias.means)
    for (k in 1:n.species) {
      jk <- (p.sdm + 2) * (k - 1) + 1:(p.sdm + 1)
      coef.block <- unstandardized.coef[jk]
      unstandardized.coef[jk[1]] <- coef.block[1] - sum(coef.block[-1] * 
                                                          sdm.means[-1])
      unstandardized.coef[jk[1] + p.sdm + 1] <- unstandardized.coef[jk[1] + 
                                                                      p.sdm + 1] - gamma.adjust
    }
    unstandardized.species.coef <- matrix(unstandardized.coef[1:(n.species * 
                                                                   (p.sdm + 2))], p.sdm + 2, n.species, dimnames = list(c(colnames(sdm.margins.ab), 
                                                                                                                          "isPO"), species))
    unstandardized.bias.coef <- unstandardized.coef[-(1:(n.species * 
                                                           (p.sdm + 2)))]
    names(unstandardized.bias.coef) <- colnames(bias.BG.model.matrix)
    tr <- list(sdm.formula = sdm.formula, bias.formula = bias.formula, model = fit, 
               normalized.species.coef = species.coef, normalized.bias.coef = bias.coef, 
               normalized.all.coef = all.coef, normalized.std.errs = std.errs, 
               all.coef = unstandardized.coef, std.errs = std.errs/sd.normalizer, 
               species.coef = unstandardized.species.coef, bias.coef = unstandardized.bias.coef, 
               linear.fit.PA = linear.fit.PA, fit.PA = fit.PA, linear.bias.fit.BG = linear.bias.fit.BG, 
               bias.fit.BG = bias.fit.BG, linear.fit.BG = linear.fit.BG, 
               fit.BG = fit.BG)
    class(tr) <- c("multispeciesPP", "list")
    tr
  }
#'
#' get_predictor_values
get_predictor_values <- function(coords, x, y, yr, clim_predictors = c("previous_year", "same_year")){
  
  clim_predictors <- match.arg(clim_predictors)
  
  clim_yr <- ifelse(clim_predictors == "previous_year", yr - 1, yr)
  
  predictors <- coords
  
  predictors$temp_median <- raster::extract(temp_summary_rasters[[paste0(clim_yr, ".median")]], predictors[c(x, y)])
  predictors$temp_min <- raster::extract(temp_summary_rasters[[paste0(clim_yr, ".min")]], predictors[c(x, y)])
  predictors$salt_median <- raster::extract(salt_summary_rasters[[paste0(clim_yr, ".median")]], predictors[c(x, y)])
  predictors$h_median <- raster::extract(h_summary_rasters[[paste0(clim_yr, ".median")]], predictors[c(x, y)])
  predictors$zeta_median <- raster::extract(zeta_summary_rasters[[paste0(clim_yr, ".median")]], predictors[c(x, y)])
  predictors$zeta_max <- raster::extract(zeta_summary_rasters[[paste0(clim_yr, ".max")]], predictors[c(x, y)])
  predictors$habitat <- raster::extract(west_coast_ESI_raster, predictors[c(x, y)])
  predictors$habitat <- factor(predictors$habitat, levels = levels(west_coast_ESI_raster_fct)[[1]]$ID)
  predictors$num_observations <- raster::extract(west_coast_effort_rasters[[as.character(yr)]]$num_observations, predictors[c(x, y)])
  predictors$num_observers <- raster::extract(west_coast_effort_rasters[[as.character(yr)]]$num_observers, predictors[c(x, y)])
  predictors$num_high_effort_visits <- raster::extract(west_coast_effort_rasters[[as.character(yr)]]$num_high_effort_visits, predictors[c(x, y)])
  predictors <- predictors %>% dplyr::filter(complete.cases(.))

  return(predictors)
}
#'
#' run_sdm
run_sdm <- function(species_set = NULL, 
                    lumping_codes = NULL,
                    focal_species_name = NULL, 
                    num_species = 10, 
                    survey_dat = survey_data_biodiversity_presabs,
                    by_year = TRUE
                    ){
  
  if (!is.null(focal_species_name)){
    species_set <- get_associated_species(focal_species = focal_species_name, num_species = num_species, association_matrix = cal_coast_species_associations) %>% 
      names()
    species_set <- gsub("_", " ", species_set)
  }
  
  background_byYear <- purrr::map(2012:2019, function(yr){
    yr_bg <- rasterToPoints(temp_summary_rasters[[paste0(yr, ".median")]], spatial = FALSE) %>% as.data.frame() %>% dplyr::select(x, y)
    yr_bg <- get_predictor_values(yr_bg, x = "x", y = "y", yr = yr, clim_predictors = "previous_year") %>%
      dplyr::mutate(year = yr)
    })
  background <- do.call("rbind", background_byYear)
    
  PO_byYear <- purrr::map(2012:2019, function(yr){
      
    model_observations <- west_coast_observations_filtered_df %>%
      dplyr::filter(scientificName %in% species_set, year == yr) %>%
      dplyr::select(scientificName, longitude, latitude, coordinateUncertaintyInMeters, year, id) %>%
      dplyr::mutate(scientificName = as.factor(as.character(scientificName))) %>% 
      dplyr::mutate(cellID = cellFromXY(temp_summary_rasters$`2019.median`, data.frame(longitude, latitude)))
    
    model_observations <- model_observations %>% dplyr::distinct(cellID, .keep_all = TRUE)
      
    model_observations <- get_predictor_values(model_observations, x = "longitude", y = "latitude", yr = yr, clim_predictors = "previous_year")
    
    model_observations
      
    })
  
  PO_data <- do.call("rbind", PO_byYear)
  PO_data <- split(PO_data, PO_data$scientificName)
  
  ### Add lumped species at the end of the PO_data list
  PO_data <- c(PO_data,
               purrr::map(lumping_codes, function(lump){
                 do.call("rbind", PO_data[lump]) 
                 })
               )
  names(PO_data) <- gsub(" ", "_", names(PO_data))
  
   
  PA_byYear <- purrr::map(2012:2019, function(yr){
      
    survey_pa <- survey_dat %>%
      dplyr::filter(species_code %in% c(species_set, names(lumping_codes)), year == yr) %>%
      dplyr::select(species_code, longitude, latitude, year, presence_absence) %>%
      dplyr::ungroup() %>% 
      dplyr::mutate(species_code = gsub(" ", "_", species_code)) %>% 
      tidyr::spread(key = species_code, value = presence_absence)
      
    survey_pa <- get_predictor_values(survey_pa, x = "longitude", y = "latitude", yr = yr, clim_predictors = "previous_year")
      
    survey_pa
      
  })
    
  PA_data <- do.call("rbind", PA_byYear)
  
  ##### Final set of taxa to model
  species_to_model <- intersect(names(PA_data), names(PO_data)) %>% setdiff(names(background))
  
  ##### Set up model formulas
  #### Bias formula
  ### Generate spline bases for bias model variables
  num_observations_knots <- quantile(background$num_observations, (1:3)/4)
  num_observations_basis <- ns(background$num_observations, knots = num_observations_knots)
  ## Write out bias formula
  bias_formula <- formula(paste("~predict(num_observations_basis, newx = num_observations)",
                          "+ num_observers"))
  #### Sdm model formula
  ### Generate spline bases for sdm model variables
  if (isTRUE(by_year)){
    sdm_spline_vars <- c("temp_min", "salt_median", "zeta_max", "year")
  } else {
    sdm_spline_vars <- c("temp_min", "salt_median", "zeta_max")
  }
  # temp_min_knots <- quantile(background$temp_min, (1:3)/4)
  # temp_min_basis <- ns(background$temp_min, knots = temp_min_knots)
  # salt_median_knots <- quantile(background$salt_median, (1:3)/4)
  # salt_median_basis <- ns(background$salt_median, knots = salt_median_knots)
  # zeta_max_knots <- quantile(background$zeta_max, (1:3)/4)
  # zeta_max_basis <- ns(background$zeta_max, knots = zeta_max_knots)
  

  
  ### Write out sdm formula
  sdm_formula <- formula(paste("~", paste("bs(", sdm_spline_vars, ", df = 3) ",
                                        sep = "", collapse = "+"),
                             "+", paste(c("h_median", "habitat"), collapse = "+")
                             )
                         )
  
  sdm <- multispeciesPP_edit(sdm.formula = sdm_formula, 
                           bias.formula = bias_formula,
                           quadrat.size = 0.1, 
                           region.size = nrow(background)*10, # in square km
                           PA = PA_data %>% as.data.frame(), 
                           PO = PO_data, 
                           BG = background %>% as.data.frame(),
                           species = species_to_model,
                           control=list(trace = TRUE, maxit = 100)
                           )
  
  sdm_output <- list(sdm = sdm, data = list(background = background, PO_data = PO_data, PA_data = PA_data))

  return(sdm_output)

}

#' ### get_prediction_rasters
get_prediction_rasters <- function(species = "Pisaster_ochraceus", sdm = inverts_sdm){
  prediction_values <- data.frame(x = sdm$data$background$x, y = sdm$data$background$y, year = sdm$data$background$year, sdm$sdm$fit.BG[, species])
  prediction_rasters <- purrr::map(c(as.list(2012:2019), list(2012:2019)), function(yr){
    yr_prediction_values <- prediction_values %>% dplyr::filter(year %in% yr)
    yr_prediction_raster <- rasterFromXYZ(yr_prediction_values %>% dplyr::select(-year), crs = CRS("+init=epsg:4326")) 
  })
  names(prediction_rasters) <- c(as.character(2012:2019), "overall")
  return(prediction_rasters)
}

#' ### get_predicted_limits
get_predicted_limits <- function(species_predictions, min_threshold = 0.1){
  predicted_limits <- purrr::map(species_predictions, function(pred){
    line <- Line(cbind(c(-125.15, -116.65), rep(xyFromCell(pred, which(pred[] >= min_threshold)) %>% as.data.frame() %>% pull(y) %>% max(), 2)))
    sp_line <- SpatialLines(list(Lines(list(line), ID="line")))  
    sp_line
  })
  return(predicted_limits)
}

#' plot_response_curve
plot_response_curve <- function(sdm, focal_species_name, predictor = "temp_min", interval = 0.1, xlabel = "Mean temperature"){
  
  bg <- sdm$data$background
  new_var <- seq(round(min(bg[, predictor])), round(max(bg[, predictor])), .1)
  new_var_seq <- rep(new_var, times = length(levels(bg$habitat)))
  new_habitat_seq <- rep(levels(bg$habitat), length(new_var)) %>% as.factor()
  bg_means <- sdm$data$background %>% dplyr::select(-x, -y, -habitat, -predictor) %>% summarise_all(mean)
  new_bg <- data.frame(new_var_seq, new_habitat_seq) %>% data.frame(bg_means[rep(1, length(new_var_seq)), ]) 
  names(new_bg)[1:2] <- c(predictor, "habitat")
  new_bg <- rbind(data.frame(new_bg, year = 2013), 
                  data.frame(new_bg, year = 2014), 
                  data.frame(new_bg, year = 2015), 
                  data.frame(new_bg, year = 2016), 
                  data.frame(new_bg, year = 2017), 
                  data.frame(new_bg, year = 2018), 
                  data.frame(new_bg, year = 2019)
  )
  sdm_predictions_new <- predict.multispeciesPP(sdm$sdm, newdata = new_bg, type = "response") %>% exp()
  sdm_predictions_new <- data.frame(cbind(new_bg %>% dplyr::select(predictor), sdm_predictions_new[, focal_species_name]))
  names(sdm_predictions_new)[2] <- "species_density"
  mean_response <- sdm_predictions_new %>% group_by_at(vars(one_of(predictor))) %>% dplyr::summarise(species_density = mean(species_density)) %>% as.data.frame()
  p <- ggplot(mean_response, aes_string(x = predictor, y = "species_density")) +
       geom_line(lwd = 1.5) + 
       theme_bw() +
       xlab(xlabel) +
       ylab("Species density (1/10km2)")
  
  print(p)
  return(p)
  
}

#' ### get_associated_species
#' #### Function to extract the n species most associated (i.e. most often recorded together) with a focal species
get_associated_species <- function(focal_species = "Pisaster_ochraceus", ### Select a focal species
                                   num_species = 10, ### Select the total number of associated species in output
                                   association_matrix = cal_coast_species_associations ### Point to the pairwise association matrix with all species
){
  ### Extract the focal species' column from the association matrix 
  association_species_ordered <- association_matrix[focal_species, ]
  ### Order the associated species vector by decreasing association 
  association_species_ordered <- association_species_ordered[, order(association_species_ordered[focal_species, ], decreasing = TRUE)]
  ### Extract the most associated species
  most_associated_species <- association_species_ordered[1:num_species]
  ### Return vector of most associated species
  return(most_associated_species)
  
}

##### -- Generate count data -- #####
#### get_count_data()
get_count_data <- function(observations = cal_coast_observations_df,
                           focal_species_name = NULL,
                           focal_species_set = NULL,
                           start_year = 2012,
                           end_year = 2019,
                           spatial_boundary = NULL,
                           target_group_size = 50){
  
  observations <- observations %>% dplyr::filter(year >= start_year, year <= end_year)
  
  if (!is.null(spatial_boundary)) {
    #' #### Isolate coordinates and turn into spatial points
    observation_points <- observations %>% dplyr::select(longitude, latitude) %>% SpatialPoints(CRS("+init=epsg:4326"))
    #' #### Identify spatial points overlapping accessible area
    observations_over_boundary <- sp::over(observation_points, spatial_boundary)
    #' #### Filter iNat data overlapping accessible area
    observations <- observations[!is.na(observations_over_boundary), ]
  }
  
  if (!is.null(focal_species_name)){
    target_group_species <- get_associated_species(focal_species = focal_species_name, num_species = target_group_size) %>% names()
  }
  
  if (!is.null(focal_species_set)){
    target_group_species <- purrr::map(focal_species_set[[1]], function(sp){
      get_associated_species(focal_species = sp, num_species = 50) %>% names()
    }) %>% unlist() %>% unique()
  }
  
  ### Filter observations from benchmark higher taxon 
  focal_taxon_observations <- observations %>% 
    dplyr::filter(scientificName %in% gsub("_", " ", target_group_species))
  
  ### Get recording data from observations across all species
  count_data <- focal_taxon_observations %>%
    dplyr::group_by(visitID, scientificName) %>%
    dplyr::count(scientificName) %>%
    dplyr::mutate(species_name = gsub(" ", "_", scientificName)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-scientificName) %>%
    tidyr::spread(key = species_name, value = n, fill = 0) 
    
  count_data <- count_data %>% dplyr::mutate(year = unlist(lapply(strsplit(visitID, "_"), "[[", 2)))
  count_data <- count_data %>% dplyr::select(visitID, year, count_data %>% dplyr::select(-visitID, -year) %>% names())
  
  if (!is.null(focal_species_name)){
    if (!(focal_species_name %in% names(count_data))){
    count_data <- count_data %>% dplyr::mutate(focal_species_name = 0) 
    names(count_data)[ncol(count_data)] <- focal_species_name
    }
  }
  
  return(count_data)
  
}

#### calculate_reporting_rate()
calculate_yearly_reporting_rate <- function(count_data, focal_species_name = NULL, focal_species_set = NULL){
  
  detection_data <- count_data %>% dplyr::select(-visitID, -year) %>% dplyr::mutate_all(~as.numeric(. > 0))
  detection_data <- count_data %>% dplyr::select(year) %>% data.frame(detection_data)
  
  reporting_rate_byYear <- detection_data %>%
    dplyr::group_by(year) %>%
    dplyr::summarise_all(sum) %>%
    dplyr::ungroup()
  
  if (!is.null(focal_species_name)){
    reporting_rate_byYear <- reporting_rate_byYear %>% 
      dplyr::mutate(focal_species_detection = reporting_rate_byYear[focal_species_name] %>% pull())
  }
  
  if (!is.null(focal_species_set)){
    reporting_rate_byYear <- reporting_rate_byYear %>% 
      dplyr::mutate(focal_species_detection = reporting_rate_byYear[focal_species_set] %>% rowSums())
  }
  
  reporting_rate_byYear <- reporting_rate_byYear %>%
    dplyr::mutate(n_visits = reporting_rate_byYear %>% dplyr::select(-year, -focal_species_detection) %>% rowSums(),
                  reporting_rate = focal_species_detection/n_visits
                  ) %>% 
    dplyr::select(year, focal_species_detection, n_visits, reporting_rate)
  
  return(reporting_rate_byYear)
  
}

#### plot_reporting_rate()
plot_reporting_rate <- function(reporting_rate_df, 
                                species_label = "Pisaster ochraceus",
                                all_plots = TRUE
                                ){
  
  rr <- ggplot(data = reporting_rate_df, aes(x = year %>% as.numeric(), y = reporting_rate)) +
    geom_line(linetype = 1, size = 1.5) +
    theme_classic() +
    theme(panel.border = element_blank(),
          panel.background = element_blank(), # bg of the panel
          plot.background = element_blank(), # bg of the plot
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 12, angle = 90, hjust = 0.5),
          axis.title.y = element_text(size = 15),
          legend.text = element_text(size = 13),
          legend.title = element_blank(),
          legend.position = "top"
    ) +
    ylab(paste0("Reporting rate for ", species_label)) +
    xlab("")
  
  if(isTRUE(all_plots)){
    
    f <- ggplot(data = reporting_rate_df, aes(x = year %>% as.numeric(), y = focal_species_detection)) +
      geom_line(linetype = 1, size = 1.5, col = "tomato") +
      theme_classic() +
      theme(panel.border = element_blank(),
            panel.background = element_blank(), # bg of the panel
            plot.background = element_blank(), # bg of the plot
            panel.grid.major = element_blank(), # get rid of major grid
            panel.grid.minor = element_blank(), # get rid of minor grid
            axis.text.x = element_text(size = 15),
            axis.text.y = element_text(size = 12, angle = 90, hjust = 0.5),
            axis.title.y = element_text(size = 15),
            legend.text = element_text(size = 13),
            legend.title = element_blank(),
            legend.position = "top"
      ) +
      ylab(paste0("Number of observations of \n ", species_label)) +
      xlab("") +
      ylim(c(0, max(reporting_rate_df$n_visits)))
      
    print(f)
    
    v <- ggplot(data = reporting_rate_df, aes(x = year %>% as.numeric(), y = focal_species_detection)) +
      geom_line(linetype = 1, size = 1.5, col = "tomato") +
      geom_line(aes(x = year %>% as.numeric(), y = n_visits), linetype = 1, size = 1.5, col = "deepskyblue3") +
      theme_classic() +
      theme(panel.border = element_blank(),
            panel.background = element_blank(), # bg of the panel
            plot.background = element_blank(), # bg of the plot
            panel.grid.major = element_blank(), # get rid of major grid
            panel.grid.minor = element_blank(), # get rid of minor grid
            axis.text.x = element_text(size = 15),
            axis.text.y = element_text(size = 12, angle = 90, hjust = 0.5),
            axis.title.y = element_text(size = 15),
            legend.text = element_text(size = 13),
            legend.title = element_blank(),
            legend.position = "top"
      ) +
      ylab(paste0("Number of observations of \n species associated with ", species_label)) +
      xlab("") +
      ylim(c(0, max(reporting_rate_df$n_visits)))
    
    print(v)
    
    print(rr)
  } else {
    print(rr)
  }
}

get_trends_data <- function(species_name, byRegion = TRUE, observations = cal_coast_observations_df, target_group_size = 50){

  trend_overall <- observations %>% 
    get_count_data(focal_species_name = species_name, target_group_size = target_group_size) %>%
    calculate_yearly_reporting_rate(focal_species_name = species_name)
  
  if (isTRUE(byRegion)){
    
    trend_byRegion <- purrr::map(c("North", "Central", "South"), function(reg){
      
      reg_poly <- region_polygons[region_polygons$Region == reg, ]
      
      rr <- observations %>% 
        get_count_data(focal_species_name = species_name, spatial_boundary = reg_poly) %>%
        calculate_yearly_reporting_rate(focal_species_name = species_name)

    }) 
  } else trend_byRegion <- list(NULL, NULL, NULL)
  
  trends_data <- c(list(trend_overall), trend_byRegion)
  names(trends_data) <- c("Statewide", "North", "Central", "South")
  
  return(trends_data)
  
}

##### -- calculate_auc() -- #####
calculate_auc <- function(species_name, intensity_predictions, validation_occurrences, niter = 100){
  
  species_points <- validation_occurrences %>% dplyr::select(longitude, latitude) %>% SpatialPoints(CRS("+init=epsg:4326"))
  
  lambda <- rasterFromXYZ(intensity_predictions[c("x", "y", species_name)]) %>% trim()
  
  proj4string(lambda) <- CRS("+init=epsg:4326")
  
  Fl <- ecdf(lambda[])
  
  lamX <- raster::extract(lambda, species_points)
  
  auc_observed <- mean(Fl(lamX))
  
  auc_random <- NULL
  
  for (i in 1:niter){
    intensity_predictions[c("x", "y")]
    random_pts <- sample_n(intensity_predictions[c("x", "y")], length(species_points)) %>% SpatialPoints(CRS("+init=epsg:4326"))
    random_lamX <- raster::extract(lambda, random_pts)
    auc_random <- c(auc_random, mean(Fl(random_lamX)))
  }
  
  # lambda <- as.im(lambda)
  # atheo <- mean(lambda[] * Fl(lambda[]))/mean(lambda)
  
  result <- c(auc_observed = auc_observed, auc_random_mean = mean(auc_random), auc_random_lower = quantile(auc_random, .025), auc_random_upper = quantile(auc_random, .975))
  
  return(result)
  
}

get.checkerboard3 <- function (occ, env, bg.coords, aggregation.factor) 
{
  occ <- as.data.frame(occ)
  rownames(occ) <- 1:nrow(occ)
  bg.coords <- as.data.frame(bg.coords)
  rownames(bg.coords) <- 1:nrow(bg.coords)
  if (length(aggregation.factor) == 1) 
    aggregation.factor <- rep(aggregation.factor, 3)
  grid <- aggregate(env[[1]], fact = aggregation.factor[1])
  grid2 <- aggregate(grid, aggregation.factor[2])
  grid3 <- aggregate(grid2, aggregation.factor[3])
  w <- gridSample(occ, grid, n = 10000, chess = "white")
  b <- gridSample(occ, grid, n = 10000, chess = "black")
  ww <- gridSample(w, grid2, n = 10000, chess = "white")
  wb <- gridSample(w, grid2, n = 10000, chess = "black")
  bw <- gridSample(b, grid2, n = 10000, chess = "white")
  bb <- gridSample(b, grid2, n = 10000, chess = "black")
  www <- gridSample(ww, grid3, n = 10000, chess = "white")
  wwb <- gridSample(ww, grid3, n = 10000, chess = "black")
  wbw <- gridSample(wb, grid3, n = 10000, chess = "white")
  wbb <- gridSample(wb, grid3, n = 10000, chess = "black")
  bww <- gridSample(bw, grid3, n = 10000, chess = "white")
  bwb <- gridSample(bw, grid3, n = 10000, chess = "black")
  bbw <- gridSample(bb, grid3, n = 10000, chess = "white")
  bbb <- gridSample(bb, grid3, n = 10000, chess = "black")
  bgw <- gridSample(bg.coords, grid, n = 10000, chess = "white")
  bgb <- gridSample(bg.coords, grid, n = 10000, chess = "black")
  bgww <- gridSample(bgw, grid2, n = 10000, chess = "white")
  bgwb <- gridSample(bgw, grid2, n = 10000, chess = "black")
  bgbw <- gridSample(bgb, grid2, n = 10000, chess = "white")
  bgbb <- gridSample(bgb, grid2, n = 10000, chess = "black")
  bgwww <- gridSample(bgww, grid3, n = 10000, chess = "white")
  bgwwb <- gridSample(bgww, grid3, n = 10000, chess = "black")
  bgwbw <- gridSample(bgwb, grid3, n = 10000, chess = "white")
  bgwbb <- gridSample(bgwb, grid3, n = 10000, chess = "black")
  bgbww <- gridSample(bgbw, grid3, n = 10000, chess = "white")
  bgbwb <- gridSample(bgbw, grid3, n = 10000, chess = "black")
  bgbbw <- gridSample(bgbb, grid3, n = 10000, chess = "white")
  bgbbb <- gridSample(bgbb, grid3, n = 10000, chess = "black")
  r <- data.frame()
  if (nrow(www) > 0) www$grp <- 1
  r <- rbind(r, www)
  if (nrow(wwb) > 0) wwb$grp <- 2
  r <- rbind(r, wwb)
  if (nrow(wbw) > 0) wbw$grp <- 3
  r <- rbind(r, wbw)
  if (nrow(wbb) > 0) wbb$grp <- 4
  r <- rbind(r, wbb)
  if (nrow(bww) > 0) bww$grp <- 5
  r <- rbind(r, bww)
  if (nrow(bwb) > 0) bwb$grp <- 6
  r <- rbind(r, bwb)
  if (nrow(bbw) > 0) bbw$grp <- 7
  r <- rbind(r, bbw)
  if (nrow(bbb) > 0) bbb$grp <- 8
  r <- rbind(r, bbb)
  occ.grp <- r[order(as.numeric(rownames(r))), ]$grp
  bgr <- data.frame()
  if (nrow(bgwww) > 0) bgwww$grp <- 1
  bgr <- rbind(bgr, bgwww)
  if (nrow(bgwwb) > 0) bgwwb$grp <- 2
  bgr <- rbind(bgr, bgwwb)
  if (nrow(bgwbw) > 0) bgwbw$grp <- 3
  bgr <- rbind(bgr, bgwbw)
  if (nrow(bgwbb) > 0) bgwbb$grp <- 4
  bgr <- rbind(bgr, bgwbb)
  if (nrow(bgbww) > 0) bgbww$grp <- 5
  bgr <- rbind(bgr, bgbww)
  if (nrow(bgbwb) > 0) bgbwb$grp <- 6
  bgr <- rbind(bgr, bgbwb)
  if (nrow(bgbbw) > 0) bgbbw$grp <- 7
  bgr <- rbind(bgr, bgbbw)
  if (nrow(bgbbb) > 0) bgbbb$grp <- 8
  bgr <- rbind(bgr, bgbbb)
  bg.grp <- bgr[order(as.numeric(rownames(bgr))), ]$grp
  noccgrp <- length(unique(occ.grp))
  nbggrp <- length(unique(bg.grp))
  out <- list(occ.grp = occ.grp, bg.grp = bg.grp)
  return(out)
}

run_sdm_validation <- function(sdm = inverts_sdm){
  
  bg_validation_sets <- purrr::map(1:8, function(block_number){
    bg_blocks <- get.checkerboard3(occ = do.call("rbind", sdm$data$PO_data)[c("longitude", "latitude")] %>% unique(), env = temp_summary_rasters$`2019.median`, bg.coords = sdm$data$background[c("x", "y")], aggregation.factor = c(3, 3, 3))$bg.grp
    sdm$data$background %>% dplyr::mutate(validation = bg_blocks == block_number)
  }) %>% set_names(paste0("validation_set", 1:8))
  
  bg_validation_set_rasters <- purrr::map(bg_validation_sets, function(x){
    rasterFromXYZ(x[c("x", "y", "validation")], crs = CRS("+init=epsg:4326"))
  }) %>% set_names(names(bg_validation_sets))
  
  bg_validation_set_proportions <- purrr::map(bg_validation_set_rasters, function(r) length(which(r[] == 1))/length(which(!is.na(r[])))) %>% unlist()
  bg_validation_sets <- bg_validation_sets[which(bg_validation_set_proportions >= .1 & bg_validation_set_proportions <= .15)]
  bg_validation_set_rasters <- bg_validation_set_rasters[which(bg_validation_set_proportions >= .1 & bg_validation_set_proportions <= .15)]
  
  #### For presence-only data
  PO_data_validation_sets <- purrr::map(names(sdm$data$PO_data), function(species_name){
    PO_occs <- sdm$data$PO_data[[species_name]]
    PO_sets <- purrr::map(bg_validation_set_rasters, function(validation_r){
      PO_occs$validation <- as.logical(raster::extract(validation_r, PO_occs[c("longitude", "latitude")]))
      PO_occs
    }) %>% set_names(names(bg_validation_sets))
    PO_sets
  }) %>% set_names(names(sdm$data$PO_data))
  
  #### For presence-absence data
  PA_data_validation_sets <- purrr::map(bg_validation_set_rasters, function(validation_r){
    sdm$data$PA_data$validation <- as.logical(raster::extract(validation_r, sdm$data$PA_data[c("longitude", "latitude")]))
    sdm$data$PA_data
  }) %>% set_names(names(bg_validation_sets))
  
  species_group_sdm_sets <-
    purrr::map(setdiff(names(bg_validation_sets), c("validation_set6", "validation_set8")), function(validation_set){
      print(validation_set)
          multispeciesPP_edit(sdm$sdm$sdm.formula, 
                              sdm$sdm$bias.formula,
                              quadrat.size = res(temp_summary_rasters$`2019.median`)[1], 
                              region.size = nrow(bg_validation_sets[[validation_set]] %>% dplyr::filter(validation == FALSE)), # in square km
                              PA = PA_data_validation_sets[[validation_set]] %>% dplyr::filter(validation == FALSE) %>% as.data.frame(), 
                              PO = lapply(PO_data_validation_sets, function(spec) spec[[validation_set]] %>% dplyr::filter(validation == FALSE) %>% as.data.frame()), 
                              BG = bg_validation_sets[[validation_set]] %>% dplyr::filter(validation == FALSE) %>% as.data.frame(),
                              species = dimnames(sdm$sdm$fit.BG)[[2]],
                              control=list(trace = TRUE, maxit = 500) # species data and background data
          )  
    }) %>% set_names(setdiff(names(bg_validation_sets), c("validation_set6", "validation_set8")))       
  
  validation <- purrr::map(names(species_group_sdm_sets), function(validation_set){
    preds <- predict.multispeciesPP(species_group_sdm_sets[[validation_set]], newdata = bg_validation_sets[[validation_set]] %>% dplyr::filter(validation == TRUE), type = "response") %>% exp()
    preds <- data.frame(cbind(bg_validation_sets[[validation_set]] %>% dplyr::filter(validation == TRUE) %>% dplyr::select(x, y), preds))
    names(preds) <- c("x", "y", dimnames(sdm$sdm$fit.BG)[[2]])
    valid <- purrr::map(names(preds)[-c(1:2)], function(sp){
      sp_preds <- preds %>% dplyr::select(x, y, all_of(sp)) %>% dplyr::mutate(cellID = cellFromXY(bg_validation_set_rasters[[validation_set]], preds %>% dplyr::select(x, y)))
      names(sp_preds)[3] <- "predicted"
      pa <- PA_data_validation_sets[[validation_set]] %>% dplyr::filter(validation == TRUE) %>% dplyr::select(longitude, latitude, all_of(sp))
      po <- PO_data_validation_sets[[sp]][[validation_set]] %>% dplyr::filter(validation == TRUE) %>% dplyr::mutate(sp = 1) %>% dplyr::select(longitude, latitude, sp)
      names(po)[3] <- sp
      sp_obs <- rbind(pa, po) 
      sp_obs <- sp_obs %>% 
        dplyr::mutate(cellID = cellFromXY(bg_validation_set_rasters[[validation_set]], sp_obs[c("longitude", "latitude")] %>% as.data.frame())) %>% 
        dplyr::select(longitude, latitude, all_of(sp), cellID)
      names(sp_obs)[3] <- "observed"
      sp_val <- left_join(sp_obs, sp_preds, by = "cellID") %>% as.data.frame()
      if (n_distinct(sp_val$observed) == 2){
        pa_val <- sdm::evaluates(x = sp_val$observed, p = sp_val$predicted)
        pa_val <- pa_val@statistics
      } else pa_val <- NA
      if (nrow(sp_val) >= 5){
        po_val <- mean(ecdf(sp_preds$predicted)(sp_val$predicted[sp_val$observed == 1])) 
      } else po_val <- NA
      
      list(PO_validation = po_val, PA_validation = pa_val)
    }) %>% set_names(names(preds)[-c(1:2)])
    valid
  }) %>% set_names(names(species_group_sdm_sets))
  
  return(validation)
}

run_sdm_validation_byYear <- function(sdm = inverts_sdm_byYear){
  
  bg_validation_sets <- purrr::map(2012:2019, function(yr){
    sdm$data$background %>% dplyr::mutate(validation = year == yr)
  })
  names(bg_validation_sets) <- 2012:2019
  
  #### For presence-only data
  PO_data_validation_sets <- purrr::map(names(sdm$data$PO_data), function(species_name){
    PO_occs <- sdm$data$PO_data[[species_name]]
    val <- purrr::map(2012:2019, function(yr){
      PO_occs %>% dplyr::mutate(validation = year == yr)
    })
    names(val) <- 2012:2019
    val
  })
  names(PO_data_validation_sets) <- names(sdm$data$PO_data)
  
  #### For presence-absence data
  PA_data_validation_sets <- purrr::map(2012:2019, function(yr){
    sdm$data$PA_data %>% dplyr::mutate(validation = year == yr)
  })
  names(PA_data_validation_sets) <- 2012:2019
  
  species_group_sdm_sets <-
    purrr::map(as.character(2012:2019), function(validation_set){
      multispeciesPP_edit(sdm$sdm$sdm.formula, 
                          sdm$sdm$bias.formula,
                          quadrat.size = res(temp_summary_rasters$`2019.median`)[1], 
                          region.size = nrow(bg_validation_sets[[validation_set]] %>% dplyr::filter(validation == FALSE)), # in square km
                          PA = PA_data_validation_sets[[validation_set]] %>% dplyr::filter(validation == FALSE) %>% as.data.frame(), 
                          PO = lapply(PO_data_validation_sets, function(spec) spec[[validation_set]] %>% dplyr::filter(validation == FALSE) %>% as.data.frame()), 
                          BG = bg_validation_sets[[validation_set]] %>% dplyr::filter(validation == FALSE) %>% as.data.frame(),
                          species = dimnames(sdm$sdm$fit.BG)[[2]],
                          control=list(trace = TRUE, maxit = 100) # species data and background data
      )  
    })
  names(species_group_sdm_sets) <- as.character(2012:2019)
  
  validation <- purrr::map(names(species_group_sdm_sets), function(validation_set){
      preds <- predict.multispeciesPP(species_group_sdm_sets[[validation_set]], newdata = bg_validation_sets[[validation_set]] %>% dplyr::filter(validation == TRUE), type = "response") %>% exp()
      preds <- data.frame(cbind(bg_validation_sets[[validation_set]] %>% dplyr::filter(validation == TRUE) %>% dplyr::select(x, y), preds))
      names(preds) <- c("x", "y", dimnames(sdm$sdm$fit.BG)[[2]])
      valid <- purrr::map(names(preds)[-c(1:2)], function(sp){
        sp_preds <- preds %>% dplyr::select(x, y, all_of(sp)) %>% dplyr::mutate(cellID = cellFromXY(temp_summary_rasters$`2019.median`, preds %>% dplyr::select(x, y)))
        names(sp_preds)[3] <- "predicted"
        pa <- PA_data_validation_sets[[validation_set]] %>% dplyr::filter(validation == TRUE) %>% dplyr::select(longitude, latitude, all_of(sp))
        po <- PO_data_validation_sets[[sp]][[validation_set]] %>% dplyr::filter(validation == TRUE) %>% dplyr::mutate(sp = 1) %>% dplyr::select(longitude, latitude, sp)
        names(po)[3] <- sp
        sp_obs <- rbind(pa, po) 
        sp_obs <- sp_obs %>% 
          dplyr::mutate(cellID = cellFromXY(temp_summary_rasters$`2019.median`, sp_obs[c("longitude", "latitude")] %>% as.data.frame())) %>% 
          dplyr::select(longitude, latitude, all_of(sp), cellID)
        names(sp_obs)[3] <- "observed"
        sp_val <- left_join(sp_obs, sp_preds, by = "cellID") %>% as.data.frame()
        if (n_distinct(sp_val$observed) == 2){
          pa_val <- sdm::evaluates(x = sp_val$observed, p = sp_val$predicted)
          pa_val <- pa_val@statistics
        } else pa_val <- NA
        if (nrow(sp_val) >= 5){
          po_val <- mean(ecdf(sp_preds$predicted)(sp_val$predicted[sp_val$observed == 1])) 
        } else po_val <- NA
        
        list(PO_validation = po_val, PA_validation = pa_val)
      }) %>% set_names(names(preds)[-c(1:2)])
      valid
    }) %>% set_names(names(species_group_sdm_sets))

}

get_trends_validation_data2 <- function(species_set = c("Pisaster_ochraceus", 
                                                       "Lottia_gigantea", 
                                                       "Tetraclita_rubescens", 
                                                       "Mytilus_californianus", 
                                                       "Pollicipes_polymerus",
                                                       "Semibalanus_cariosus",
                                                       "Dermasterias_imbricata",
                                                       "Egregia_menziesii", 
                                                       "Pelvetiopsis_limitata",
                                                       "Sargassum_muticum",
                                                       "Silvetia_compressa"),
                                       species_groups = list("Anthopleura elegantissima/sola" = c("Anthopleura_elegantissima", "Anthopleura_sola", "Anthopleura_xanthogrammica"),
                                                             "Phyllospadix spp" = c("Phyllospadix_scouleri", "Phyllospadix_torreyi")
                                                             ),
                                       byRegion = TRUE
                                       )
{
  
  rr_trends <- c(
    purrr::map(species_set, function(sp){
    calculate_yearly_reporting_rate2(focal_species_name = sp)
  }) %>% set_names(species_set),
  purrr::map(species_groups, function(sp){
    calculate_yearly_reporting_rate2(focal_species_set = sp)
  }) %>% set_names(names(species_groups))
  )
  
  survey_trends <- survey_data_lt %>% 
    dplyr::filter(species_code %in% c(gsub("_", " ", species_set), names(species_groups))) %>%
    dplyr::select(year, species_code, total) %>% 
    dplyr::group_by(year, species_code) %>%
    dplyr::summarise(total = mean(total, na.rm = TRUE)) %>%
    dplyr::ungroup() 
  
  survey_trends <- split(survey_trends, as.factor(survey_trends$species_code))
  names(survey_trends)[names(survey_trends) %in% gsub("_", " ", species_set)] <- gsub(" ", "_", names(survey_trends)[names(survey_trends) %in% gsub("_", " ", species_set)])
  survey_trends <- survey_trends[names(rr_trends)]
  
  trends_overall <- purrr::map(names(rr_trends), function(sp){
    out <- left_join(rr_trends[[sp]], survey_trends[[sp]], by = "year") %>% 
      dplyr::mutate(species_name = sp, region = "statewide") 
    out %>% dplyr::select(species_name, region, out %>% dplyr::select(-species_name, -region, -species_code) %>% names())  
  }) %>% set_names(names(rr_trends))
  
  if (isTRUE(byRegion)){
    
    trends_byRegion <- purrr::map(c("North", "Central", "South"), function(reg){
      
      reg_poly <- region_polygons[region_polygons$Region == reg, ]
      
      rr_trends <- c(
        purrr::map(species_set, function(sp){
          calculate_yearly_reporting_rate2(focal_species_name = sp, spatial_boundary = reg_poly)
        }) %>% set_names(species_set),
        purrr::map(species_groups, function(sp){
          calculate_yearly_reporting_rate2(focal_species_set = sp, spatial_boundary = reg_poly)
        }) %>% set_names(names(species_groups))
      )
      
      survey_trends <- survey_data_lt %>% 
        dplyr::filter(species_code %in% c(gsub("_", " ", species_set), names(species_groups))) %>%
        dplyr::select(latitude, longitude, year, species_code, total) %>% 
        dplyr::filter(latitude >= ymin(reg_poly) & latitude <= ymax(reg_poly)) %>% 
        dplyr::group_by(year, species_code) %>%
        dplyr::summarise(total = mean(total, na.rm = TRUE)) %>%
        dplyr::ungroup() 
      
      survey_trends <- split(survey_trends, as.factor(survey_trends$species_code))
      names(survey_trends)[names(survey_trends) %in% gsub("_", " ", species_set)] <- gsub(" ", "_", names(survey_trends)[names(survey_trends) %in% gsub("_", " ", species_set)])
      survey_trends <- survey_trends[names(rr_trends)]
      
      trends_overall <- purrr::map(names(rr_trends), function(sp){
        out <- left_join(rr_trends[[sp]], survey_trends[[sp]], by = "year") %>% 
          dplyr::mutate(species_name = sp, region = reg) 
        out %>% dplyr::select(species_name, region, out %>% dplyr::select(-species_name, -region, -species_code) %>% names())  
      }) %>% set_names(names(rr_trends))
      
      trends_overall
      
    }) %>% set_names(c("North", "Central", "South"))
  }
  
  trends_data <- c(list(trends_overall), trends_byRegion)
  names(trends_data) <- c("Statewide", "North", "Central", "South")
  
  return(trends_data)
  
}

calculate_yearly_reporting_rate2 <- function(observations = cal_coast_observations_df,
                                             focal_species_name = NULL,
                                             focal_species_set = NULL,
                                             start_year = 2012,
                                             end_year = 2019,
                                             spatial_boundary = NULL,
                                             target_group_size = 50){
  
  
  observations <- observations %>% dplyr::filter(year >= start_year, year <= end_year)
  
  if (!is.null(spatial_boundary)) {
    #' #### Isolate coordinates and turn into spatial points
    observation_points <- observations %>% dplyr::select(longitude, latitude) %>% SpatialPoints(CRS("+init=epsg:4326"))
    #' #### Identify spatial points overlapping accessible area
    observations_over_boundary <- sp::over(observation_points, spatial_boundary)
    #' #### Filter iNat data overlapping accessible area
    observations <- observations[!is.na(observations_over_boundary), ]
  }
  
  if (!is.null(focal_species_name)){
    target_group_species <- get_associated_species(focal_species = focal_species_name, num_species = 50) %>% names()
    
    ### Filter observations from benchmark higher taxon 
    focal_taxon_observations <- observations %>% 
      dplyr::filter(scientificName %in% gsub("_", " ", target_group_species))
    
    reporting_rate_byYear <- focal_taxon_observations %>% 
      dplyr::filter(scientificName == gsub("_", " ", focal_species_name)) %>% dplyr::group_by(year) %>% dplyr::summarise(focal_species_detection = n_distinct(visitID)) %>% 
      left_join(focal_taxon_observations %>% dplyr::group_by(year) %>% dplyr::summarise(target_group_detection = n_distinct(visitID)), by = "year") %>% 
      dplyr::mutate(reporting_rate = focal_species_detection/target_group_detection)
    
  }
  
  if (!is.null(focal_species_set)){
    target_group_species <- purrr::map(focal_species_set, function(sp){
      get_associated_species(focal_species = sp, num_species = 50) %>% names()
    }) %>% unlist() %>% unique()
    
    ### Filter observations from benchmark higher taxon 
    focal_taxon_observations <- observations %>% 
      dplyr::filter(scientificName %in% gsub("_", " ", target_group_species))
    
    reporting_rate_byYear <- focal_taxon_observations %>% 
      dplyr::filter(scientificName %in% gsub("_", " ", focal_species_set)) %>% dplyr::group_by(year) %>% dplyr::summarise(focal_species_detection = n_distinct(visitID)) %>% 
      left_join(focal_taxon_observations %>% dplyr::group_by(year) %>% dplyr::summarise(target_group_detection = n_distinct(visitID)), by = "year") %>% 
      dplyr::mutate(reporting_rate = focal_species_detection/target_group_detection)
    
  }
  
  return(reporting_rate_byYear)
  
}

get_trends_validation_data1 <- function(species_set = c("Pisaster_ochraceus", 
                                                       "Lottia_gigantea", 
                                                       "Tetraclita_rubescens", 
                                                       "Mytilus_californianus", 
                                                       "Pollicipes_polymerus",
                                                       "Semibalanus_cariosus",
                                                       "Dermasterias_imbricata",
                                                       "Egregia_menziesii", 
                                                       "Pelvetiopsis_limitata",
                                                       "Sargassum_muticum",
                                                       "Silvetia_compressa"),
                                       species_groups = list("Anthopleura elegantissima/sola" = c("Anthopleura_elegantissima", "Anthopleura_sola", "Anthopleura_xanthogrammica"),
                                                             "Phyllospadix spp" = c("Phyllospadix_scouleri", "Phyllospadix_torreyi")
                                       ),
                                       byRegion = TRUE,
                                       target_group_size = 50
)
{
  
  rr_trends <- c(
    purrr::map(species_set, function(sp){
      get_count_data(focal_species_name = sp, target_group_size = target_group_size) %>%
      calculate_yearly_reporting_rate(focal_species_name = sp)
    }) %>% set_names(species_set),
    purrr::map(species_groups, function(sp){
      get_count_data(focal_species_set = sp, target_group_size = target_group_size) %>% 
      calculate_yearly_reporting_rate(focal_species_set = sp)
    }) %>% set_names(names(species_groups))
  )
  
  survey_trends <- survey_data_lt %>% 
    dplyr::filter(species_code %in% c(gsub("_", " ", species_set), names(species_groups))) %>%
    dplyr::select(year, species_code, total) %>% 
    dplyr::group_by(year, species_code) %>%
    dplyr::summarise(total = mean(total, na.rm = TRUE)) %>%
    dplyr::ungroup() %>% 
    dplyr::mutate(year = as.character(year))
  
  survey_trends <- split(survey_trends, as.factor(survey_trends$species_code))
  names(survey_trends)[names(survey_trends) %in% gsub("_", " ", species_set)] <- gsub(" ", "_", names(survey_trends)[names(survey_trends) %in% gsub("_", " ", species_set)])
  survey_trends <- survey_trends[names(rr_trends)]
  
  trends_overall <- purrr::map(names(rr_trends), function(sp){
    out <- left_join(rr_trends[[sp]], survey_trends[[sp]], by = "year") %>% 
      dplyr::mutate(species_name = sp, region = "statewide") 
    out %>% dplyr::select(species_name, region, out %>% dplyr::select(-species_name, -region, -species_code) %>% names())  
  }) %>% set_names(names(rr_trends))
  
  if (isTRUE(byRegion)){
    
    trends_byRegion <- purrr::map(c("North", "Central", "South"), function(reg){
      
      reg_poly <- region_polygons[region_polygons$Region == reg, ]
      
      rr_trends <- c(
        purrr::map(species_set, function(sp){
          get_count_data(focal_species_name = sp, spatial_boundary = reg_poly, target_group_size = target_group_size) %>%
          calculate_yearly_reporting_rate(focal_species_name = sp)
        }) %>% set_names(species_set),
        purrr::map(species_groups, function(sp){
          get_count_data(focal_species_set = sp, spatial_boundary = reg_poly, target_group_size = target_group_size) %>%
          calculate_yearly_reporting_rate(focal_species_set = sp)
        }) %>% set_names(names(species_groups))
      )
      
      survey_trends <- survey_data_lt %>% 
        dplyr::filter(species_code %in% c(gsub("_", " ", species_set), names(species_groups))) %>%
        dplyr::select(latitude, longitude, year, species_code, total) %>% 
        dplyr::filter(latitude >= ymin(reg_poly) & latitude <= ymax(reg_poly)) %>% 
        dplyr::group_by(year, species_code) %>%
        dplyr::summarise(total = mean(total, na.rm = TRUE)) %>%
        dplyr::ungroup() %>% 
        dplyr::mutate(year = as.character(year))
      
      survey_trends <- split(survey_trends, as.factor(survey_trends$species_code))
      names(survey_trends)[names(survey_trends) %in% gsub("_", " ", species_set)] <- gsub(" ", "_", names(survey_trends)[names(survey_trends) %in% gsub("_", " ", species_set)])
      survey_trends <- survey_trends[names(rr_trends)]
      
      trends_overall <- purrr::map(names(rr_trends), function(sp){
        out <- left_join(rr_trends[[sp]], survey_trends[[sp]], by = "year") %>% 
          dplyr::mutate(species_name = sp, region = reg) 
        out %>% dplyr::select(species_name, region, out %>% dplyr::select(-species_name, -region, -species_code) %>% names())  
      }) %>% set_names(names(rr_trends))
      
      trends_overall
      
    }) %>% set_names(c("North", "Central", "South"))
  }
  
  trends_data <- c(list(trends_overall), trends_byRegion)
  names(trends_data) <- c("Statewide", "North", "Central", "South")
  
  return(trends_data)
  
}

plot_trends_data <- function(species_name = "Pisaster_ochraceus", region = "Statewide", validation_data = species_trends_validation_data){
  
  plot_dat <- validation_data[[region]][[species_name]]
  scaling_parameter <- max(plot_dat$total, na.rm = TRUE)/max(plot_dat$reporting_rate, na.rm = TRUE)
  scaling_parameter <- ifelse(!is.na(scaling_parameter), scaling_parameter, 1)
  
  if(nrow(plot_dat %>% dplyr::filter(complete.cases(.))) >= 3){
    trend_cor <- cor.test(plot_dat$total, plot_dat$reporting_rate)
  } else {
    trend_cor <- list(estimate = NA, p.value = NA) 
  }
  
  p <- plot_dat %>% ggplot(aes(x = year %>% as.integer())) + 
    geom_line(aes(y = reporting_rate, colour = "iNaturalist"), size = 1.5) +
    geom_line(aes(y = total/(scaling_parameter), colour = "MARINe"), size = 1.5) +
    scale_y_continuous(sec.axis = sec_axis(~.*scaling_parameter, name = "Mean count/percentage cover (MARINe)")) +
    scale_colour_manual(values = c("deepskyblue3", "tomato")) + 
    labs(y = "Reporting rate (iNaturalist)", x = "") +
    theme_classic() +
    theme(panel.border = element_blank(),
          panel.background = element_blank(), # bg of the panel
          plot.background = element_blank(), # bg of the plot
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          axis.text.x = element_text(size = 15, family = "Noto Mono"),
          axis.text.y = element_text(size = 12, angle = 90, hjust = 0.5, family = "Noto Mono"),
          axis.title.y = element_text(size = 15, family = "Noto Mono"),
          legend.text = element_text(size = 13, family = "Noto Mono"),
          legend.title = element_blank(),
          legend.position = "bottom"
    ) + 
    annotate("text", family = "Noto Mono", x = 2015, y = 0.001, label = paste("r = ", ifelse(is.numeric(trend_cor$estimate), round(trend_cor$estimate, 3), NA), ", p = ", ifelse(is.numeric(trend_cor$p.value), round(trend_cor$p.value, 3), NA), sep = ""), size = 5)
    #ggtitle(species_name)
  
  p
}
