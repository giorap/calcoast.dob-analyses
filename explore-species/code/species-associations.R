#' ---
#' title: DOB Cal Coast <br> <strong> Explore Species - Species Associations
#' ---
#' 
#' ## Identify species associations
#' ### Identify species with the most observations across the accessible area
#' #### Identify "common" species - species with at least 100 observations across all visits where more than two species were observed
cal_coast_species_association_matrix <- get_species_association_matrix()
#' ### Identify the degree of association between pairs of species, that is the proportion of visits where both species in a pair have been recorded together
#' #### This could be done in a number of ways. For consistency with the ecological literature and for simplicity, we use the Jaccard index, which is appropriate for binary data.
#' #### Calculate the Jaccard index for each pair of species
cal_coast_species_associations <- vegdist(cal_coast_species_association_matrix, method = "jaccard") %>% as.matrix() %>% as.data.frame()
#' #### Since vegdist calculates the Jaccard index as a dissimilarity, we have to subtract all calculated values from 1 in order to obtain associations 
cal_coast_species_associations <- 1 - cal_coast_species_associations
