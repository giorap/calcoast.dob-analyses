#' ---
#' title: DOB Cal Coast <br> <strong> Explore Species - Process Outputs
#' ---
#' 
#' ## Output useful objects
#' ### Species distribution models
#' #### Overall
#' #### Inverts
saveRDS(inverts_sdm, "output/inverts_sdm.rds")
#' #### Algae
saveRDS(algae_sdm, "output/algae_sdm.rds")
#' ### Species distribution model predictions
saveRDS(sdm_predictions, "output/sdm_predictions.rds")
#' #### By year
#' #### Inverts
saveRDS(inverts_sdm_byYear, "output/inverts_sdm_byYear.rds")
#' #### Algae
saveRDS(algae_sdm_byYear, "output/algae_sdm_byYear.rds")
#' ### Species distribution model predictions
saveRDS(sdm_predictions_byYear, "output/sdm_predictions_byYear.rds")
#'
#' ### Temporal trends
saveRDS(species_trends, "output/species_trends.rds")
#'
#' ### Physical summaries
saveRDS(temp_means_byRegion, "output/temp_means_byRegion.rds")
saveRDS(salt_means_byRegion, "output/salt_means_byRegion.rds")
saveRDS(zeta_means_byRegion, "output/zeta_means_byRegion.rds")

