#' ---
#' title: DOB Cal Coast <br> <strong> Explore Species - Species Distribution Models
#' ---
#' Overall distribution
#' ## Run species distribution models
#' ### For invertebrates
inverts_sdm <- run_sdm(species_set = c("Pisaster ochraceus", 
                                              "Lottia gigantea", 
                                              "Tetraclita rubescens", 
                                              "Mytilus californianus", 
                                              "Pollicipes polymerus",
                                              "Semibalanus cariosus",
                                              "Anthopleura elegantissima",
                                              "Anthopleura sola",
                                              "Anthopleura xanthogrammica",
                                              "Chthamalus dalli", 
                                              "Chthamalus fissus",
                                              "Balanus glandula"
                                              ),
                              lumping_codes = list("Anthopleura elegantissima/sola" = c("Anthopleura elegantissima", "Anthopleura sola", "Anthopleura xanthogrammica"),
                                                   "Chthamalus spp/Balanus glandula" = c("Chthamalus dalli", "Chthamalus fissus", "Balanus glandula")
                                                   ),
                              survey_dat = survey_data_lt_presabs,
                              by_year = FALSE
                              )
#' ### For kelp and seaweed
algae_sdm <- run_sdm(species_set = c("Egregia menziesii", 
                                     "Eisenia arborea", 
                                     "Pelvetiopsis limitata",
                                     "Phyllospadix scouleri",
                                     "Phyllospadix torreyi",
                                     "Sargassum muticum",
                                     "Silvetia compressa"
                                     ),
                     lumping_codes = list("Phyllospadix spp" = c("Phyllospadix scouleri",
                                            "Phyllospadix torreyi")
                                          ),
                     survey_dat = survey_data_lt_presabs,
                     by_year = FALSE
                     )
#' ## Generate predictions
sdm_predictions <- purrr::map(c(dimnames(inverts_sdm$sdm$species.coef)[[2]], 
                                       dimnames(algae_sdm$sdm$species.coef)[[2]]),
                                     function(spec){
                                       if (spec %in% dimnames(inverts_sdm$sdm$species.coef)[[2]]){
                                         focal_sdm <- inverts_sdm
                                       }
                                       if (spec %in% dimnames(algae_sdm$sdm$species.coef)[[2]]){
                                         focal_sdm <- algae_sdm
                                       } 
                                       get_prediction_rasters(species = spec, sdm = focal_sdm)
                                     }
)
names(sdm_predictions) <- c(dimnames(inverts_sdm$sdm$species.coef)[[2]], dimnames(algae_sdm$sdm$species.coef)[[2]])
#' ## Validate species distribution models
inverts_sdm_validation <- run_sdm_validation(sdm = inverts_sdm)
algae_sdm_validation <- run_sdm_validation(sdm = algae_sdm)
#' 
#' 
#' ## Yearly change
#' ## Run species distribution models
#' ### For invertebrates
inverts_sdm_byYear <- run_sdm(species_set = c("Pisaster ochraceus", 
                                              "Lottia gigantea", 
                                              "Tetraclita rubescens", 
                                              "Mytilus californianus", 
                                              "Pollicipes polymerus",
                                              "Semibalanus cariosus",
                                              "Anthopleura elegantissima",
                                              "Anthopleura sola",
                                              "Anthopleura xanthogrammica",
                                              "Chthamalus dalli", 
                                              "Chthamalus fissus",
                                              "Balanus glandula"
                                              ),
                              lumping_codes = list("Anthopleura elegantissima/sola" = c("Anthopleura elegantissima", "Anthopleura sola", "Anthopleura xanthogrammica"),
                                                   "Chthamalus spp/Balanus glandula" = c("Chthamalus dalli", "Chthamalus fissus", "Balanus glandula")
                                                   ),
                              survey_dat = survey_data_lt_presabs,
                              by_year = TRUE
                              )
#' ### For kelp and seaweed
algae_sdm_byYear <- run_sdm(species_set = c("Egregia menziesii", 
                                            "Eisenia arborea", 
                                            "Pelvetiopsis limitata",
                                            "Phyllospadix scouleri",
                                            "Phyllospadix torreyi",
                                            "Sargassum muticum",
                                            "Silvetia compressa"
                                            ),
                            lumping_codes = list("Phyllospadix spp" = c("Phyllospadix scouleri",
                                                                        "Phyllospadix torreyi")
                                                 ),
                            survey_dat = survey_data_lt_presabs,
                            by_year = TRUE
                            )
#' ## Generate predictions
sdm_predictions_byYear <- purrr::map(c(dimnames(inverts_sdm_byYear$sdm$species.coef)[[2]], 
                                       dimnames(algae_sdm_byYear$sdm$species.coef)[[2]]),
                                     function(spec){
                                       if (spec %in% dimnames(inverts_sdm_byYear$sdm$species.coef)[[2]]){
                                           focal_sdm <- inverts_sdm_byYear
                                           }
                                       if (spec %in% dimnames(algae_sdm_byYear$sdm$species.coef)[[2]]){
                                         focal_sdm <- algae_sdm_byYear
                                         } 
                                       get_prediction_rasters(species = spec, sdm = focal_sdm)
                                       }
                                     )
names(sdm_predictions_byYear) <- c(dimnames(inverts_sdm_byYear$sdm$species.coef)[[2]], dimnames(algae_sdm_byYear$sdm$species.coef)[[2]])
#' ## Validate species distribution models
inverts_sdm_validation_byYear <- run_sdm_validation_byYear(sdm = inverts_sdm_byYear)
algae_sdm_validation_byYear <- run_sdm_validation_byYear(sdm = algae_sdm_byYear)

#' ### Summarise validation results
inverts_sdm_validation_df <- data.frame(po_AUC = purrr::map(inverts_sdm_validation, function(x) purrr::map(x, "PO_validation")) %>% unlist(),
                                        pa_AUC = purrr::map(inverts_sdm_validation, function(x) purrr::map(x, function(y) y$PA_validation[2])) %>% unlist(),
                                        species = rep(names(inverts_sdm_validation[[1]]), times = 3),
                                        set = rep(1:3, each = length(names(inverts_sdm_validation[[1]])))
) %>% group_by(species) %>% dplyr::summarise(mean_po_AUC = mean(po_AUC, na.rm = TRUE), mean_pa_AUC = mean(pa_AUC, na.rm = TRUE))

algae_sdm_validation_df <- data.frame(po_AUC = purrr::map(algae_sdm_validation, function(x) purrr::map(x, "PO_validation")) %>% unlist(),
                                        pa_AUC = purrr::map(algae_sdm_validation, function(x) purrr::map(x, function(y) y$PA_validation[2])) %>% unlist(),
                                        species = rep(names(algae_sdm_validation[[1]]), times = 3),
                                        set = rep(1:3, each = length(names(algae_sdm_validation[[1]])))
) %>% group_by(species) %>% dplyr::summarise(mean_po_AUC = mean(po_AUC, na.rm = TRUE), mean_pa_AUC = mean(pa_AUC, na.rm = TRUE))

inverts_sdm_validation_time_df <- data.frame(po_AUC = purrr::map(inverts_sdm_validation_byYear, function(x) purrr::map(x, "PO_validation")) %>% unlist(),
                                        pa_AUC = purrr::map(inverts_sdm_validation_byYear, function(x) purrr::map(x, function(y) y$PA_validation[2])) %>% unlist(),
                                        species = rep(names(inverts_sdm_validation_byYear[[1]]), times = 8),
                                        set = rep(1:8, each = length(names(inverts_sdm_validation_byYear[[1]])))
) %>% group_by(species) %>% dplyr::summarise(mean_po_AUC = mean(po_AUC, na.rm = TRUE), mean_pa_AUC = mean(pa_AUC, na.rm = TRUE))

algae_sdm_validation_time_df <- data.frame(po_AUC = purrr::map(algae_sdm_validation_byYear, function(x) purrr::map(x, "PO_validation")) %>% unlist(),
                                      pa_AUC = purrr::map(algae_sdm_validation_byYear, function(x) purrr::map(x, function(y) y$PA_validation[2])) %>% unlist(),
                                      species = rep(names(algae_sdm_validation_byYear[[1]]), times = 8),
                                      set = rep(1:8, each = length(names(algae_sdm_validation_byYear[[1]])))
) %>% group_by(species) %>% dplyr::summarise(mean_po_AUC = mean(po_AUC, na.rm = TRUE), mean_pa_AUC = mean(pa_AUC, na.rm = TRUE))
