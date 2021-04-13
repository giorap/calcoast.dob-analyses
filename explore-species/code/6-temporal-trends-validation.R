#' ---
#' title: DOB Cal Coast <br> <strong> Explore Species - Validate Trends
#' ---
#'
#' ### Bind iNaturalist- and MARINe-derived trends for a validation set of species
species_trends_validation_data1 <- get_trends_validation_data1(target_group_size = 50)

#' ### Calculate spearman's rank correlations between iNaturalist- and MARINe-derived trends for a validation set of species
trend_cor1 <- purrr::map(species_trends_validation_data1, function(x){
  purrr::map(x, function(y){
    ifelse(nrow(y) > 0, cor(y$reporting_rate, y$total, method = "spearman", use = "complete.obs"), NA)
    }) 
}) %>% unlist()

trend_cor20 <- purrr::map(species_trends_validation_data_20, function(x){
  purrr::map(x, function(y){
    ifelse(nrow(y) > 0, cor(y$reporting_rate, y$total, method = "spearman", use = "complete.obs"), NA)
  }) 
}) %>% unlist()
#' ### Calculate predictors of correlation between iNaturalist- and MARINe-derived trends
#' #### First, identify the IDs of high effort visits (i.e. visits with more than 10 observations made)
high_effort_visits <- cal_coast_observations %>% group_by(visitID) %>% dplyr::summarize(count = n()) %>% dplyr::filter(count >= 10) %>% pull(visitID)
#' #### Calculate predictors
species_trends_validation_df1 <- data.frame(trend_cor = trend_cor1 %>% unlist(), 
                                           region = rep(c("Statewide", "North", "Central", "South"), each = length(species_trends_validation_data1$Statewide)),
                                           species = rep(names(species_trends_validation_data1$Statewide), times = 4),
                                           number_visits_where_focal_species_detected = purrr::map(species_trends_validation_data1, function(x){
                                             purrr::map(x, function(y){
                                               sum(y$focal_species_detection, na.rm = TRUE)
                                             }) %>% unlist()
                                           }) %>% unlist(),
                                           number_relevant_visits = purrr::map(species_trends_validation_data1, function(x){
                                             purrr::map(x, function(y){
                                               sum(y$n_visits, na.rm = TRUE)
                                             }) %>% unlist()
                                           }) %>% unlist(),
                                           # number_observers = cal_coast_observations %>% 
                                           #   dplyr::filter(scientificName %in% gsub("_", " ", names(species_trends_validation_data$Statewide))) %>% 
                                           #   dplyr::group_by(scientificName, region_ID) %>% 
                                           #   dplyr::summarize(n = n_distinct(recordedBy)) %>% 
                                           #   rbind(cal_coast_observations %>% 
                                           #           dplyr::filter(scientificName %in% gsub("_", " ", names(species_trends_validation_data$Statewide))) %>% 
                                           #           dplyr::group_by(scientificName) %>% 
                                           #           dplyr::summarize(n = n_distinct(recordedBy)) %>% 
                                           #           dplyr::mutate(region_ID = "Statewide")) %>% 
                                           #   dplyr::ungroup() %>% 
                                           #   as.data.frame() %>% 
                                           #   add_row(scientificName = c("Pelvetiopsis limitata", "Semibalanus cariosus"), region_ID = c("South", "South"), n = 0) %>% 
                                           #   dplyr::arrange(factor(scientificName, levels = gsub("_", " ", species_set)), factor(region_ID, levels = c("Statewide", "North", "Central", "South"))) %>% 
                                           #   pull(n),
                                           # number_high_effort_visits = cal_coast_observations %>% 
                                           #   dplyr::filter(scientificName %in% gsub("_", " ", names(species_trends_validation_data$Statewide)), visitID %in% high_effort_visits) %>% 
                                           #   distinct(.keep_all = TRUE) %>% 
                                           #   dplyr::count(scientificName, region_ID, .drop = FALSE) %>% 
                                           #   rbind(cal_coast_observations %>% 
                                           #           dplyr::filter(scientificName %in% gsub("_", " ", names(species_trends_validation_data$Statewide)), visitID %in% high_effort_visits) %>% 
                                           #           distinct(.keep_all = TRUE) %>% 
                                           #           dplyr::count(scientificName, .drop = FALSE) %>% 
                                           #           dplyr::mutate(region_ID = "Statewide")) %>% 
                                           #   dplyr::ungroup() %>% 
                                           #   as.data.frame() %>% 
                                           #   add_row(scientificName = c("Dermasterias imbricata", "Pelvetiopsis limitata", "Semibalanus cariosus"), region_ID = c("South", "South", "South"), n = 0) %>% 
                                           #   dplyr::arrange(factor(scientificName, levels = gsub("_", " ", species_set)), factor(region_ID, levels = c("Statewide", "North", "Central", "South"))) %>% 
                                           #   pull(n),
                                           largest_population_fluctuation = purrr::map(species_trends_validation_data1, function(x){
                                             purrr::map(x, function(y){
                                               y <- y %>% dplyr::mutate(total = ifelse(total > 0, total, 0.000001))
                                               out <- purrr::map(1:(length(y$total)-1), function(z){
                                                 ((y$total[z + 1] - y$total[z])/y$total[z]) * 100
                                               }) %>% unlist()
                                               max(abs(out), na.rm = TRUE)
                                             }) 
                                           }) %>% unlist(),
                                           mean_population_size = purrr::map(gsub("_", " ", names(species_trends_validation_data1$Statewide)), function(sp){
                                             survey_dat <- rbind(survey_data_lt_plots, survey_data_lt_counts %>% dplyr::filter(species_code %in% setdiff(species_code, survey_data_lt_plots %>% distinct(species_code) %>% pull())))
                                             survey_dat %>% dplyr::filter(species_code == sp) %>% group_by(species_code) %>% dplyr::summarise(mean = mean(total[which(total > 0)], na.rm = TRUE)) %>% ungroup() %>% pull(mean)
                                           }) %>% unlist(),
                                           largest_reporting_fluctuation = purrr::map(species_trends_validation_data1, function(x){
                                             purrr::map(x, function(y){
                                               y <- y %>% dplyr::mutate(reporting_rate = ifelse(reporting_rate > 0, reporting_rate, 0.000001))
                                               out <- purrr::map(1:(length(y$reporting_rate)-1), function(z){
                                                 ((y$reporting_rate[z + 1] - y$reporting_rate[z])/y$reporting_rate[z]) * 100
                                               }) %>% unlist()
                                               max(abs(out), na.rm = TRUE)
                                             }) 
                                           }) %>% unlist(),
                                           mean_reporting_rate = purrr::map(species_trends_validation_data1, function(x){
                                             purrr::map(x, function(y){
                                               mean(y$reporting_rate, na.rm = TRUE)
                                             }) %>% unlist()
                                           }) %>% unlist()
                                           
)

species_trends_validation_df2 <- data.frame(trend_cor = trend_cor2 %>% unlist(), 
                                            region = rep(c("Statewide", "North", "Central", "South"), each = length(species_trends_validation_data2$Statewide)),
                                            species = rep(names(species_trends_validation_data2$Statewide), times = 4),
                                            number_visits_where_focal_species_detected = purrr::map(species_trends_validation_data2, function(x){
                                              purrr::map(x, function(y){
                                                sum(y$focal_species_detection, na.rm = TRUE)
                                              }) %>% unlist()
                                            }) %>% unlist(),
                                            number_relevant_visits = purrr::map(species_trends_validation_data2, function(x){
                                              purrr::map(x, function(y){
                                                sum(y$target_group_detection, na.rm = TRUE)
                                              }) %>% unlist()
                                            }) %>% unlist(),
                                            largest_population_fluctuation = purrr::map(species_trends_validation_data2, function(x){
                                              purrr::map(x, function(y){
                                                y <- y %>% dplyr::mutate(total = ifelse(total > 0, total, 0.000001))
                                                out <- purrr::map(1:(length(y$total)-1), function(z){
                                                  ((y$total[z + 1] - y$total[z])/y$total[z]) * 100
                                                }) %>% unlist()
                                                max(abs(out), na.rm = TRUE)
                                              }) 
                                            }) %>% unlist(),
                                            median_population_size = purrr::map(species_trends_validation_data2, function(x){
                                              purrr::map(x, function(y){
                                                median(y$total, na.rm = TRUE)
                                              }) %>% unlist()
                                            }) %>% unlist(),
                                            largest_reporting_fluctuation = purrr::map(species_trends_validation_data2, function(x){
                                              purrr::map(x, function(y){
                                                y <- y %>% dplyr::mutate(reporting_rate = ifelse(reporting_rate > 0, reporting_rate, 0.000001))
                                                out <- purrr::map(1:(length(y$reporting_rate)-1), function(z){
                                                  ((y$reporting_rate[z + 1] - y$reporting_rate[z])/y$reporting_rate[z]) * 100
                                                }) %>% unlist()
                                                max(abs(out), na.rm = TRUE)
                                              }) 
                                            }) %>% unlist(),
                                            mean_reporting_rate = purrr::map(species_trends_validation_data2, function(x){
                                              purrr::map(x, function(y){
                                                mean(y$reporting_rate, na.rm = TRUE)
                                              }) %>% unlist()
                                            }) %>% unlist()
                                            
)
#' ### Liner mixed effects model of preditors of correlation between iNaturalist- and MARINe-derived trends
#' #### Scale predictor variables
species_trends_validation_model1 <- purrr::map(unique(species_trends_validation_df1$region), function(reg){
  species_trends_validation_df <- species_trends_validation_df1 %>% dplyr::filter(species != "Dermasterias_imbricata", region == reg) %>% 
    dplyr::filter(complete.cases(.)) %>%
    dplyr::mutate_at(vars(-trend_cor, -region, -species), scale)
  species_trends_validation_lm <- lm(trend_cor ~ mean_population_size + number_visits_where_focal_species_detected + number_relevant_visits + mean_reporting_rate, data = species_trends_validation_df)
  options(na.action=na.fail)
  species_trends_validation_lm_dredge <- dredge(species_trends_validation_lm, rank = "AIC")
  species_trends_validation_lm_models <- get.models(species_trends_validation_lm_dredge, subset = which(weight/weight[1] >= 0.05))
  if (length(species_trends_validation_lm_models) > 1) model_averaging <- model.avg(species_trends_validation_lm_models, rank = "AIC") else model_averaging <- species_trends_validation_lm_models[[1]]
  list(full_model = species_trends_validation_lm, model_selection = species_trends_validation_lm_dredge, models = species_trends_validation_lm_models, model_averaging = model_averaging)
}) %>% set_names(unique(species_trends_validation_df1$region))

ggplot(species_trends_validation_df1 %>% dplyr::filter(species != "Dermasterias_imbricata"), aes(x = mean_population_size, y = trend_cor, color = region)) +
  geom_point() +
  # geom_smooth(se = FALSE) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.position = c(0.8, 0.8),
        legend.title = element_blank()) +
  ylab("iNaturalist-survey trend correlation ") +
  xlab("Mean population size")



species_trends_validation_model2 <- purrr::map(unique(species_trends_validation_df2$region), function(reg){
  species_trends_validation_df <- species_trends_validation_df2 %>% dplyr::filter(region == reg) %>% 
    dplyr::filter(complete.cases(.)) %>% 
    dplyr::mutate_at(vars(-trend_cor, -region, -species), scale)
  species_trends_validation_lm <- lm(trend_cor ~ mean_population_size + largest_population_fluctuation + number_visits_where_focal_species_detected + number_relevant_visits + largest_reporting_fluctuation + mean_reporting_rate, data = species_trends_validation_df)
  options(na.action=na.fail)
  species_trends_validation_lm_dredge <- dredge(species_trends_validation_lm, rank = "AIC")
  species_trends_validation_lm_models <- get.models(species_trends_validation_lm_dredge, subset = which(weight/weight[1] >= 0.05))
  if (length(species_trends_validation_lm_models) > 1) model_averaging <- model.avg(species_trends_validation_lm_models, rank = "AIC") else model_averaging <- species_trends_validation_lm_models[[1]]
  list(full_model = species_trends_validation_lm, model_selection = species_trends_validation_lm_dredge, models = species_trends_validation_lm_models, model_averaging = model_averaging)
}) %>% set_names(unique(species_trends_validation_df2$region))

purrr::map(species_trends_validation_model2, function(reg) reg$model_averaging$coefficients[2, ])
