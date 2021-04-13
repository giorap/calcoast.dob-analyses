#' ---
#' title: DOB Cal Coast <br> <strong> Explore Species - Temporal Trends
#' ---
#' ## Extract temporal trends
#' ### All observations
species_trends_20 <- purrr::map(names(cal_coast_species_associations), get_trends_data, byRegion = TRUE, observations = cal_coast_observations_df, target_group_size = 20) %>% set_names(names(cal_coast_species_associations))
#' 
#' ### MPA observations only
sp_names <- cal_coast_observations %>% as.data.frame() %>% group_by(within_MPA) %>% dplyr::count(scientificName) %>% ungroup() %>% dplyr::filter(n >= 50) %>% dplyr::count(scientificName) %>% dplyr::filter(n == 2) %>% pull(scientificName)
sp_names <- gsub(" ", "_", sp_names)
species_trends_insideMPAs <- purrr::map(sp_names, get_trends_data, byRegion = TRUE, observations = cal_coast_observations_df %>% dplyr::filter(within_MPA == "within_MPA")) %>% set_names(sp_names)
#' ### Outside MPA observations only
species_trends_outsideMPAs <- purrr::map(sp_names, get_trends_data, byRegion = TRUE, observations = cal_coast_observations_df %>% dplyr::filter(within_MPA == "outside_MPA")) %>% set_names(sp_names)

species_trends_insideMPAs_df <- purrr::map(names(species_trends_insideMPAs), function(sp){
  trend <- species_trends_insideMPAs[[sp]]$Central %>% 
    dplyr::select(year, reporting_rate)
  names(trend)[2] <- sp
  trend
}) %>% join_all(by = "year") %>% 
  dplyr::mutate(overall_trend = apply(.[, -1], 1, median), within_MPA = "within_MPA")

species_trends_outsideMPAs_df <- purrr::map(names(species_trends_outsideMPAs), function(sp){
  trend <- species_trends_outsideMPAs[[sp]]$Central %>% 
    dplyr::select(year, reporting_rate)
  names(trend)[2] <- sp
  trend
}) %>% join_all(by = "year") %>% 
  dplyr::mutate(overall_trend = apply(.[, -1], 1, median), within_MPA = "outside_MPA")

species_trends_inside_vs_outside_MPAs_df <- rbind(species_trends_insideMPAs_df %>% dplyr::select(year, overall_trend, within_MPA),
                           species_trends_outsideMPAs_df %>% dplyr::select(year, overall_trend, within_MPA)
)

species_trends_inside_vs_outside_MPAs_df %>% ggplot(aes(x = year %>% as.integer(), y = overall_trend, color = within_MPA)) + 
  geom_line(size = 0.7) +
  geom_point(size = 1.3) +
  scale_y_continuous(labels = scales::percent_format(accuracy = .01))

