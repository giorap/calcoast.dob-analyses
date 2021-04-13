#' ---
#' title: DOB Cal Coast <br> <strong> Explore Places - Temporal Beta Diversity
#' ---
#' ### Calculate beta diversity across time periods using the Morisita-Horn index
#' #### All species
#' #### MPAs
temporal_stability_byMPA <- estimate_temporal_stability(boundary_ID = "MPA_ID", temporal_resolution = "year", num_iterations = 100) %>%
  clean_temporal_stability()
#' #### Watersheds
temporal_stability_byWatershed <- estimate_temporal_stability(boundary_ID = "watershed_ID", temporal_resolution = "year", num_iterations = 100) %>%
  clean_temporal_stability()
#' #### Counties
temporal_stability_byCounty <- estimate_temporal_stability(boundary_ID = "county_ID", temporal_resolution = "year", num_iterations = 100) %>% 
  clean_temporal_stability()
#' #### Hexagons
temporal_stability_byHexagon <- estimate_temporal_stability(boundary_ID = "hexagon_ID", temporal_resolution = "year", num_iterations = 100) %>% 
  clean_temporal_stability()
#' #### Regions
temporal_stability_byRegion <- estimate_temporal_stability(boundary_ID = "region_ID", temporal_resolution = "year", num_iterations = 100) %>% 
  clean_temporal_stability()
#' #### Within vs Outside MPA
temporal_stability_within_vs_outside_MPAs <- estimate_temporal_stability(boundary_ID = "within_MPA_byRegion", temporal_resolution = "year", num_iterations = 100) %>% 
  clean_temporal_stability()

change_data <- temporal_stability_within_vs_outside_MPAs
plot_dat <- change_data$place_stability %>%
  dplyr::filter(time_period >= 2013 & time_period < 2020, place %in% c("within_MPA_North", "within_MPA_Central", "within_MPA_South")) %>% 
  dplyr::select(time_period, place, upper, mean, lower)

plot_dat[c("upper", "mean", "lower")] <- round(-1 * (100 * (plot_dat[c("upper", "mean", "lower")] - change_data$global_stability["mean"])/change_data$global_stability["mean"]), 6)
change_data$global_stability <- round(-1 * (100 * (change_data$global_stability - change_data$global_stability["mean"])/change_data$global_stability["mean"]), 6)

plot_dat %>%
  ggplot(aes(x = time_period, y = mean, color = place)) +
  geom_line(size = 1) +
  geom_point(size = 1.3) +
  geom_hline(yintercept = change_data$global_stability["mean"], linetype = "dashed") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.position = c(0.8, 0.8),
        legend.title = element_blank()) +
  scale_x_continuous(breaks = plot_dat$time_period) +
  ylab("% change") +
  xlab("")

