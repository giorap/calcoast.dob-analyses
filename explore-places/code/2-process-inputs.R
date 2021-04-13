#' ---
#' title: DOB Cal Coast <br> <strong> Explore Places - Process Inputs
#' ---
#' #### This code loads and processes inputs for place-based analyses
#' ### Generate species count data over boundaries
#' #### Overall
cal_coast_count_data_byMPA <- cal_coast_observations %>% as.data.frame() %>% get_count_data(group_field = "MPA_ID") 
cal_coast_count_data_byWatershed <- cal_coast_observations %>% as.data.frame() %>% get_count_data(group_field = "watershed_ID") 
cal_coast_count_data_byCounty <- cal_coast_observations %>% as.data.frame() %>% get_count_data(group_field = "county_ID") 
cal_coast_count_data_byHexagon <- cal_coast_observations %>% as.data.frame() %>% get_count_data(group_field = "hexagon_ID")
#' #### Over time
#' #### By year
cal_coast_count_data_byMPA_byYear <- cal_coast_observations %>% as.data.frame() %>% get_count_data(group_field = "MPA_ID", extra_grouping_var = "year") %>% dplyr::filter(complete.cases(.))
cal_coast_count_data_byWatershed_byYear <- cal_coast_observations %>% as.data.frame() %>% get_count_data(group_field = "watershed_ID", extra_grouping_var = "year") %>% dplyr::filter(complete.cases(.))
cal_coast_count_data_byCounty_byYear <- cal_coast_observations %>% as.data.frame() %>% get_count_data(group_field = "county_ID", extra_grouping_var = "year") %>% dplyr::filter(complete.cases(.))
cal_coast_count_data_byHexagon_byYear <- cal_coast_observations %>% as.data.frame() %>% get_count_data(group_field = "hexagon_ID", extra_grouping_var = "year") %>% dplyr::filter(complete.cases(.))
#' #' #### By quarter
#' cal_coast_count_data_byMPA_byQuarter <- cal_coast_observations %>% as.data.frame() %>% get_count_data(group_field = "MPA_ID", extra_grouping_var = "quarter") %>% dplyr::filter(complete.cases(.))
#' cal_coast_count_data_byWatershed_byQuarter <- cal_coast_observations %>% as.data.frame() %>% get_count_data(group_field = "watershed_ID", extra_grouping_var = "quarter") %>% dplyr::filter(complete.cases(.))
#' cal_coast_count_data_byCounty_byQuarter <- cal_coast_observations %>% as.data.frame() %>% get_count_data(group_field = "county_ID", extra_grouping_var = "quarter") %>% dplyr::filter(complete.cases(.))
#' cal_coast_count_data_byHexagon_byQuarter <- cal_coast_observations %>% as.data.frame() %>% get_count_data(group_field = "hexagon_ID", extra_grouping_var = "quarter") %>% dplyr::filter(complete.cases(.))
