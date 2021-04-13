#' ---
#' title: DOB Cal Coast <br> <strong> Explore Places - Spatial Beta Diversity
#' ---
#' ### Calculate beta diversity across spatial boundary units using the Morisita-Horn index
#' #### All species
#' #### MPAs
uniqueness_estimates_byMPA <- estimate_beta_diversity_from_counts(count_data = cal_coast_count_data_byMPA, min_species = 10, suffix = "all")
#' #### Watersheds
uniqueness_estimates_byWatershed <- estimate_beta_diversity_from_counts(count_data = cal_coast_count_data_byWatershed, min_species = 10, suffix = "all")
#" #### Counties
uniqueness_estimates_byCounty <- estimate_beta_diversity_from_counts(count_data = cal_coast_count_data_byCounty, min_species = 10, suffix = "all")
#' #### Hexagons
uniqueness_estimates_byHexagon <- estimate_beta_diversity_from_counts(count_data = cal_coast_count_data_byHexagon, min_species = 10, suffix = "all")
#'
#' #### Individual taxonomic groups
#' #### MPAs
uniqueness_estimates_byMPA_byTaxon <- purrr::map(taxonomic_groups, function(taxon) {
  estimate_beta_diversity_from_counts(count_data = cal_coast_count_data_byMPA, min_species = 5, focal_taxon = taxon, suffix = taxon)
})
#' #### Watersheds
uniqueness_estimates_byWatershed_byTaxon <- purrr::map(taxonomic_groups, function(taxon) {
  estimate_beta_diversity_from_counts(count_data = cal_coast_count_data_byWatershed, min_species = 5, focal_taxon = taxon, suffix = taxon)
})
#' #### Counties
uniqueness_estimates_byCounty_byTaxon <- purrr::map(taxonomic_groups, function(taxon) {
  estimate_beta_diversity_from_counts(count_data = cal_coast_count_data_byCounty, min_species = 5, focal_taxon = taxon, suffix = taxon)
})
#' #### Hexagons
uniqueness_estimates_byHexagon_byTaxon <- purrr::map(taxonomic_groups, function(taxon) {
  estimate_beta_diversity_from_counts(count_data = cal_coast_count_data_byHexagon, min_species = 5, focal_taxon = taxon, suffix = taxon)
})