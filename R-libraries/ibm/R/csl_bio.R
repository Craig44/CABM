#' @title csl2_bio
#'
#' @description
#' A function that writes a @observation block that is compatible with Casal2 for use in simulations, for observation type = proportions_at_age
#'
#' @author Craig Marsh
#' @param obs_file <sting> Complete path of the file to write the observation in
#' @param bio_expected <data.frame> an observation data frame from the IBM output
#' @param error_value <scalar,vector> an error value for the observation can be a single value (same for all years) or one for each year of the observation
#' @param catchability_label <string> label for the casal2 @catchibility block
#' @param years <vector<int>> years to apply the observation in
#' @param label <string> A label for the observation.
#' @param cells <vector<string>> cells to create observations for, this will create an observation for each cell
#' @param selectivity_labels <vector<string>> selectivity labels for each category
#' @param time_step <string> The time step in Casal2 to apply the observation in 
#' @param mortality_proportion <scalar> The proportion through the mortality block to execute the observation for
#' @param categories <list$vector<string>> A list of categories, must be equal to number of cells, You may have to use Casal2 shorthand syntax as well e.g. EN.HG+EN.EN to mean the observation is over a few categories
#' @param append <bool> are we appending to the observation file or overwriting?
#' @param likelihood <string> The likelihood for Casal2 to apply for this observation
#' @return void it creates a file
#' @export
csl2_bio = function(obs_file, bio_expected, label, selectivity_labels, error_value, catchability_label, mortality_proportion, years, cells, time_step ,categories, append = T,likelihood) {
  if (length(error_value) == 1)
    error_value = rep(error_value, length(years))
  if (length(categories) != length(cells))
    stop(paste0("expected the same number of categories for each cell, you supplied '",length(categories), "' categories but '", length(cells),"' cells please sort this out"))
  for (i in 1:length(cells)) {
    write(paste0("@observation ", label,"_",cells[i],"_bio" ),obs_file, append = T)
    write("type biomass",obs_file, append = T)
    write(paste0("likelihood ", likelihood),obs_file, append = T)
    write(paste0("selectivities ", paste(rep(selectivity_labels, length(categories[[i]])), collapse = " ")),obs_file, append = append)
    write(paste0("time_step ",time_step),obs_file, append = T)
    write(paste0("catchability ",catchability_label),obs_file, append = T)
    write(paste0("categories ",paste(categories[[i]], collapse = " ")),obs_file, append = T)
    write(paste0("years ",paste(years, collapse = " ")),obs_file, append = T)
    write(paste0("time_step_proportion ",mortality_proportion),obs_file, append = T)
    this_cell_values = bio_expected[bio_expected$year %in%  years & bio_expected$cell == cells[i] ,"simulated"]
    write(paste0("obs ",paste(round(this_cell_values,1), collapse = " ")),obs_file, append = T)
    write(paste0("error_value ",paste(error_value,collapse = " "), "\n\n"),obs_file, append = T)
  }
}