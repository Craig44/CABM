#' @title csl2_fishery_age
#'
#' @description
#' A function that writes a @observation block that is compatible with Casal2 for use in simulations, for observation type = process_removals_by_age
#'
#' @author Craig Marsh
#' @param obs_file <sting> Complete path of the file to write the observation in
#' @param obs_df <data.frame> an observation data frame from the IBM output
#' @param error_value <scalar,vector> an error value for the observation can be a single value (same for all years) or one for each year of the observation
#' @param age_range <vector> a vector of the form c(min_age, max_age)
#' @param plus_group <bool> is the max age a plus group
#' @param fishery <vector<string>> A vector of fisheries labels can take up to two
#' @param years <vector<int?>> years to apply the observation in
#' @param cells <vector<string>> cells to create observations for, this will create an observation for each cell
#' @param time_step <string> The time step in Casal2 to apply the observation in 
#' @param mort_process <string> the mortality process that this observation come from in Casal2
#' @param categories <list$vector<string>> A list of categories, must be equal to number of cells, You may have to use Casal2 shorthand syntax as well e.g. EN.HG+EN.EN to mean the observation is over a few categories
#' @param append <bool> are we appending to the observation file or overwriting?
#' @param likelihood <string> The likelihood for Casal2 to apply for this observation
#' @return void it creates a file
#' @export

## age settings for process removals by age
csl2_fishery_age = function(obs_file, obs_df, error_value, age_range, plus_group, fishery, years, cells, time_step, mort_process,categories, append = T,likelihood) {
  if (length(age_range) != 2)
    stop("expected age_range to have two values one for min_age and one for max_age")
  min_age = age_range[1]
  max_age = age_range[2]
  if (length(error_value) == 1)
    error_value = rep(error_value, length(years))
  if (length(categories) != length(cells))
    stop(paste0("expected the same number of categories for each cell, you supplied '",length(categories), "' categories but '", length(cells),"' cells please sort this out"))
  
  for (i in 1:length(cells)) {
    write(paste0("@observation ", fishery[i],"_",cells[i], "_age"),obs_file, append = append)
    write("type process_removals_by_age",obs_file, append = T)
    write(paste0("likelihood ", likelihood),obs_file, append = T)
    write(paste0("plus_group ", plus_group) ,obs_file, append = T)
    write(paste0("method_of_removal ", fishery[i]),obs_file, append = T)
    write(paste0("mortality_instantaneous_process ",mort_process),obs_file, append = T)
    write(paste0("min_age ",min_age),obs_file, append = T)
    write(paste0("max_age ",max_age),obs_file, append = T)
    write(paste0("time_step ",time_step),obs_file, append = T)
    write(paste0("categories ",paste(categories[[i]], collapse = " ")),obs_file, append = T)
    write(paste0("years ",paste(years, collapse = " ")),obs_file, append = T)
    write("table obs",obs_file, append = T)
    ## obs
    for (y in 1:length(years)) {
      this_year = obs_df[obs_df$year == years[y] & obs_df$cell == cells[i],]
      ## truncate data
      min_age_ndx = which(this_year$age == min_age) ## minus plus group? maybe not
      max_age_ndx = which(this_year$age == max_age) ## minus plus group? maybe not
      
      plus_ndx = this_year$age >= max_age 
      plus_group = sum(this_year$simulated[plus_ndx])
      new_age_freq = this_year$simulated[min_age_ndx:max_age_ndx]
      new_age_freq[length(new_age_freq)] = plus_group
      new_age_freq = round(new_age_freq / sum(new_age_freq),4)
      write(paste0(years[y], " ",paste(new_age_freq, collapse = " ")),obs_file, append = T)
    }
    write("end_table\n",obs_file, append = T)
    ## error
    
    write("table error_values",obs_file, append = T)
    for (y in 1:length(years)) {
      write(paste0(years[y], " ", error_value[y]),obs_file, append = T)
    }
    write("end_table\n\n\n",obs_file, append = T)
  }
}