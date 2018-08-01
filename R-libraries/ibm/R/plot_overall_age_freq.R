#' @title plot_overall_age_freq
#'
#' @description
#' A plotting function to plot age frequency at a point in time of the world vs teh fishery
#'
#' @author Craig Marsh
#' @param model <ibm_output> object that are generated from one of the extract() functions.
#' @param world_age_frequency_label <string>
#' @param fishery_age_frequency_label <string> whether numbers or scaled by an initialisation value
#' @param proportion <bool> plot frequency or proportion
#' @param year <int> the year to plot for
#' @return A plot of of age frequency
#' @rdname plot_overall_age_freq
#' @export plot_overall_age_freq
#' @examples
#' library(ibm)


plot_overall_age_freq = function(model, world_age_frequency_label, fishery_age_frequency_label, proportion = TRUE, year, add_legend = TRUE) {
  ## check report label exists
  if (!world_age_frequency_label %in% names(model))
    stop(Paste("In model the report label '", world_age_frequency_label, "' could not be found. The report labels available are ", paste(names(model),collapse = ", ")))
  if (!fishery_age_frequency_label %in% names(model))
    stop(Paste("In model the report label '", fishery_age_frequency_label, "' could not be found. The report labels available are ", paste(names(model),collapse = ", ")))
  ## get the report out
  world_report = get(world_age_frequency_label, model)
  fishery_report = get(fishery_age_frequency_label, model)
  if(fishery_report$`1`$type != "process")
    stop(Paste("report " , fishery_age_frequency_label, " needs to be a process report, sort it out"))
  
  if (world_report$`1`$type != "world_age_frequency") {
    stop(Paste("report " , world_age_frequency_label, " needs to be a world_age_frequency report, sort it out"))
  }
  # check year exists for both
  if (!year %in% fishery_report$'1'$years)
    stop(Paste("couldn't find year ", year , " in the fishing process make sure this process happens in this year"))
  fishery_year_index = which(fishery_report$'1'$years %in% year)
  year_index = 0;
  for (i in 1:length(world_report)) {
    if (world_report[[i]]$year == year) {
      year_index = i;
    }
  }
  age_spread = length(fishery_report$`1`$'age_frequency'[fishery_year_index,]) - 1
  if (year_index == 0)
    stop(Paste("couldn't find year ", year , " in the world age frequency, something is not right"))
  data = NULL;
  if (proportion) {
    data = rbind(world_report[[year_index]]$values / sum(world_report[[year_index]]$values), fishery_report$`1`$'age_frequency'[fishery_year_index,2:(age_spread + 1)] / sum(fishery_report$`1`$'age_frequency'[fishery_year_index,2:(age_spread + 1)]))
  } else {
    data = rbind(world_report[[year_index]]$values, fishery_report$`1`$'age_frequency'[fishery_year_index,2:(age_spread + 1)])
  }
  data$type = c("world","fishery")
  merged_data = melt(data)
  
  P = NULL
  if (add_legend) {
    P = ggplot(merged_data, aes(x = variable, y = value, color = type)) + geom_point() +
      scale_fill_discrete(name="type",
                          labels=c("World", "Fishery"))+ 
      xlab("Age")+ylab("fequency")
  } else {
    P = ggplot(merged_data, aes(x = variable, y = value, color = type)) + geom_point() +
      xlab("Age")+ylab("type") +
      theme(legend.position="none")
  }
    
  print(P)
  ## consider returning the data too if people want to modify the plots
}
