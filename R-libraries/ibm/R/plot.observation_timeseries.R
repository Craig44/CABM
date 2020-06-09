#' @title plot.observation_timeseries
#'
#' @description
#' A plotting to plot which observation occur when over the model cycle
#'
#' @author Craig Marsh
#' @param ibm_model <ibm_output> object that are generated from one of the extract.run() functions.
#' @param model_years <integers> years to create the plot over, x axis
#' @rdname plot.observation_timeseries
#' @export plot.observation_timeseries
plot.observation_timeseries = function(model_years, ibm_model) {
  n_reports = length(ibm_model);
  report_labs = names(ibm_model)
  DF = data.frame(obs = as.character(), years = as.numeric())
  labs_ = vector();
  obs_years = vector()
  for(i in 1:n_reports) {
    if (report_labs[i] != "model_run_time") {
      if (ibm_model[[i]]$type == "observation") {
        print(paste0("found obs = ", report_labs[i]))
        years = as.numeric(unique(ibm_model[[i]]$Values$year))
        labs_ = c(labs_, rep(report_labs[i], length(years)))
        obs_years = c(obs_years, years)
      }
    }
  }
  DF = data.frame(obs = labs_, years = obs_years)
  DF$years = DF$years
  
  gplot = ggplot(data = DF, aes(x = years, y = obs, col = obs, group = obs)) + 
    geom_line(size = 1.4) + 
    xlim(min(model_years), max(model_years)) + 
    theme(legend.position = "none") + 
    xlab("Observations")
  print(gplot)
}
