#' @title ssb_multiple_runs
#'
#' @description
#' A plotting function to plot Derived Quantities for the 'ibm_output' object.
#'
#' @author Craig Marsh
#' @param SSB_outputs_df a data frame of derived quantities
#' @param labels labels for the legend
#' @param Ylim y axis limits
#' @return A ggplot to compare multiple ibm run's and also Casal2 runs.
#' @rdname ssb_multiple_runs
#' @export ssb_multiple_runs
#' @examples
#' library(ibm)
#' # plotting Standard Output
#' data <- extract.run(file = system.file("extdata", "run.log", package="ibm"))
#' data2 <- extract.run(file = system.file("extdata", "run2.log", package="ibm"))
#' names(data)
#' SSB_df = cbind(data$derived_quants$'1'$SSB, data2$derived_quants$'1'$SSB)
#' ssb_multiple_runs(SSB_df, labels = c("high M", "low M"), Ylim = c(0,90000))
## lets look at 

ssb_multiple_runs = function(SSB_outputs_df, labels, Ylim = NULL) {
  colnames(SSB_outputs_df) = labels
  melt_ssb = melt(SSB_outputs_df)
  colnames(melt_ssb) = c("year", "run", "SSB")
  
  p = ggplot(melt_ssb, aes(x = year, y = SSB, group = run, color = run)) + 
    geom_line(size = 1.3)
  if (!is.null(Ylim)) {
    p = p + ylim(Ylim)
  }
  print(p)
  return(p)
}