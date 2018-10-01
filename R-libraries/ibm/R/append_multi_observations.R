#' @title append_multi_observations
#'
#' @description
#' A function that takes multi run output and creates a single data frame that is compatible with ggplot for observations
#'
#' @author Craig Marsh
#' @param model <ibm_output> object that are generated from one of the extract() functions.
#' @param report_label <string>
#' @return data.frame
#' @export 


append_multi_observations = function(model, report_label) {
  ## check report label exists
  if (!report_label %in% names(model))
    stop(Paste("In model the report label '", report_label, "' could not be found. The report labels available are ", paste(names(model),collapse = ", ")))
  ## get the report out
  this_ob = get(report_label, model)
  n_copies = length(names(this_ob))
  final_df = NULL;
  for (i in 1:n_copies) {
    this_copy = this_ob[[i]]$Values
    this_copy$simulation = i
    final_df = rbind(final_df, this_copy)
  }
  return(final_df)
}
