#' Utility function
#'
#' @author Craig Marsh
#' @keywords internal
#'
"get.unique_subcommands_list" <- function() {
  command = c("table", "end_table")
  type = c("table_label", "end_table")
  ibm_list = list(command = command, type = type)
  return(ibm_list)
}

