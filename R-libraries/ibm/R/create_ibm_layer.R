#' @title create_ibm_layer
#'
#' @description
#' A function that creates a file consists of @layer block.
#'
#' @author Craig Marsh
#' @param label <sting> label for the @layer command
#' @param type <string> what type of layer is it
#' @param filename <string> the full filename of the file it will overwrite any file with the same address
#' @param Matrix the matrix to fill the Layer with
#' @param proporstion <bool> if true then it tells, the ibm this matrix should sum to 1
#' @return void it creates a file
#' @export

create_ibm_layer = function(label, type, filename, Matrix, proportion = FALSE) {
  # check type is allowd
  allowed_types = c("numeric","numeric_meta","categorical","integer", "biomass")
  
  if (!type %in% allowed_types)
    stop(paste0("type '", type, "' is not an allowed layer, you need to either check what you are doing or change the R functions"));
  
  cat(paste0("@layer ", label,"\n"),file = filename, append = FALSE);
  cat(paste0("type ", type,"\n"),file = filename, append = TRUE);
  if (proportion)
    cat("proportions true\n",file = filename, append = TRUE);
  cat("table layer\n",file = filename, append = TRUE);
  write(t(Matrix),file = filename, append = TRUE, ncolumns  = ncol(Matrix));
  cat("end_table\n",file = filename, append = TRUE);
  return(NULL);
}
