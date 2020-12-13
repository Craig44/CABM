#' @title extract function for readin in IBM output that has been generated from a '-r' & '-s' run mode.  
#'
#' @description
#' An extract function that reads IBM output that are produced from a '-r' & '-s' model run. This funciton
#' also create a 'ibm_output' class which can be used in plotting and summary functions. See the casal2 manual for more information.
#'
#' @author Dan Fu & C.Marsh
#' @param file the name of the input file containing model output to extract
#' @param path Optionally, the path to the file
#' @param fileEncoding Optional, allows the R-library to read in files that have been encoded in alternative UTF formats, see the manual for the error message that would indicate when to use this switch.
#' @return a 'ibm_output' object which is essentially a list, that can be integrated using the str() function.
#' @export
#' @examples
#' library(ibm)
#' data <- extract.run(file = system.file("extdata", "output.log", package="ibm"))
#' class(data)
"extract.run" <-
function (file, path = "",fileEncoding = "") {
  set.class <- function(object,new.class){
    # use in the form
    #  object <- set.class(object,"new class")
    attributes(object)$class <- c(new.class,attributes(object)$class[attributes(object)$class!=new.class])
    object
  }

  filename = make.filename(path = path, file = file)
  file = convert.to.lines(filename, fileEncoding = fileEncoding)
    
  temp = get.lines(file, starts.with = "\\*",fixed=F)
  if (length(temp) != 0) {
    if(!is.even(length(temp))) {
        ## find the report which doesn't have a *end
        nd_element = seq(2, length(temp), by = 2)
        ndx = which(temp[nd_element] != "*end")[1]
        stop (Paste("Each report section must begin with '*' and end with '* end'. The report beginning with ", temp[ndx], " did not have a trailing *end"))
    }
    temp = temp[is.odd(1:length(temp))]
    counter = length(temp)
    ## iterate over all reports and see if this is a multi input run, this is identified by checking if -i in the header &
    ## if ALL reports are duplicated. Some reports will be duplicated because they are year based reports
    mult_input_run = FALSE;
    if (grepl(pattern = "-i", x = file[2]) | grepl(pattern = "-s", x = file[2])) {
      if (sum(!duplicated(temp)) == length(unique(temp))) {
        mult_input_run = TRUE;
        cat("loading a run from -i or -s format\n")
      }
    }
    get_time = substring(get.lines(file, starts.with = "Total elapsed time:",fixed=F), first = 21)
    
    multi_year_reports = c("summarise_agents", "world_age_frequency","numeric_layer","age_frequency_by_cell","age_length_matrix_by_cell")
    
    result = list()
    for (i in 1:counter) {
        header = split.header(temp[i])
        label = header[1]
        type = header[2]
        report = get.lines(file, clip.to = temp[i])
        report = get.lines(report,clip.from = "*end")
        report = tryCatch({
	   make.list(report)
	}, warning = function(warning_condition) {
	    cat(paste0("For report label ",label ," found the following warning ", warning_condition))
	}, error = function(error_condition) {
	    cat(paste0("For report label ",label ," found the following error ", error_condition))
	})
        report$type = type
        if (mult_input_run) {
          ## deal with the multi input run.
          if (!is.in(label, names(result))) {
            ## if first appearence of report register it as 1
            result[[label]][["1"]] = list();
            if (type %in% multi_year_reports) {
              result[[label]][["1"]][[as.character(report$year)]] = report
              result[[label]][["1"]][[as.character(report$year)]][["type"]] = NULL
              result[[label]][["1"]][[as.character(report$year)]][["year"]] = NULL
              
              result[[label]][["1"]][["type"]] = type
            } else {
              result[[label]][["1"]] = report
            }
          } else if (is.in(label, names(result)) & type %in% multi_year_reports) {  
            ## if we have seen this report find out if we are adding to -i or year
            if (report$year %in% names(result[[label]][[as.character(length(result[[label]]))]])) {
              result[[label]][[as.character(length(result[[label]])+1)]][[as.character(report$year)]] = report
              result[[label]][[as.character(length(result[[label]]))]][[as.character(report$year)]][["type"]] = NULL
              result[[label]][[as.character(length(result[[label]]))]][[as.character(report$year)]][["year"]] = NULL
              result[[label]][[as.character(length(result[[label]]))]][["type"]] = type

            } else {
              result[[label]][[as.character(length(result[[label]]))]][[as.character(report$year)]] = report
              result[[label]][[as.character(length(result[[label]]))]][[as.character(report$year)]][["type"]] = NULL
              result[[label]][[as.character(length(result[[label]]))]][[as.character(report$year)]][["year"]] = NULL
              
            }
          } else {
            ## a simple report
            result[[label]][[as.character(length(result[[label]])+1)]] = report
          }
          file = get.lines(file, clip.to = "*end")
        } else {
          ## deal with the single input run.
          if (type %in% multi_year_reports) {       
            result[[label]][[as.character(report$year)]] = report
            result[[label]][[as.character(report$year)]][["type"]] = NULL;
            result[[label]][["type"]] = type
          } else {
            result[[label]] = report
          }
          file = get.lines(file, clip.to = "*end")
        }
    }
    result[["model_run_time"]] = get_time
    result<-set.class(result,"ibm_output")
    result
  } else {
    warning("File is empty, no reports found")
  }
}


  
