#' @title extract MPD function for readin in IBM output that has been generated from a -rrun mode.  
#'
#' @description
#' An extract function that reads IBM output that are produced from a '-r' model run. This funciton
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

    result = list()
    for (i in 1:counter) {
        header = split.header(temp[i])
        label = header[1]
        type = header[2]
        report = get.lines(file, clip.to = temp[i])
        report = get.lines(report,clip.from = "*end")
        report = make.list(report)
        report$type = type
        if (!is.in(label, names(result))) {       
           result[[label]] = list()
           result[[label]][["1"]] = report
        } else {
             result[[label]][[as.character(length(result[[label]])+1)]] = report
        }
        file = get.lines(file, clip.to = "*end")

    } 
    result<-set.class(result,"ibm_output")
    result
  } else {
    warning("File is empty, no reports found")
  }
}


  
