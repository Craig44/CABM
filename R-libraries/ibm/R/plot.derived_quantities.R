#' @title plot.derived_quantities default
#'
#' @description
#' A plotting function to plot Derived Quantities for the 'ibm_output' object.
#'
#' @author Craig Marsh
#' @param model <ibm_output> object that are generated from one of the extract() functions.
#' @param report_label <string>
#' @param type <string> whether numbers or scaled by an initialisation value
#' @param plot.it Whether to generate a default plot or return the values as a matrix.
#' @param ... remaining plotting options
#' @return A plot of derived quantities over time if plot.it = T, if plot.it = F it will return a matrix of derived quantities.
#' @rdname plot.derived_quantities
#' @export plot.derived_quantities
#' @examples
#' library(ibm)
#' # plotting Standard Output
#' data <- extract.run(file = system.file("extdata", "estimate.log", package="casal2"))
#' names(data)
#' plot.derived_quantities(model = data, report_label = "biomass")
#' # if you are unhappy with the default plotting you can use plot.it = FALSE and create a plot of your own.
#' SSB = plot.derived_quantities(model = data, report_label = "biomass", plot.it = FALSE)
#' plot(rownames(SSB),SSB, main = "My SSB", type = "l", ylim = c(0,90000))

"plot.derived_quantities"<-
function(model, report_label="", type = "number", xlim, ylim, xlab, ylab, main, col,plot.it = T, ...){
  UseMethod("plot.derived_quantities",model)
}

#' @return \code{NULL}
#'
#' @rdname plot.derived_quantities
#' @method plot.derived_quantities ibm_output
#' @export
"plot.derived_quantities.ibm_output" = function(model, report_label="", type = "number", xlim, ylim, xlab, ylab, main, col,plot.it = T, ...) {
  if (!type %in% c("number", "percent")) {
    stop ("the parameter type must be: 'number' or 'percent'")
  }
  muliple_iterations_in_a_report = FALSE;
  N_runs = 1;
  temp_DF = NULL;

  ## check report label exists
  if (!report_label %in% names(model))
    stop(Paste("In model the report label '", report_label, "' could not be found. The report labels available are ", paste(names(model),collapse = ", ")))
  ## get the report out
  this_report = get(report_label, model)
  ## check that the report label is of type derived_quantity
  if (this_report$'1'$type != "derived_quantity") {
    stop(Paste("The report label ", report_label, " in model is not a derived quantity plz Check you have specified the correct report_label."))     
  }
  if (length(this_report) > 1) {
      print("multi iteration report found")
      muliple_iterations_in_a_report = TRUE;
      N_runs = length(this_report);
  }
  if (!muliple_iterations_in_a_report) {
    ## only a single trajectory
    derived_quantities = names(this_report$'1')[names(this_report$'1') != "type"]
    ## create a multi-plot panel
    if (plot.it)
      par(mfrow = c(1,length(derived_quantities)))
    for (i in 1:length(derived_quantities)) {
      this_derived_q = get(derived_quantities[i], this_report$'1')
      values = this_derived_q$values
      zero_ndx = values == 0;
      if (any(zero_ndx)) {
        ## I am going to assume that these are in the projection phase and am going to remove them
        values = values[!zero_ndx]      
      }
      years = as.numeric(names(values))
      ## does the user want it plotted as percent B0
      if (type == "percent") {
        B0 = "hello"
        values = values / this_derived_q$B0 * 100
      }  
      if(missing(ylim)) {
        ymax = max(values) + quantile(values, 0.05) 
        ylim = c(0, ymax)
      }    
      if(missing(xlim)) 
        xlim = c(min(years) - 1, max(years) + 1)    
      if(missing(ylab)) {
        if (type == "percent")
          ylab = "%B0"
        else
          ylab = "Biomass (t)"
      }
      if(missing(xlab)) 
        xlab = "Years"
      if(missing(col)) 
        col = "black"      
      if(missing(main)) 
        main = derived_quantities[i]  
      if (plot.it == TRUE) {
        plot(years, values, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main, type = "o", ...)
      } else {
        temp_DF = cbind(temp_DF,values);
      }
    }  
    if (plot.it == FALSE)
      colnames(temp_DF) = derived_quantities
  } else {
    ## Multiple trajectory's
    derived_quantities = names(this_report$'1')[names(this_report$'1') != "type"]
    
    #print(Paste("found - ",length(this_report), " iterations"))
    if(missing(col)) {
      palette(gray(seq(0.,.90,len = N_runs)))
      Cols = palette()
    } else {
      Cols = col;
    }
    temp_DF = list();
    if (plot.it)
      par(mfrow = c(1,length(derived_quantities)))
    for (k in 1:length(derived_quantities)) {
    temp_df = c()
      first_run = TRUE;
      Legend_txt = c();
      for (i in 1:length(this_report)) {
        this_derived_quantity = this_report[[i]]
        this_derived_quantity = get(derived_quantities[k], this_derived_quantity)
        values = this_derived_quantity$values
        zero_ndx = values == 0;
        if (any(zero_ndx)) {
          ## I am going to assume that these are in the projection phase and am going to remove them
          values = values[!zero_ndx]      
        }
        years = as.numeric(names(values))
        ## does the user want it plotted as percent B0
        if (type == "percent")
          values = values / this_derived_quantity$B0 * 100

        Legend_txt = c(Legend_txt,this_derived_quantity$B0)  

        if (first_run) {      
          if(plot.it == FALSE) {
            temp_df = values
          }      
          if(missing(ylim)) {
            ymax = max(values) + quantile(values, 0.15) 
            ylim = c(0, ymax)
          }    
          if(missing(xlim)) 
            xlim = c(min(years) - 1, max(years) + 1)    
          if(missing(ylab)) {
            if (type == "percent")
              ylab = "%B0"
            else
              ylab = "Biomass (t)"
          }
          if(missing(xlab)) 
            xlab = "Years"    
          if(missing(main)) 
            main = derived_quantities[k]      
          if (plot.it == TRUE) {  
            plot(years, values , xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main, type = "o", lwd = 3, col = Cols[i] ,...)        
          }
        } else {
          if (plot.it == TRUE) {
            lines(years, values, type = "o", lwd = 3, col = Cols[i])          
          } else {
            temp_df = rbind(temp_df,values)
          }
        }
        first_run = FALSE;
      }  
      if (plot.it == TRUE) {
        legend('bottomleft',legend = Legend_txt, col = Cols, lty = 2, lwd = 3)
      } else {
        rownames(temp_df) = as.character(1:length(this_report))
        temp_DF[[derived_quantities[k]]] = temp_df
      }
    }
  }
  if (plot.it == FALSE)
    return(temp_DF)
    
  invisible();  
}
