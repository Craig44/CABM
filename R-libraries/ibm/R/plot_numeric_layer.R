#' @title plot_numeric_layer
#'
#' @description
#' A plotting function to generate a snapshot of heatmaps that can be used in an animation, to show how the 
#' the layer changes over time
#'
#' @author Craig Marsh
#' @param model <ibm_output> object that are generated from one of the extract() functions.
#' @param report_label <string>
#' @return A plot of of age frequency
#' @rdname plot_numeric_layer
#' @export plot_overall_age_freq
#' @examples
#' library(ibm)


plot_numeric_layer = function(model, report_label) {
  ## check report label exists
  if (!report_label %in% names(model))
    stop(Paste("In model the report label '", report_label, "' could not be found. The report labels available are ", paste(names(model),collapse = ", ")))
  ## get the report out
  report = get(report_label, model)
  if(report$`1`$type != "numeric_layer")
    stop(Paste("report " , report_label, " needs to be a numeric_layer report, sort it out"))
  
  ## go through data and find, min, max and some other data characteristics that might help with aesthetics
  vec_data = vector()
  for (i in 1:length(report)) {
    vec_data = c(vec_data,as.numeric(report[[i]]$values))
  }
  breaks = c(min(vec_data) - 1, quantile(vec_data, c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)), max(vec_data) + 1)
  cols = colorRampPalette(brewer.pal(9,"YlOrRd"))(length(breaks) - 1)
  legend_lab = vector();
  for (i in 1:(length(breaks) - 1)) {
    legend_lab[i] = round(breaks[i] + ((breaks[i + 1] - breaks[i])/2))
  }
  
  dir.create("Figures")
  setwd("Figures")
  for (i in 1:length(report)) {
    this_data = melt(report[[i]]$values)
    brks <- cut(this_data$value,breaks)
    brks <- gsub(","," - ",brks,fixed=TRUE)
    this_data$brks <- sapply(brks, function(x){strsplit(x, "[( ]")[[1]][[2]]})
    if (i == 1) {
      # add legend otherwise skip it
      P = ggplot(this_data, aes(y = rev(Var1), x = Var2))  + geom_tile(aes(fill=brks)) +
        scale_fill_manual("Biomass",values = cols, labels = legend_lab, guide = guide_legend(reverse=TRUE)) +
        xlab("")+ylab("")
    } else {
      P = ggplot(this_data, aes(y = rev(Var1), x = Var2))  + geom_tile(aes(fill=brks)) +
        scale_fill_manual("Biomass",values = cols, labels = legend_lab, guide = guide_legend(reverse=TRUE)) +
        theme(legend.position="none") +
      xlab("")+ylab("")
    }
    ggsave(P, filename = paste0(report[[i]]$year,"_", report_label,".png"))
  }
}
  