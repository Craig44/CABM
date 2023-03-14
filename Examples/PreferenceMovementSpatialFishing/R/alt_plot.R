plot_mulitple_mats = function(model, report_label, directory = "", pattern = "removals_", file_name = "", Title = "", breaks = NULL, width = NA, height = NA, plot_catch = T) {
  if (!report_label %in% names(model)) 
    stop(Paste("In model the report label '", report_label, 
               "' could not be found. The report labels available are ", 
               paste(names(model), collapse = ", ")))
  report = get(report_label, model)
  if (any(names(report) %in% "type")) {
    if (report$type != "process") {
      stop(Paste("The report label ", report_label, 
                 " in model is not a 'process' plz Check you have specified the correct report_label."))
    }
    if (is.null(report$starting_value_for_lambda)) {
      stop(Paste("The report label ", report_label, 
                 " could not find starting_value_for_lambda, this means this process isn't of type 'mortality_effort_based' which this function expects. Please check"))
    }
  } else {
    stop("multi iteration report found can't deal with this type of output yet, plz modify R code")
    muliple_iterations_in_a_report = TRUE
  }
  if (directory == "") 
    directory = getwd()
  rep_ndx = NULL
  if (plot_catch) {
    rep_ndx = which(grepl(names(report), pattern = pattern))
  } else {
    rep_ndx = which(grepl(names(report), pattern = pattern))
  }
  vec_data = vector()
  for (i in 1:length(rep_ndx)) {
    vec_data = c(vec_data, as.numeric(report[[rep_ndx[i]]]))
  }
  if (is.null(breaks)) 
    breaks = c(min(vec_data) - 1, quantile(vec_data, c(0.1, 
                                                       0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)), max(vec_data) + 
                 1)
  dup_ndx = which(duplicated(breaks))
  first_break = 0
  for (i in 1:length(dup_ndx)) {
    if (i == 1) 
      first_break = breaks[dup_ndx[i]]
    first_break = first_break + 0.1
    breaks[dup_ndx[i]] = first_break
  }
  breaks = unique(breaks)
  cols = colorRampPalette(brewer.pal(9, "YlOrRd"))(length(breaks) - 
                                                     1)
  legend_lab = vector()
  for (i in 1:(length(breaks) - 1)) {
    legend_lab[i] = round(breaks[i] + ((breaks[i + 1] - breaks[i])/2))
  }
  if (file_name == "") 
    file_name = report_label
  setwd(directory)
  for (i in 1:length(rep_ndx)) {
    this_data = melt(report[[rep_ndx[i]]])
    this_year = as.numeric(gsub(".*?([0-9]+).*", "\\1", 
                                names(report[rep_ndx[i]])))
    brks <- cut(this_data$value, breaks)
    brks = factor(brks, levels = levels(brks))
    this_data$brks <- brks
    P = ggplot(this_data, aes(y = rev(Var1), x = Var2, fill = brks)) + 
      geom_tile() + geom_tile(colour = "white", size = 0.25, 
                              show.legend = FALSE) + ggtitle(paste0(this_year, 
                                                                    " ", Title)) + labs(x = "", y = "") + 
      scale_y_discrete(expand = c(0, 0)) + scale_x_discrete(expand = c(0, 
                                                                       0)) + scale_fill_manual("Key", values = cols, 
                                                                                               na.value = "grey90", labels = round(breaks), 
                                                                                               guide = guide_legend(reverse = T), drop = FALSE)
    ggsave(P, filename = paste0(this_year, "_", file_name, 
                                ".png"), width = width, height = height)
  }
  return(NULL)
}