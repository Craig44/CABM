#' plot_categorical_layer
#' @param categorical_matrix matrix that you want to visulise
#' @param fill_color boolean whether to include the fill color for each category group
#' @return ggplot
#' @export
plot_categorical_layer = function(categorical_matrix, fill_color = T) {
  molten_cat_mat = reshape2::melt(categorical_matrix)
  colnames(molten_cat_mat) = c("Row", "Col", "value")
  gplt = NULL
  if(fill_color) {
    gplt = ggplot(molten_cat_mat, aes(x = Col, y = Row, fill = value)) +
      geom_tile(col = "black") +
      geom_text(aes(x = Col, y = Row, label = value), color = "black", size = 4) +
      labs(x = "", y = "", fill = "") +
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            legend.position = "none",
            panel.border = element_blank(),
            panel.background = element_blank(),
            plot.background = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank())
  } else {
    gplt = ggplot(molten_cat_mat, aes(x = Col, y = Row)) +
      geom_tile(col = "black", fill = NA) +
      geom_text(aes(x = Col, y = Row, label = value), color = "black", size = 4) +
      labs(x = "", y = "", fill = "") +
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            legend.position = "none",
            panel.border = element_blank(),
            panel.background = element_blank(),
            plot.background = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) 
  }
  return(gplt)
}