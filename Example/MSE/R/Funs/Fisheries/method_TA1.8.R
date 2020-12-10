#' @title Method.TA1.8
#' @description
#' This function is useful for deciding on the data weights of one or more at-age or at-length data sets with assumed multinomial error structure in a stock assessment. Can produce a diagnostic plot if the analysis is for a single data set
#'
#' @author Chris Francis
#' @param bin_lab vector of bin labels, either ages or length mid-point
#' @param observed matrix of composition (internally resacaled in each year to be proportions) with rows = bins (ages or lengths) and cols years
#' @param expected matrix of composition (internally resacaled in each year to be proportions) with rows = bins (ages or lengths) and cols years
#' @param error_value vector of effective sample sizes
#' @param plot.it If TRUE, plot the index and the smoothed fit. Otherwise, return a dataframe of the year, index, smoothed fitted value, and cv)
#'
#' @return Outputs a mutiplier, w, so that N2y = w x N1y, where N1y and N2y are the stage-1 and stage-2 multinomial sample sizes for the data set in year y.
#'
#' @note Method TA1.8 is described in Appendix A of the following paper Francis, R.I.C.C. (2011). Data weighting in statistical fisheries stock assessment models.
#' Canadian Journal of Fisheries and Aquatic Sciences 68: 1124-1138. (With corrections to the equation in Francis R.I.C.C. (2011) Corrigendum: Data weighting in statistical fisheries stock assessment models.
#' @export

Method.TA1.8 <- function (bin_lab, observed, expected, error_value, plot.it = FALSE) {
    ## rescale to sum = 1
    Obs <- sweep(observed, 2, colSums(observed), "/")
    Exp <- sweep(expected, 2, colSums(expected), "/")
    Nassumed <- error_value
    Ry = Sy = c()
    
    ##aa <- as.numeric(strsplit(dimnames(Obs)[[2]],split = "\\["), nchar(dimnames(Obs)[[2]])))
    My <- cbind(Obs = apply(Obs, 2, function(x) sum(bin_lab * x)),
                Exp = apply(Exp, 2, function(x) sum(bin_lab * x)))
    Ry <- c(Ry, My[, "Obs"] - My[, "Exp"])
    Sy <- c(Sy, sqrt(apply(Exp, 2, function(x) sum(x * bin_lab^2)) - My[, "Exp"]^2))
  
  
  wj <- 1/var(Ry * sqrt(Nassumed)/Sy, na.rm = T)
  
  if (plot.it) {
    ses <- Sy/sqrt(Nassumed)
    Obs.bnds <- My[, "Obs"] + cbind(-2 * ses, 2 * ses)
    
    if (is.null(ylim)) {
      ylim <- range(Obs.bnds)
    }
    if (is.null(xlim)) {
      xlim <- range(years)
    }
    
    plot(years, My[, "Obs"], type = "n", ylab = "", xlab = "", xlim = xlim, ylim = ylim, las = 1)
    points(years, My[, "Obs"], pch = "x", col = 3)
    segments(years, Obs.bnds[, 1], years, Obs.bnds[, 2], col = 3)
    lines(years, My[, "Exp"], col = "red")
  }
  
  return (wj)
}
