#' BH recruitmen
#'
#' @param SSB Spawning stock biomass
#' @param B0 equilibrium SSB
#' @param h steepness parameter as defined by Doonan and Mace 1981, represents 
#'
#' @return Number of recruits acording to the Beverton holt relationship
#' @export
#'
BH = function(SSB,B0,h) {
  
  ssb_ratio = SSB / B0
  part_2 = (1 - ((5*h - 1) / (4*h)) * ( 1 - ssb_ratio))
  val = ssb_ratio / part_2
  return(val)
}

#' BH_non_equil recruitmen
#'
#' @param SSB Spawning stock biomass
#' @param alpha classic params alpha BH param
#' @param beta classic params beta BH param
#'
#' @return Number of recruits acording to the Beverton holt relationship
#' @export
#'
BH_non_equil = function(SSB, alpha, beta) {
  R = SSB / (alpha + beta * SSB)
  return(R)
}

