#'
#' Selectivities
#'
#'

#' Double normal, bell curve with different standard deviations either side of the mean (not symetric)
#' returns values that have max(X) = 1 where x = mu and min(X) = 0
#' @param x <scalar/vector> values to calculate age for
#' @param mu <scalar> meain of double normal
#' @param sig_l <scalar> standard deviation for the left hand side of distribution
#' @param sig_r <scalar> standard deviation for the right hand side of distribution
#' @return vector of double normal values
d_norm<- function(X,mu, sig_l,sig_r) {
  store<- vector()
  for( i in 1:length(X)) {
    if( X[i] <= mu) {
      store[i]<- 2^-((X[i]-mu)/sig_l)^2
    } else {
      store[i]<- 2^-((X[i]-mu)/sig_r)^2
    }
  }
  return(store)
}

## Double Exponential, U shape described by two expontial functions where the converge at (x_0,y_0) (midpoint)
#' values that have max(X) = 1 where x = mu and min(X) = 0
#' @param x <scalar/vector> values to calculate age for
#' @param x_1 <scalar> x value for y1 
#' @param x_2 <scalar> x value for y1
#' @param x_0 <scalar> x value for y0
#' @param y_1 <scalar> y value at x_1
#' @param y_2 <scalar> y value at x_2
#' @param y_0 <scalar> y value at x_0
#' @return vector of double exponential values
d_exp <- function(X, x_1, x_2, x_0, y_0, y_1, y_2) {
  store <- vector()
  for( i in 1:length(X)) {
    if( X[i] <= x_0) {
      store[i]<- min(1,y_0*(y_1/y_0)^((X[i]-x_0)/(x_1-x_0)))
    } else {
      store[i]<- min(1, y_0*(y_2/y_0)^((X[i]-x_0)/(x_2-x_0)))
    }
  }
  return(store)
}

#' Constant selectivity
#' @param C <scalar> consant
#' @return vector of double normal values
cons<- function(C) {
  if(C < 0| C > 1) {
    warning("C must be bound between 0 and 1")
  }
  return(C)
}

#' Knife Edge, selectivity X < E = 0, else = 1
#' @param x <scalar/vector> values to calculate age for
#' @param E <scalar> meain of double normal
#' @return vector of double normal values

knife<- function(X, E) {
  store<- vector()
  for( i in 1:length(X)) {
    if( X[i] < E)  {
      store[i] <- 0
    } else {
      store[i] <- 1
    }
  }
  return(store)
}

## 
#' Logistic selectivity
#' @param x <scalar/vector> values to calculate age for
#' @param a50 <scalar> value of x that has y  = 0.5
#' @param a95 <scalar> difference in x betwenn logis(X_95) = 0.95 logis(X_50) = 0.5. i.e. a95 = X_95 - X_50
#' @return vector of Logistic values
logis<- function(X, a50, a95) {
  1/(1+19^((a50-X)/a95)) 
}

#' Inverse Logistic selectivity
#' @param x <scalar/vector> values to calculate age for
#' @param a50 <scalar> value of x that has y  = 0.5
#' @param a95 <scalar> difference in x betwenn logis(X_95) = 0.95 logis(X_50) = 0.5. i.e. a95 = X_95 - X_50
#' @return vector of Logistic values
inv_logis<- function(X,a50,a95) {
  1- 1/(1+19^((a50-X)/a95)) 
}

#' Exponentialselectivity
#' @param x <scalar/vector> values to calculate age for
#' @param lambda <scalar> rate parameter
#' @return vector of exponetial values
Exp<- function(X,lambda) {
  exp(-X*lambda)
}




