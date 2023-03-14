#' constant selectivity
#' @param C values
#' @return constant selectivty
#' @export
cons<- function(C)
{
  if(C<0| C>1)
  {
    warning("C must be bound between 0 and 1")
  }
  return(C)
}

#' Knife Edge selectivity
#' @export

knife<- function(X,E)
{
  store<- vector()
  for( i in 1:length(X))
  {
  if( X[i] < E)
  {store[i]<- 0} else {
   store[i]<- 1}
  }
  return(store)
  }

#' Logistic selectivity
#' @export
logis<- function(X,a50,a95)
{
 1/(1+19^((a50-X)/a95)) 
}

#' Logistic selectivity
#' @export
logis<- function(X,a50,a95)
{
 1/(1+19^((a50-X)/a95)) 
}

#'  Maturity calculated in SS
#' @export
ss_mat = function (X, a50, slope) {
  1/(1 + exp(slope * (X - a50)))
}

#' Inverse Logisitic selectivity
#' @export

inv_logis<- function(X,a50,a95)
{
  1- 1/(1+19^((a50-X)/a95)) 
}

#' Exponential selectivity
#' @export

Exp<- function(X,lambda)
{
  exp(-X*lambda)
}


#' Double normal selectivity
#' @export

d_norm<- function(X,mu, sig_l,sig_r)
{
  store<- vector()
  for( i in 1:length(X))
  {
    if( X[i] <= mu)
    {store[i]<- 2^-((X[i]-mu)/sig_l)^2} else {
      store[i]<- 2^-((X[i]-mu)/sig_r)^2}
  }
  return(store)
}


#' Double Exponential selectivity
#' @export

d_exp <- function(X,x_1,x_2,x_0,y_0,y_1,y_2)
{
  store<- vector()
  for( i in 1:length(X))
  {
    if( X[i] <= x_0)
    {store[i]<- min(1,y_0*(y_1/y_0)^((X[i]-x_0)/(x_1-x_0)))} else {
      store[i]<- min(1, y_0*(y_2/y_0)^((X[i]-x_0)/(x_2-x_0)))}
  }
  return(store)
}


