## test IBM's Logistic normal random number generator
## with Chris's code as a validation step
## need to source Chris's original logistic normal R-code if you fully want to run this
## can be found here 'https://github.com/Craig44/LogisticNormalLL/tree/master/Rstudio'
library(mvtnorm)

source("Chris_original_code.R")

### Casal2 & IBM logistic normal simulator
rlogisNorm = function(n = 100, exp_prop, covar) {
  if (ncol(covar) != length(exp_prop))
    stop("dimensions wrong")
  store_mat = matrix(0.0,nrow = length(exp_prop), ncol = n)
  cho_covar = t(chol(covar)) ## this is what IBM produces the DoCholeskyDecomposition is the transpose of R's chol() function
  for( k in 1:n) {
    for (i in 1:length(exp_prop)) {
      rng_norm = rnorm(n = nrow(covar))
      row_sum = 0.0;
      for(j in 1:nrow(covar)) 
        row_sum = row_sum + cho_covar[i, j] * rng_norm[j]
      store_mat[i,k] = exp(row_sum + log(exp_prop[i]))
    }
  }
  sweep(store_mat,2,apply(store_mat,2,sum),'/')
}

## one with a multivariate normal from an independent package
rlogisNorm_alt = function(n = 100, exp_prop, covar) {
  if (ncol(covar) != length(exp_prop))
    stop("dimensions wrong")
  store_mat = exp(rmvnorm(n = n, mean = log(exp_prop), sigma = covar))
  t(sweep(store_mat,1,apply(store_mat,1,sum),'/'))
}

Sigma = 0.5
phi = c(0.8,-0.2)
exp_prop = c(0.026988805, 0.045033722, 0.055701479, 0.045370913, 0.120156357, 0.101589149, 0.057220146, 0.027703137, 0.017341034, 0.027337446, 0.024559262, 0.033917239, 0.031879276, 0.052872193, 0.046389372, 0.034637439, 0.035200240, 0.031583224, 0.044784866, 0.020459028, 0.022896556, 0.017558524, 0.018343697, 0.008215704, 0.005650193, 0.007288457, 0.007903854, 0.002525200, 0.028893487)
covar = covmat.logistnorm(sigma = Sigma, phi = phi, binnam = paste0("X",2:30), ARMA = F)

N = 100000
alt_sim = rlogisNorm_alt(n = N, exp_prop = exp_prop, covar = covar)
sim = rlogisNorm(n = N, exp_prop = exp_prop, covar = covar)
## compare with Chris Francis code to double triple check
orig = (rlogistnorm(n = N, expprop = exp_prop, covmat = covar))

apply(orig, 1, mean)
apply(alt_sim, 1, mean)
apply(sim, 1, mean)

apply(orig, 1, sd)
apply(alt_sim, 1, sd)
apply(sim, 1, sd)

## same expectation
ages = 2:30
plot(ages, apply(orig, 1, mean), type = "l",lwd = 2, col = "red")
lines(ages, apply(alt_sim, 1, mean), col = "blue",lwd =2, lty = 2)
lines(ages, apply(sim, 1, mean), col = "black",lwd =2, lty = 2)

plot(ages, apply(orig, 1, sd), type = "l",lwd = 2, col = "red")
lines(ages, apply(alt_sim, 1, sd), col = "blue",lwd =2, lty = 2)
lines(ages, apply(sim, 1, sd), col = "black",lwd =2, lty = 2)

## All looks good

