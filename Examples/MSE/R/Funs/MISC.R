# every time a function is added here it needs to be documneted, don't be a fool and skip over it?
# list of functions
#'rmvnorm_prec simualte parameters from the joint precision matrix derived from a TMB objects
#' @param mu vector of MLE both fixed and random effect parameters
#' @param prec precision matrix, derived from sdreport(obj, getJointPrecision = T)
#' @param n.sims integer number of simulations
#' @param random_seed integer seed
rmvnorm_prec <- function(mu, prec, n.sims, random_seed ) {
  set.seed( random_seed )
  z = matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
  L = Matrix::Cholesky(prec, super=TRUE)
  z = Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
  z = Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
  z = as.matrix(z)
  return(mu + z)
}
#' ALN transformation
#' sum(x) = 1 & length(x) = n -> y length(y) = n - 1
#' The last value is made to be the 'reference' element 
#' @param x vector of compositions doesn't actually have to sum = 1
#' @return additive log ratio transformed variable
alr = function(x) {
  y = log(x[-length(x)] / x[length(x)])
  return(y)
}
#' ALN transformation inverse
#' sum(x) = 1 & length(x) = n -> y length(y) = n - 1
#' @param y vector of compositions that have been alr
#' @return additive log ratio transformed variable
alrinv = function(y) {
  x = c(exp(y), 1)
  x1 = x / sum(x)
  return(x1)
}
#' get_par_vals, from a named list, pull out the element named = param_label
#' 
#' 
get_par_vals = function(x, param_label) {get(x = param_label, pos = x)}
## some code to summarise key outputs
#' Root mean square error of RELATIVE-ERRORS
RMSE = function(lst, param_label = "SSB", true_val) {
  estimated = Reduce(rbind, lapply(lst, FUN = get_par_vals, param_label = param_label))
  Error = sweep(estimated, MARGIN = 2, FUN = "-", true_val)
  R_E = sweep(Error, MARGIN = 2, FUN = "/", true_val)
  Error_sqrd = R_E^2
  sqrt(apply(Error_sqrd, FUN = mean, MARGIN = 2))
}
#' relative error
RE = function(lst, param_label = "SSB", true_val) {
  estimated = Reduce(rbind, lapply(lst, FUN = get_par_vals, param_label = param_label))
  Error = sweep(estimated, MARGIN = 2, FUN = "-", true_val)
  R_E = sweep(Error, MARGIN = 2, FUN = "/", true_val)
  return(R_E)
}
#' MARE Median Absolute Relative error
MARE = function(lst, param_label = "SSB", true_val) {
  estimated = Reduce(rbind, lapply(lst, FUN = get_par_vals, param_label = param_label))
  Error = sweep(estimated, MARGIN = 2, FUN = "-", true_val)
  R_E = sweep(Error, MARGIN = 2, FUN = "/", true_val)
  MA_RE = apply(abs(R_E), FUN = median, MARGIN = 2)
  return(MA_RE)
}
#'
#' Demonstrate different standardisation methods for
#' comparing variables on different scales
#'
#'
# Convert Variables to have the same Lower (0) and Upper Limits (1)
convert_var = function(x) {
  return((x - min(x)) / diff(range(x)))
}

#' MARE Median  Relative error
MRE = function(lst, param_label = "SSB", true_val) {
  estimated = Reduce(rbind, lapply(lst, FUN = get_par_vals, param_label = param_label))
  Error = sweep(estimated, MARGIN = 2, FUN = "-", true_val)
  R_E = sweep(Error, MARGIN = 2, FUN = "/", true_val)
  MA_RE = apply((R_E), FUN = median, MARGIN = 2)
  return(MA_RE)
}

#' MARE Median Relative error
quant_ARE = function(lst, param_label = "SSB", true_val, quantile = c(0.025,0.975)) {
  estimated = Reduce(rbind, lapply(lst, FUN = get_par_vals, param_label = param_label))
  Error = sweep(estimated, MARGIN = 2, FUN = "-", true_val)
  R_E = sweep(Error, MARGIN = 2, FUN = "/", true_val)
  MA_RE = apply((R_E), FUN = quantile, quantile, MARGIN = 2)
  return(MA_RE)
}

#' gmean: geometric mean
#' @param x vector of positive values
#' @return the geometric mean
gmean <- function(x)prod(x)^(1/length(x)) ## geometric mean, ot the best formulation as you can get under and overflow
## an alternative to the above, accoutns for NA's and is more stable for under/overflow
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

#'
#' flip_rows return a matrix x, the has been had its rows reversed i.e. x[1:5,] we return x[5:1,].
#' This is used for plotting image() functions, as they often reformat the matrix
#'
flip_cols = function(x) {
  return(x[,ncol(x):1])
}
#'
#' flip_rows return a matrix x, the has been had its rows reversed i.e. x[1:5,] we return x[5:1,] 
#' This is used for plotting image() functions, as they often reformat the matrix
#'
flip_rows = function(x) {
  return(x[nrow(x):1,])
}


#' return MLE estaimtes for fixed effects
#' @param obj An optimised list that has been build by MakeAdFun
#' 
get_tmb_fixed_effects = function(obj) {
  if (length(obj$env$random) == 0) {
    return(obj$env$last.par.best)
  }
  return(obj$env$last.par.best[-obj$env$random])
}

#' check_tmb_convergence
#' @author C.Marsh
#' @date 
#' @description use TMB object to check gradients and if the hessian is definite positive
#' @param obj An optimised list that has been build by MakeAdFun
#' @param delta Gradient threshold for defining a converged model
#' @depends TMB
#' 
check_tmb_convergence = function(obj, delta = 0.001) {
  best_fixed_eff_pars = get_tmb_fixed_effects(obj)
  grads = tryCatch( expr = obj$gr(best_fixed_eff_pars));
  if (inherits(grads, "error"))
    return("Could not get gradients with an error. generally happens from unconverged models please check")
  if (any(is.na(grads)))
    return("Found gradients with NaN, conclusion = unconverged model")
  labs = names(obj$par)
  if(max(abs(grads)) > delta)
    return(paste0("param with label = ", labs[which.max(abs(grads) > delta)], " has gradient > ", delta, " |grad| = ", round(max(abs(grads)), 5)))
  hess = optimHess(best_fixed_eff_pars, fn = obj$fn, gr = obj$gr)
  sdr = sdreport(obj, getJointPrecision = TRUE)
  if (is.character(try(chol(sdr$cov.fixed), silent = TRUE)))
    return("Covariance of fixed effect no positive Definitive")
  if ("jointPrecision" %in% names(sdr)) {
    if (is.character(try(chol(sdr$jointPrecision), silent = TRUE)))
      return("Joint Precision not positive Definitive")
  }
  return("No evidence of non convergence")
}
  
#' TMB helper function
#' @author C.Marsh
#' @date 5/9/2018
#' @description this function returns a list of factors used in the map argument of the MakeADFun function
#' @param par_list a named list that you give to the par argument in the MakeADFun
#' @param pars_to_exclude a vector of strings with names of parmeters you want to FIX in the objective object.
#' @param vec_pars_to_adjust a vector string of parameter labels that we want to exclude certain elements.
#' @param vec_elements_to_exclude a list with number of elements = length(vec_pars_to_adjust). each list element 
#' contains a vector of elements that we want to exclude from estimation.
#' @return a list of factors used in the MakeADFun function
fix_pars = function(par_list, pars_to_exclude, vec_pars_to_adjust = NULL, vec_elements_to_exclude = NULL) {
  if (!any(pars_to_exclude %in% names(par_list))) {
    stop(paste0("The parameters ", paste(pars_to_exclude[!pars_to_exclude %in% names(par_list)],collapse = " ")," in exclusion parameters could not be found in the 'par_list', please sort this out"))
  }
  pars = names(par_list)
  mapped_pars = list();
  tailor_vectors = FALSE
  if (!is.null(vec_pars_to_adjust)) {
    tailor_vectors = TRUE
    if (!all(vec_pars_to_adjust %in% pars_to_exclude))
      stop("parmaeters noted in vec_pars_to_adjust, need to also be in pars_to_exclude")
  }
  param_factor = 1;
  for(i in 1:length(pars)) {
    if (pars[i] %in% pars_to_exclude) {
      params_in_this_par = par_list[[pars[i]]];
      if (tailor_vectors & (pars[i] %in% vec_pars_to_adjust)) {
        include_element_index = c(1:length(params_in_this_par))[-vec_elements_to_exclude]
        params_vals = factor(rep(NA, length(params_in_this_par)), levels = factor(param_factor:(param_factor + length(include_element_index) - 1)))
        params_vals[include_element_index] = factor(param_factor:(param_factor + length(include_element_index) - 1))#, levels = factor(include_element_index))
        param_factor = param_factor + length(include_element_index)
        mapped_pars[[pars[i]]] = params_vals;
      } else {
        mapped_pars[[pars[i]]] = rep(factor(NA),length(params_in_this_par));
      }
    } else {
      params_in_this_par = par_list[[pars[i]]];
      params_vals = factor(param_factor:(param_factor + length(params_in_this_par) - 1))
      param_factor = param_factor + length(params_in_this_par)
      mapped_pars[[pars[i]]] = params_vals
    }
  }
  return(mapped_pars);
}

#' eigen_decomp_covariance Do an eigen decomposition to look at poorly estimated parameters from MLE fit
#' @param covariance_matrix symetric covariance matrix
#' @param param_labels vector of param labels (optional)
#' @param delta a cut off value for 'poorly' defined parameters.
#' @return: data frame of eiegen values for the matrix and index of good and bad parameters based on delta
#'
eigen_decomp_covariance = function(covariance_matrix, param_labels = NULL, delta = .Machine$double.eps) {
  ## check covariance is invertable
  if (!isSymmetric(covariance_matrix))
    stop("covariance matrix is not symetric something is wrong here.")
  ## check positive semi defintie matrix
  if(class(try(solve(covariance_matrix),silent=T)) != "matrix")
    stop("covariance not invertible")
  ## calculate hessian
  hess = solve(covariance_matrix)
  ## eigen decomposition
  Eig = eigen(hess)
  WhichBad = which(Eig$values < sqrt(delta))
  df = NULL;
  if (is.null(param_labels)) 
    param_labels = as.character(1:ncol(covariance_matrix))
  
   
  if (length(WhichBad) == 0) {    
    message( "All parameters are identifiable" )
  } else {
    # Check for parameters
    RowMax = apply( Eig$vectors[,WhichBad,drop=FALSE], MARGIN=1, FUN=function(vec){max(abs(vec))} )
    df = data.frame("Param"=param_labels, "eigenvalues", Eig$values ,"Param_check"=ifelse(RowMax>0.1, "Bad","OK"))
  }
  return(df)
}

## get_time_var get a time-series parameter from mcmc data frame for plotting return as DF that is ggplot applicable
#' @param DF an object that is created from join_chains(as.array(stanfit)) 
#' @param long_format return data from an long-format = T, 
#' @param var param var to be in the format of 'var[timestep]' in the colnames. e.g. var = "expected_observations"  colnames(DF) = ... "expected_observations[1]", "expected_observations[2]" etc.
#' @import(reshape2)
#' @return dataframe in long format of the variable with time-ndx
get_time_var =  function(DF, var, long_format = TRUE) {
  if (all(!grepl(var, colnames(DF)))) {
    stop("could not find the var pattern in the colnames in the data frame provided please sort this out")
  }
  df_of_interest = DF[,grepl(var, colnames(DF))]
  if (!long_format)
    return(df_of_interest)
  # strip out [] or ()
  colnames(df_of_interest) = sub("[^(]+\\[([^)]+)\\].*", "\\1", colnames(df_of_interest), perl=TRUE)
  long_df = melt(df_of_interest)
  colnames(long_df) = c("time","Var");
  return(long_df)
}

# multi_ci Utility funciton for calculating upper and lower bounds for a multinomial value "Normal Approximation Method" of the Binomial Confidence Interval:
#' @param P vector of proportions
#' @param N = sample size int
#' @param alpha significance level
#' @return upper and lower confidence limits
multi_ci = function(P, N, alpha = 0.05) {
  SE = sqrt((P * (1 - P)) / N);
  half_z = abs(qnorm(alpha / 2)); ## symmetric so abs is fine
  return(list = c(U_CI = P + half_z * SE, L_CI = P - half_z * SE))
}
## useful function for STAN
## borrowed from https://github.com/betanalpha/knitr_case_studies/blob/master/qr_regression/stan_utility.R
# Check transitions that ended with a divergence
check_div <- function(fit, sample_tmb = FALSE) {
  sampler_params <- NULL
  if(sample_tmb) {
    sampler_params <- fit[["sampler_params"]]
    for(i in 1:length(sampler_params))
      sampler_params[[i]] = sampler_params[[i]][(fit$warmup + 1):nrow(sampler_params[[i]]), ]
    
  } else {
    sampler_params <- rstan::get_sampler_params(fit, inc_warmup=FALSE)
  }
  divergent <- do.call(rbind, sampler_params)[,'divergent__']
  n = sum(divergent)
  N = length(divergent)
  
  return(sprintf('%s of %s iterations ended with a divergence (%s%%)',
                 n, N, 100 * n / N))
  if (n > 0)
    return('  Try running with larger adapt_delta to remove the divergences')
}

# Check transitions that ended prematurely due to maximum tree depth limit
check_treedepth <- function(fit, max_depth = 10) {
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
  treedepths <- do.call(rbind, sampler_params)[,'treedepth__']
  n = length(treedepths[sapply(treedepths, function(x) x == max_depth)])
  N = length(treedepths)
  
  return(sprintf('%s of %s iterations saturated the maximum tree depth of %s (%s%%)',
                 n, N, max_depth, 100 * n / N))
  if (n > 0)
    return('  Run again with max_depth set to a larger value to avoid saturation')
}

# Checks the energy Bayesian fraction of missing information (E-BFMI)
check_energy <- function(fit) {
  sampler_params <- rstan::get_sampler_params(fit, inc_warmup=FALSE)
  no_warning <- TRUE
  for (n in 1:length(sampler_params)) {
    energies = sampler_params[n][[1]][,'energy__']
    numer = sum(diff(energies)**2) / length(energies)
    denom = var(energies)
    if (numer / denom < 0.2) {
      print(sprintf('Chain %s: E-BFMI = %s', n, numer / denom))
      no_warning <- FALSE
    }
  }
  if (no_warning)
    return('E-BFMI indicated no pathological behavior')
  else
    return('  E-BFMI below 0.2 indicates you may need to reparameterize your model')
}

# Checks the effective sample size per iteration
check_n_eff <- function(fit, pars) {
  fit_summary <- summary(fit, probs = c(0.025, 0.25, 0.50, 0.75, 0.975))$summary
  
  row_ndx = vector()
  for(i in 1:length(pars))
    row_ndx = c(row_ndx, which(grepl(rownames(fit_summary), pattern = pars[i])))
  
  ## number of params
  N <- dim(fit_summary[row_ndx,])[[1]]
  
  iter <- dim(rstan::extract(fit)[[1]])[[1]]
  
  no_warning <- TRUE
  for (n in 1:N) {
    ratio <- fit_summary[row_ndx[n],"n_eff"] / iter
    if (ratio < 0.2) {
      return(sprintf('n_eff / iter for parameter %s is %s!',
                     rownames(fit_summary[row_ndx,])[n], ratio))
      no_warning <- FALSE
    }
  }
  if (no_warning)
    return('n_eff / iter looks reasonable for all parameters')
  else
    return('  n_eff / iter below 0.2 indicates that the effective sample size has likely been overestimated')
}

# Checks the potential scale reduction factors
check_rhat <- function(fit, pars) {
  fit_summary <- rstan::summary(fit, probs = c(0.025, 0.25, 0.50, 0.75, 0.975))$summary
  
  row_ndx = vector()
  for(i in 1:length(pars))
    row_ndx = c(row_ndx, which(grepl(rownames(fit_summary), pattern = pars[i])))
  
  N <- dim(fit_summary[row_ndx,])[[1]]
  
  no_warning <- TRUE
  for (n in 1:N) {
    rhat <- fit_summary[row_ndx[n],"Rhat"]
    if (rhat > 1.05 || is.infinite(rhat) || is.nan(rhat)) {  ## changed based on recent research found at https://arxiv.org/pdf/1903.08008.pdf
      return(sprintf('Rhat for parameter %s is %s!',
                     rownames(fit_summary[row_ndx,])[n], rhat))
      no_warning <- FALSE
    }
  }
  if (no_warning)
    return('Rhat looks reasonable for all parameters')
  else
    return('  Rhat above 1.05  indicates that the chains very likely have not mixed')
}

check_all_diagnostics <- function(fit, pars) {
  if (check_n_eff(fit, pars) != "n_eff / iter looks reasonable for all parameters")
    return("failed: check_n_eff");
  if (check_rhat(fit, pars) != "Rhat looks reasonable for all parameters")
    return("failed: check_rhat");
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
  divergent <- do.call(rbind, sampler_params)[,'divergent__']
  n = sum(divergent)
  N = length(divergent)
  if (check_div(fit) != sprintf('%s of %s iterations ended with a divergence (%s%%)',n, N, 100 * n / N)) 
    return("failed: check_div");
  max_depth = 10
  treedepths <- do.call(rbind, sampler_params)[,'treedepth__']
  n = length(treedepths[sapply(treedepths, function(x) x == max_depth)])
  N = length(treedepths)
  if (check_treedepth(fit) != sprintf('%s of %s iterations saturated the maximum tree depth of %s (%s%%)', n, N, max_depth, 100 * n / N))
    return("failed: check_treedepth");
  if (check_energy(fit) != "E-BFMI indicated no pathological behavior")
    return("failed: check_energy");
  return("passed")
}

# Returns parameter arrays separated into divergent and non-divergent transitions
partition_div <- function(fit) {
  nom_params <- extract(fit, permuted=FALSE)
  n_chains <- dim(nom_params)[2]
  params <- as.data.frame(do.call(rbind, lapply(1:n_chains, function(n) nom_params[,n,])))
  
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
  divergent <- do.call(rbind, sampler_params)[,'divergent__']
  params$divergent <- divergent
  
  div_params <- params[params$divergent == 1,]
  nondiv_params <- params[params$divergent == 0,]
  
  return(list(div_params, nondiv_params))
}

#' @description this function takes all chains and returns a single data frame
#' @param params_array is an array that is generated from as.array of a 'stanfit' object
join_chains = function(params_array) {
  n_chains = dim(params_array)[2]
  merged_df = as.data.frame(apply(params_array, MARGIN = c(3),FUN = rbind))
  merged_df$chain_id = NA
  merged_df$sample_id = NA
  
  row_counter = 1;
  for(i in 1:n_chains) {
    merged_df$chain_id[row_counter:(row_counter + nrow(params_array[,i,]) - 1)] = i;
    merged_df$sample_id[row_counter:(row_counter + nrow(params_array[,i,]) - 1)] = 1:nrow(params_array[,i,])
    
    row_counter = row_counter + nrow(params_array[,i,])
    
  }
  return(merged_df)
}

#' @description this function takes all chains and returns a single data frame for a JAGs model
#' @param params_array is an mcmc.list that is generated from coda.samples or jags.samples
merge_jags = function(params_array) {
  n_chains = length(params_array)
  n_samples = nrow(params_array[[1]])
  
  merged_df = NULL;
  
  row_counter = 1;
  for(i in 1:n_chains) {
    temp_df = as.data.frame(params_array[[i]])
    temp_df$chain_id = i
    temp_df$sample_id = 1:n_samples;
    merged_df = rbind(merged_df, temp_df)
  }
  return(merged_df)
}
#' utiltiy functions to pull color hues gray scale I am guessing
#' @param n number of color shades
color_vector <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=50, c=50)[1:n]
}
#' utiltiy functions to pull color hues gray scale I am guessing
#' @param n number of color shades
color_vector_chain <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=80, c=50)[1:n]
}

#' my own trace plot using ggplot, so we can add the true value during a simulation study
#' A trace plot that I tweaked from rstan that allows me to add the true value on 
#' the plot, helps diagnose use stan_trace unless you want to view the true value on the trace plot
#' @param DF that is a consequence of join_chains() function
#' @param pars string vector of parameters you want to plot
#' @param true_vals vector of true values, one to one relationship with pars
#' @param FUN if you want the mean, median or some other function as well
#' @param nrow rows in final plots same as cols
#' @param true_col cols of true value.
gg_trace = function(DF, pars, true_vals = NULL, FUN = median, nrow = 1, ncol = length(pars), true_col = "blue", true_size = 1) {
  if(!is.null(true_vals)) {
    if (length(pars) != length(true_vals))
      stop("need to supply a vector of the same length for pars and true vals.");
  }
  ## remove vector brackets 
  ndx = strsplit(colnames(DF),"\\[")
  labs = Reduce(c, lapply(ndx, function(x) {x[1]}))
  if (any(!pars %in% colnames(DF))) {
    ## check that it isn't a vector parameter
    if (any(!pars %in% labs)) {
      stop(paste0("could not find 'pars' ", pars[pars %in% labs], " in column name of DF please sort out"))
    } else {
      cat("if you are supplying vector or matrix, you need to supply dimensions in 'pars' parameters, otherwise I am ignoring it")
    }
  }
  df_of_interest =  DF[,colnames(DF) %in% c(pars,"chain_id","sample_id")]
  df_melted = melt(df_of_interest, id.vars = c("chain_id","sample_id"), value.name = "value",variable.name = "parameter")
  
  ## keep colors the same as stan
  cahin_colors = rgb(matrix(c(230, 97, 1,
                              153, 142, 195,
                              84, 39, 136,
                              241, 163, 64,
                              216, 218, 235,
                              254, 224, 182),
                            byrow = TRUE, ncol = 3),
                     names = paste(1:6), maxColorValue = 255)
  clrs = cahin_colors[1:max(df_melted$chain_id)]
  
  base = ggplot(data = df_melted,
                aes(x = sample_id, y = value, col = as.factor(chain_id)))
  
  
  graph <-
    base +
    geom_path() +
    scale_color_manual(name = "chain", values = clrs) +
    labs(x="",y="");
  
  if(!is.null(true_vals)) {
    plot_lines = data.frame("parameter" = pars, truth = true_vals)
    fun_value = tapply(df_melted$value, df_melted$parameter, FUN = FUN)
    fun_value = c(fun_value[3], fun_value[1],fun_value[2])
    
    plot_lines$fun_value = fun_value[match(pars,names(fun_value))]
    max_sample = max(df_melted$sample_id)
    
    ## add to graph
    graph = graph +
      scale_linetype_manual(values = c("True value" = "solid", "Point estimate" = "dashed")) +
      geom_segment(data = plot_lines, aes(x = 1, xend = max_sample,y = truth, yend = truth), size = true_size, col = true_col, linetype = "solid") + 
      geom_segment(data = plot_lines, aes(x = 1, xend = max_sample,y = fun_value, yend = fun_value), size = true_size, col = "black", linetype = "dashed")
  }
  
  if (length(pars) == 1)
    graph <- graph + ylab(pars)
  else
    graph <- graph + facet_wrap(~parameter, nrow = nrow, ncol = ncol, scales = "free")
  
  print(graph);
}
#' General function "rqr" for calculating Randomized Quantile Residual (RQRs)
#' from the thesis https://math.usask.ca/~longhai/researchteam/theses/BAI-THESIS-2018.pdf
#' @param object GLM functions glmmTMB
#' 
rqr <- function(object) {
  family<-family(object)$family
  mu<-predict(object,zitype="conditional")
  size<-sigma(object)
  p <- predict(object,zitype="zprob")
  n <- object$modelInfo$nobs
  y <- object$frame[,1]
  dzpois <- function(x,lambda,p) {
    return((1-p)*dpois(x,lambda)+p*(x==0))
  }
  pzpois <- function(x,lambda,p) {
    return((1-p)*ppois(x,lambda)+p*(x>=0))
  }
  dznbinom <- function(x,size,mu,p){
    return((1-p)*dnbinom(x,size = size, mu = mu)+p*(x==0))
  }
  pznbinom <- function(x,size,mu,p) {
    return((1-p)*pnbinom(x,size = size, mu = mu)+p*(x>=0))
  }
  if(1*(object$modelInfo$allForm$ziformula==~0)==1) {##no zero-inflation
    if(family[1]=="gaussian") {
      pvalue=pnorm(y,mu,size)
    }
    else if(family[1]=="poisson") {
      pvalue = ppois(y-1,mu) + dpois(y,mu) * runif(n)
    }
    else if(family[1]=="nbinom2") {
      pvalue=pnbinom (y-1,size=size, mu=mu) + dnbinom (y,size = size, mu = mu) * runif(n)
    }
  } else if(1*(object$modelInfo$allForm$ziformula==~0)==0) {##zero-inflation
    if(family[1]=="poisson") {
      pvalue=pzpois(y-1,mu,p) + dzpois (y, mu,p) * runif(n)
    } else if(family[1]=="nbinom2") {
      pvalue=pznbinom (y-1,size,mu,p) + dznbinom (y,size,mu,p) * runif(n)
    }
  }
  pvalue<-pmin(pmax(pvalue,10^{-10}),1-10^{-10})
  qnorm=qnorm(pvalue)
  test=shapiro.test(qnorm)$p.value
  list(pvalue=pvalue,qnorm=qnorm,test=test)
}
###General function "rqrhurdle" for calculating RQRs
rqrhurdle<-function(model_count, model_zero, data) {
  library(actuar)
  name.y <- names(model_count$frame)[[1]]
  y <- data[,name.y]
  m<-model_zero$modelInfo$nobs
  family<-family(model_count)$family
  mu<-predict(model_count,newdata=data,zitype="conditional")
  size<-sigma(model_count)
  ## success rate in negnative binomial
  prob_nb<-size/(size+mu)
  ## probability of 0
  #pi<- predict(model_zero,zitype="zprob")
  pi<- 1-predict(model_zero,newdata = data, type="response")
  if(family[1]=="truncated_poisson") {
    pvalue=pzmpois(y-1,mu,pi)+dzmpois(y,mu,pi)* runif(m)
  } else if(family[1]=="truncated_nbinom2") {
    pvalue=pzmnbinom(y-1,size,prob_nb,pi)+dzmnbinom(y,size,prob_nb,pi)*runif(m)
  }
  pvalue<-pmin(pmax(pvalue,10^{-10}),1-10^{-10})
  qnorm=qnorm(pvalue)
  test=shapiro.test(qnorm)$p.value
  list(pvalue=pvalue,qnorm=qnorm,test=test)
}
##functions for calculating pearson residuals
pearson<-function(object) {
  family<-family(object)$family
  mu<-predict(object,zitype="conditional")
  size<-sigma(object)
  p0<-predict(object,zitype="zprob")
  y<-object$frame[,1]
  if(1*(object$modelInfo$allForm$ziformula==~0)==1) {##no zero-inflation
    pearson<-residuals(object,type="pearson")
  } else if(1*(object$modelInfo$allForm$ziformula==~0)==0) {##zero-inflation 
    if(family[1]=="poisson") {
      mu_hat<-mu*(1-p0)
      var_hat<-mu_hat*(1+p0*mu)
      pearson<-(y-mu_hat)/sqrt(var_hat)
    } else if(family[1]=="nbinom2") {
      mu_hat<-mu*(1-p0)
      var_hat<-mu_hat*(1+mu/size)+mu^2/(p0^2+p0)
      pearson<-(y-mu_hat)/sqrt(var_hat)
    }
  }
  shapiro<-shapiro.test(pearson)$p.value
  list(pearson=pearson,shapiro=shapiro)
}
pearsonhurdle<-function(model_count, model_zero, data) {
  name.y <- names(model_count$frame)[[1]]
  y <- data[,name.y]
  m<-model_zero$modelInfo$nobs
  family<-family(model_count)$family
  mu<-predict(model_count,newdata=data,zitype="conditional")
  size<-sigma(model_count)
  p0<-1-predict(model_zero,newdata = data, type="zprob")
  prob_nb<-size/(size+mu)
  if(family[1]=="truncated_poisson")
  {
    mu_hat<-(1-p0)*mu/(1-exp(-mu))
    var_hat<-mu_hat*(1+mu)-mu_hat^2
    pearson<-(y-mu_hat)/sqrt(var_hat)
  } else if(family[1]=="truncated_nbinom2") {
    mu_hat<-(1-p0)*mu/(1-prob_nb^size)
    var_hat<-mu_hat*(mu+1+mu/size)-mu_hat^2
    pearson<-(y-mu_hat)/sqrt(var_hat)
  }
  shapiro<-shapiro.test(pearson)$p.value
  list(pearson=pearson,shapiro=shapiro)
}

#' DistanceLongLat: calculate distance between a set of coordinates (as the crow flies)
#' @param long1 longitude first point decimal degrees > 180 are WESTINGS
#' @param long2 longitude second point decimal degrees > 180 are WESTINGS
#' @param lat1 lattitude first point decimal degrees < 0 Southings
#' @param lat2 lattitude second point decimal degrees < 0 Southings
#' @param metres bool, F = nautical miles, T = metres
#' @return a distance
DistanceLongLat<-function(long1, long2, lat1, lat2, metres = F) {
# Inputs start=(long1,lat1) and end=(long2,lat2) in decimal degrees
# OR assumes that locator is used to define exactly 2 points
#
# Assumes longitude numbers are positive and that numbers > 180 are WESTINGS
# Assumes lattitude negative numbers are SOUTHINGS
# Outputs nautical miles if metres=F, else distance in metres
	if(missing(long1) | missing(long2) | missing(lat1) | missing(lat2)) {
		cat("Using function \"locator(2)\" to locate end points\n")
		x <- nz.locator(2)
		long1 <- x$x[1]
		long2 <- x$x[2]
		lat1 <- x$y[1]
		lat2 <- x$y[2]
		print(unlist(c(long1, long2, lat1, lat2)))
	}
	long1 <- (long1 * pi)/180
	long2 <- (long2 * pi)/180
	lat1 <- ifelse(lat1 > 180, 360 - lat1,  - lat1)
	lat2 <- ifelse(lat2 > 180, 360 - lat2,  - lat2)
	lat1 <- (lat1 * pi)/180
	lat2 <- (lat2 * pi)/180
	d <- 2 * asin(sqrt((sin((lat1 - lat2)/2))^2 + cos(lat1) * cos(lat2) * (sin((long1 - long2)/2))^2))
	nm <- (d * 180 * 60)/pi
	if(metres)
		nm <- nm * 1.852 * 1000
	return(nm)
}
#'
#' Calculate the cononical index advocated for in francis 1999
#' theory is in the CPUE folder, basically given GLM coeffecients it gives geometrically scaled index termed the canonical
#' @param GLM fitted glm object
#' @param years years to calcuate index
#' @param base <integer> An index (1:number of factors in year.name) which indicates the reference year (generally the intercept)
#' @param year.name lable for the year coeffecient (factor)
#'
canonical.index = function(GLM, years, base = 1, year.name = "year") {
  res <- as.data.frame(GLM$coefficients)
  index <- regexpr(year.name, row.names(res)) > 0
  X <- res[index, 1]
  V <- summary(GLM)$cov.unscaled[index, , drop = FALSE][, index]
  n <- length(X) + 1
  A <- matrix(-1/n, n, n - 1)
  # a trick to get the geometric mean, looks like they are dropping a level but not.
  A[-base, ] <- A[-base, ] + diag(rep(1, n - 1))
  CPUE <- A %*% X
  COV <- sqrt(diag((A %*% V) %*% t(A)))
  index <- exp(CPUE)
  lower.CI <- exp(CPUE - 2 * COV)
  upper.CI <- exp(CPUE + 2 * COV)
  ans <- data.frame(years, index, lower.CI, upper.CI)
  ans$se <- COV
  ans$cv <- sqrt(exp(COV^2) - 1)
  return(ans)
} 


#' Von Bert mean length at age relationship
#' @param age ages
#' @param K rate parameter
#' @param length at age 0
#' @param L_inf asympototic mean length at age
#' @return mean length at age
VB<- function(age,K,L_inf,t0) {
 return(L_inf * (1-exp(-K*(age -t0))))
}
#' schnute mean length at age relationship
#' @param age ages
#' @param y1 mean length at reference age t1
#' @param y2 mean length at reference age t2
#' @param t1 reference age for y1
#' @param t2 reference age for y2
#' @param a schnute growth params
#' @param b schnute growth params
#' @return mean length at age
schnute <- function(y1,y2,t1,t2,a,b,age) {
  if(!(a == 0) & !(b == 0)) {
    (y1^b + (y2^b - y1^b) * ((1-exp(-a*(age-t1)))/(1-exp(-a*(t2-t1)))))^(1/b)
  } else if(!(a == 0) & b == 0){
    (y1*exp(log(y2/y1) * ((1-exp(-a*(age-t1)))/(1-exp(-a*(t2-t1))))))
  } else if(a == 0 & !(b == 0)) {
    ((y1^b + (y2^b - y1^b)* (age-t1)/(t2-t1))^(1/b))
  } else {
    (y1*exp(log(y2/y1)*(age-t1)/(t2-t1)))
  }
}

#' bound_unit constrains Y from -inf -> inf to be between -1 -> 1 
#' @param Y scalar range [-inf, inf]
#' @return X to be between [-1,1] 
bound_unit = function(Y) {
  return(Y / sqrt(1.0 + Y * Y))
}

#' inv_bound_unit constrains Y from -inf -> inf to be between -1 -> 1 
#' @param X scalar range [-1,1] 
#' @return Y to be between [-inf, inf]
inv_bound_unit = function(X) {
  return(sqrt((X*X) / (1 - X*X)) * ifelse(X < 0, -1, 1))
}
#' logit bounds X which is between 0-1 to -inf -> inf based on the logit transformation
#' equivalent to qlogis(X)
#' @param X scalar range [0,1]
#' @return Y to be between [-inf, inf]
logit = function(X) {
  log(X / (1 - X)) 
}
#' invlogit Inverse logit transformation, equivalent to plogis(Y)
#' @param Y scalar between [-inf, inf] 
#' @return X between [0,1]
invlogit<- function(Y) {
  1/(1 + exp(-Y))
}

#' logit_general bounds X which is between [lb,ub] to -inf -> inf based on the logit transformation
#' @param X scalar range [lb,ub]
#' @param ub upper bound for X
#' @param lb lower bound for X
#' @return Y to be between [-inf, inf]
logit_general = function(X, lb, ub) {
  X1 = (X - lb) / (ub - lb)
  log(X1/(1 - X1))
}
#' invlogit_general bounds X which is between -inf -> inf to [lb,ub] based on the logit transformation
#' @param Y scalar range [-inf, inf]
#' @param ub upper bound for X
#' @param lb lower bound for X
#' @return X to be between [lb,ub]
invlogit_general = function(Y, lb, ub) {
  Y1 = 1 / (1 + exp(-Y))
  lb + (ub - lb)*Y1
}
#' create a function for simulating starting values for model robustness
#' @param n <integer> number of draws you want
#' @param dist <string> distribution to use to draw starting values, allowed "unif", "norm" 
#' @param LB <scalar> Lower bound of starting value, shouldn't use lower bound of parameter
#' @param UB <scalar> Upper Bound of starting balues.
#' @param dist_pars <vector> has distributional parameters, i.e for norm = c(mu, sigma)
#' @return vector of starting values
ran_start = function(n = 100, dist = "unif", LB = -Inf, UB = Inf, dist_pars = NULL) {
  if (!any(dist %in% c("unif", "norm")))
    stop("parameter 'dist', needs to be either unif or norm")
  start_vals = vector();
  if (dist == "unif") {
    start_vals = runif(n, LB, UB)
  } else if (dist == "norm") {
    start_vals = rnorm(n, dist_pars[1],dist_pars[2])
    start_vals[start_vals < LB] = LB
    start_vals[start_vals > UB] = UB
    
  } else {
    stop("something went wrong")
  }
  return(start_vals)
}

#' deg2decdeg convert from degrees degrees:minutes:seconds -> decimal degrees
#' @param x <scalar> format degrees:minutes:seconds
#' @return dec_degrees
#' @example
#' deg2decdeg(paste(-43,50,23,sep = ":"))
deg2decdeg = function(x) {
  sp_x = strsplit(x, split = ":")
  class(sp_x[[1]]) = "numeric"
  return( (abs(sp_x[[1]][1]) + (sp_x[[1]][2] / 60) + (sp_x[[1]][3] / 3600)) * ifelse(sp_x[[1]][1] < 0, -1, 1))
}
## vectorise the function
deg2decdeg_v = Vectorize(FUN = deg2decdeg)

#' decdeg2deg convert from  decimal degrees -> degrees degrees:minutes:seconds
#' @param x <scalar> format decimal degrees
#' @return degrees:minutes:seconds
#' @example
#' deg2decdeg(paste(-43,50,23,sep = ":"))
decdeg2deg = function(x) {
  x1 = abs(x)
  deg = floor(x1)
  min =  floor((x1 - deg) * 60)
  sec = floor(((x1 - deg) * 60 - min) * 60)
  return(paste(deg * ifelse(x < 0, -1, 1), min, sec, sep = ":"))
}
## vectorise the function
decdeg2deg_v = Vectorize(FUN = decdeg2deg)

#' long2UTM returns a zone for UTM based on longitude. Expects long can be either [-180, 180] or [0,360]
#' @param long <scalar> longitude, to get UTM zone
#' @return UTM zone
long2UTM <- function(long) {
  ## Function to get the UTM zone for a given longitude
  (floor((long + 180)/6) %% 60) + 1
}

#' LongLatToUTM takes a data frame whith lat and long in decimal degree format. Then
#' returns a spatial dataframe that is in UTM format.
#' @param df <data.frame> that has at least two columns labelled long lat
#' @param zone <integer> zone for UTM.
#' @return SpatialDataFrame
#' @depends library(sp)
LongLatToUTM <- function(df, zone = NULL){
  if(!all(c("lat","long") %in% colnames(df)))
    stop("LongLatToUTM() df, requires the columns 'lat' & 'long' to be present.")
  ## Args: df, data frame must have x and y columns. Should be from same UTM zone.
  ## Create a spatial dataframe
  coordinates(df) <- ~ long + lat
  proj4string(df) <- CRS("+proj=longlat +datum=WGS84")  
  
  ## Get zones for from a range of longitudes of data.
  ## Take the most represenative zone of the data.
  if (is.null(zone)) {
    longs = seq(from = min(df$long, na.rm = T), to = max(df$long, na.rm = T), by = 0.1)
    zones <- long2UTM(longs)
    tab_zones = table(zones)
    zone = as.numeric(names(tab_zones)[which.max(tab_zones)])
  }
  
  if (length(unique(zone)) > 1) {
    ## take the zone with the most 
    stop("values from different UTM zones")
  }

  ## Change CRS of the spatial data frame and convert to data frame
  res <- spTransform(df, CRS(paste0("+proj=utm +zone=", zone, "+datum=WGS84")))
  return(res)
}
