#####################
## Update Data
#####################
# load TMB data
load(file = file.path(savedobjs_dir, "First_tmb_dir.RData"))
if(debug)
  cat("year = ", current_year, "\n")

years_diff = abs(max(TMB_data$years) - current_year)
TMB_data$years = min(TMB_data$years):current_year
TMB_data$nyears = length(TMB_data$years)
TMB_data$ycs_estimated = rep(1, TMB_data$nyears)
TMB_data$ycs_estimated[(TMB_data$nyears - 1):TMB_data$nyears] = 0 ## don't estimate the last two ycs
if(years_diff > 0) {
  for(ndx in 1:years_diff) {
    TMB_data$catchMeanLength = cbind(TMB_data$catchMeanLength, TMB_data$catchMeanLength[,ncol(TMB_data$catchMeanLength)])
    TMB_data$stockMeanLength = cbind(TMB_data$stockMeanLength, TMB_data$stockMeanLength[,ncol(TMB_data$stockMeanLength)])
    TMB_data$propMat = cbind(TMB_data$propMat, TMB_data$propMat[,ncol(TMB_data$propMat)])
    TMB_data$propZ_ssb = c(TMB_data$propZ_ssb, TMB_data$propZ_ssb[length(TMB_data$propZ_ssb)])
  }
}
#cat("extension = ", extension,"\n")
#cat("current year = ", current_year,"\n")
fishingproc = extract.run(file = file.path(abm_sim_dir, paste0("FishingProc", extension)))
sim_survey_bio = extract.ibm.file(file = file.path(abm_sim_dir, paste0("chatTANbiomass", extension)), quiet  = T)
sim_survey_age = extract.ibm.file(file = file.path(abm_sim_dir, paste0("chatTANage", extension)), quiet  = T)
sim_fishery_age = extract.ibm.file(file = file.path(abm_sim_dir, paste0("FishingAge", extension)), quiet  = T)

# year indicators, I bit niggly as years are inputs so won't match simulated obs
TMB_data$survey_year_indicator = as.integer(TMB_data$years %in% as.numeric(sim_survey_bio[[1]]$years$value[1:length(sim_survey_bio[[1]]$obs$value)]))
TMB_data$survey_sample_time = rep(0.5, sum(TMB_data$survey_year_indicator))
TMB_data$propZ_survey = rep(0.5, sum(TMB_data$survey_year_indicator))
## implement obs
TMB_data$survey_obs = as.numeric(sim_survey_bio[[1]]$obs$value)
TMB_data$survey_cv = as.numeric(sim_survey_bio[[1]]$error_value$value)
## implement obs
SAge_sim_ibm = Reduce(rbind, sim_survey_age[[1]]$Table$obs)
SAge_N_eff = Reduce(rbind, sim_survey_age[[1]]$Table$error_values)
class(SAge_sim_ibm) <- "numeric"
class(SAge_N_eff) <- "numeric"
## convert to numbers? for dirichlet
SAge_prop = sweep(SAge_sim_ibm, MARGIN = 1, STATS = rowSums(SAge_sim_ibm), "/")

FAge_sim_ibm = Reduce(rbind, sim_fishery_age[[1]]$Table$obs)
class(FAge_sim_ibm) <- "numeric"
FAge_N_eff = Reduce(rbind, sim_fishery_age[[1]]$Table$error_values)
class(FAge_N_eff) <- "numeric"
TMB_data$fishery_year_indicator = as.integer(TMB_data$years %in% as.numeric(sim_fishery_age[[1]]$years$value[1:nrow(FAge_N_eff)]))

## truncate
FAge_trunc = FAge_sim_ibm[,TMB_data$fishery_min_age:TMB_data$fishery_max_age]
FAge_trunc[,ncol(FAge_trunc)] = FAge_trunc[,ncol(FAge_trunc)] + rowSums(FAge_sim_ibm[,TMB_data$fishery_max_age:ncol(FAge_sim_ibm)])
FAge_prop = sweep(FAge_trunc, MARGIN = 1, STATS = rowSums(FAge_trunc), "/")

## fishing obs
obs_catch_ndx = grepl(names(fishingproc$Fishing), pattern = "actual_catches")
obs_catch = Reduce(c, lapply(fishingproc$Fishing[obs_catch_ndx], FUN = function(x){as.numeric(x[1])}))
TMB_data$catches = obs_catch[1:TMB_data$nyears]

TMB_data$survey_trans_comp = t(SAge_prop * init_age_comp_sample_size) #SAge_N_eff)
TMB_data$fishery_trans_comp = t(FAge_prop * init_age_comp_sample_size) #FAge_N_eff[,1])
## Update random start fucntion
## params like ycs and F's will have to be adapted with time
# start_vals
start_vals = function() {
  vals = 
    list(
      logit_f_a50 = logit_general(runif(1,3,9), TMB_data$sel_a50_bounds[1], TMB_data$sel_a50_bounds[2]),
      logit_f_ato95 = logit_general(runif(1,3,9), TMB_data$sel_ato95_bounds[1], TMB_data$sel_ato95_bounds[2]),
      logit_survey_a50 = logit_general(runif(1,3,9), TMB_data$sel_a50_bounds[1], TMB_data$sel_a50_bounds[2]),
      logit_survey_ato95 = logit_general(runif(1,3,9), TMB_data$sel_ato95_bounds[1], TMB_data$sel_ato95_bounds[2]),
      ln_extra_survey_cv = log(0.001),
      trans_survey_error = log(runif(1,0.1,0.4)),
      trans_fishery_error = log(runif(1,0.1,0.4)),
      logit_surveyQ = qlogis(runif(1,0.05,0.2)),
      ln_F = log(runif(length(TMB_data$catches),0.02,0.2)),
      ln_R0 = runif(1,14.5,16),
      ln_catch_sd = log(0.02),
      ln_sigma_r = log(runif(1,0.2,0.7)),
      ln_ycs_est = log(exp(rnorm(n = sum(TMB_data$ycs_estimated), 0, 0.2)))
    )
  return(vals)
}

####################
## Now estimate
####################
init_vals = start_vals()
map_fix_pars = fix_pars(par_list = init_vals, pars_to_exclude = c("ln_catch_sd", "trans_survey_error","trans_fishery_error","ln_extra_survey_cv"))
est_obj <- MakeADFun(TMB_data, init_vals, DLL="ASMycs", random = c("ln_ycs_est"), map = map_fix_pars, silent = T)
re_estimate_pars = est_obj$par
pre_rep = est_obj$report()
#est_obj$env$tracepar = TRUE # track params
#options(warn = 2)
#options()$warn
opt_ = nlminb(est_obj$par, est_obj$fn, est_obj$gr, control = list(iter.max = 10000, eval.max = 10000, abs.tol = 1e-20, rel.tol = 1e-10))
## check convergence and data-weighting
if(DATA_WEIGHTING) {
  Rep_est = est_obj$report(est_obj$env$last.par.best)
  ## if multinomial then do Francis re-weighting
  for(k in 1:3) {
    Re_estiamte = F
    if(TMB_data$survey_comp_likelihood == 1) {
      wj = Method.TA1.8(bin_lab = survey_ages, observed = TMB_data$survey_trans_comp, expected = Rep_est$survey_comp_fitted, error_value = colSums(TMB_data$survey_trans_comp))
      if(abs(wj - 1.0) > 0.1) {
        TMB_data$survey_trans_comp = matrix(as.integer(TMB_data$survey_trans_comp * wj + 0.5), nrow = nrow(TMB_data$survey_trans_comp), ncol = ncol(TMB_data$survey_trans_comp))## needs to be an integer for internal TMB simualtion code to work
        Re_estiamte = TRUE
      }
    }
    if(TMB_data$fishery_comp_likelihood == 1) {
      wj = Method.TA1.8(bin_lab = fishery_ages, observed = TMB_data$fishery_trans_comp, expected = Rep_est$fishery_comp_fitted, error_value = colSums(TMB_data$fishery_trans_comp))
      if(abs(wj - 1.0) > 0.1) {
        TMB_data$fishery_trans_comp = matrix(as.integer(TMB_data$fishery_trans_comp * wj + 0.5), nrow = nrow(TMB_data$fishery_trans_comp), ncol = ncol(TMB_data$fishery_trans_comp)) ## needs to be an integer for internal TMB simualtion code to work
        Re_estiamte = TRUE
      }
    }  
    if(Re_estiamte) {
      est_obj <- MakeADFun(TMB_data, start_vals, DLL= "ASMycs", random = c("ln_ycs_est"), map = map_fix_pars, silent = T)
      opt_ = nlminb(start = est_obj$par, objective = est_obj$fn, gradient = est_obj$gr, control = list(eval.max = 10000, iter.max = 10000))#, x.tol = 1e-15, rel.tol = 1e-12))
      ## check convergence
      ##conv_msg = check_tmb_convergence(obj = est_obj, delta = grad_tol)
      if(opt_$convergence != 0) {
        next;
      }
      Rep_est = est_obj$report(est_obj$env$last.par.best)
    } else {
      break; ## I don't need
    }
  }
}

#check_converg = check_tmb_convergence(est_obj, converge_tol);
#if (check_converg != "No evidence of non convergence") {
#  opt_$convergence = 1
#}
#if (opt_$convergence != 0) {
#  loopnum = 4
#  # Re-run to further decrease final gradient
#  for(k in seq(2,loopnum,length=max(0,loopnum-1)) ){
#    Temp = opt_[c('iterations','evaluations')]
#    these_start_val = start_vals()
#    for(n_p in 1:length(unique(names(re_estimate_pars)))) {
#      this_par = get(x = unique(names(re_estimate_pars))[n_p], these_start_val)
#      re_estimate_pars[names(re_estimate_pars) %in% unique(names(re_estimate_pars))[n_p]] = this_par
#    }
#    opt_ = nlminb(re_estimate_pars, objective = est_obj$fn, gradient = est_obj$gr, control = list(eval.max = 10000, iter.max = 10000))
#    opt_[['iterations']] = opt_[['iterations']] + Temp[['iterations']]
#    opt_[['evaluations']] = opt_[['evaluations']] + Temp[['evaluations']]
#    check_converg = check_tmb_convergence(est_obj, converge_tol)
#    if (check_converg == "No evidence of non convergence") {
#      opt_$convergence = 0
#      break;
#    }
#  }
#}
#cat("opt - ", opt_$convergence, " - convg ", check_converg, "\n")
check_converg = check_tmb_convergence(est_obj, converge_tol);
if(debug)
  cat("convergence - opt - ", opt_$convergence, " - convg ", check_converg, "\n")

MLE_est = est_obj$report(est_obj$env$last.par.best)
## save MLE estimates and data
save(MLE_est, TMB_data, file = file.path(est_out_dir, paste0("mse_year_",current_year, extension, ".RData")))
#####################
## HCR
#####################
## find reference B0 and stock status
if(is_final_year == 0) {
  if(debug)
    cat("Calculate projections\n")
  proj_pars = NULL
  mle_covar = sdreport(est_obj, getJointPrecision = T)
  if(length(est_obj$env$random) > 0) {
    proj_pars = t(rmvnorm_prec(est_obj$env$last.par.best, mle_covar$jointPrecision, n.sims = n_proj_sims, random_seed = trunc(runif(1,1,1e5))))
  } else {
    proj_pars = mvrnorm(n = n_proj_sims, est_obj$env$last.par.best, Sigma = mle_covar$cov.fixed)
  }
  proj_data = TMB_data
  proj_data$n_proj_years = reference_time
  proj_data$proj_method = 1 # Catch not F
  proj_data$future_vals = 0.0 # zero catch
  #proj_data$future_F = 0.0
  ## Build the projection model
  proj_fun = MakeADFun(proj_data, init_vals, DLL= "ASMycs_proj", random = c("ln_ycs_est"), map = map_fix_pars, silent = T, type = "Fun")
  ref_bio = matrix(NA,nrow = n_proj_sims, ncol = length(50:101) ) ## use the last projected 50 years with F = 0 as the reference Biomass
  rel_status = vector();
  for(i in 1:n_proj_sims) {
    set.seed(trunc(runif(1,1, 1e5)))
    proj_rep = proj_fun$report(proj_pars[i, ])
    ref_bio[i, ] = proj_rep$ssb_proj[50:101]
    rel_status[i] = MLE_est$ssb[length(MLE_est$ssb)] / mean(ref_bio[i, ]) * 100
  }
  ## now check if we need to rebuild or keep the status quo?
  fishing_ls = list()
  current_stock_status = mean(rel_status)
  reference_bio = mean(ref_bio)
  future_catch = NULL;
  if(debug)
    cat("stock status = ", current_stock_status, "\n")
  if (management_state == "rebuild") {
    if(debug)
      cat("We are in a rebuild phase, check on track\n")
    ## check we are on track to successful rebuild
    years_left = rebuild_time - (current_year - start_rebuild_year)
    proj_data$n_proj_years = years_left;
    ## Keep catch same as avg last 5 years
    target_reference = 0.3
    proj_fun$env$data$future_vals = TMB_data$catches[TMB_data$nyears]
    ssb_above_target = 0;
    for(i in 1:n_proj_sims) {
      set.seed(trunc(runif(1,1, 1e5)))
      proj_rep = proj_fun$report(proj_pars[i,])
      if (proj_rep$ssb_proj[rebuild_time + 1] > mean(ref_bio[i, ]) * target_reference)
        ssb_above_target = ssb_above_target + 1
    }
    if (ssb_above_target < 0.5) {
      if(debug)
        cat("We are not on track update rebuild plan\n")
      
      ## adapt catch
      catches_to_search = seq(from = max(signif(trunc(min(TMB_data$catches)),1),catch_tonnes_to_search), to = signif(trunc(max(TMB_data$catches)),1), by = 100)
      prob_above_status = vector()
      for(f_ndx in 1:length(catches_to_search)) {
        proj_fun$env$data$future_vals = catches_to_search[f_ndx]
        ssb_above_target = 0;
        for(i in 1:n_proj_sims) {
          set.seed(trunc(runif(1,1, 1e5)))
          proj_rep = proj_fun$report(proj_pars[i,])
          if (proj_rep$ssb_proj[length(proj_rep$ssb_proj)] > mean(ref_bio[i, ]) * target_reference)
            ssb_above_target = ssb_above_target + 1
        }
        prob_above_status[f_ndx] = ssb_above_target / n_proj_sims
      }
      future_catch = catches_to_search[tail(which(prob_above_status > 0.5),1)]
    } else {
      future_catch = TMB_data$catches[TMB_data$nyears]
      if(debug)
        cat("Rebuild plan on track keep going\n")
      
    }
  } else {
    if(current_stock_status > 50) {
      if(debug)
        cat("Stock in good shape above 50% Bref\n")
      ## Change Catch most likely to be higher than previous catch
      target_reference = 0.4
      ## increse Catch 
      catches_to_search = seq(from = max(signif(trunc(min(TMB_data$catches)),1),catch_tonnes_to_search), to = signif(trunc(max(TMB_data$catches)),1), by = 100)
      prob_above_status = vector()
      for(f_ndx in 1:length(catches_to_search)) {
        proj_fun$env$data$future_vals = catches_to_search[f_ndx]
        ssb_above_target = 0;
        for(i in 1:n_proj_sims) {
          set.seed(trunc(runif(1,1, 1e5)))
          proj_rep = proj_fun$report(proj_pars[i,])
          if (proj_rep$ssb_proj[length(proj_rep$ssb_proj)] > mean(ref_bio[i, ]) * target_reference)
            ssb_above_target = ssb_above_target + 1
        }
        prob_above_status[f_ndx] = ssb_above_target / n_proj_sims
      }
      future_catch = catches_to_search[tail(which(prob_above_status > 0.5),1)]
      if(debug)
        cat("Updated catch to ", future_catch, "\n")
      
    } else if(current_stock_status > 30) {  
      if(debug)
        cat("above 30% check we are on track with recent catch\n")
      
      ## Keep catch same as avg last 5 years
      target_reference = 0.3
      proj_fun$env$data$future_vals = TMB_data$catches[TMB_data$nyears]
      ssb_above_target = 0;
      for(i in 1:n_proj_sims) {
        set.seed(trunc(runif(1,1, 1e5)))
        proj_rep = proj_fun$report(proj_pars[i,])
        if (proj_rep$ssb_proj[length(proj_rep$ssb_proj)] > mean(ref_bio[i, ]) * target_reference)
          ssb_above_target = ssb_above_target + 1
      }
      if (ssb_above_target < 0.5) {
        if(debug)
          cat("Not going to make target so update\n")
        ## increse Catch 
        catches_to_search = seq(from = max(signif(trunc(min(TMB_data$catches)),1),catch_tonnes_to_search), to = signif(trunc(max(TMB_data$catches)),1), by = 100)
        prob_above_status = vector()
        for(f_ndx in 1:length(catches_to_search)) {
          proj_fun$env$data$future_vals = catches_to_search[f_ndx]
          ssb_above_target = 0;
          for(i in 1:n_proj_sims) {
            set.seed(trunc(runif(1,1, 1e5)))
            proj_rep = proj_fun$report(proj_pars[i,])
            if (proj_rep$ssb_proj[length(proj_rep$ssb_proj)] > mean(ref_bio[i, ]) * target_reference)
              ssb_above_target = ssb_above_target + 1
          }
          prob_above_status[f_ndx] = ssb_above_target / n_proj_sims
        }
        future_catch = catches_to_search[tail(which(prob_above_status > 0.5),1)]
      } else {
        if(debug)
          cat("We are on track keep going\n")
        future_catch = TMB_data$catches[TMB_data$nyears]
      }
    } else if(current_stock_status > 15 & current_stock_status <= 20) {
      if(debug)
        cat("Entering rebuild phase\n")
      ## Rebuild plan
      management_state = "rebuild"
      start_rebuild_year = current_year;
      target_reference = 0.4
      ## increse Catch 
      catches_to_search = seq(from = max(signif(trunc(min(TMB_data$catches)),1),catch_tonnes_to_search), to = signif(trunc(max(TMB_data$catches)),1), by = 100)
      prob_above_status = vector()
      for(f_ndx in 1:length(catches_to_search)) {
        proj_fun$env$data$future_vals = catches_to_search[f_ndx]
        ssb_above_target = 0;
        for(i in 1:n_proj_sims) {
          set.seed(trunc(runif(1,1, 1e5)))
          proj_rep = proj_fun$report(proj_pars[i,])
          if (proj_rep$ssb_proj[length(proj_rep$ssb_proj)] > mean(ref_bio[i, ]) * target_reference)
            ssb_above_target = ssb_above_target + 1
        }
        prob_above_status[f_ndx] = ssb_above_target / n_proj_sims
      }
      future_catch = catches_to_search[tail(which(prob_above_status > 0.5),1)]
    } else {
      if(debug)
        cat("Closing the fishery < 15% reference stock status\n")
      future_catch = 0.01 ## can't be zero as we have a lognormal dist in estimator
    }
  }
  ## save projections with this catch, will be interesting to compare with actual trajectory
  proj_fun$env$data$future_vals = future_catch
  proj_ssbs = matrix(NA, nrow = n_proj_sims, ncol = length(proj_rep$ssb_proj))
  for(i in 1:n_proj_sims) {
    set.seed(trunc(runif(1,1, 1e5)))
    proj_rep = proj_fun$report(proj_pars[i,])
    proj_ssbs[i, ] = proj_rep$ssb_proj
  }
  save(proj_ssbs, future_catch, file = file.path(est_out_dir, paste0("proj_year_",current_year, extension, ".RData")))

  ## create object that will be returned to the ABM C++ interface
  Fishing_ls = list();
  fishing_vector = future_catch
  names(fishing_vector) = "Fishing"
  for(year_ndx in (current_year + 1):next_ass_year)
    Fishing_ls[[as.character(year_ndx)]] = fishing_vector
  Fishing_ls
} else {
  if(debug)
    cat("in final year of run no projections needed\n")
  ## final year just want to save estimated state not care about projections.
  Fishing_ls = list();
  fishing_vector = TMB_data$catches[TMB_data$nyears]
  names(fishing_vector) = "Fishing"
  Fishing_ls[[as.character(current_year)]] = fishing_vector
  Fishing_ls
}