#cat("correctly read me =)\n") # this will be printed into the > outtxt file, which will break the R-library read functions
library(TMB)
library(ibm)
library(mvtnorm)
## if you change the directory i.e. setwd()
## this will chnage the directory in the C++ and will effect the reporting.
debug = TRUE; ## if true it will print to the abm report. And ruin the R extract code, but useful when debugging MSE initally
model_label = "M1" # used for saving objects
## everything will be relative to this
abm_dir = getwd() ## assuming you are running the ibm model in the ibm directory 
abm_sim_dir = file.path(abm_dir, "sim")## assuming you are running the ibm model in the ibm directory 
r_dir = file.path(abm_dir, "..","R")
#setwd(r_dir)
tmb_dir = file.path(abm_dir, "..","TMB")
abm_sim_dir = file.path(abm_dir, "sim")
savedobjs_dir = file.path(r_dir, "savedOBJs")
est_out_dir = file.path(r_dir, "est_outputs", model_label)
if(debug)
  cat("tmb_dir ", tmb_dir, "\n")
if(debug)
  cat("abm_dir ", abm_dir, "\n")
if(debug)
  cat("r_dir ", r_dir, "\n")

## source auxillary functions
ls_files = list.files(file.path(r_dir, "Funs"), pattern = "\\.R$")
for(i in 1:length(ls_files))
  source(file.path(r_dir,"Funs", ls_files[i]))
ls_files = list.files(file.path(r_dir, "Funs", "Fisheries"), pattern = "\\.R$")
for(i in 1:length(ls_files))
  source(file.path(r_dir, "Funs","Fisheries", ls_files[i]))
## MLE settings
init_age_comp_sample_size = 500
DATA_WEIGHTING = F
converge_tol = 0.01
## Assumes you have already Compiled TMB model
#file.create("sink_info.txt",showWarnings  = F)
#compile(file = file.path(tmb_dir, "ASMycs.cpp"),trace = F, tracesweep = F)
dyn.load(dynlib(file.path(tmb_dir,"ASMycs"))) ## load tmb
#dyn.unload(dynlib(file.path(tmb_dir, "ASMycs_proj")))
#compile(file.path(tmb_dir,"ASMycs_proj.cpp"), flags = "",DLLFLAGS="");
dyn.load(dynlib(file.path(tmb_dir, "ASMycs_proj")))
## Build TMB-data and save, which will have be reloaded in each MSE step 
## i.e. every RunHCR.R call
#abm = extract.run(file = "run.out", path = abm_dir)
#save(abm, file = file.path(savedobjs_dir, "abm_output.RData"), version = 2)
load(file = file.path(savedobjs_dir, "abm_output.RData"))
abm_pop = extract.ibm.file(file = "Population.ibm", path = abm_dir, quiet = T)
abm_obs = extract.ibm.file(file = "Observation.ibm", path = abm_dir, quiet = T)
years = abm$model_attributes$model_years
## change years so it starts at the first assesment year
years = min(years):min(as.numeric(abm_pop$model$assessment_years$value))
ages = abm$model_attributes$ages
years_with_init = c(min(years) - 1, years)
mean_length_at_age = schnute(as.numeric(abm_pop$`process[Growth]`$y1$value), as.numeric(abm_pop$`process[Growth]`$y2$value), as.numeric(abm_pop$`process[Growth]`$tau1$value), as.numeric(abm_pop$`process[Growth]`$tau2$value), as.numeric(abm_pop$`process[Growth]`$alpha$value), as.numeric(abm_pop$`process[Growth]`$beta$value), ages)
abm_F = abm$Fishing$F_by_fishery_year["Fishing",]
abm_F = abm_F[names(abm_F) %in% years]
# - Actual Catches
abm_catch_ndx = grepl(names(abm$Fishing), pattern = "actual_catches")
abm_catch = Reduce(c, lapply(abm$Fishing[abm_catch_ndx], FUN = function(x){as.numeric(x[1])}))
abm_catch = abm_catch[abm$model_attributes$model_years %in% years]
# - Recruitment Variability
recruit_cv = 0.6
Fish_exp = reformat.compositional.data(model = abm, report_label = "FishingAge", quiet = T)$`r1-c1`$fit
survey_exp = ibm::reformat.compositional.data(model = abm, report_label = "chatTANage", quiet = T)$`r1-c1`$fit
survey_age_ndx = unique(abm$chatTANage$Values$year) %in% years
fishery_ages_ndx = unique(abm$FishingAge$Values$year) %in% years
Fish_exp = Fish_exp[fishery_ages_ndx, ]
survey_exp = survey_exp[survey_age_ndx, ]

####################
## Build Biology
####################
TMB_data = list()
TMB_data$years = years
TMB_data$ages = ages
TMB_data$nyears = length(TMB_data$years)
TMB_data$nages = length(TMB_data$ages)
TMB_data$maxAgePlusGroup = 1
TMB_data$sel_a50_bounds = c(0.1,50)
TMB_data$sel_ato95_bounds = c(0.1,50)
TMB_data$mean_weight_a = as.numeric(abm_pop$`process[Growth]`$a$value)
TMB_data$mean_weight_b = as.numeric(abm_pop$`process[Growth]`$b$value)
TMB_data$natMor = as.numeric(abm_pop$`process[Natural_Mortality]`$m$value)
TMB_data$sigma_r =  sqrt(log(recruit_cv^2 + 1))
TMB_data$steepness = as.numeric(abm_pop$`process[Recruitment]`$steepness$value)
TMB_data$ycs_estimated = rep(1,TMB_data$nyears)
TMB_data$ycs_estimated[(length(TMB_data$ycs_estimated) - 1):length(TMB_data$ycs_estimated)] = 0
TMB_data$stockRecruitmentModelCode = 2 ## BH
TMB_data$standardise_ycs = 1
prop_mat = rep(0,TMB_data$nages)
prop_mat[as.numeric(abm_pop$`selectivity[MaturationSel]`$l$value):as.numeric(abm_pop$`selectivity[MaturationSel]`$h$value)] = as.numeric(abm_pop$`selectivity[MaturationSel]`$v$value) #hak_pop$maturity_props$all[4:length(hak_pop$maturity_props$all)]
prop_mat[as.numeric(abm_pop$`selectivity[MaturationSel]`$h$value):length(prop_mat)] = 1 
class(prop_mat) = "numeric"
TMB_data$propMat = matrix(prop_mat, ncol = TMB_data$nyears, nrow = TMB_data$nages, byrow = F)
TMB_data$catchMeanLength = TMB_data$stockMeanLength = matrix(mean_length_at_age, ncol = TMB_data$nyears, nrow = TMB_data$nages, byrow = F)
TMB_data$propZ_survey = rep(0.5, TMB_data$nyears)
TMB_data$propZ_ssb = rep(0.5, TMB_data$nyears)
TMB_data$catches = as.numeric(abm_catch)
#####################
## observation stuff
#####################
TMB_data$ageing_error_matrix = matrix(0, nrow = TMB_data$nages, ncol = TMB_data$nages)
diag(TMB_data$ageing_error_matrix) = 1;

TMB_data$survey_year_indicator = as.integer(TMB_data$years %in% abm$chatTANbiomass$Values$year)
TMB_data$survey_obs = abm$chatTANbiomass$Values$expected
TMB_data$survey_cv = rep(0.1,sum(TMB_data$survey_year_indicator))
TMB_data$survey_min_age = 3
TMB_data$survey_max_age = 25
TMB_data$survey_age_plus_group = 1
survey_ages = TMB_data$survey_min_age:TMB_data$survey_max_age
TMB_data$survey_comp_type = 0 ## 0 = prop at age, 1 = numbers
TMB_data$survey_covar_structure = 1 # iid
TMB_data$survey_comp_likelihood = 1
n_surv_ages = length(TMB_data$survey_min_age:TMB_data$survey_max_age)
if(TMB_data$survey_comp_type == 0 & TMB_data$survey_comp_likelihood == 0)
  n_surv_ages = n_surv_ages - 1
TMB_data$survey_trans_comp = matrix(1000 * t(survey_exp), nrow = n_surv_ages, ncol = sum(TMB_data$survey_year_indicator))

TMB_data$fishery_year_indicator = as.integer(TMB_data$years %in% rownames(Fish_exp));
TMB_data$fishery_min_age = 1
TMB_data$fishery_max_age = 30
fishery_ages = TMB_data$fishery_min_age:TMB_data$fishery_max_age
TMB_data$fishery_comp_type = 0 ## 0 = prop at age, 1 = numbers
TMB_data$fishery_age_plus_group = 1
TMB_data$fishery_covar_structure = 1 # iid
TMB_data$fishery_comp_likelihood = 1 ## 0 = MVN, 1 = multinomial, 2 = dirichlet-multinomial
n_fishery_ages = length(TMB_data$fishery_min_age:TMB_data$fishery_max_age)
if(TMB_data$fishery_comp_type == 0 & TMB_data$fishery_comp_likelihood == 0)
  n_fishery_ages = n_fishery_ages - 1
## Truncate the fishery abm doesn't truncate for some reason
Fish_exp_trunc = Fish_exp[,fishery_ages]
Fish_exp_trunc[,ncol(Fish_exp_trunc)] = rowSums(Fish_exp[,c((TMB_data$fishery_max_age + 1):ncol(Fish_exp))])
if(n_fishery_ages != nrow(t(Fish_exp_trunc)))
  print("found incompatible issue")
TMB_data$fishery_trans_comp = matrix(1000 * t(Fish_exp_trunc), nrow = n_fishery_ages, ncol = sum(TMB_data$fishery_year_indicator))
F_a50 =  abm$FSel$a50
F_ato95 =  abm$FSel$ato95

# parameters
params = list(
  logit_f_a50 = logit_general(F_a50, TMB_data$sel_a50_bounds[1], TMB_data$sel_a50_bounds[2]),
  logit_f_ato95 = logit_general(F_ato95, TMB_data$sel_ato95_bounds[1], TMB_data$sel_ato95_bounds[2]),
  logit_survey_a50 = logit_general(as.numeric(abm_pop$`selectivity[chatTANsel]`$a50$value), TMB_data$sel_a50_bounds[1], TMB_data$sel_a50_bounds[2]),
  logit_survey_ato95 = logit_general(as.numeric(abm_pop$`selectivity[chatTANsel]`$ato95$value), TMB_data$sel_a50_bounds[1], TMB_data$sel_a50_bounds[2]),
  ln_extra_survey_cv = log(0.001),
  trans_survey_error = log(0.1),
  trans_fishery_error = log(0.1),
  logit_surveyQ = qlogis(as.numeric(abm_obs$`observation[chatTANbiomass]`$catchability$value)),
  ln_F = log(as.numeric(abm_F)),
  ln_R0 = log(abm$model_attributes$Recruitment * abm$Rec$initial_recruits),
  ln_catch_sd = log(0.02),
  ln_sigma_r = log(0.7),
  #recruit_devs = log(abm$Rec$ycs_values)
  ln_ycs_est = log(abm$Rec$ycs_values[TMB_data$ycs_estimated == 1])
)

#map_fix_pars = fix_pars(par_list = params, pars_to_exclude = c("ln_catch_sd", "trans_survey_error","trans_fishery_error","ln_extra_survey_cv", "ln_sigma_r"))
#test_obj <- MakeADFun(TMB_data, params, DLL="ASMycs", map = map_fix_pars, silent = T)
#test_obj$fn()
save(TMB_data, params, file = file.path(savedobjs_dir, "First_tmb_dir.RData"))


