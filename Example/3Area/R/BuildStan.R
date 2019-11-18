#
detach("package:ibm", unload=TRUE)

# Set up Stan
library(rstan)
library(ggplot2)
library(reshape2)
library(ibm)
library(abind)
## change some settings
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores() - 2)
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')

stan_dir = file.path(getwd(),"../Stan")
list.files(stan_dir)

## Single stock Recruit function
stan_sim_model = stan_model(file.path(stan_dir, "SpatialModelSim.stan"))
## IBM
ibm = extract.run(file = "../base_ibm/single.out")

years = 1990:2013
stan_data = list()
stan_data$Y = length(years)
stan_data$A = 30
stan_data$R = 3
stan_data$T = 2*stan_data$R

stan_data$proportion_mortality_spawning = 0.0
stan_data$proportion_mortality_biomass = 1.0
stan_data$u_max = 0.8
stan_data$catch_penalty = 1
stan_data$penalty_in_log_space = 1

## observations
# Fishery age comp
# CPUE biomass
# no survey age data
stan_data$fishery_obs_indicator = matrix(1,nrow = stan_data$R, ncol = stan_data$Y)
stan_data$survey_age_indicator = matrix(0,nrow = stan_data$R, ncol = stan_data$Y)
stan_data$biomass_indicator = matrix(1,nrow = stan_data$R, ncol = stan_data$Y)

stan_data$M = 0.2
stan_data$a = 4.467e-08
stan_data$b = 2.793
stan_data$L_inf = 61.78748773
stan_data$k = 0.134
stan_data$t0 = -1.2
stan_data$h = 0.75
stan_data$mat_a50 = 3.9
stan_data$mat_ato95 = 4.2
stan_data$apply_prior = 0 
stan_data$ycs_prior_applies_to_standardised = 1;

#stan_data$proportion_recruitment_by_region = c(0.3,0.3,0.4)
stan_data$proportion_recruitment_by_region = c(0.33333,0.33333,0.33333)

stan_data$prob_move = matrix(0, nrow = 3, ncol = 3, byrow = T)
stan_data$prob_move[1,] = as.numeric(ibm$Jump_One$`year-area-1990_1-1`$destination_values / sum(ibm$Jump_One$`year-area-1990_1-1`$destination_values))
stan_data$prob_move[2,] = as.numeric(ibm$Jump_One$`year-area-1990_2-1`$destination_values / sum(ibm$Jump_One$`year-area-1990_2-1`$destination_values))
stan_data$prob_move[3,] = as.numeric(ibm$Jump_One$`year-area-1990_3-1`$destination_values / sum(ibm$Jump_One$`year-area-1990_3-1`$destination_values))

stan_data$debug = 0

set.seed(15)
x = rlnorm(n = stan_data$Y, log(1), 0.6)
round(x / mean(x),2)

stan_data$fishery_at_age_samples = matrix(rbinom(3*24, size = 500, prob = 0.7), nrow = 3, ncol = 24)

## parameters 
init_pars = list()
init_pars$p = 0.5
set.seed(15)
x = rlnorm(n = stan_data$Y, log(1), 0.6)
round(x / mean(x),2)
init_pars$unity_YCS = rep(1, stan_data$Y) / stan_data$Y
init_pars$ln_R0 = 16
init_pars$f_a50 = 2
init_pars$f_ato95 = 3
init_pars$Q = 0.2832
init_pars$overdispersion_tag_ll = 0.4
init_pars$p = 0.5
init_pars$rho = 0.5
init_pars$sigma = 1.3
init_pars$adj_biomass_cv = 0.03
## use observations from IBM quicker
ibm$Jump_One$`year-area-1990_1-1`$initial_numbers_in_cell
sum(ibm$Jump_One$`year-area-1990_1-1`$destination_values)

ibm_ssb = ibm::plot.derived_quantities(ibm, "derived_quants", plot.it = F)
ibm_ssb[,"SSB"]
# R0
ibm$model_attributes$Recruitment * ibm$Rec$initial_recruits
stan_sim$par$R0

years = unique(ibm$EN_age_sample$Values$year)
ages = 1:30
stan_data$Y_f = length(as.numeric(rownames(ibm$fishing$catches)))
stan_data$fish_age_years = as.numeric(rownames(ibm$fishing$catches))
fish_obs_A1 = matrix(ibm$EN_age_sample$Values$expected[ibm$EN_age_sample$Values$age != 0], ncol = length(years), nrow = length(ages))
fish_obs_A2 = matrix(ibm$HG_age_sample$Values$expected[ibm$HG_age_sample$Values$age != 0], ncol = length(years), nrow = length(ages))
fish_obs_A3 = matrix(ibm$BP_age_sample$Values$expected[ibm$BP_age_sample$Values$age != 0], ncol = length(years), nrow = length(ages))

stan_data$fishery_at_age_obs = round(abind(fish_obs_A1,fish_obs_A2,fish_obs_A3, rev.along = 3))
## numbers at age
#stan_data$fishery_at_age_obs = round(stan_data$fishery_at_age_obs * 250)
#stan_data$survey_at_age_obs[1,,]
#fish_obs_A1 
stan_data$biomass_obs = rbind(ibm$CPUE_EN$Values$expected, ibm$CPUE_HG$Values$expected, ibm$CPUE_BP$Values$expected)
stan_data$biomass_error = matrix(rnorm(stan_data$Y * stan_data$R, mean = 0.13, 0.05), nrow = 3, ncol = stan_data$Y)
stan_data$biomass_error[stan_data$biomass_error<=0.05] = 0.1

process_error_sqrd = 0.04^2
obs_cv = sqrt(stan_data$biomass_error^2 - process_error_sqrd)

R0 = 40535890 #(ibm$Rec$initial_recruits) * ibm$model_attributes$Recruitment
catches = NULL;
for(i in 1:length(years)) {
  this_catch = eval(expr = parse(text = paste0("ibm$fishing$`actual_catches-", years[i],"`")))
  catches = cbind(catches, as.numeric(this_catch[,1]))
}

stan_data$catches = catches
stan_data$years = years
stan_data$tag_years = c(1995,1995,1995,2005,2005,2005)


stan_data$unity_YCS = ibm$Rec$ycs_values / sum(ibm$Rec$ycs_values)
stan_data$ln_R0 = log(R0)
stan_data$f_a50 = 2.4
stan_data$f_ato95 = 7
stan_data$Q = 0.3832

stan_data$fishery_obs_indicator = stan_data$biomass_indicator
stan_data$tag_release_by_age = array(0,dim= c(stan_data$R,stan_data$A, stan_data$T))
stan_data$tag_release_by_age[1,,1] = as.matrix(ibm$Tagging$`tag_release_by_age-1995`)[1,]
stan_data$tag_release_by_age[2,,2] = as.matrix(ibm$Tagging$`tag_release_by_age-1995`)[2,]
stan_data$tag_release_by_age[3,,3] = as.matrix(ibm$Tagging$`tag_release_by_age-1995`)[3,]
stan_data$tag_release_by_age[1,,4] = as.matrix(ibm$Tagging$`tag_release_by_age-2005`)[1,]
stan_data$tag_release_by_age[2,,5] = as.matrix(ibm$Tagging$`tag_release_by_age-2005`)[2,]
stan_data$tag_release_by_age[3,,6] = as.matrix(ibm$Tagging$`tag_release_by_age-2005`)[3,]
# check total number of tags
sum(stan_data$tag_release_by_age)
sum(ibm$Tagging$`tag_release_by_age-1995`) + sum(ibm$Tagging$`tag_release_by_age-2005`)

stan_data$tag_recapture_indicator = array(0,dim= c(stan_data$R, stan_data$T,stan_data$Y))


release_years = which(years %in% stan_data$tag_years)
## years after that we want to record recaptures 2 years post release
recapture_years = c(release_years, release_years + 1, release_years + 2, release_years + 3)
stan_data$tag_recapture_indicator[,1,sort(recapture_years)[1:4]] = 1
stan_data$tag_recapture_indicator[,2,sort(recapture_years)[1:4]] = 1
stan_data$tag_recapture_indicator[,3,sort(recapture_years)[1:4]] = 1
stan_data$tag_recapture_indicator[,4,sort(recapture_years)[5:8]] = 1
stan_data$tag_recapture_indicator[,5,sort(recapture_years)[5:8]] = 1
stan_data$tag_recapture_indicator[,6,sort(recapture_years)[5:8]] = 1


ibm$fishing$`tag_recapture_info-1995-1-1`$scanned_fish * ibm$fishing$`tag_recapture_info-1995-1-1`$prob_tagged_fish
dim(ibm$fishing$`tag_recapture_info-1995-1-1`$values)
ibm$fishing$`tag_recapture_info-1995-1-1`$agents_sampled * ibm$fishing$`tag_recapture_info-1995-1-1`$prob_tagged_fish


years[recapture_years -1]
recapture_obs = array(0.0,dim = c(stan_data$T, stan_data$R, stan_data$Y))
recapture_obs[1,1,6] = ncol(ibm$fishing$`tag_recapture_info-1995-1-1`$values)
recapture_obs[1,2,6] = ncol(ibm$fishing$`tag_recapture_info-1995-2-1`$values)
recapture_obs[1,3,6] = ncol(ibm$fishing$`tag_recapture_info-1995-3-1`$values)

recapture_obs[1,1,7] = ncol(ibm$fishing$`tag_recapture_info-1996-1-1`$values)
recapture_obs[1,2,7] = ncol(ibm$fishing$`tag_recapture_info-1996-2-1`$values)
recapture_obs[1,3,7] = ncol(ibm$fishing$`tag_recapture_info-1996-3-1`$values)

recapture_obs[1,1,8] = ncol(ibm$fishing$`tag_recapture_info-1997-1-1`$values)
recapture_obs[1,2,8] = ncol(ibm$fishing$`tag_recapture_info-1997-2-1`$values)
recapture_obs[1,3,8] = ncol(ibm$fishing$`tag_recapture_info-1997-3-1`$values)

recapture_obs[1,1,9] = ncol(ibm$fishing$`tag_recapture_info-1998-1-1`$values)
recapture_obs[1,2,9] = ncol(ibm$fishing$`tag_recapture_info-1998-2-1`$values)
recapture_obs[1,3,9] = ncol(ibm$fishing$`tag_recapture_info-1998-3-1`$values)

recapture_obs[2,1,16] = ncol(ibm$fishing$`tag_recapture_info-2005-1-1`$values)
recapture_obs[2,2,16] = ncol(ibm$fishing$`tag_recapture_info-2005-2-1`$values)
recapture_obs[2,3,16] = ncol(ibm$fishing$`tag_recapture_info-2005-3-1`$values)

recapture_obs[2,1,17] = ncol(ibm$fishing$`tag_recapture_info-2006-1-1`$values)
recapture_obs[2,2,17] = ncol(ibm$fishing$`tag_recapture_info-2006-2-1`$values)
recapture_obs[2,3,17] = ncol(ibm$fishing$`tag_recapture_info-2006-3-1`$values)

recapture_obs[2,1,18] = ncol(ibm$fishing$`tag_recapture_info-2007-1-1`$values)
recapture_obs[2,2,18] = ncol(ibm$fishing$`tag_recapture_info-2007-2-1`$values)
recapture_obs[2,3,18] = ncol(ibm$fishing$`tag_recapture_info-2007-3-1`$values)

recapture_obs[2,1,19] = ncol(ibm$fishing$`tag_recapture_info-2008-1-1`$values)
recapture_obs[2,2,19] = ncol(ibm$fishing$`tag_recapture_info-2008-2-1`$values)
recapture_obs[2,3,19] = ncol(ibm$fishing$`tag_recapture_info-2008-3-1`$values)
stan_data$tag_recapture_obs = recapture_obs
# stan_data$tag_recapture_obs[2,,] = recapture_obs[2,3,19]
# recapture_obs[2,,] = recapture_obs[2,,] / sum(stan_data$tag_release_by_age[,,2])

stan_data$rho = 0.5
stan_data$sigma = 0.4
stan_data$N_sims = 100
stan_data$tag_neg_bin_overdispersion = 17.5

stan_data$tag_recapture_obs = array(100, dim = c(stan_data$R, stan_data$T, stan_data$Y))

stan_sim = optimizing(stan_sim_model,data=stan_data,init=init_pars,
                      iter=1,algorithm="LBFGS",
                      verbose=F,as_vector = FALSE)

par(mfrow = c(1,1))
plot(rownames(ibm_ssb),ibm_ssb[,"SSB"], type = "l", lwd = 2, lty = 2, ylim = c(0,300000))
lines(rownames(ibm_ssb),stan_sim$par$SSB, col = "red", lwd = 2)
legend('bottomleft', col = c("black","red"), legend = c("ABM","Stan"), lty = c(2,1), lwd = 2)
box()
stan_sim$par$B0
ibm$Rec$b0

stan_sim$par$standardised_ycs
ibm$Rec$ycs_values
# read in simulated data
sim_obs_location= "../base_ibm/sim_obs"
ibm_sim = extract.run("../base_ibm/simulated_obs.log")
ibm_n_runs = length(ibm_sim$derived_quants)
for(i in 1:ibm_n_runs)
  lines(rownames(ibm_ssb),ibm_sim$derived_quants[[i]]$SSB$values, col = adjustcolor(col = "black", alpha.f = 0.3), lwd = 2)


stan_data$N_sims = 1
stan_en_tag_95 = stan_hg_tag_95 = stan_bp_tag_95 = array(NA, dim = c(3, ibm_n_runs, 4))
stan_en_tag_00 = stan_hg_tag_00 = stan_bp_tag_00 = array(NA, dim = c(3, ibm_n_runs, 4))

for(i in 1:ibm_n_runs) {
  catches = NULL;
  for(k in 1:length(years)) {
    this_catch = eval(expr = parse(text = paste0("ibm_sim$fishing[[i]]$`actual_catches-", years[k],"`")))
    catches = cbind(catches, as.numeric(this_catch[,1]))
  }
  print(sum(catches))
  stan_data$catches = catches
  stan_data$tag_release_by_age = array(0,dim= c(stan_data$R,stan_data$A, stan_data$T))
  stan_data$tag_release_by_age[1,,1] = as.matrix(ibm_sim$Tagging[[i]]$`tag_release_by_age-1995`)[1,]
  stan_data$tag_release_by_age[2,,2] = as.matrix(ibm_sim$Tagging[[i]]$`tag_release_by_age-1995`)[2,]
  stan_data$tag_release_by_age[3,,3] = as.matrix(ibm_sim$Tagging[[i]]$`tag_release_by_age-1995`)[3,]
  stan_data$tag_release_by_age[1,,4] = as.matrix(ibm_sim$Tagging[[i]]$`tag_release_by_age-2005`)[1,]
  stan_data$tag_release_by_age[2,,5] = as.matrix(ibm_sim$Tagging[[i]]$`tag_release_by_age-2005`)[2,]
  stan_data$tag_release_by_age[3,,6] = as.matrix(ibm_sim$Tagging[[i]]$`tag_release_by_age-2005`)[3,]
  #(sum(stan_data$tag_release_by_age))
  
  stan_sim = optimizing(stan_sim_model,data=stan_data,init=init_pars,
                      iter=1,algorithm="LBFGS",
                      verbose=F,as_vector = FALSE)
  stan_en_tag_95[1,i,1] = sum(stan_sim$par$recapture_expectations[1,1,,stan_data$years == 1995])
  stan_en_tag_95[1,i,2] = sum(stan_sim$par$recapture_expectations[1,1,,stan_data$years == 1996])
  stan_en_tag_95[1,i,3] = sum(stan_sim$par$recapture_expectations[1,1,,stan_data$years == 1997])
  stan_en_tag_95[1,i,4] = sum(stan_sim$par$recapture_expectations[1,1,,stan_data$years == 1998])
  
  print(stan_en_tag_95[1,i,1])
  
  stan_en_tag_95[2,i,2] = sum(stan_sim$par$recapture_expectations[2,1,,stan_data$years == 1996])
  stan_en_tag_95[2,i,3] = sum(stan_sim$par$recapture_expectations[2,1,,stan_data$years == 1997])
  stan_en_tag_95[2,i,4] = sum(stan_sim$par$recapture_expectations[2,1,,stan_data$years == 1998])
  
  stan_en_tag_95[3,i,1] = sum(stan_sim$par$recapture_expectations[3,1,,stan_data$years == 1995])
  stan_en_tag_95[3,i,2] = sum(stan_sim$par$recapture_expectations[3,1,,stan_data$years == 1996])
  stan_en_tag_95[3,i,3] = sum(stan_sim$par$recapture_expectations[3,1,,stan_data$years == 1997])
  stan_en_tag_95[3,i,4] = sum(stan_sim$par$recapture_expectations[3,1,,stan_data$years == 1998])
}

sim_file_names = unique(sapply(strsplit(list.files(sim_obs_location), split = "\\."), "[", 1))
sim_file_names
extensions = unique(sapply(strsplit(list.files(sim_obs_location), split = "\\."), "[", 2))

en_tag_95 = hg_tag_95 = bp_tag_95 = array(NA, dim = c(3, length(extensions), 4))
en_tag_00 = hg_tag_00 = bp_tag_00 = array(NA, dim = c(3, length(extensions), 4))

for (i in 1:length(extensions)) {
  temp_en_en_tag_95 = extract.ibm.file(file = file.path(sim_obs_location, paste0("/tag_recapture_1995_EN_EN.",extensions[i])))
  temp_en_hg_tag_95 = extract.ibm.file(file = file.path(sim_obs_location, paste0("/tag_recapture_1995_EN_HG.",extensions[i])))
  temp_en_bp_tag_95 = extract.ibm.file(file = file.path(sim_obs_location, paste0("/tag_recapture_1995_EN_BP.",extensions[i])))
  
  en_tag_95[1,i,1] = sum(as.numeric(temp_en_en_tag_95$`observation[tag_recapture_1995_EN_EN_sim]`$Table$obs$`1995`))
  en_tag_95[1,i,2] = sum(as.numeric(temp_en_en_tag_95$`observation[tag_recapture_1995_EN_EN_sim]`$Table$obs$`1996`))
  en_tag_95[1,i,3] = sum(as.numeric(temp_en_en_tag_95$`observation[tag_recapture_1995_EN_EN_sim]`$Table$obs$`1997`))
  en_tag_95[1,i,4] = sum(as.numeric(temp_en_en_tag_95$`observation[tag_recapture_1995_EN_EN_sim]`$Table$obs$`1998`))
  
  en_tag_95[2,i,1] = sum(as.numeric(temp_en_hg_tag_95$`observation[tag_recapture_1995_EN_HG_sim]`$Table$obs$`1995`))
  en_tag_95[2,i,2] = sum(as.numeric(temp_en_hg_tag_95$`observation[tag_recapture_1995_EN_HG_sim]`$Table$obs$`1996`))
  en_tag_95[2,i,3] = sum(as.numeric(temp_en_hg_tag_95$`observation[tag_recapture_1995_EN_HG_sim]`$Table$obs$`1997`))
  en_tag_95[2,i,4] = sum(as.numeric(temp_en_hg_tag_95$`observation[tag_recapture_1995_EN_HG_sim]`$Table$obs$`1998`))
  
  en_tag_95[3,i,1] = sum(as.numeric(temp_en_bp_tag_95$`observation[tag_recapture_1995_EN_BP_sim]`$Table$obs$`1995`))
  en_tag_95[3,i,2] = sum(as.numeric(temp_en_bp_tag_95$`observation[tag_recapture_1995_EN_BP_sim]`$Table$obs$`1996`))
  en_tag_95[3,i,3] = sum(as.numeric(temp_en_bp_tag_95$`observation[tag_recapture_1995_EN_BP_sim]`$Table$obs$`1997`))
  en_tag_95[3,i,4] = sum(as.numeric(temp_en_bp_tag_95$`observation[tag_recapture_1995_EN_BP_sim]`$Table$obs$`1998`))
}

sum(ibm$tag_recapture_1995_EN$Values$expected[ibm$tag_recapture_1995_EN$Values$year == 1995 & ibm$tag_recapture_1995_EN$Values$cell == "EN"])
sum(stan_sim$par$recapture_expectations[1,1,,stan_data$years == 1995])
## plot it
## Release at 1995 in "EN"
jpeg(filename = "../Figures/tag_recaptures_EN_1995.jpeg", units = "in", res = 300, width = 6, height = 14)
par(mfrow = c(4,3), mar = c(3,3,1,1), oma = c(2,2,2,0))
## 1995
boxplot(cbind(stan_en_tag_95[1,,1],en_tag_95[1,,1]), main = "recapture EN", xaxt = "n")
abline(h = sum(ibm$tag_recapture_1995_EN$Values[ibm$tag_recapture_1995_EN$Values$year == 1995 & ibm$tag_recapture_1995_EN$Values$cell == "EN","expected"]), col = "red", lty = 2, lwd = 2)

boxplot(cbind(stan_en_tag_95[2,,1],en_tag_95[2,,1]), main = "recapture HG", xaxt = "n")
abline(h = sum(ibm$tag_recapture_1995_EN$Values[ibm$tag_recapture_1995_EN$Values$year == 1995 & ibm$tag_recapture_1995_EN$Values$cell == "HG","expected"]), col = "red", lty = 2, lwd = 2)

boxplot(cbind(stan_en_tag_95[3,,1],en_tag_95[3,,1]), main = "recapture BP", xaxt = "n")
abline(h = sum(ibm$tag_recapture_1995_EN$Values[ibm$tag_recapture_1995_EN$Values$year == 1995 & ibm$tag_recapture_1995_EN$Values$cell == "BP","expected"]), col = "red", lty = 2, lwd = 2)
## 1996
boxplot(cbind(stan_en_tag_95[1,,2],en_tag_95[1,,2]), main = "", xaxt = "n")
abline(h = sum(ibm$tag_recapture_1995_EN$Values[ibm$tag_recapture_1995_EN$Values$year == 1996 & ibm$tag_recapture_1995_EN$Values$cell == "EN","expected"]), col = "red", lty = 2, lwd = 2)

boxplot(cbind(stan_en_tag_95[2,,2],en_tag_95[2,,2]), main = "", xaxt = "n")
abline(h = sum(ibm$tag_recapture_1995_EN$Values[ibm$tag_recapture_1995_EN$Values$year == 1996 & ibm$tag_recapture_1995_EN$Values$cell == "HG","expected"]), col = "red", lty = 2, lwd = 2)

boxplot(cbind(stan_en_tag_95[3,,2],en_tag_95[3,,2]), main = "", xaxt = "n")
abline(h = sum(ibm$tag_recapture_1995_EN$Values[ibm$tag_recapture_1995_EN$Values$year == 1996 & ibm$tag_recapture_1995_EN$Values$cell == "BP","expected"]), col = "red", lty = 2, lwd = 2)

## 1997
boxplot(cbind(stan_en_tag_95[1,,3],en_tag_95[1,,3]), main = "")
abline(h = sum(ibm$tag_recapture_1995_EN$Values[ibm$tag_recapture_1995_EN$Values$year == 1997 & ibm$tag_recapture_1995_EN$Values$cell == "EN","expected"]), col = "red", lty = 2, lwd = 2)

boxplot(cbind(stan_en_tag_95[2,,3],en_tag_95[2,,3]), main = "")
abline(h = sum(ibm$tag_recapture_1995_EN$Values[ibm$tag_recapture_1995_EN$Values$year == 1997 & ibm$tag_recapture_1995_EN$Values$cell == "HG","expected"]), col = "red", lty = 2, lwd = 2)

boxplot(cbind(stan_en_tag_95[3,,3],en_tag_95[3,,3]), main = "")
abline(h = sum(ibm$tag_recapture_1995_EN$Values[ibm$tag_recapture_1995_EN$Values$year == 1997 & ibm$tag_recapture_1995_EN$Values$cell == "BP","expected"]), col = "red", lty = 2, lwd = 2)
## 1998
boxplot(cbind(stan_en_tag_95[1,,4],en_tag_95[1,,4]), main = "")
abline(h = sum(ibm$tag_recapture_1995_EN$Values[ibm$tag_recapture_1995_EN$Values$year == 1998 & ibm$tag_recapture_1995_EN$Values$cell == "EN","expected"]), col = "red", lty = 2, lwd = 2)

boxplot(cbind(stan_en_tag_95[2,,4],en_tag_95[2,,4]), main = "")
abline(h = sum(ibm$tag_recapture_1995_EN$Values[ibm$tag_recapture_1995_EN$Values$year == 1998 & ibm$tag_recapture_1995_EN$Values$cell == "HG","expected"]), col = "red", lty = 2, lwd = 2)

boxplot(cbind(stan_en_tag_95[3,,4],en_tag_95[3,,4]), main = "")
abline(h = sum(ibm$tag_recapture_1995_EN$Values[ibm$tag_recapture_1995_EN$Values$year == 1998 & ibm$tag_recapture_1995_EN$Values$cell == "BP","expected"]), col = "red", lty = 2, lwd = 2)
mtext(outer = T, line = -0.2, font = 2, text = "1998", las = 3, adj = 0.125, side = 2)
mtext(outer = T, line = -0.2, font = 2, text = "1997", las = 3, adj = 0.40, side = 2)
mtext(outer = T, line = -0.2, font = 2, text = "1996", las = 3, adj = 0.67, side = 2)
mtext(outer = T, line = -0.2, font = 2, text = "1995", las = 3, adj = 0.94, side = 2)
dev.off()

###
sum(stan_sim$par$recapture_expectations[1,1,,6])
sum(stan_sim$par$recapture_expectations[1,1,,7])

stan_data$tag_release_by_age
