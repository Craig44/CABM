#'
#' Summarise/Debug MSE output
#'
# assumes getwd() is in the R directory
model_label = "M1" # used for saving objects
## everything will be relative to this
r_dir = getwd() 
abm_dir = file.path(r_dir, "..","ibm") ## assuming you are running the ibm model in the ibm directory 
abm_sim_dir = file.path(abm_dir, "sim")## assuming you are running the ibm model in the ibm directory 
#setwd(r_dir)
tmb_dir = file.path(abm_dir, "..","TMB")
abm_sim_dir = file.path(abm_dir, "sim")
savedobjs_dir = file.path(r_dir, "savedOBJs")
est_out_dir = file.path(r_dir, "est_outputs", model_label)
## source auxillary functions
ls_files = list.files(file.path(r_dir, "Funs"), pattern = "\\.R$")
for(i in 1:length(ls_files))
  source(file.path(r_dir,"Funs", ls_files[i]))
ls_files = list.files(file.path(r_dir, "Funs", "Fisheries"), pattern = "\\.R$")
for(i in 1:length(ls_files))
  source(file.path(r_dir, "Funs","Fisheries", ls_files[i]))

## neccessary libraries
library(ibm)
library(mvtnorm)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
## abm output
extension = ".1"
abm = extract.run(file = "mse.out",path = file.path("..","ibm"))
abm$Rec$`1`$actual_YCS
fishingproc = extract.run(file = file.path(abm_sim_dir, paste0("FishingProc", extension)))
sim_survey_bio = extract.ibm.file(file = file.path(abm_sim_dir, paste0("chatTANbiomass", extension)), quiet  = T)
sim_survey_age = extract.ibm.file(file = file.path(abm_sim_dir, paste0("chatTANage", extension)), quiet  = T)
sim_fishery_age = extract.ibm.file(file = file.path(abm_sim_dir, paste0("FishingAge", extension)), quiet  = T)

## read TMB info
mse_years = unique(substring(list.files(est_out_dir, pattern = "mse_year"), first = 10, last = 13))
mse_extensions = unique(Reduce(c,lapply(strsplit(list.files(est_out_dir, pattern = "mse_year"), split = "\\."), FUN = function(x){x[2]})))
mse_est = list()  ## dims first is extensions then its years with estimated TMB objects with in
mse_data = list()  ## dims first is extensions then its years with estimated TMB objects with in
SSB_df = catch_df = F_df = NULL

for(i in 1:length(mse_extensions)) {
  this_extension = mse_extensions[i]
  mse_est[[i]] = list()
  mse_data[[i]] = list()
  for(j in 1:length(mse_years)) {
    current_year = mse_years[j]
    this_obj = load(file = file.path(est_out_dir, paste0("mse_year_",current_year, ".", this_extension, ".RData")))
    mse_est[[i]][[j]] = MLE_est
    mse_data[[i]][[j]] = TMB_data
    this_ssb = data.frame(run = rep(i, length(MLE_est$ssb)), ssb = MLE_est$ssb, ass_year = rep(mse_years[j], length(MLE_est$ssb)), year = c(min(TMB_data$years),TMB_data$years))
    this_catch = data.frame(run = rep(i, length(TMB_data$catches)), catch = TMB_data$catches, ass_year = rep(mse_years[j], length(TMB_data$catches)), year = TMB_data$years)
    this_F = data.frame(run = rep(i, length(MLE_est$annual_F)), est_F = MLE_est$annual_F, ass_year = rep(mse_years[j], length(MLE_est$annual_F)), year = TMB_data$years)
    SSB_df = rbind(SSB_df, this_ssb)
    catch_df = rbind(catch_df, this_catch)
    F_df = rbind(F_df, this_F)
  }
}

SSB_df$run = factor(SSB_df$run)
SSB_df$ass_year = factor(SSB_df$ass_year)
catch_df$run = factor(catch_df$run)
catch_df$ass_year = factor(catch_df$ass_year)
F_df$run = factor(F_df$run)
F_df$ass_year = factor(F_df$ass_year)

## check simulated data isn't changing in the past.
ggplot(SSB_df, aes(x = year, y = ssb, col = run, linetype = ass_year)) +
  geom_line(size = 1.4) +
  ylim(0,50000)
ggplot(catch_df, aes(x = year, y = catch, col = run, linetype = ass_year)) +
  geom_line(size = 1.4) +
  ylim(0,6000)
ggplot(F_df, aes(x = year, y = est_F, col = run, linetype = ass_year)) +
  geom_line(size = 1.4) +
  ylim(0,0.5)
## look at estimated SSB
init_year = min(abm$model_attributes$`1`$model_years) - 1
plot(init_year:mse_years[1], mse_est[[1]][[1]]$ssb, type = "l", lwd =3, col = "black", xlab = "Time", ylab = "SSB", ylim = c(0,40000), xlim = c(init_year, as.numeric(max(mse_years))))

for(i in 1:length(mse_extensions)) {
  for(j in 1:length(mse_est[[1]])) {
    lines(init_year:mse_years[i], mse_est[[1]][[i]]$ssb, lwd = 3, col = "blue", lty = 2)
  }
}
## get the truth
true_ssb = abm$derived_quants$`1`$SSB$values
lines(names(true_ssb), true_ssb, lwd = 3, col = "red", lty = 3)
names(mse_data[[1]][[1]])

## Check that we aren't re-simulating over old observations.
mse_data[[1]][[1]]$survey_obs
mse_data[[1]][[2]]$survey_obs
mse_data[[1]][[3]]$survey_obs
mse_data[[1]][[4]]$survey_obs

cbind(mse_data[[1]][[1]]$survey_trans_comp[,1],
mse_data[[1]][[2]]$survey_trans_comp[,1],
mse_data[[1]][[3]]$survey_trans_comp[,1],
mse_data[[1]][[4]]$survey_trans_comp[,1])

cbind(mse_data[[1]][[1]]$survey_trans_comp[,42],
      mse_data[[1]][[2]]$survey_trans_comp[,42],
      mse_data[[1]][[3]]$survey_trans_comp[,42],
      mse_data[[1]][[4]]$survey_trans_comp[,42])
cbind(mse_data[[1]][[1]]$fishery_trans_comp[,42],
      mse_data[[1]][[2]]$fishery_trans_comp[,42],
      mse_data[[1]][[3]]$fishery_trans_comp[,42],
      mse_data[[1]][[4]]$fishery_trans_comp[,42])

## look at the Bias if any in the estimation
mated




