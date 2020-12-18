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
## abm output
extension = ".1"
abm = extract.run(file = "mse_test.out",path = file.path("..","ibm"))
fishingproc = extract.run(file = file.path(abm_sim_dir, paste0("FishingProc", extension)))
sim_survey_bio = extract.ibm.file(file = file.path(abm_sim_dir, paste0("chatTANbiomass", extension)), quiet  = T)
sim_survey_age = extract.ibm.file(file = file.path(abm_sim_dir, paste0("chatTANage", extension)), quiet  = T)
sim_fishery_age = extract.ibm.file(file = file.path(abm_sim_dir, paste0("FishingAge", extension)), quiet  = T)

## read TMB info
mse_years = unique(substring(list.files(est_out_dir), first = 10, last = 13))
mse_extensions = unique(Reduce(c,lapply(strsplit(list.files(est_out_dir), split = "\\."), FUN = function(x){x[2]})))
mse_est = list()  ## dims first is extensions then its years with estimated TMB objects with in
mse_data = list()  ## dims first is extensions then its years with estimated TMB objects with in

for(i in 1:length(mse_extensions)) {
  this_extension = mse_extensions[i]
  mse_est[[i]] = list()
  mse_data[[i]] = list()
  for(j in 1:length(mse_years)) {
    current_year = mse_years[j]
    this_obj = load(file = file.path(est_out_dir, paste0("mse_year_",current_year, ".", this_extension, ".RData")))
    mse_est[[i]][[j]] = MLE_est
    mse_data[[i]][[j]] = TMB_data
  }
}

## check simulated data isn't changing in the past.

## look at estimated SSB
init_year = min(abm$model_attributes$model_years) - 1
plot(init_year:mse_years[1], mse_est[[1]][[1]]$ssb, type = "l", lwd =3, col = "black", xlab = "Time", ylab = "SSB", ylim = c(0,40000), xlim = c(init_year, as.numeric(max(mse_years))))
for(i in 2:length(mse_est[[1]])) {
  lines(init_year:mse_years[i], mse_est[[1]][[i]]$ssb, lwd = 3, col = "blue", lty = 2)
}
## get the truth
true_ssb = abm$derived_quants$SSB$values
lines(names(true_ssb), true_ssb, lwd = 3, col = "red", lty = 3)

