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
abm = extract.run(file = "mse_test.out",path = file.path("..","ibm"))
fishingproc = extract.run(file = file.path(abm_sim_dir, paste0("FishingProc", extension)))
sim_survey_bio = extract.ibm.file(file = file.path(abm_sim_dir, paste0("chatTANbiomass", extension)), quiet  = T)
sim_survey_age = extract.ibm.file(file = file.path(abm_sim_dir, paste0("chatTANage", extension)), quiet  = T)
sim_fishery_age = extract.ibm.file(file = file.path(abm_sim_dir, paste0("FishingAge", extension)), quiet  = T)

## read TMB info
extension = ".1"
mse_years = unique(substring(list.files(est_out_dir), first = 10, last = 13))
mse_extensions = unique(Reduce(c,lapply(strsplit(list.files(est_out_dir), split = "\\."), FUN = function(x){x[2]})))
mse_est = list()  ## dims first is extensions then its years with estimated TMB objects with in
for(i in 1:length(mse_extensions)) {
  this_extension = mse_extensions[i]
  mse_est[[i]] = list()
  for(j in 1:length(mse_years)) {
    current_year = mse_years[j]
    this_obj = load(file = file.path(est_out_dir, paste0("mse_year_",current_year, ".", this_extension, ".RData")))
    mse_est[[i]][[j]] = this_obj
  }
}

current_year = 2016
load(file = file.path(est_out_dir, paste0("mse_year_",current_year, extension, ".RData")))
est_2016 = MLE_est
current_year = 2020
load(file = file.path(est_out_dir, paste0("mse_year_",current_year, extension, ".RData")))
est_2020 = MLE_est
current_year = 2024
load(file = file.path(est_out_dir, paste0("mse_year_",current_year, extension, ".RData")))
est_2024 = MLE_est

init_year = 1974
plot(init_year:2016, est_2016$ssb, type = "l", lwd =3, col = "black", xlab = "Time", ylab = "SSB", ylim = c(0,30000), xlim = c(init_year, 2036))
lines(init_year:2020, est_2020$ssb, lwd = 3, col = "blue", lty = 2)


