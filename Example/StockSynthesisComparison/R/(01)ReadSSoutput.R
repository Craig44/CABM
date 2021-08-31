#'
#'
#' (01)ReadSSoutput.R
#'
#'
#'
library(ggplot2)
library(ibm)
library(r4ss) ## remotes::install_github("r4ss/r4ss")
library(dplyr)
library(reshape2)
library(stockassessmenthelper) # devtools::install_github("Craig44/stockassessmenthelper", build_vignettes  = TRUE)
ss_col ="red"
abm_col = "black"
abm_dir = file.path("..","abm_w_length_age_sel")
ss_dir = file.path("..","SS")
fig_path = file.path("..","Figures")
#abm = extract.run(file = "run.out", path = abm_dir)
abm = extract.run(file = "run.out", path = abm_dir)
fishery_labs = colnames(abm$fishing$`actual_catches-1931`)


ss_run = SS_output(dir = ss_dir)
ss_ages = 0:80
abm_ages = abm$model_attributes$ages

ss_run$parameters$Active_Cnt
ss_run$Pstar_sigma
ss_run$Spawn_timing_in_season
ss_run$SpawnOutputUnits
ss_run$recruit

ssb_run = ss_run$recruit

## plot SSB
png(file.path(fig_path,"SSB.png"), units = "in", width = 6, height = 6, res = 250)
par(mfrow =c (1,1))
plot(ssb_run$Yr, ssb_run$SpawnBio, xlab ="Year", ylab = "SSB (t)", type = "l", lwd = 4, lty= 1, ylim = c(0, 2e5), xlim = c(1931, 2020), col = ss_col)
#lines(ssb_est$Yr, ssb_est$SpawnBio, lwd = 4, lty = 2, col = 'blue')
lines(names(abm$derived_quants$SSB$values), abm$derived_quants$SSB$values, lwd = 3, lty = 2, col = abm_col)
legend('topright', legend = c("SS", "ABM"), lwd = 3, col = c(ss_col, abm_col))
dev.off()

growdat_seas1 <- ss_run$endgrowth[ss_run$endgrowth$Seas == 1, ]
growdat_seas2 <- ss_run$endgrowth[ss_run$endgrowth$Seas == 2, ]
head(ss_run$endgrowth)

## age-length
abm_step1 = abm$agents_step1$'1990'$values
abm_step2 = abm$agents_step2$'1990'$values
abm_step3 = abm$agents_step3$'1990'$values
ndx = abm_step1$age == 1
abm_step1$length[ndx]


## length at age
png(file.path(fig_path, "growth.png"), units = "in", width = 6, height = 8, res = 250)
par(mfrow = c(2,1), mar = c(3,4,1,1))
plot(growdat_seas1$Age_Beg, growdat_seas1$Len_Beg, type = "l", lwd = 3, xlab = "", ylab = "", main = "Beginning", col = ss_col, lty = 1)
points(abm_step1$age, abm_step1$length, pch = 16, col = alpha(abm_col, alpha = 0.3))
lines(growdat_seas1$Age_Beg, vonbert(growdat_seas1$Age_Beg, K = 0.146478, L_inf = 57.4842, t0 = -1), col = "blue", lwd = 3, lty =3)
plot(growdat_seas1$Age_Mid, growdat_seas1$Len_Mid, type = "l", lwd = 3, xlab = "", ylab = "", main = "Middle", col = ss_col, lty = 1)
points(abm_step2$age + 0.5, abm_step2$length, pch = 16, col =  alpha(abm_col, alpha = 0.3))
lines(growdat_seas1$Age_Beg + 0.5, vonbert(growdat_seas1$Age_Beg + 0.5, K = 0.146478, L_inf = 57.4842, t0 = -1), col = "blue", lwd = 3, lty =3)
mtext(side = 2, text = "Length (cm)", adj = 0.5, outer = T, line = -2, cex =1.5)
mtext(side = 1, text = "Age", adj = 0.56, outer = T, line = -2, cex =1.5)
dev.off()


png(file.path(fig_path, "weight_at_age.png"), units = "in", width = 6, height = 8, res = 250)
par(mfrow = c(2,1), mar = c(3,4,1,1))
plot(growdat_seas1$Age_Beg, growdat_seas1$Wt_Beg, type = "l", lwd = 3, xlab = "", ylab = "", main = "Beginning", col = ss_col, lty = 1)
points(abm_step1$age, abm_step1$weight * 1000, pch = 16, col = alpha(abm_col, alpha = 0.3))
plot(growdat_seas1$Age_Mid, growdat_seas1$Wt_Mid, type = "l", lwd = 3, xlab = "", ylab = "", main = "Middle", col = ss_col, lty = 1)
points(abm_step2$age + 0.5, abm_step2$weight * 1000, pch = 16, col =  alpha(abm_col, alpha = 0.3))
mtext(side = 2, text = "Weight (kgs)", adj = 0.5, outer = T, line = -2, cex =1.5)
mtext(side = 1, text = "Age", adj = 0.56, outer = T, line = -2, cex =1.5)
dev.off()

cbind(growdat_seas1$Len_Mid, growdat_seas1$Len_Beg)

weight_at_age = tapply(abm_step2$weight, INDEX = list(abm_step2$age), FUN = mean)
ndx = match(table = names(weight_at_age), x = ss_ages)
weight_at_age

## compare weight at age for many fisheries and proceses
weights = cbind(growdat_seas1$`Mat*Fecund`, growdat_seas1$Wt_Beg, growdat_seas1$Wt_Mid, growdat_seas1$'SelWt:_1', growdat_seas1$'SelWt:_2')
colnames(weights) = c("SSB weights","Beginning","Mid","Fishery1","Fishery2")
head(weights)
## weight at age
plot(ss_ages, ss_run$wtatage[1,7:ncol(ss_run$wtatage)], col = ss_col, xlab = "Age", ylab = "Weight (kgs)")
plot(ss_ages, ss_run$wtatage[1,7:ncol(ss_run$wtatage)], col = ss_col, xlab = "Age", ylab = "Weight (kgs)")

points(abm$agents_step2$'1931'$values$age, abm$agents_step2$'1931'$values$weight * 1000, pch = 16, col = abm_col)
points(abm$agents_step1$'1931'$values$age, abm$agents_step1$'1931'$values$weight * 1000, pch = 16, col = "blue")
points(abm$agents_step3$'1931'$values$age, abm$agents_step3$'1931'$values$weight * 1000, pch = 16, col = "purple")

## weight at age by fishery
plot(ss_ages, growdat_seas1$'SelWt:_1', type = "l", lwd = 3, col = ss_col, xlab = "", ylab = "")
points(abm$agents_step2$'1931'$values$age, abm$agents_step2$'1931'$values$weight * 1000, col = abm_col, pch = 16)
growdat_seas1$'SelWt:_2'

## redo equation 1.6 length in plus group
A = max(ss_ages)
length_plus_group = 0
L_inf = 57.4842
lmin = 1
a3 = 0
k = 0.146478
init_length_at_age = vector()
init_length_at_age[1] = lmin + 0


####################
## Selectivities
####################
##
len_sel = ss_run$sizeselex

ndx = len_sel$Fleet %in% which (ss_run$FleetNames %in% c("FISHERY_Rec1","FISHERY_Rec2")) & len_sel$Yr == 1931
len_of_interest = len_sel %>% filter(ndx, Factor == "Lsel")
len_of_interest = len_of_interest[,6:ncol(len_of_interest)]
png(file.path(fig_path, "Selectivities.png"), units = "in", width = 8, height = 8, res = 250)
par(mfrow = c(2,3))
plot(colnames(len_of_interest), len_of_interest[1,], type = "l", lwd = 3, xlab = "", ylab = "", main = "RECO", col = ss_col, lty = 1)
lines(as.numeric(names(abm$sel_RECO$Values)), abm$sel_RECO$Values, col = abm_col, lwd =3, lty =2 )
#lines(colnames(len_of_interest_est), len_of_interest_est[1,], type = "l", lwd = 3, lty = 3, col = "blue")

plot(colnames(len_of_interest), len_of_interest[2,], type = "l", lwd = 3, xlab = "", ylab = "", col = ss_col, main ="RECI")
lines(names(abm$sel_RECO$Values), abm$sel_RECI$Values, col = abm_col, lwd =3, lty =2 )
legend('topright', legend = c("SS", "ABM"), lwd = 3, col = c(ss_col, abm_col))

age_sel = ss_run$ageselex
ndx = age_sel$Fleet %in% which (ss_run$FleetNames %in% c("FISHERY_BT","FISHERY_BPT","FISHERY_JP")) & age_sel$Yr == 1931
age_of_interest = age_sel %>% filter(ndx, Factor == "Asel")
age_of_interest = age_of_interest[,8:ncol(age_of_interest)]

plot(colnames(age_of_interest), age_of_interest[1,], type = "l", lwd = 3, xlab = "", ylab = "", col = ss_col, main = "BT")
lines(names(abm$sel_BT$Values), abm$sel_BT$Values, col = abm_col, lwd =3, lty =2 )
#lines(colnames(age_of_interest_est), age_of_interest_est[1,], col = "blue", lwd =3, lty =3 )

plot(colnames(age_of_interest), age_of_interest[2,], type = "l", lwd = 3, xlab = "", ylab = "", col = ss_col, main = "BPT")
lines(names(abm$sel_BPT$Values), abm$sel_BPT$Values, col = abm_col, lwd =3, lty =2 )

plot(colnames(age_of_interest), age_of_interest[3,], type = "l", lwd = 3, xlab = "", ylab = "", col = ss_col, main = "JPL")
lines(names(abm$sel_JPL$Values), abm$sel_JPL$Values, col = abm_col, lwd =3, lty =2 )

plot(growdat_seas1$Real_Age, growdat_seas1$Age_Mat, type = "l", lwd = 3, xlab = "", ylab = "", col = ss_col, main = "Maturity")
lines(names(abm$mature_ogive$Values), abm$mature_ogive$Values, col = abm_col, lwd =3, lty =2 )
dev.off()
write.table(x = len_of_interest, file = file.path(abm_dir, "length_sel.txt"))

########## Recruit devs
ss_run$MGparmAdj$Eggs_scalar_Fem_GP_1
ss_run$MGparmAdj$Eggs_scalar_Fem_GP_1

devs = log(c(abm$Rec$ycs_values[2:length(abm$Rec$ycs_values)])) + 0.5*0.6^2
cbind(devs, ssb_run$dev[1:length(devs)])
names(devs) = paste0(1931:2020,"R")



new_dev = vector()
counter = 1;
for(i in 1:length(devs)) {
  new_dev[counter] = names(devs)[i]
  counter = counter + 1
  new_dev[counter] = devs[i]
  counter = counter + 1
}
write.table(t(new_dev), file = file.path(abm_dir, "devs.txt"), quote = F, col.names = F, row.names = F)

ssb_run$biasadjuster
cbind(ssb_run$raw_dev, ssb_run$dev)
ssb_run$exp_recr * exp(ssb_run$dev - 0.5 * 0.6^2)
ssb_run$pred_recr
## check the bias is being taken into account
cbind(ssb_run$exp_recr * exp(ssb_run$dev - 0.5 * 0.6^2), ssb_run$pred_recr)
## compare with abm 0 note abm starts at age 1

ss_recruit_years = 1931:2020 ## ss recruits at age 0
abm_recruit_years = 1930:2020 ## ss recruits at age 0


png(file.path(fig_path, "Recruitment.png"), units = "in", width = 6, height = 6, res = 250)
par(mfrow = c(1,1))
plot(ss_recruit_years, (ssb_run$pred_recr * exp(-0.075))[1:length(ss_recruit_years)], xlab = "time", col= ss_col, ylab = "recruits (000's)", type = "l", lwd = 3)
lines(abm_recruit_years, (abm$Rec$recruits * abm$model_attributes$Recruitment) / 1000, lty = 2, lwd = 3, col= abm_col)
legend('topleft', legend = c("SS", "ABM"), lwd = 3, col = c(ss_col, abm_col))
dev.off()
########## SSB



####################
## Catch
###################
fishery_lab = substring(names(abm$fishing)[grepl(pattern = "age_freq-", names(abm$fishing))], first = 10)
ss_dat$FleetNames
fishery_lab

ss_catch = ss_run$catch
ss_catch$Fleet_Name
ss_catch$Fishery = fishery_lab[ss_catch$Fleet]
ss_catch$model = "SS"
# ABM
catch_ndx = grepl(names(abm$fishing), pattern = "actual_catches")
fishery_labs = colnames(abm$fishing$`actual_catches-1931`)
abm_catches = Reduce(rbind, lapply(abm$fishing[which(catch_ndx)], FUN = 
                                     function(x) { as.vector(t(x)) }))
rownames(abm_catches) = substring(names(abm$fishing)[which(catch_ndx)], first = 16)
colnames(abm_catches) = fishery_labs
abm_catch_df = melt(t(abm_catches))
colnames(abm_catch_df) = c("Fishery", "Yr", "Exp")
abm_catch_df$model = "ABM"
abm_catch_df$Fishery = as.character(abm_catch_df$Fishery)


catch_df = rbind(ss_catch[,c("Fishery", "Yr", "Exp", "model")], abm_catch_df)

catch_df$Year = as.integer(catch_df$Yr)
ggplot(catch_df, aes(x = Year, y = Exp, col = model, linetype = model)) +
  geom_line(size = 1.5) + facet_grid(~ Fishery) +
  ylim(0,5000) +
  ylab("Expected Catch")
ggsave(filename = file.path(fig_path, "catch.png"),width = 9, height = 6)

names(ss_run)
head(ss_run$derived_quants)
head(ss_run$fatage)
max(ss_run$fatage[,8:ncol(ss_run$fatage)])

#############################
## Numbers at age 
##############################
Nage = ss_run$natage[,13:ncol(ss_run$natage)]
table(ss_run$natage_annual_1_no_fishery$Yr == 1931)

## recalculate ssb
# B0
sum(Nage[1,] * growdat_seas1$Age_Mat *  4.467e-05 * growdat_seas1$Len_Beg^2.793) 
sum(Nage[1,] * ss_run$endgrowth$Mat_F_wtatage)
ssb_run$SpawnBio[1]
cbind(growdat_seas1$Age_Mat *  4.467e-05 * ss_run$endgrowth$Len_Beg^2.793, ss_run$endgrowth$Mat_F_wtatage)
head(ss_run$endgrowth)

abm$Rec$B0

## plot SSB
par(mfrow =c (1,1))
plot(ssb_run$Yr, ssb_run$SpawnBio, xlab ="Year", ylab = "SSB (t)", type = "l", lwd = 4, lty= 1, ylim = c(0, 2e5), xlim = c(1931, 2020), col = ss_col)
#lines(ssb_est$Yr, ssb_est$SpawnBio, lwd = 4, lty = 2, col = 'blue')
lines(names(abm$derived_quants$SSB$values), abm$derived_quants$SSB$values, lwd = 3, lty = 2, col = abm_col)

abline(h = sum(Nage[1,] * growdat_seas1$Age_Mat *  4.467e-05 * growdat_seas1$Len_Beg^2.793) , col = "blue", lwd =3, lty = 2)
legend('topright', legend = c("SS", "ABM"), lwd = 3, col = c(ss_col, abm_col))


## look at relative AF
par(mfrow = c(1,1))
plot(names(abm$age_freq_ssb$`1931`$values), abm$age_freq_ssb$`1931`$values / 1000, type = "l", lwd = 3, ylab = "Absolute AF")
lines(names(Nage[1,2:ncol(Nage)]), Nage[1,2:ncol(Nage)], lty = 2, col = "red", lwd = 3)
lines(1:80, ss_run$natage_annual_1_no_fishery[1,5:ncol(ss_run$natage_annual_1_no_fishery)] , col = "blue", lwd = 3, lty = 3)
lines(1:80, ss_run$natage_annual_2_with_fishery[1,5:ncol(ss_run$natage_annual_2_with_fishery)], col = "purple", lwd = 3, lty = 3)

plot(names(abm$age_freq_ssb$`1931`$values), abm$age_freq_ssb$`1931`$values / sum(abm$age_freq_ssb$`1931`$values), type = "l", lwd = 3, ylab = "relative AF")
lines(names(Nage[1,2:ncol(Nage)]), Nage[1,2:ncol(Nage)] / sum(Nage[1,2:ncol(Nage)]), lty = 2, col = "red", lwd = 3)
lines(1:80, ss_run$natage_annual_1_no_fishery[1,5:ncol(ss_run$natage_annual_1_no_fishery)] / sum(ss_run$natage_annual_1_no_fishery[1,5:ncol(ss_run$natage_annual_1_no_fishery)]), col = "blue", lwd = 3, lty = 3)
lines(1:80, ss_run$natage_annual_2_with_fishery[1,5:ncol(ss_run$natage_annual_2_with_fishery)] / sum(ss_run$natage_annual_2_with_fishery[1,5:ncol(ss_run$natage_annual_2_with_fishery)]), col = "purple", lwd = 3, lty = 3)


#################
## Age comp
##############
head(ss_run$catage)
ss_catch_at = ss_run$catage
## remove uneccessary cols
ss_catch_at = ss_catch_at[,c("Fleet","Yr", 0:80)]
head(ss_catch_at)
ss_catch_at$fishery = fishery_labs[ss_catch_at$Fleet]
ss_lng = melt(ss_catch_at, id.vars = c("Fleet", "Yr","fishery"))
drops <- c("Fleet")
ss_lng = ss_lng[,!(colnames(ss_lng) %in% drops)]
head(ss_lng)
colnames(ss_lng)= c("year", "fishery", "age", "numbers")
ss_lng$numbers = ss_lng$numbers * 1000
ss_lng$model = "SS"
# reformat abm removals by age and gear
full_abm_AF = NULL
for(f in 1:length(fishery_labs)) {
  this_fish_AF = get(paste0("age_freq-",fishery_labs[f]),abm$fishing)
  tot_ndx = this_fish_AF$cell == "Total"
  this_fish_AF = this_fish_AF[tot_ndx, ]
  ## drop cell column
  drops <- c("cell")
  this_fish_AF = this_fish_AF[,!(colnames(this_fish_AF) %in% drops)]
  this_AF_long = melt(this_fish_AF,  id.vars = c("year"))
  colnames(this_AF_long) = c("year", "age", "numbers")
  ## truncate now to have the same plus group as SS
  this_AF_long$fishery = fishery_labs[f]
  full_abm_AF= rbind(full_abm_AF, this_AF_long)
}
dim(full_abm_AF)
full_abm_AF$model = "abm"
##
temp_AF = abm$age_freq_ssb$'1931'$values * abm$sel_BT$Values

dim(abm$fishing$`census_info-JPL-1965-1-1`)
abm$fishing$`census_info-JPL-1965-1-1`[,1:5]
table(abm$fishing$`census_info-JPL-1965-1-1`[1,])

## merge them for plotting
full_catch_at_age = rbind(ss_lng, full_abm_AF)
full_catch_at_age$year = factor(full_catch_at_age$year, ordered = T)
ggplot(full_catch_at_age %>% filter(year %in% 1975:1975, fishery == "BT"), aes(x = age, y = numbers, col = model, group = interaction(year,model), fill = model)) +
  geom_line(size = 2) +
  facet_wrap(~year)

ggplot(full_catch_at_age %>% filter(year %in% 1975:1975, fishery == "JPL"), aes(x = age, y = numbers, col = model, group = interaction(year,model), fill = model)) +
  geom_line(size = 2) +
  facet_wrap(~year)
