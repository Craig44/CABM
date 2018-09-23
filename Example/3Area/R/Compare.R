## Read in Casal2 output
detach("package:casal2", unload=TRUE)
detach("package:ibm", unload=TRUE)

## this script is just comparing, steady state models just constant recruitment and M
library(ggplot2)
library(reshape2)
library(cyrils)
library(casal)
setwd("../CASAL")

cas = extract.mpd("run_b0.log")
cas_dq = cas$quantities$SSBs

setwd("../Casal2")
library(casal2)

cas2 = extract.mpd("run_original.log")

cas2_dq = plot.derived_quantities(cas2, report_label = "derived_quants", plot.it = F)


## compare CASAL with Casal2
years = rownames(cas2_dq)
par(mfrow = c(1,3))
plot(years, cas2_dq[,"SSB_EN"], xlab = "years", ylab = "SSB (t)", main = "EN", type = "l", lwd = 2, ylim = c(0,200000))
lines(years, cas_dq$EN, lwd = 2, col = "red",lty = 2)

plot(years, cas2_dq[,"SSB_HG"], xlab = "years", ylab = "SSB (t)", main = "HG", type = "l", lwd = 2, ylim = c(0,150000))
lines(years, cas_dq$HAGU, lwd = 2, col = "red",lty = 2)

plot(years, cas2_dq[,"SSB_BP"], xlab = "years", ylab = "SSB (t)", main = "BP", type = "l", lwd = 2, ylim = c(0,60000))
lines(years, cas_dq$BOP, lwd = 2, col = "red",lty = 2)

## lets look at B0
cas$quantities$B0
c(cas2$Rec_EN$`1`$b0, cas2$Rec_HG$`1`$b0,cas2$Rec_BP$`1`$b0)

## lets look at R0
cas$quantities$R0
c(cas2$Rec_EN$`1`$r0, cas2$Rec_HG$`1`$r0,cas2$Rec_BP$`1`$r0)

## lets look at F
rbind(cas$quantities$fishing_pressures$BP_LL,cas2$instantaneous_mort$`1`$`fishing_pressure[BP_LL]`)
rbind(cas$quantities$actual_catches$BP_LL,cas2$instantaneous_mort$`1`$`actual_catch[BP_LL]`)

rbind(cas$quantities$actual_catches$HG_LL,cas2$instantaneous_mort$`1`$`actual_catch[HG_LL]`)
rbind(cas$quantities$fishing_pressures$HG_LL,cas2$instantaneous_mort$`1`$`fishing_pressure[HG_LL]`)

rbind(cas$quantities$actual_catches$EN_LL,cas2$instantaneous_mort$`1`$`actual_catch[EN_LL]`)
rbind(cas$quantities$fishing_pressures$EN_LL,cas2$instantaneous_mort$`1`$`fishing_pressure[EN_LL]`)


cbind(cas$quantities$`Ogive parameter values`$`selectivity[Sel_LL].all`, cas2$Sel_LL$`1`$Values)

cas2$mean_weight_1$`1`$EN.EN$mean_weights$values
cas_en = 4.467e-008 * c(12.55  ,   16.43 ,    19.92   ,  23.07 ,    25.91  ,   28.46  ,   30.77 ,    32.85  ,   34.72 ,    36.41,     37.93,     39.31,     40.54,     41.66,     42.66,     43.57,     44.38,     45.12,     45.78,     46.38)^2.793 

cas_hg = 4.467e-008 * c(11.76,     15.36,     18.7   ,   21.79   ,  24.67    , 27.34   ,  29.81    , 32.11     ,34.25   ,  36.23   ,  38.06   ,  39.77  ,   41.35  ,   42.82  ,   44.18  ,   45.45    , 46.62   ,  47.71   ,  48.73 ,    49.66 )^2.793 
cas_bp = 4.467e-008 *  c(13.94  ,   17.68  ,   21.14  ,   24.35 ,    27.31  ,   30.06  ,   32.6 ,     34.96 ,    37.14  ,   39.16   ,  41.03 ,    42.77,     44.37,     45.85,     47.23 ,    48.5  ,    49.68  ,   50.77,     51.79,     52.72)^2.793
cbind(cas2$mean_weight_1$`1`$EN.EN$mean_weights$values, cas_en)
cbind(cas2$mean_weight_1$`1`$HG.HG$mean_weights$values, cas_hg)
cbind(cas2$mean_weight_1$`1`$BP.BP$mean_weights$values, cas_bp)



cbind(cas2$mean_weight_1$`1`$EN.EN$mean_weights$values,cas2$mean_weight_1$`1`$EN.HG$mean_weights$values,cas2$mean_weight_1$`1`$EN.BP$mean_weights$values)
cbind(cas2$mean_weight_1$`1`$BP.EN$mean_weights$values,cas2$mean_weight_1$`1`$BP.HG$mean_weights$values,cas2$mean_weight_1$`1`$BP.BP$mean_weights$values)
cbind(cas2$mean_weight_1$`1`$HG.EN$mean_weights$values,cas2$mean_weight_1$`1`$HG.HG$mean_weights$values,cas2$mean_weight_1$`1`$HG.BP$mean_weights$values)

cas2 = extract.mpd("ibm_comparison.out")
cas2_dq = plot.derived_quantities(cas2, report_label = "derived_quants", plot.it = F)

## Read in the IBM output
setwd("../ibm")
library(ibm)


ibm = extract.run("run_100.log")

ibm$Movement_home$`1`$`1950_1-1_destination`
ibm$Movement_home$`1`$`1950_2-1_destination`
ibm$Movement_home$`1`$`1950_3-1_destination`

# find proportions of agents among the three regions to short cut it
rowSums(ibm$init_2$`1`$values[,-1]) /  sum(rowSums(ibm$init_2$`1`$values[,-1]))

ibm_40 = extract.run("output_40.log")
ibm_50 = extract.run("check_100.out")

## look at initial partition
prop_init = ibm$init_2$`1`$values[,-c(1,2)] / rowSums(ibm$init_2$`1`$values[,-c(1,2)])

en_init = colSums(cas2$Init$`1`$values[1:3,-1])
en_prop = en_init / sum(en_init)

hg_init = colSums(cas2$Init$`1`$values[4:6,-1])
hg_prop = hg_init / sum(hg_init)

bp_init = colSums(cas2$Init$`1`$values[7:9,-1])
bp_prop = bp_init / sum(bp_init)


## plot
par(mfrow = c(1,3))
plot(1:20,prop_init[1,], type = "l", lwd = 2, xlab = "Ages", ylab = "Initial proportions",ylim = c(0,0.26))
lines(1:20,en_prop, col = "red", lwd = 2)
#lines(0:20,prop_40[1,], col = "blue", lwd = 2, lty = 2)
#lines(0:20,prop_50[1,], col = "orange", lwd = 2, lty = 2)

plot(1:20,prop_init[2,], type = "l", lwd = 2, xlab = "Ages", ylab = "Initial proportions",ylim = c(0,0.26))
lines(1:20,hg_prop, col = "red", lwd = 2)
#lines(0:20,prop_40[2,], col = "blue", lwd = 2, lty = 2)
#lines(0:20,prop_50[2,], col = "orange", lwd = 2, lty = 2)

plot(1:20,prop_init[3,], type = "l", lwd = 2, xlab = "Ages", ylab = "Initial proportions",ylim = c(0,0.26))
lines(1:20,bp_prop, col = "red", lwd = 2)
#lines(0:20,prop_40[3,], col = "blue", lwd = 2, lty = 2)
#lines(0:20,prop_50[3,], col = "orange", lwd = 2, lty = 2)

names(ibm)
ibm$model_attributes$`1`$global_scalar

ibm$Rec_BP$`1`$b0
ibm$Rec_EN$`1`$b0
ibm$agents$`1`$values$weight
ibm$agents$`1`$values$length
ibm$agents$`1`$values$age
ibm$agents$`1`$values[1:10,]
ibm$derived_quants$`1`$SSB_EN

ibm_dq = plot.derived_quantities(ibm, report_label = "derived_quants", plot.it = F)
years = rownames(ibm_dq)
plot(years,ibm_dq[,c("SSB_HG")], type = "l", ylim = c(120000,140000))
abline(h =ibm$Rec_HG$`1`$b0, col = "red")

plot(years,ibm_dq[,c("SSB_EN")], type = "l", ylim = c(0000,140000))
abline(h =ibm$Rec_EN$`1`$b0, col = "red")

plot(years,ibm_dq[,c("SSB_BP")], type = "l", ylim = c(0,140000))
abline(h =ibm$Rec_BP$`1`$b0, col = "red")

merged = melt(rbind(ibm_dq[,c("SSB_HG",  "SSB_BP",  "SSB_EN")], cas2_dq[,c("SSB_HG",  "SSB_BP",  "SSB_EN")]))
merged$model = rep(c(rep("ibm", nrow(ibm_dq)), rep("Casal2", nrow(cas2_dq))),3)
colnames(merged) = c("year", "region", "SSB", "model")
## melt down the matrices
#jpeg("Threading_with_movement.jpg")
ggplot(merged, aes(x=year, y=SSB, linetype = model, col = region)) + 
  geom_line(size=1) +
  xlab("years")+ylab("Biomass (t)") + 
  ylim(c(0,200000)) +
  scale_colour_discrete(name  ="Region",
                          breaks=c("SSB_HG", "SSB_BP", "SSB_EN"),
                          labels=c("HG", "BP","EN")) +
  #scale_linetype_discrete(name  ="model",
  #                      breaks=c("non-threaded", "threaded"),
  #                      labels=c("no", "yes")) +
  ggtitle("with movement")
#dev.off()



N1 = 10000
N2 = N1 * 0.148
N1 = N1 - N2
N3 = N2 * 0.36
N2 =N2 - N3

(c(N1,N2,N3))

N0 = 10000
N1a = N0 * 0.9220
N2a = N0 * 0.09472
N3a = N0 * 0.002964
c(N1a,N2a,N3a)
##########################################
## create Inputs for teh IBM especially the catch ones which such
##########################################
list.files()
dir.create("CatchLayers")
setwd("CatchLayers")
catch_mat = matrix(0,nrow = 3,ncol = 1)
for (i in 1:length(cas2$instantaneous_mort$`1`$year)) {
  year = cas2$instantaneous_mort$`1`$year[i]
  catch_mat[1,] = cas2$instantaneous_mort$`1`$`actual_catch[EN_LL]`[i]
  catch_mat[2,] = cas2$instantaneous_mort$`1`$`actual_catch[HG_LL]`[i]
  catch_mat[3,] = cas2$instantaneous_mort$`1`$`actual_catch[BP_LL]`[i]
  Filename = make.filename(file = paste0(year,"_catch.ibm"), path = getwd())
  create_ibm_layer(label = paste0(year,"_catch"), type = "numeric", filename = Filename,catch_mat)
}
## now add teh include statements in the config.ibm
setwd("../")

for (i in 1:length(cas2$instantaneous_mort$`1`$year)) {
  year = cas2$instantaneous_mort$`1`$year[i]
  Line = paste0('!include "Layers/CatchLayers/',year,'_catch.ibm"')
  write(Line,file="config.ibm",append=TRUE)
}

##########################################
# A little experiment on how to move individuals around box's (Box transfer)
##########################################
N = 1000000 ## in cell 1
prob = c(0.333333,0.033333,0.6333333)
area = c(1,2,3)
## prob = probabiliity of moving to another cell
area_freq = area_move2 = area_move = vector(length = 3, mode = "numeric");
for (i in 1:N) {
  ## Current implementation
  possible_area = sample(area,1)
  area_freq[possible_area] = area_freq[possible_area] + 1;
  if( runif(1) <= prob[possible_area]) {
    area_move[possible_area] = area_move[possible_area] + 1;
  }
  ## proposed implementation
  ndx = which(as.numeric(rmultinom(1,size = 1, prob)) > 0)
  area_move2[ndx] = area_move2[ndx] + 1;
}
area_move2 / sum(area_move2)
area_move / sum(area_move)
N - sum(area_move2)
N - sum(area_move)

# standard approach for sampling from a multinomial like this is to compare a random standard 
# uniform value to the cumulative sums of the probabilities and return the first index for which 
# the cumulative sum is greater than the random uniform



## compare SSB's
years = as.numeric(rownames(ibm_dq))
plot(years, ibm_dq[,"SSB_HG"], type = "l", lwd = 2, col = "red", xlab = "years", ylab = "SSB (t)", ylim = c(80000,156000))
lines(years, ibm_thread_dq[,"SSB_HG"], lwd = 2, col = "blue")
lines(years, ibm_thread1_dq[,"SSB_HG"], lwd = 2, col = "blue", lty = 2)
lines(years, ibm_thread2_dq[,"SSB_HG"], lwd = 2, col = "green", lty = 2)
lines(years, ibm_thread3_dq[,"SSB_HG"], lwd = 2, col = "yellow", lty = 2)

legend('bottomleft', legend = c("threaded", "not threaded"), col = c("blue", "red"), lwd = 2)



## look at age frequencies
threaded_model = as.numeric(ibm_threaded$init_2$`1`$values[2:32] / sum(ibm_threaded$init_2$`1`$values[2:32]))
nonthreaded_model = as.numeric(ibm$init_2$`1`$values[2:32] / sum(ibm$init_2$`1`$values[2:32]))
## reformat for GGPLOT
new_data = melt(data.frame(casal2_model,ibm_model), variable.name = "model")
new_data$age = c(0:30,0:30)

#jpeg("initial_age_distribution.jpg")
ggplot(new_data,aes(x=age,y=value,fill=model))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_discrete(name="Model",
  labels=c("Casal2", "IBM"))+ 
  xlab("Age")+ylab("Proportion") + 
  ggtitle("Initial age distribution")
#dev.off()

## Look at the fishing age structure
rbind(
ibm_threaded$fishing$`1`$values[5,2:41],
ibm_threaded3$fishing$`1`$values[5,2:41],
ibm$fishing$`1`$values[5,2:41])

ibm_threaded$fishing$`1`$biomass_removed
ibm_threaded3$fishing$`1`$biomass_removed
ibm$fishing$`1`$biomass_removed

## comparison of R0
tot = ibm$Rec_BP$`1`$r0 + ibm$Rec_EN$`1`$r0 + ibm$Rec_HG$`1`$r0
ibm$Rec_BP$`1`$r0 / tot
ibm$Rec_EN$`1`$r0 / tot
ibm$Rec_HG$`1`$r0 / tot


cas2_tot = cas2$Rec_EN$`1`$r0 + cas2$Rec_BP$`1`$r0 + cas2$Rec_HG$`1`$r0
cas2$Rec_EN$`1`$r0 / cas2_tot
cas2$Rec_HG$`1`$r0 / cas2_tot
cas2$Rec_BP$`1`$r0 / cas2_tot


### Estimating growth
growth_df = read.csv("growth.csv", header = T)
## estimate a loose VB for each area based on this data
fit_VB = function(pars, area = "EN") {
  k = pars[1]
  l_inf = pars[2]
  t0 = pars[3]
  
  data = growth_df[growth_df$stock == area, -c(1,2)]
  length_hat = VB(1:20,k,l_inf,t0)
  SSE = sweep(data, MARGIN = 2, length_hat, FUN = "-")
  return(sum(SSE * SSE))
}

opt_EN = nlminb(c(0.2,80,1.2), objective = fit_VB, area = "EN")
opt_EN$convergence
opt_EN$par
## look at fit
EN_fit = VB(1:20, opt_EN$par[1],opt_EN$par[2],opt_EN$par[3])
en_data = growth_df[growth_df$stock == "EN", -c(1,2)]
plot(1:20, en_data[1,], xlab = "age", ylab = "length", ylim = c(0,60), pch = 19, col = "blue")
for(i in 2:nrow(en_data))
  points(1:20, en_data[i,], pch = 19, col = "blue")
lines(1:20, EN_fit, col = "red", lwd = 2)

opt_HG = nlminb(c(0.2,80,1.2), objective = fit_VB, area = "HG")
opt_HG$convergence
opt_HG$par


opt_BP = nlminb(c(0.2,80,1.2), objective = fit_VB, area = "BP")
opt_BP$convergence
opt_BP$par


## Some R code to simulate a baranov process F + M in annual time step
N0 = 1000000
M = 0.2
age_freq = c(1)
for(i in 2:8) {
  age_freq[i] = age_freq[i - 1] * exp(-M)
}
mean_weight = c(0.2,0.4,0.6,0.8,1.2,1.6,2.3,3.0)


N1 = round(N0 * age_freq)

F_mort = 0.3#
sel = d_norm(1:8, 3,1.1,8)
F_a = F_mort * sel
Z = F_a + M
C_baranov = sum((F_a / Z) * (1 - exp(-Z)) * N1)
  
N2 = N1 * exp(-Z)
## now as a individual based model
catch_by_age = vector();
catch_by_age1 = vector();

M_by_age = M_by_age1 = vector();
N2a = vector();
individual = 1;
set.seed(127)
for (a in 1:length(N1)) {
  catch = 0;
  catch1 = 0
  M_this_age = M_this_age1 = 0;
  N_age = 0;
  M_or_F = c(M,F_a[a])
  for (i in 1:(N1[a])) {
    ## randomly check if it survives M or F
    if (runif(1) < (1 - exp(-Z[a]))) {
      if (runif(1) < F_a[a] / (F_a[a] + M)) {
        catch1 = catch1 + 1;
     } else {
        M_this_age1 = M_this_age1 + 1;
     }
    }
    
    this_m_or_f = sample(c(1,2), 1)
    if (runif(1) < (1 - exp(-M_or_F[this_m_or_f]))) {
      if (this_m_or_f == 1) {
        ## dies of M
        M_this_age = M_this_age + 1;
      } else {
        catch = catch + 1;
      }
    } else if (runif(1) < (1 - exp(-(Z[a] - M_or_F[this_m_or_f])))) {
      ## dies of M
      if (this_m_or_f == 2) {
        M_this_age = M_this_age + 1;
      } else {
        catch = catch + 1;
      }
    } else {
      ## survives
      N_age = N_age + 1;
    }
    individual = individual + 1;
  }
  N2a[a] = N_age;
  M_by_age[a] = M_this_age
  catch_by_age[a] = catch
  M_by_age1[a] = M_this_age1
  catch_by_age1[a] = catch1
}

c(sum(catch_by_age), sum(catch_by_age1))
sum(catch_by_age) - C_baranov
cbind(N2a, N2)

cbind(M_by_age, M_by_age1)
