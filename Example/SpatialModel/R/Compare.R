## this script is just comparing, steady state models just constant recruitment and M
library(ggplot2)
library(reshape2)

## Read in the IBM output
setwd("../ibm")
library(ibm)

## the first thing we want to do with complex models is check the burn-in period

ibm_30 = extract.run("output_30.log")
ibm_50 = extract.run("output_50.log")
ibm_80 = extract.run("output_80.log")
ibm_120 = extract.run("output_120.log")

names(ibm_50)
## compare the initial partition between, 30, 50, 80 and 120 years of burnin
mat = ibm_50$init_2$`1`$values

temp = rbind(ibm_30$init_2$`1`$values,ibm_50$init_2$`1`$values, ibm_80$init_2$`1`$values, ibm_120$init_2$`1`$values)
temp$burnin =  c(rep("30", nrow(mat)),rep("50", nrow(mat)), rep("80", nrow(mat)), rep("120", nrow(mat)))

merged = melt(temp)
colnames(merged) = c("cell", "burnin", "age", "frequency")

S
# Use vars() to supply variables from the dataset:
p + facet_grid(rows = vars(drv))


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

# A little experiment on how to move individuals around box's (Box transfer)
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
