## Read in Casal2 output
detach("package:Casal2", unload=TRUE)
detach("package:ibm", unload=TRUE)

## this script is just comparing, steady state models just constant recruitment and M
library(ggplot2)
library(reshape2)
library(casal2)
setwd("../Casal2")

cas2 = extract.mpd("output.log")
cas2_dq = plot.derived_quantities(cas2, report_label = "derived_quants", plot.it = F)

## Read in the IBM output
setwd("../ibm")
library(ibm)

ibm = extract.run("output.log")
ibm2 = extract.run("with_variation.out")

names(ibm)

ibm_dq = plot.derived_quantities(ibm, report_label = "derived_quants", plot.it = F)
ibm_dq2 = plot.derived_quantities(ibm2, report_label = "derived_quants", plot.it = F)

## compare SSB's
years = as.numeric(rownames(ibm_dq))
plot(years, ibm_dq[,"SSB"], type = "l", lwd = 2, col = "red", xlab = "years", ylab = "SSB (t)", ylim = c(0,36000))
lines(years, cas2_dq[,"SSB"], lwd = 2, col = "blue")
lines(years, ibm_dq2[,"SSB"], lwd = 2, col = "black")
legend('bottomleft', legend = c("casal2", "ibm"), col = c("blue", "red"), lwd = 2)

## lets look at fishing
ibm_comp = ibm$fishing$`1`$age_frequency
ibm_comp2 = ibm2$fishing$`1`$age_frequency



## look at age frequencies
casal2_model = as.numeric(cas2$Init$`1`$values[2:32] / sum(cas2$Init$`1`$values[2:32]))
ibm_model = as.numeric(ibm$init_2$`1`$values[2:32] / sum(ibm$init_2$`1`$values[2:32]))
ibm_model2 = as.numeric(ibm2$init_2$`1`$values[2:32] / sum(ibm2$init_2$`1`$values[2:32]))

## reformat for GGPLOT
new_data = melt(data.frame(casal2_model,ibm_model,ibm_model2), variable.name = "model")
new_data$age = c(0:30,0:30,0:30)

jpeg("initial_age_distribution.jpg")
ggplot(new_data,aes(x=age,y=value,fill=model))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_discrete(name="Model",
  labels=c("Casal2", "IBM","IBM latest"))+ 
  xlab("Age")+ylab("Proportion") + 
  ggtitle("Initial age distribution")
dev.off()


casal2_model = as.numeric(cas2$Annual_part$'53'$values[2:32] / sum(cas2$Annual_part$'53'$values[2:31]))
ibm_model = as.numeric(ibm$total_age_freq$`53`$values / sum(ibm$total_age_freq$`53`$values))
ibm_model2 = as.numeric(ibm2$total_age_freq$`53`$values / sum(ibm2$total_age_freq$`53`$values))

## reformat for GGPLOT
new_data = melt(data.frame(casal2_model,ibm_model,ibm_model2), variable.name = "model")
new_data$age = c(0:30,0:30,0:30)
jpeg("final_age_distribution.jpg")
ggplot(new_data,aes(x=age,y=value,fill=model))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_discrete(name="Model",
  labels=c("Casal2", "IBM","IBM latest"))+ 
  xlab("Age")+ylab("Proportion") + 
  ggtitle("Final age distribution")
dev.off()
## look at the fisher age distribution


