## Read in Casal2 output
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
ibm2 = extract.run("run.log")

names(ibm)

ibm_dq = plot.derived_quantities(ibm, report_label = "derived_quants", plot.it = F)
ibm_dq2 = plot.derived_quantities(ibm2, report_label = "derived_quants", plot.it = F)
ibm_dq3 = plot.derived_quantities(ibm3, report_label = "derived_quants", plot.it = F)

## compare SSB's
years = as.numeric(rownames(ibm_dq))
plot(years, ibm_dq[,"SSB"], type = "l", lwd = 2, col = "red", xlab = "years", ylab = "SSB (t)", ylim = c(0,36000))
lines(years, cas2_dq[,"SSB"], lwd = 2, col = "blue")
lines(years, ibm_dq2[,"SSB"], lwd = 2, col = "black")
#lines(years, ibm_dq3[,"SSB"], lwd = 2, lty = 2, col = "orange")
legend('bottomleft', legend = c("Casal2", "IBM-no variability", "IBM with variability"), col = c("red","blue","black"), lty = c(1,1,1),lwd = 2)

## lets look at fishing
ibm_comp = ibm$fishing$`1`$age_frequency
ibm_comp2 = ibm2$fishing$`1`$age_frequency

## look at age and length relationship
df = ibm3$agents$`3`$values
df1 = ibm$agents$`3`$values

plot(df$age, df$length,pch = 19, xlab = "age", ylab = "length", main = "age length of 1000 agents")
points(df1$age, df1$length,pch = 19, col = "red")

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


## plot preference functions
png(filename = "pref_plot.png", units = "in", width = 6, height = 4, res = 250)
par(mfrow = c(1,1), cex.main = 2, cex.lab = 2, mar = c(5,6,1,1))
plot(names(ibm2$normal$Values), ibm2$normal$Values, type = "l", lwd = 3, xlab = expression(x[r]), ylab = expression(P[r]*x), ylim = c(0,1))
lines(names(ibm2$double_normal$Values), ibm2$double_normal$Values, lwd = 3, lty = 2, col = "red")
lines(names(ibm2$logistic$Values), ibm2$logistic$Values, lwd = 3, lty = 2, col = "blue")
lines(names(ibm2$inverse_logistic$Values), ibm2$inverse_logistic$Values, lwd = 3, lty = 1, col = "purple")
legend('bottomright', legend = c("Normal", "Double Normal","Logistic","Inverse Logistic"), lty = c(1,2,2,1), col = c("black","red","blue","purple"), lwd = 3)
dev.off()
