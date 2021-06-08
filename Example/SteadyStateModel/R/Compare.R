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

names(ibm)

ibm_dq = plot.derived_quantities(ibm, report_label = "derived_quants", plot.it = F)

## compare SSB's
years = as.numeric(rownames(ibm_dq))
plot(years, ibm_dq[,"SSB"], type = "l", lwd = 3, col = "red", xlab = "years", ylab = "SSB (t)", ylim = c(15000,56000))
lines(years, cas2_dq[,"SSB"], lwd = 3, col = "blue", lty = 2)
legend('bottomleft', legend = c( "IBM", "Casal2"), col = c("red","blue"), lty = c(1,1,2),lwd = 2)
## look at age frequencies
casal2_model = as.numeric(cas2$Init$values[2:31])
ibm_model = as.numeric(ibm$init_2$values )

## gROWTH


## reformat for GGPLOT
new_data = melt(data.frame(casal2_model,ibm_model), variable.name = "model")
new_data$age = c(1:30,1:30)

jpeg("initial_age_distribution.jpg")
ggplot(new_data,aes(x=age,y=value,fill=model))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_discrete(name="Model",
  labels=c("Casal2", "IBM","IBM2"))+ 
  xlab("Age")+ylab("Proportion") + 
  ggtitle("Initial age distribution")
dev.off()


casal2_model = as.numeric(cas2$Annual_part$'53'$values[2:32] / sum(cas2$Annual_part$'53'$values[2:31]))
ibm_model = as.numeric(ibm$total_age_freq$`53`$values / sum(ibm$total_age_freq$`53`$values))
## reformat for GGPLOT
new_data = melt(data.frame(casal2_model,ibm_model), variable.name = "model")
new_data$age = c(0:30,0:30)
jpeg("final_age_distribution.jpg")
ggplot(new_data,aes(x=age,y=value,fill=model))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_discrete(name="Model",
  labels=c("Casal2", "IBM"))+ 
  xlab("Age")+ylab("Proportion") + 
  ggtitle("Final age distribution")
dev.off()