## Test spatial F
#' @author C.Marsh
#' @date 1/11/2018
#' @description This file helps debug a spatial F that is based on perfect knowledge of the spation distribution of the
#' population (ideal free distribution)
#' 
#' 
detach("package:ibm", unload=TRUE)

set.seed(123)
library(ibm)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(animation)
library(fields)
library(cyrils) # used for image_mat() and other handy functions
setwd("..")
BaseDir = getwd();
ibmDir = paste0(BaseDir,"/ibm/")
LayerDir = paste0(BaseDir,"/ibm/Layers/")
RDir = paste0(BaseDir,"/R/")
setwd("R")

## Draw some effort values, catchs values and YCS values for the model
#effort = rlnorm(100,-2.3,0.5)
#write.table(round(t(effort),3), file = paste0(ibmDir,"Effort.txt"),quote = FALSE, row.names = F, col.names = F)
#catch_history = c(rnorm(9,seq(1000,10000,by = 500), 300), rnorm(11, 8000, 300), rnorm(8, 2000, 300))
#plot(1991:2018,catch_history,type = "l")
#write.table(round(t(catch_history),1), file = paste0(ibmDir,"catch_history.txt"),quote = FALSE, row.names = F, col.names = F)
#ycs = r_lnorm(length(1990:2018),1,0.6)
#write.table(round(t(ycs),2), file = paste0(ibmDir,"ycs.txt"),quote = FALSE, row.names = F, col.names = F)

output = extract.run(file = "run.log", path= ibmDir)
#output_70 = extract.run(file = "output_70.log", path= ibmDir) ## 70 years initialisation
#output = output_70
years = 1990:2018
scalar = output$model_attributes$Recruitment
## An ugly plot of SSB
plot.derived_quantities(output,report_label = "derived_quants", lwd = 2, pch = 16, cex = 0.2, lty = 2)
## Recruitment
plot(years, output$Rec$recruits * scalar /1000, type = "l",lwd = 2, ylim = c(0,70000), ylab = "Recruits (000's)")
abline(h = output$Rec$initial_recruits  * scalar /1000, col = "blue", lty = 2, lwd = 2)
legend('topright', legend = c("R0"), col = c("blue"), lty = c(2), lwd = 2)

agent_info = output$Summer_agents$`1990`$values
plot(agent_info$age, agent_info$length)

output$summer_fishery$`1991_message`
## check we removed what we thought we wanted to
output$summer_fishery$catches - output$summer_fishery$actual_catch
output$summer_fishery$`2018_message`
output$summer_fishery$`2017_message`
output$summer_fishery$`2000_message`

## A crude proxy for F over time for this fishery
output$summer_fishery$lambda * mean(output$summer_fishery$effort_values)

## Look at Spatial distribution over time
dir.create(paste0(ibmDir,"Figures/Exploitation"), showWarnings = FALSE)
plot_numeric_layer(output,report_label = "Summer_total_biomass_by_cell", directory = paste0(ibmDir,"Figures/Exploitation"), file_name = "Summer_bio", Title = "Summer")
plot_numeric_layer(output,report_label = "Autumn_total_biomass_by_cell", directory = paste0(ibmDir,"Figures/Exploitation"), file_name = "Autumn_bio", Title = "Autumn")
plot_numeric_layer(output,report_label = "Winter_total_biomass_by_cell", directory = paste0(ibmDir,"Figures/Exploitation"), file_name = "Winter_bio", Title = "Winter")
plot_numeric_layer(output,report_label = "Spring_total_biomass_by_cell", directory = paste0(ibmDir,"Figures/Exploitation"), file_name = "Spring_bio", Title = "Spring")


## make an animation
## First order the file names so the animation is in a sensible chronological order.
images = list.files(paste0(ibmDir,"Figures/Exploitation"), pattern = ".png")
year = sort(unique(as.numeric(substr(images, 0,4))))
quarter_time_step = c("Sum","Aut","Win","Spr")
file_names = c();
for (year_ndx in 1:length(year)) {
  this_year = images[grepl(year[year_ndx],x = images)]
  for (quarter_ndx in 1:length(quarter_time_step)) {
    this_quarter = this_year[grepl(quarter_time_step[quarter_ndx],x = this_year)]
    file_names = c(file_names, this_quarter)
  }
}

setwd(paste0(ibmDir,"Figures/Exploitation"))
system(paste("magick convert -delay 120 ",paste(file_names, collapse = " ")," biomass.gif", sep = " "))
setwd(RDir)


output$summer_fishery$vulnerable_by_cell_1991
output$summer_fishery$removals_1991

## compare age based observations
output$summer_fishery_scaled_age_freq
output$overall_fishery_age

sum(output$summer_fishery_scaled_age_freq$Values$expected)
sum(output$overall_fishery_age$Values$expected)
sum(output$overall_fishery_length$Values$expected)


ages = output$summer_fishery_scaled_age_freq$Values$age
plot(ages, output$summer_fishery_scaled_age_freq$Values$expected / sum(output$summer_fishery_scaled_age_freq$Values$expected), type = "l", col = "red", lwd = 2, xlab  ="Age", ylab = "Expected Frequency")
lines(ages, output$overall_fishery_age$Values$expected / sum(output$overall_fishery_age$Values$expected), col = "blue", lwd = 2);


## fishery length distribution
plot(output$overall_fishery_length$Values$length,output$overall_fishery_length$Values$expected, type = "l", lwd = 2, xlab = "length", ylab = "Numbers")





