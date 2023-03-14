## Check model output so that the model is behaving as expected
#' @author C.Marsh
#' @date 30/10/2018
#' @description This file reads in some of the initial model runs to check all is working as expected
#' 
#' 
detach("package:ibm", unload=TRUE);

set.seed(123);
library(ibm);
library(ggplot2);
library(RColorBrewer);
library(reshape2);
library(animation);
setwd("..");
BaseDir = getwd();
ibmDir = paste0(BaseDir,"/ibm/")
LayerDir = paste0(BaseDir,"/ibm/Layers/")
RDir = paste0(BaseDir,"/R/")
setwd("R")


output = extract.run(file = "run.log", path= ibmDir)
#output = extract.run(file = "output.log", path= ibmDir)

names(output)
names(output$Move_offshore)
## biomass in the 4 seasons for 1990
output$Autumn_total_biomass_by_cell$`1990`
output$Spring_total_biomass_by_cell$`1990`
output$inshore_total_biomass_by_cell$`1990`
output$Winter_total_biomass_by_cell$`1990`


output$Summer_agents$`1990`$values
## spatial distribution of agents over time
plot(output$Spring_agents$`1990`$values$long,output$Spring_agents$`1990`$values$lat,  xlim = c(18,19),ylim = c(2,3), pch = 19, col ="blue")
points(output$Summer_agents$`1990`$values$long,output$Summer_agents$`1990`$values$lat,   pch = 19, col ="red")
points(output$Winter_agents$`1990`$values$long,output$Winter_agents$`1990`$values$lat,  pch = 19, col ="orange")
points(output$Autumn_agents$`1990`$values$long,output$Autumn_agents$`1990`$values$lat,  pch = 19, col ="green")

## look at SSB over time
plot.derived_quantities(output,report_label = "derived_quants", pch = 16, cex = 0.3, lwd = 2, lty = 2, col = "red")

## Look at Spatial distribution over time
plot_numeric_layer(output,report_label = "Summer_total_biomass_by_cell", directory = paste0(ibmDir,"Figures"), file_name = "Summer_bio", Title = "Summer")
plot_numeric_layer(output,report_label = "Autumn_total_biomass_by_cell", directory = paste0(ibmDir,"Figures"), file_name = "Autumn_bio", Title = "Autumn")
plot_numeric_layer(output,report_label = "Winter_total_biomass_by_cell", directory = paste0(ibmDir,"Figures"), file_name = "Winter_bio", Title = "Winter")
plot_numeric_layer(output,report_label = "Spring_total_biomass_by_cell", directory = paste0(ibmDir,"Figures"), file_name = "Spring_bio", Title = "Spring")

## make an animation
## First order the file names so the animation is in a sensible chronological order.
images = list.files(paste0(ibmDir,"Figures"), pattern = ".png")
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

setwd(paste0(ibmDir,"Figures"))
system(paste("magick convert -delay 120 ",paste(file_names, collapse = " ")," biomass.gif", sep = " "))
setwd(RDir)

## Look at movement rates
sum(output$Move_offshore$`year-area-1990_1-1`$destination_values)
output$Move_offshore$`year-area-1990_1-1`$initial_numbers_in_cell

names(output$offshore_preference_movement)
output$offshore_preference_movement$initialisation_preference
output$offshore_preference_movement$initialisation_meridional
output$offshore_preference_movement$initialisation_zonal

output$offshore_preference_movement$initialisation_preference

output$inshore_preference_movement$time_intervals
output$inshore_preference_movement$standard_dev
names(output$inshore_preference_movement)
output$inshore_preference_movement$initialisation_preference
output$inshore_preference_movement$initialisation_zonal
output$inshore_preference_movement$initialisation_meridional

library(cyrils)
library(fields)
par(mfrow = c(2,2))
image.plot(image_mat(output$inshore_preference_movement$meridonal_1990),main = "Meridonal gradient")
image.plot(image_mat(output$inshore_preference_movement$zonal_1990),main = "Zonal gradient")
#image.plot(image_mat(output$inshore_preference_movement$`1`$initialisation_preference),main = "Preference")

image.plot(image_mat(output$inshore_preference_movement$average_meridional_jump_1990),main = "average merid jump")
image.plot(image_mat(output$inshore_preference_movement$average_zonal_jump_1990),main = "average zonal jump")

par(mfrow = c(2,2))
image.plot(image_mat(output$offshore_preference_movement$meridonal_1990),main = "Meridonal gradient")
image.plot(image_mat(output$offshore_preference_movement$zonal_1990),main = "Zonal gradient")
#image.plot(image_mat(output$inshore_preference_movement$`1`$initialisation_preference),main = "Preference")

image.plot(image_mat(output$offshore_preference_movement$average_meridional_jump_1990),main = "average merid jump")
image.plot(image_mat(output$offshore_preference_movement$average_zonal_jump_1990),main = "average zonal jump")


################
## look at scaled age freq from this model
################
output = extract.run(file = "run.log", path= ibmDir)
output1 = extract.run(file = "run_ageing_error.log", path= ibmDir)
output_e = extract.run(file = "run_equal.log", path= ibmDir)
output_p = extract.run(file = "run_prop.log", path= ibmDir)

length(output$model_attributes$length_mid_points)


len = colnames(output$summer_fishery_scaled_age_freq$`1`$length_freq_by_year_stratum[,-1])
plot(len, apply(output$summer_fishery_scaled_age_freq$`1`$ALK_1991, MARGIN = 2, FUN = sum, na.rm = T)/sum(output$summer_fishery_scaled_age_freq$`1`$ALK_1991) , xlab = "length bin", ylim = c(0,0.11), ylab = "LF", type = "l", lwd = 2)
lines(len, apply(output_e$summer_fishery_scaled_age_freq$`1`$ALK_1991, MARGIN = 2, FUN = sum, na.rm = T)/sum(output_e$summer_fishery_scaled_age_freq$`1`$ALK_1991), lwd = 2, col = "red")
lines(len, apply(output_p$summer_fishery_scaled_age_freq$`1`$ALK_1991, MARGIN = 2, FUN = sum, na.rm = T)/sum(output_p$summer_fishery_scaled_age_freq$`1`$ALK_1991), lwd = 2, col = "blue", lty = 3)
lines(len, output$summer_fishery_scaled_age_freq$`1`$length_freq_by_year_stratum[,-1]/sum(output$summer_fishery_scaled_age_freq$`1`$length_freq_by_year_stratum[,-1]), lwd = 2, col = "green", lty = 2)

## plot comparison
ages = output$summer_fishery_scaled_age_freq$Values$age
plot(ages, output$summer_fishery_scaled_age_freq$Values$expected / sum(output$summer_fishery_scaled_age_freq$Values$expected), type = "l", col = "red", lwd = 3, xlab  ="Age", ylab = "Expected Frequency", ylim= c(0,0.15))
lines(ages, output$overall_fishery_age$Values$expected / sum(output$overall_fishery_age$Values$expected), col = "blue", lwd = 2);
lines(ages, output_e$summer_fishery_scaled_age_freq$Values$expected / sum(output_e$summer_fishery_scaled_age_freq$Values$expected), col = "orange", lwd = 2);
lines(ages, output_p$summer_fishery_scaled_age_freq$Values$expected / sum(output_p$summer_fishery_scaled_age_freq$Values$expected), col = "green", lwd = 2);
lines(ages, output1$summer_fishery_scaled_age_freq$Values$expected / sum(output1$summer_fishery_scaled_age_freq$Values$expected), col = "gray60", lwd = 2);

alk_age = output$summer_fishery_scaled_age_freq$Values$expected
dir_age = output$overall_fishery_age$Values$expected 

cbind(alk_age, dir_age)

plot(output$Summer_agents$`1990`$values$age, output$Summer_agents$`1990`$values$length)












