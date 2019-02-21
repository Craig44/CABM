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


output_50 = extract.run(file = "output_50_burnin.log", path= ibmDir)
output = extract.run(file = "run.log", path= ibmDir)
#output = extract.run(file = "output.log", path= ibmDir)

names(output)
output$move_offshore$`1`$values
output$Autumn_total_biomass_by_cell$`1`$values
output$Spring_total_biomass_by_cell$`1`$values
output$inshore_total_biomass_by_cell$`1`$values
output$Winter_total_biomass_by_cell$`1`$values


output$Summer_agents$`1`$values

plot(output$Spring_agents$`1`$values$long,output$Spring_agents$`1`$values$lat,  xlim = c(18,19),ylim = c(2,3), pch = 19, col ="blue")
points(output$Summer_agents$`1`$values$long,output$Summer_agents$`1`$values$lat,  xlim = c(18,19),ylim = c(2,3), pch = 19, col ="red")

## look at SSB over time
plot.derived_quantities(output,report_label = "derived_quants")

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

image.plot(image_mat(report[[i]]$values), breaks = breaks, col = colorRampPalette(brewer.pal(9,"YlOrRd"))(length(breaks) - 1))

## Look at movement rates
sum(output$Move_offshore$`1`$`year-area-1990_1-1`$destination_values)
output$Move_offshore$`1`$`year-area-1990_1-1`$initial_numbers_in_cell

names(output$offshore_preference_movement$`1`)
output$offshore_preference_movement$`1`$initialisation_preference
output$offshore_preference_movement$`1`$initialisation_meridional
output$offshore_preference_movement$`1`$initialisation_zonal

output$offshore_preference_movement$`1`$initialisation_preference

output$inshore_preference_movement$`1`$time_intervals
output$inshore_preference_movement$`1`$standard_dev
names(output$inshore_preference_movement$`1`)
output$inshore_preference_movement$`1`$initialisation_preference
output$inshore_preference_movement$`1`$initialisation_zonal
output$inshore_preference_movement$`1`$initialisation_meridional

library(cyrils)
library(fields)
par(mfrow = c(2,2))
image.plot(image_mat(output$inshore_preference_movement$`1`$meridonal_1990),main = "Meridonal gradient")
image.plot(image_mat(output$inshore_preference_movement$`1`$zonal_1990),main = "Zonal gradient")
#image.plot(image_mat(output$inshore_preference_movement$`1`$initialisation_preference),main = "Preference")

image.plot(image_mat(output$inshore_preference_movement$`1`$average_meridional_jump_1990),main = "average merid jump")
image.plot(image_mat(output$inshore_preference_movement$`1`$average_zonal_jump_1990),main = "average zonal jump")

par(mfrow = c(2,2))
image.plot(image_mat(output$offshore_preference_movement$`1`$meridonal_1990),main = "Meridonal gradient")
image.plot(image_mat(output$offshore_preference_movement$`1`$zonal_1990),main = "Zonal gradient")
#image.plot(image_mat(output$inshore_preference_movement$`1`$initialisation_preference),main = "Preference")

image.plot(image_mat(output$offshore_preference_movement$`1`$average_meridional_jump_1990),main = "average merid jump")
image.plot(image_mat(output$offshore_preference_movement$`1`$average_zonal_jump_1990),main = "average zonal jump")


################
## look at scaled age freq
################
output = extract.run(file = "run.log", path= ibmDir)
output1 = extract.run(file = "run_ageing_error.log", path= ibmDir)
output_e = extract.run(file = "run_equal.log", path= ibmDir)
output_p = extract.run(file = "run_prop.log", path= ibmDir)

length(output$model_attributes$`1`$length_mid_points)

sum(apply(output$summer_fishery_scaled_age_freq$`1`$ALK_1991, MARGIN = 1, FUN = sum, na.rm = T))
sum(apply(output1$summer_fishery_scaled_age_freq$`1`$ALK_1991, MARGIN = 1, FUN = sum, na.rm = T))

sum(output$summer_fishery_scaled_age_freq$`1`$length_freq_by_year_stratum[,-1])
sum(output1$summer_fishery_scaled_age_freq$`1`$length_freq_by_year_stratum[,-1])

output1$summer_fishery_scaled_age_freq$`1`$length_freq_by_year_stratum[,-1]/sum(output1$summer_fishery_scaled_age_freq$`1`$length_freq_by_year_stratum[,-1]) - 
t(apply(output1$summer_fishery_scaled_age_freq$`1`$ALK_1991, MARGIN = 2, FUN = sum, na.rm = T)) / sum(apply(output1$summer_fishery_scaled_age_freq$`1`$ALK_1991, MARGIN = 2, FUN = sum, na.rm = T))


dim(output$summer_fishery_scaled_age_freq$`1`$ALK_1991)
apply(output$summer_fishery_scaled_age_freq$`1`$ALK_1991, MARGIN = 1, FUN = sum, na.rm = T)
apply(output1$summer_fishery_scaled_age_freq$`1`$ALK_1991, MARGIN =1, FUN = sum, na.rm = T)

sum(output$summer_fishery_scaled_age_freq$`1`$ALK_1991)
sum(output_e$summer_fishery_scaled_age_freq$`1`$ALK_1991)
sum(output_p$summer_fishery_scaled_age_freq$`1`$ALK_1991)

len = colnames(output$summer_fishery_scaled_age_freq$`1`$length_freq_by_year_stratum[,-1])
plot(len, apply(output$summer_fishery_scaled_age_freq$`1`$ALK_1991, MARGIN = 2, FUN = sum, na.rm = T)/sum(output$summer_fishery_scaled_age_freq$`1`$ALK_1991) , xlab = "length bin", ylim = c(0,0.11), ylab = "LF", type = "l", lwd = 2)
lines(len, apply(output_e$summer_fishery_scaled_age_freq$`1`$ALK_1991, MARGIN = 2, FUN = sum, na.rm = T)/sum(output_e$summer_fishery_scaled_age_freq$`1`$ALK_1991), lwd = 2, col = "red")
lines(len, apply(output_p$summer_fishery_scaled_age_freq$`1`$ALK_1991, MARGIN = 2, FUN = sum, na.rm = T)/sum(output_p$summer_fishery_scaled_age_freq$`1`$ALK_1991), lwd = 2, col = "blue", lty = 3)
lines(len, output$summer_fishery_scaled_age_freq$`1`$length_freq_by_year_stratum[,-1]/sum(output$summer_fishery_scaled_age_freq$`1`$length_freq_by_year_stratum[,-1]), lwd = 2, col = "green", lty = 2)

## plot comparison
ages = output$summer_fishery_scaled_age_freq$`1`$Values$age
plot(ages, output$summer_fishery_scaled_age_freq$`1`$Values$expected / sum(output$summer_fishery_scaled_age_freq$`1`$Values$expected), type = "l", col = "red", lwd = 3, xlab  ="Age", ylab = "Expected Frequency", ylim= c(0,0.15))
lines(ages, output$overall_fishery_age$`1`$Values$expected / sum(output$overall_fishery_age$`1`$Values$expected), col = "blue", lwd = 2);
lines(ages, output_e$summer_fishery_scaled_age_freq$`1`$Values$expected / sum(output_e$summer_fishery_scaled_age_freq$`1`$Values$expected), col = "orange", lwd = 2);
lines(ages, output_p$summer_fishery_scaled_age_freq$`1`$Values$expected / sum(output_p$summer_fishery_scaled_age_freq$`1`$Values$expected), col = "green", lwd = 2);
lines(ages, output1$summer_fishery_scaled_age_freq$`1`$Values$expected / sum(output1$summer_fishery_scaled_age_freq$`1`$Values$expected), col = "gray60", lwd = 2);

alk_age = output$summer_fishery_scaled_age_freq$`1`$Values$expected
dir_age = output$overall_fishery_age$`1`$Values$expected 

cbind(alk_age, dir_age)

plot(output$Summer_agents$`1`$values$age, output$Summer_agents$`1`$values$length)

## look at the ALK
output$summer_fishery_scaled_age_freq$`1`$ALK_1991
output$model_attributes$`1`$Recruitment













