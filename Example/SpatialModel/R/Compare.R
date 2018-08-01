## this script is just comparing, steady state models just constant recruitment and M
library(ggplot2)
library(gridExtra)
library(reshape2)
library(RColorBrewer) 
library(fields) ## Used for heatmap legend
detach("package:ibm", unload=TRUE)
library(ibm)

## Read in the IBM output
setwd("../ibm")
initialised = TRUE;
## the first thing we want to do with complex models is check the burn-in period
if (!initialised) {
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
  
  ## plot age frequency for each row of spatial gird
  for (i in 1:5) {
    temp_data = merged[substring(merged$`cell`,0,1) == as.character(i),]
    p <- ggplot(temp_data, aes(x = age, y = frequency, color = burnin)) + geom_point()
    p + facet_grid(cell ~ ., scales="free_y")
    print(p)
    Sys.sleep(2);
  }
}

##########################
##  Debugging the model ##
##########################
ibm_box = extract.run("output_box_transfer_w_fishing.log")
ibm_box_thread = extract.run("output_box_transfer_w_fishing_thread.log")
ibm_box1 = extract.run("output_box_transfer_w_fishing1.log")
ibm_box_thread1 = extract.run("output_box_transfer_w_fishing_thread1.log")

age_spread = length(0:30)

## check the seed is working
names(ibm_box)

p = plot_overall_age_freq(ibm_box, world_age_frequency_label = "world_age_freq", fishery_age_frequency_label = "fishing", year = 1990)
p1 = plot_overall_age_freq(ibm_box, world_age_frequency_label = "world_age_freq", fishery_age_frequency_label = "fishing", year = 2000)

grid.arrange(p$plot, p1$plot, ncol = 2)

## create a biomass grid over time of the cells


########################
##  look at movement  ##
########################
ibm_pref = extract.run("movment_preference_with_fishing.out")
ibm_box = extract.run("output_box_transfer_w_fishing.log")

names(ibm_pref)
names(ibm_box)

## create heat maps of changing abundance over time with a legend


