#'
#'
#' Summarise
#'

library(ibm)
library(ggplot2)
library(reshape2)
library(fields)
library(RColorBrewer)

flip_rows = function(x) {
  return(x[nrow(x):1,])
}
r_dir = getwd()
ibm = extract.run(file = file.path("..","ibm", "run.out"))
## these two are for validation
ibm_old = extract.run(file = file.path("..","ibm", "run_old_method.out"))
ibm_new = extract.run(file = file.path("..","ibm", "run_new_method.out"))
ibm_old$Fishing$actual_catch
ibm_new$Fishing$actual_catch
ibm_old$Fishing$seconds_to_minimise
ibm_new$Fishing$seconds_to_minimise




n_row = 20
n_col = 50
n_cells = n_row * n_col

spatial_grid = expand.grid(1:n_row, 1:n_col)

# lat_mid,long_mid, lat_bins,long_bins,  bath_mat1
xlim= c(174, 186.5)
ylim= c(-42.5, -45)

x_mid = seq(from = xlim[1], to = xlim[2], length = n_col)
x_res = unique(round(diff(x_mid) / 2,2))
y_mid = seq(from =  ylim[1], to = ylim[2], length = n_row)
y_res = unique(round(diff(y_mid) / 2,2))

x_upper = c(min(x_mid) - x_res, x_mid + x_res)
y_upper = c(max(y_mid) - y_res, y_mid + y_res)


## compare old with new

image.plot(x_mid, rev(y_mid), t(flip_rows(ibm_old$Fishing$effort_by_cell_1972)), main = "effort old", xlab = "", ylab = "")
image.plot(x_mid, rev(y_mid), t(flip_rows(ibm_new$Fishing$effort_by_cell_1972)), main = "effort new", xlab = "", ylab = "")
image.plot(x_mid, rev(y_mid), t(flip_rows(ibm_new$Fishing$preference_by_cell_1972)), main = "pref new", xlab = "", ylab = "")
image.plot(x_mid, rev(y_mid), t(flip_rows(ibm_old$Fishing$vulnerable_by_cell_1972)), main = "pref new", xlab = "", ylab = "")

plot(years, ibm_new$Fishing$catches, type = "l", lwd = 3)
lines(years, ibm_new$Fishing$actual_catch, lwd = 3, lty = 2, col = "red")



dir.create(file.path("..","ibm", "Animation"))
dir.create(file.path("..","ibm", "Animation","shallow"))
dir.create(file.path("..","ibm", "Animation","deep"))
dir.create(file.path("..","ibm", "Animation","vulnerable"))
dir.create(file.path("..","ibm", "Animation","catch"))
dir.create(file.path("..","ibm", "Animation","preference"))
dir.create(file.path("..","ibm", "Animation","effort"))

years = names(ibm$derived_quants$total_SSB$values)
ages = ibm$model_attributes$ages


## plot spatial maps
vals = vector()
years_to_sim = 1972:2016
for(j in 1:length(years_to_sim)) {
  vals = c(vals, as.vector(ibm:::evalit(paste0("ibm$'shallow_bio'$'",years_to_sim[j],"'$values"))))
}
old_vals = vector()
years_to_sim = 1972:2016
for(j in 1:length(years_to_sim)) {
  old_vals = c(old_vals, as.vector(ibm:::evalit(paste0("ibm$'deep_bio'$'",years_to_sim[j],"'$values"))))
}

breaks = c(min(vals) - 1, quantile(vals, c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)), max(vals) + 1)
breaks_old = c(min(old_vals) - 1, quantile(old_vals, c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)), max(old_vals) + 1)

plot_numeric_layer(ibm, report_label = "shallow_bio", directory = file.path("..","ibm", "Animation","shallow"), file_name = "shallow_biomass", Title = "Shallow Biomass", breaks = breaks, width = 10, height = 7)
setwd(r_dir)
plot_numeric_layer(ibm, report_label = "deep_bio", directory = file.path("..","ibm", "Animation","deep"), file_name = "deep_biomass", Title = "Deep Biomass", breaks = breaks_old, width = 10, height = 7)
setwd(r_dir)
plot_numeric_layer(ibm, report_label = "Depth", directory = file.path("..","ibm"), file_name = "Depth", Title = "Depth", breaks = NULL, width = 10, height = 7)
setwd(r_dir)
plot_spatial_catch(ibm, report_label = "Fishing", directory = file.path("..","ibm", "Animation","catch"), file_name = "catch", Title = "catch", width = 10, height = 7, plot_catch = T)
setwd(r_dir)
plot_spatial_catch(ibm, report_label = "Fishing", directory = file.path("..","ibm", "Animation","vulnerable"), file_name = "vulnerable", Title = "vulnerable", width = 10, height = 7, plot_catch = F)
setwd(r_dir)
plot_spatial_catch(ibm, report_label = "Fishing", directory = file.path("..","ibm", "Animation","vulnerable"), file_name = "vulnerable", Title = "vulnerable", width = 10, height = 7, plot_catch = F)
setwd(r_dir)


plot_mulitple_mats(model= ibm, report_label = "Fishing",pattern = "F_by_cell_", directory = file.path("..","ibm", "Animation","effort"), file_name = "effort", Title = "effort", width = 10, height = 7)
setwd(r_dir)
plot_mulitple_mats(model= ibm, report_label = "Fishing",pattern = "preference_by_cell_", directory = file.path("..","ibm", "Animation","preference"), file_name = "preference", Title = "preference", width = 10, height = 7)
setwd(r_dir)

## First order the file names so the animation is in a sensible chronological order.
images = list.files(file.path("..","ibm", "Animation","shallow"), pattern = ".png")
year = sort(unique(as.numeric(substr(images, 0,4))))
file_names = c();
for (year_ndx in 1:length(years_to_sim)) {
  this_year = images[grepl(years_to_sim[year_ndx],x = images)]
  file_names = c(file_names, this_year)
}

images = list.files(file.path("..","ibm", "Animation","deep"), pattern = ".png")
year = sort(unique(as.numeric(substr(images, 0,4))))
file_names_old = c();
for (year_ndx in 1:length(years_to_sim)) {
  this_year = images[grepl(years_to_sim[year_ndx],x = images)]
  file_names_old = c(file_names_old, this_year)
}
images = list.files(file.path("..","ibm", "Animation","catch"), pattern = ".png")
year = sort(unique(as.numeric(substr(images, 0,4))))
file_names_catch = c();
for (year_ndx in 1:length(years_to_sim)) {
  this_year = images[grepl(years_to_sim[year_ndx],x = images)]
  file_names_catch = c(file_names_catch, this_year)
}
images = list.files(file.path("..","ibm", "Animation","vulnerable"), pattern = ".png")
year = sort(unique(as.numeric(substr(images, 0,4))))
file_names_vuln = c();
for (year_ndx in 1:length(years_to_sim)) {
  this_year = images[grepl(years_to_sim[year_ndx],x = images)]
  file_names_vuln = c(file_names_vuln, this_year)
}

setwd(file.path("..","ibm", "Animation","shallow"))
system(paste("magick convert -delay 120 ",paste(file_names, collapse = " ")," biomass.gif", sep = " "))
setwd(r_dir)
setwd(file.path("..","ibm", "Animation","deep"))
system(paste("magick convert -delay 120 ",paste(file_names_old, collapse = " ")," biomass.gif", sep = " "))
setwd(r_dir)
setwd(file.path("..","ibm", "Animation","catch"))
system(paste("magick convert -delay 120 ",paste(file_names_catch, collapse = " ")," biomass.gif", sep = " "))
setwd(r_dir)
setwd(file.path("..","ibm", "Animation","vulnerable"))
system(paste("magick convert -delay 120 ",paste(file_names_vuln, collapse = " ")," biomass.gif", sep = " "))
setwd(r_dir)

preference = ibm$Fishing$preference_by_cell_1990
vulnerable = ibm$Fishing$vulnerable_by_cell_1990
effort = ibm$Fishing$effort_by_cell_1990

#vulnerable = ibm$Fishing$effort_by_cell_1982
base_pref = ibm$Fishing$base_preference
par(mfrow = c(3,1))
image.plot(x_mid, rev(y_mid), t(flip_rows(preference)), main = "Preference", xlab = "", ylab = "")
image.plot(x_mid, rev(y_mid), t(flip_rows(vulnerable)), main = "vulnerable", xlab = "", ylab = "")
image.plot(x_mid, rev(y_mid), t(flip_rows(base_pref)), main = "base_pref", xlab = "", ylab = "")
image.plot(x_mid, rev(y_mid), t(flip_rows(effort)), main = "effort", xlab = "", ylab = "")

vulnerable_pref = vulnerable / max(vulnerable)
image.plot(x_mid, rev(y_mid), t(flip_rows(vulnerable_pref)), main = "vulnerable_pref", xlab = "", ylab = "")
pref  = (base_pref + vulnerable_pref) / 2
image.plot(x_mid, rev(y_mid), t(flip_rows(pref)), main = "pref", xlab = "", ylab = "")
