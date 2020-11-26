#'
#' Demonstrate how preference movement works.
#'
library(ibm)
library(fields)
library(RColorBrewer)
depth_file = extract.ibm.file("Depth.ibm")
#run = extract.run("run.out")
run = extract.run("run_alt.out")

y_lim = -1.0 * c(42.44, 42.56, 42.69, 42.81, 42.94, 43.07, 43.19, 43.32, 43.44, 43.57, 43.70, 43.82, 43.95, 44.08, 44.20, 44.33, 44.45, 44.58, 44.71, 44.83, 44.96)
x_lim = c(173.37, 173.63, 173.90, 174.17, 174.43, 174.70, 174.97, 175.24, 175.51, 175.78, 176.04, 176.31, 176.58, 176.85, 177.12, 177.38, 177.65, 177.92, 178.19, 178.46, 178.73, 178.99, 179.26, 179.53, 179.80, 180.07, 180.33, 180.60, 180.87, 181.14, 181.41, 181.68, 181.94, 182.21, 182.48, 182.75, 183.02, 183.28, 183.55, 183.82, 184.09, 184.36, 184.63, 184.89, 185.16, 185.43)
-1 * run$model_attributes$latitude_mid_points
y_mid = diff(y_lim) + y_lim[1:(length(y_lim) - 1)]
x_mid = diff(x_lim) + x_lim[1:(length(x_lim) - 1)]

#x_mid = run$model_attributes$longitude_mid_points
#y_mid = run$model_attributes$latitude_mid_points
n_col = length(x_mid)
n_row = length(y_mid)
dim(depth_file$`layer[Depth]`$Table$layer)
depth = depth_file$`layer[Depth]`$Table$layer
class(depth) = "numeric"
depth_plot = depth
depth_plot[depth_plot > quantile(depth, 0.8)] = quantile(depth, 0.8) + 5
brks = c(min(depth_plot) - 2, seq(from = quantile(depth, 0.2), to = quantile(depth, 0.8), by = 20), max(depth_plot) + 5)
## visualise covariate in this case death
image.plot(x_mid, rev(y_mid), t(depth_plot[nrow(depth_plot):1,]), breaks = brks, col = tim.colors(length(brks) - 1))
par(mfrow = c(2,1))
image(t(depth[nrow(depth):1,]), xaxt = "n", yaxt = "n")
image(xaxt = "n", yaxt = "n", t(run$depth$`1990`$values[nrow(run$depth$`1990`$values):1,]))

# Calculate gradient
merid = run$old_pref$initialisation_meridional
zonal = run$old_pref$initialisation_zonal
merid_diff = run$old_pref$meridonal_1970
summary(as.vector(merid - merid_diff))


run$old_pref$d_max
run$old_pref$standard_dev
diffuse =  run$old_pref$d_max * (1 - (run$old_pref$initialisation_preference / (run$old_pref$zeta + run$old_pref$initialisation_preference)))
std_dev = sqrt(2 * diffuse * run$old_pref$time_intervals)
## visualise it
image.plot(xaxt = "n", yaxt = "n", t(std_dev[nrow(std_dev):1,]), main = "Standard Deviation", xlab = "x", ylab = "y")
par(mfrow = c(3,1))
image.plot(xaxt = "n", yaxt = "n", t(zonal[nrow(zonal):1,]), main = "Zonal Gradient (E-W)", xlab = "x", ylab = "y")
image.plot(xaxt = "n", yaxt = "n", t(merid[nrow(merid):1,]), main = "Meridional Gradient (N-S)", xlab = "x", ylab = "y")
image.plot(xaxt = "n", yaxt = "n", t(run$old_pref$initialisation_preference[nrow(run$old_pref$initialisation_preference):1,]), main = "Preference", xlab = "x", ylab = "y")

this_pref = run$old_pref$initialisation_preference
squares = 5
my_zonal = (this_pref[,squares:n_col] - this_pref[,1:(n_col - squares + 1)]) / (squares - 1)
my_merid = (this_pref[1:(n_row - squares + 1),] - this_pref[squares:n_row, ]) / (squares - 1)

my_merid[1,]
(this_pref[1,] - this_pref[5,]) / 4
round(merid[2,], 2)
round((this_pref[3,] - this_pref[1,] ) / 2, 2)
merid[1, ]
(this_pref[2,] - this_pref[1,] )

dim(my_zonal)
dim(my_merid)
dim(this_pref)
if(squares == 5) {
  my_zonal = cbind((this_pref[, 3] - this_pref[, 1]) / 2, (this_pref[, 4] - this_pref[, 1]) / 3, my_zonal)
  my_zonal = cbind(my_zonal, (this_pref[, n_col] - this_pref[, n_col - 3]) / 3, (this_pref[, n_col] - this_pref[, n_col - 2]) / 2)
  
  my_merid = rbind((this_pref[1, ] - this_pref[squares - 2, ]) / 2,(this_pref[1,] - this_pref[squares - 1, ]) / 3, my_merid)
  my_merid = rbind(my_merid, (this_pref[n_row - 3, ] - this_pref[n_row, ]) / 3, (this_pref[n_row - 2,] - this_pref[n_row, ]) / 2)
} else if (squares == 3) {
  my_zonal = cbind((this_pref[, 2] - this_pref[, 1]), my_zonal)
  my_zonal = cbind(my_zonal, (this_pref[, n_col] - this_pref[, n_col - 1]))
  
  my_merid = rbind((this_pref[1, ] - this_pref[squares - 1, ]), my_merid)
  my_merid = rbind(my_merid, (this_pref[n_row - 1, ] - this_pref[n_row, ]))
}
dim(my_zonal)
dim(my_merid)

round(zonal[1,], 3)
round(my_zonal[1,], 3)

round(zonal[,2], 3)
round(my_zonal[,2], 3)

round(zonal[,n_col], 2)
round(my_zonal[,n_col], 2)

round(zonal[,n_col - 1], 2)
round(my_zonal[,n_col - 1], 2)

summary(as.vector(merid - my_merid))
summary(round(as.vector(zonal - my_zonal),3))
which(zonal - my_zonal != 0, arr.ind = T)
which(merid - my_merid != 0, arr.ind = T)
table(as.vector(zonal - my_zonal) == 0)
table(merid - my_merid == 0)


par(mfrow = c(2,2))
image.plot(xaxt = "n", yaxt = "n", t(zonal[nrow(zonal):1,]), main = "Zonal Gradient (E-W)", xlab = "x", ylab = "y")
#par(mfrow = c(2,1))
image.plot(xaxt = "n", yaxt = "n", t(merid[nrow(merid):1,]), main = "Meridional Gradient (N-S)", xlab = "x", ylab = "y")
image.plot(xaxt = "n", yaxt = "n", t(my_zonal[nrow(my_zonal):1,]), main = "Zonal Gradient (E-W)", xlab = "x", ylab = "y")
image.plot(xaxt = "n", yaxt = "n", t(my_merid[nrow(my_merid):1,]), main = "Meridional Gradient (N-S)", xlab = "x", ylab = "y")



## re-write the C++ code
square_offset = (squares - 1) / 2
my_zonal_alt = my_merid_alt = matrix(NA, nrow = n_row, ncol = n_col)
for(row in 1:n_row) {
  for(col in 1:n_col) {
    # meredional_gradient
    if(row <= square_offset) {
      offset = square_offset - (row - 1)
      my_merid_alt[row, col] = (this_pref[1, col] - this_pref[squares - offset, col]) / (row + 1);
    } else if (row > (n_row - square_offset)) {
      offset = (n_row - (row - 1))
      my_merid_alt[row, col] = (this_pref[n_row - (1 +  offset), col] - this_pref[n_row, col]) / (offset + 1);
    } else {
      my_merid_alt[row, col] = (this_pref[row - square_offset, col] - this_pref[row + square_offset, col]) / (squares - 1);
    }
    
    #zonal
    if(col <= square_offset) {
      offset = square_offset - (col - 1)
      my_zonal_alt[row, col] = (this_pref[row, squares - offset] - this_pref[row, 1]) / (col + 1);
    } else if (col > (n_col - square_offset)) {
      offset = (n_col - (col - 1))
      my_zonal_alt[row, col] = (this_pref[row, n_col] - this_pref[row, n_col - (1 +  offset)]) / (offset + 1);
    } else {
      my_zonal_alt[row, col] = (this_pref[row, col + square_offset] - this_pref[row, col - square_offset]) / (squares - 1);
    }
  }
}

## re-write the C++ code
my_zonal_alt_c = my_merid_alt_c = matrix(NA, nrow = n_row, ncol = n_col)
for(row in 0:(n_row - 1)) {
  for(col in 0:(n_col - 1)) {
    # meredional_gradient
    if(row < square_offset) {
      offset = square_offset - (row)
      my_merid_alt_c[row + 1, col + 1] = (this_pref[1, col + 1] - this_pref[squares - offset, col + 1]) / (row + 2);
    } else if (row >= (n_row - square_offset)) {
      offset = (n_row - (row ))
      my_merid_alt_c[row + 1, col + 1] = (this_pref[n_row - (1 +  offset), col + 1] - this_pref[n_row, col + 1]) / (offset + 1);
    } else {
      my_merid_alt_c[row + 1, col + 1] = (this_pref[row - square_offset + 1, col + 1] - this_pref[row + square_offset + 1, col + 1]) / (squares - 1);
    }
    
    #zonal
    if(col < square_offset) {
      offset = square_offset - (col)
      my_zonal_alt_c[row + 1, col + 1] = (this_pref[row + 1, squares - offset] - this_pref[row + 1, 1]) / (col + 2);
    } else if (col >= (n_col - square_offset)) {
      offset = (n_col - (col))
      my_zonal_alt_c[row + 1, col + 1] = (this_pref[row + 1, n_col] - this_pref[row + 1, n_col - (1 +  offset)]) / (offset + 1);
    } else {
      my_zonal_alt_c[row + 1, col + 1] = (this_pref[row + 1, col + square_offset + 1] - this_pref[row + 1, col - square_offset + 1]) / (squares - 1);
    }
  }
}
which(zonal - my_zonal_alt_c != 0, arr.ind = T)
which(merid - my_merid_alt_c != 0, arr.ind = T)
table(round(zonal - my_zonal_alt_c, 4) == 0)
table(round(merid - my_merid_alt_c,4) == 0)

merid[row + 1, col + 1]
round(my_zonal_alt_c[1,],3)
round(zonal[1,],3)

round(zonal - my_zonal_alt_c, 2)
round(merid - my_merid_alt_c, 2)

## numeric layer
par(mfrow = c(2,1))
image.plot(xaxt = "n", yaxt = "n", t(run$old_bio$`1990`$values[nrow(run$old_bio$`1990`$values):1,]), main = "Mature Biomass", xlab = "x", ylab = "y")
image.plot(xaxt = "n", yaxt = "n", t(run$young_bio$`1990`$values[nrow(run$young_bio$`1990`$values):1,]), main = "Immature Biomass", xlab = "x", ylab = "y")

par(mfrow = c(1,1))
plot(names(run$pref_fun_young$Values), run$pref_fun_young$Values, type = "l", lwd = 3, xlab = "biomass", ylab = "Preference")
vals = seq(from = 10, to = 100, by = 5)
pref = dnorm(vals, 50, 10)
pref = pref / max(pref)
plot(vals, pref, type = "l", lwd = 3)
## calcualte the gradients
grad = (pref[3:(length(pref))] - pref[1:(length(pref) - 2)]) / 2
round(grad, 2)

par(mfrow = c(2,2))
plot(x_mid, apply(run$old_pref$initialisation_preference, MARGIN = 2, mean), type = "l", lwd = 3, ylab = "Average Preference")
plot(x_mid, apply(zonal, MARGIN = 2, mean), type = "l", lwd = 3, ylab = "Average Zonal")
lines(x_mid, apply(my_zonal, MARGIN = 2, mean), lwd = 3, lty = 2, col = "red")
abline(h = 0, lty = 2, col = "blue", lwd = 3)
plot(y_mid, apply(run$old_pref$initialisation_preference, MARGIN = 1, mean), type = "l", lwd = 3, ylab = "Average Preference")
plot(y_mid, apply(merid, MARGIN = 1, mean), type = "l", lwd = 3, ylab = "Average Zonal")
abline(h = 0, lty = 2, col = "blue", lwd = 3)
##################################################
## Simulate some movements from different locations.
##################################################
d_max = run$old_pref$d_max
zeta = run$old_pref$zeta
d_max_s = c(0.01, 0.05, 0.1, 0.15,0.2)
zeta_s = c(0.01,0.05,0.1, 0.15,0.2)

for(d_ndx in 1:length(d_max_s)) {
  for(z_ndx in 1:length(zeta_s)) {
    d_max = d_max_s[d_ndx]
    zeta = zeta_s[z_ndx]
    diffuse =  d_max * (1 - (run$old_pref$initialisation_preference / (zeta + run$old_pref$initialisation_preference)))
    std_dev = sqrt(2 * diffuse * run$old_pref$time_intervals)
    
    start_mat = expand.grid(x = 1:length(x_mid), y = 1:length(y_mid))
    n_steps = 100
    cols = colorRampPalette(brewer.pal(n = 7, name = "Greys"))(n_steps)
    final_loc = NULL
    start_loc = NULL
    n_sim_agents = 20; # from each starting point
    dev.off()
    image.plot(x_mid, abs(y_mid), t(run$old_pref$initialisation_preference[nrow(run$old_pref$initialisation_preference):1,]))
    plot(1,1,type = "n", ylim = c(min(y_lim), max(y_lim)), xlim = c(min(x_lim), max(x_lim)), xlab = "x", ylab = "y")
    coords_for_each_start = list()
    for(k in 1:nrow(start_mat)) {
      set.seed(123)
      for(agent_ndx in 1:n_sim_agents) {
        x_cell = start_mat[k,1]
        y_cell = start_mat[k,2]
        lat = y_mid[y_cell]
        lon = x_mid[x_cell]
        coords = NULL;
        coords = rbind(coords, c(lon, lat))
        start_loc = rbind(start_loc, c(lon, lat))
        #points(lon, lat, pch = 16, col = "red")
        for(i in 2:n_steps) {
          #points(lon, lat, pch = 16, col = "gray60")
          #update_lon = lon + rnorm(1, zonal[y_cell, x_cell], std_dev[y_cell, x_cell])
          #update_lat = lat + rnorm(1, merid[y_cell, x_cell], std_dev[y_cell, x_cell])
          
          update_lon = lon + (zonal[y_cell, x_cell] + rnorm(1) * std_dev[y_cell, x_cell])
          update_lat = lat + (merid[y_cell, x_cell] + rnorm(1)* std_dev[y_cell, x_cell])
          #update_lon = lon + rnorm(1, my_zonal[y_cell, x_cell], std_dev[y_cell, x_cell])
          #update_lat = lat + rnorm(1, my_merid[y_cell, x_cell], std_dev[y_cell, x_cell])
          old_lon = lon
          old_lat = lat
          # check boundaries
          if(update_lon < max(x_lim) & update_lon > min(x_lim)) {
            lon = update_lon
            # update cell coordinates
            for(j in 2:length(x_lim)) {
              if (lon < x_lim[j]) {
                x_cell = j - 1
                break
              }
            }
          }
          if(update_lat < max(y_lim) & update_lat > min(y_lim)) {
            lat = update_lat
            for(j in 2:length(y_lim)) {
              if (lat > y_lim[j]) {
                y_cell = j - 1
                break
              }
            }
          }
          coords = rbind(coords, c(lon, lat))
          
          #points(c(lon,old_lon), c(lat, old_lat), type = "o", pch = 16, col = adjustcolor(col = cols[i], alpha.f = 0.5), lwd = 3)
          #Sys.sleep(0.3)
        }
        #points(lon, lat, pch = 16, col = "green")
        final_loc = rbind(final_loc, c(lon, lat))
        if(agent_ndx == 1) {
          coords_for_each_start[[k]] = coords
        }
      }
    }
    #hist(rnorm(10000, 2, 1.3))
    #hist(2 + rnorm(10000) * 1.3)
    
    png(filename = paste0("std_settings_mature_nsteps_",n_steps,"_zeta_",zeta, "_dmax_", d_max, ".png"), height = 8, width = 6, res = 250, units = "in")
    par(mfrow = c(2,1))
    plot(1,1,type = "n", ylim = c(min(y_lim), max(y_lim)), xlim = c(min(x_lim), max(x_lim)), xlab = "x", ylab = "y", main = paste0("zeta = ", zeta, " d_max = ", d_max))
    image(x_mid, rev(y_mid), t(run$old_pref$initialisation_preference[nrow(run$old_pref$initialisation_preference):1,]), add = T)
    #points(start_loc[,1], start_loc[,2], pch = 16, col = "blue")
    
    plot(1,1,type = "n", ylim = c(min(y_lim), max(y_lim)), xlim = c(min(x_lim), max(x_lim)), xlab = "x", ylab = "y")
    image(x_mid, rev(y_mid), t(run$old_pref$initialisation_preference[nrow(run$old_pref$initialisation_preference):1,]), add = T)
    points(final_loc[,1], final_loc[,2], pch = 16, col = adjustcolor(col = "purple", alpha.f = 0.2))
    dev.off()
    #image.plot(xaxt = "n", yaxt = "n", t(run$old_bio$`1990`$values[nrow(run$old_bio$`1990`$values):1,]), main = "Biomass", xlab = "x", ylab = "y")
  }
}

## look at the amount of movement
diff(coords_for_each_start[[1]][,1])
par(mfrow = c(1,1))
plot(1,1,type = "n", ylim = c(min(y_lim), max(y_lim)), xlim = c(min(x_lim), max(x_lim)), xlab = "x", ylab = "y")
image(x_mid, rev(y_mid), t(run$old_pref$initialisation_preference[nrow(run$old_pref$initialisation_preference):1,]), add = T)
for(i in 1:length(coords_for_each_start)) {
  points(coords_for_each_start[[i]][,1], coords_for_each_start[[i]][,2], type = "o", pch = 16, col =cols, lwd = 3)
  points(coords_for_each_start[[i]][1,1],coords_for_each_start[[i]][1,2], col = "red", pch = 16)
  points(coords_for_each_start[[i]][nrow(coords_for_each_start[[i]]),1], coords_for_each_start[[i]][nrow(coords_for_each_start[[i]]),2], col = "black", pch = 16)
  Sys.sleep(0.4)
  points(coords_for_each_start[[i]][,1], coords_for_each_start[[i]][,2], type = "o", pch = 16, col ="white", lwd = 3)
}


par(mfrow = c(3,1))
image.plot(x_mid, y_mid, t(zonal[nrow(zonal):1,]), main = "Zonal Gradient (E-W)", xlab = "x", ylab = "y")
points(final_loc[,1], final_loc[,2], pch = 16, col = "purple")
image.plot(x_mid, y_mid, t(merid[nrow(merid):1,]), main = "Meridional Gradient (N-S)", xlab = "x", ylab = "y")
points(final_loc[,1], final_loc[,2], pch = 16, col = "purple")
image.plot(x_mid, rev(y_mid), t(run$old_pref$initialisation_preference[nrow(run$old_pref$initialisation_preference):1,]), main = "Preference", xlab = "x", ylab = "y")
points(final_loc[,1], final_loc[,2], pch = 16, col = "purple")

#########
## Try a different D_max
#########
diffuse =  0.1 * (1 - (run$old_pref$initialisation_preference / (run$old_pref$zeta + run$old_pref$initialisation_preference)))
std_dev = sqrt(2 * diffuse * run$old_pref$time_intervals)
summary(std_dev)
lat = y_mid[start_y]
lon = x_mid[start_x]
dev.off()
image.plot(x_mid, y_mid, t(depth_plot[nrow(depth_plot):1,]), breaks = brks, col = tim.colors(length(brks) - 1))
points(lon, lat, pch = 16, col = cols[1])
x_cell = start_x
y_cell = start_y
coords = NULL;
coords = rbind(c(lon, lat))
set.seed(123)
for(i in 1:n_steps) {
  cat("i = ", i,"\n")
  #points(lon, lat, pch = 16, col = "gray60")
  update_lon = lon + rnorm(1, zonal[y_cell, x_cell], std_dev[y_cell, x_cell])
  update_lat = lat + rnorm(1, merid[y_cell, x_cell], std_dev[y_cell, x_cell])
  # check boundaries
  if(update_lon < max(x_mid) & update_lon > min(x_mid))
    lon = update_lon
  if(update_lat < max(y_mid) & update_lat > min(y_mid))
    lat = update_lat
  ## else stay where you are
  # update cell coordinates
  for(j in 2:length(x_lim)) {
    if (lon < x_lim[j]) {
      x_cell = j - 1
      break
    }
  }
  for(j in 2:length(y_lim)) {
    if (lat < y_lim[j]) {
      y_cell = j - 1
      break
    }
  }
  coords = rbind(coords, c(lon, lat))
  points(lon, lat, pch = 16, col = cols[i])
  #Sys.sleep(0.5)
}

#########
## Try a different Zeta
#########
diffuse =  0.01 * (1 - (run$old_pref$initialisation_preference / (run$old_pref$zeta + run$old_pref$initialisation_preference)))
std_dev = sqrt(2 * diffuse * run$old_pref$time_intervals)
summary(std_dev)
lat = y_mid[start_y]
lon = x_mid[start_x]
dev.off()
image.plot(x_mid, y_mid, t(run$old_pref$initialisation_preference[nrow(run$old_pref$initialisation_preference):1,]))
points(lon, lat, pch = 16, col = "black")
x_cell = start_x
y_cell = start_y
coords = NULL;
coords = rbind(c(lon, lat))
set.seed(123)
for(i in 1:n_steps) {
  cat("i = ", i,"\n")
  points(lon, lat, pch = 16, col = "gray60")
  update_lon = lon + rnorm(1, zonal[y_cell, x_cell], std_dev[y_cell, x_cell])
  update_lat = lat + rnorm(1, merid[y_cell, x_cell], std_dev[y_cell, x_cell])
  # check boundaries
  if(update_lon < max(x_mid) & update_lon > min(x_mid))
    lon = update_lon
  if(update_lat < max(y_mid) & update_lat > min(y_mid))
    lat = update_lat
  ## else stay where you are
  # update cell coordinates
  for(j in 2:length(x_lim)) {
    if (lon < x_lim[j]) {
      x_cell = j - 1
      break
    }
  }
  for(j in 2:length(y_lim)) {
    if (lat < y_lim[j]) {
      y_cell = j - 1
      break
    }
  }
  coords = rbind(coords, c(lon, lat))
  points(lon, lat, pch = 16, col = "black")
  #Sys.sleep(0.5)
}



#########
## Try a different Zeta
#########
pref = seq(from = 0.01, to =0.99, by = 0.01)
par(mfrow= c(1,1))
plot(pref, 1 - pref / (0.2 + pref), type = "l", lwd = 3, xlab = "Preference", ylab = "multiplier on std", ylim = c(0,1))
lines(pref, 1 - pref / (0.1 + pref), type = "l", lwd = 3, lty = 2, col = "red")
lines(pref, 1 - pref / (0.05 + pref), type = "l", lwd = 3, lty = 2, col = "blue")
lines(pref, 1 - pref / (0.5 + pref), type = "l", lwd = 3, lty = 2, col = "darkgreen")
legend('topright', legend = c(0.05,0.1, 0.2, 0.5), col = c("blue", "red","black","darkgreen"), lwd = 3)
diffuse =  0.05 * (1 - (run$old_pref$initialisation_preference / (0.2 + run$old_pref$initialisation_preference)))
std_dev = sqrt(2 * diffuse * run$old_pref$time_intervals)
summary(std_dev)
lat = y_mid[start_y]
lon = x_mid[start_x]
dev.off()
image.plot(x_mid, y_mid, t(depth_plot[nrow(depth_plot):1,]), breaks = brks, col = tim.colors(length(brks) - 1))
points(lon, lat, pch = 16, col = "black")
x_cell = start_x
y_cell = start_y
coords = NULL;
coords = rbind(c(lon, lat))
set.seed(123)
for(i in 1:n_steps) {
  cat("i = ", i,"\n")
  points(lon, lat, pch = 16, col = "gray60")
  update_lon = lon + rnorm(1, zonal[y_cell, x_cell], std_dev[y_cell, x_cell])
  update_lat = lat + rnorm(1, merid[y_cell, x_cell], std_dev[y_cell, x_cell])
  # check boundaries
  if(update_lon < max(x_mid) & update_lon > min(x_mid))
    lon = update_lon
  if(update_lat < max(y_mid) & update_lat > min(y_mid))
    lat = update_lat
  ## else stay where you are
  # update cell coordinates
  for(j in 2:length(x_lim)) {
    if (lon < x_lim[j]) {
      x_cell = j - 1
      break
    }
  }
  for(j in 2:length(y_lim)) {
    if (lat < y_lim[j]) {
      y_cell = j - 1
      break
    }
  }
  coords = rbind(coords, c(lon, lat))
  points(lon, lat, pch = 16, col = "black")
  Sys.sleep(0.5)
}
