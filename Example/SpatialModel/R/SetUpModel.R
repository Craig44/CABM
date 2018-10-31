## Help Set up a hypothetical spatial model for demonstration
#' @author C.Marsh
#' @date 30/10/2018
#' @description This file just helps build layers mainly based on a hypothetical life history and spatial distribution
#' 
set.seed(123)
library(ibm)
setwd("..")
BaseDir = getwd();
LayerDir = paste0(BaseDir,"/ibm/Layers/")
RDir = paste0(BaseDir,"/R/")
setwd("R")
nrow = 10
ncol = 10

## Base layer
Base_layer = matrix(1,nrow = nrow, ncol = ncol) ## could also fill this in with are of cell
layers_to_add_to_config = c("base_layer.ibm")
create_ibm_layer(label = "Base", type = "numeric", filename = paste0(LayerDir,"base_layer.ibm"), Matrix = Base_layer)

## Cell layer categorical
cell_mat = matrix(paste(paste0("r",rep(1:10,10)),paste0("C",sort(rep(1:10,10))),sep = "-"), nrow = nrow, ncol =ncol)
layers_to_add_to_config = c(layers_to_add_to_config,"cell_layer.ibm")
create_ibm_layer(label = "Cells", type = "numeric", filename = paste0(LayerDir,"cell_layer.ibm"), Matrix = Base_layer)


## Annual Cycle
# - Quartely time steps
# - Summer, Autumn, Winter,Spring
# - Summer
# -- Ageing (implicit)
# -- Growth (0.35)
# -- Spawning
# -- Recruitment
# -- Summer Fishery + M
# - Autumn
# -- Move off shore
# -- Growth 0.25
# -- preference movement
# - Winter
# -- Preference movement
# -- Growth 0.1
# -- Winter fishery and M
# - Spring
# -- Growth 0.3
# -- Move inshore
# -- preference movement


## Box transfer movement
## Offshore W->E
line = paste(rep(1:10,10),sort(rep(1:10,10)),sep = "-")
line = c("origin_cell", line)
file = write(line,paste0(LayerDir,"offshore_cell_movements.ibm"),ncolumns = length(line))
line = paste("offshore",paste(rep(1:10,10),sort(rep(1:10,10)),sep = "-"), sep = "_")
line = c("probability_layers", line)
file = write(line,paste0(LayerDir,"offshore_cell_movements.ibm"),ncolumns = length(line), append = T)

## Create layers for this offshore movement, a purely horizontal shift
for (i in 1:nrow) {
  for (j in 1:ncol) {
    Mat = matrix(0,nrow = nrow, ncol = ncol)
    if (j <= 5) {
      Mat[i,6:ncol] = rep(1/5,5)
    } else {
      Mat[i,j] = 1
    }
    create_ibm_layer(label = paste0("offshore_",i, "-", j), type = "numeric", filename = paste0(LayerDir,"Offshore/offshore_",i, "-", j,".ibm"), Matrix = Mat, proportion = T)
    layers_to_add_to_config = c(layers_to_add_to_config, paste0("Offshore/offshore_",i, "-", j,".ibm"))
  }
}

## Inshore E->W
line = paste(rep(1:10,10),sort(rep(1:10,10)),sep = "-")
line = c("origin_cell", line)
file = write(line,paste0(LayerDir,"inshore_cell_movements.ibm"),ncolumns = length(line))
line = paste("inshore",paste(rep(1:10,10),sort(rep(1:10,10)),sep = "-"), sep = "_")
line = c("probability_layers", line)
file = write(line,paste0(LayerDir,"inshore_cell_movements.ibm"),ncolumns = length(line), append = T)

## Create layers for this offshore movement, a purely horizontal shift
for (i in 1:nrow) {
  for (j in 1:ncol) {
    Mat = matrix(0,nrow = nrow, ncol = ncol)
    if (j >= 5) {
      Mat[i,1:5] = rep(1/5,5)
    } else {
      Mat[i,j] = 1
    }
    create_ibm_layer(label = paste0("inshore_",i, "-", j), type = "numeric", filename = paste0(LayerDir,"Inshore/inshore_",i, "-", j,".ibm"), Matrix = Mat, proportion = T)
    layers_to_add_to_config = c(layers_to_add_to_config, paste0("Inshore/inshore_",i, "-", j,".ibm"))
  }
}

## Preference movement
## Offshore 
Depth = matrix(0,nrow = nrow, ncol = ncol)
## make a shelf like effect
Depth[,1] = 16
Depth[,2] = 15
Depth[,3] = 16
Depth[,4] = 18
Depth[,5] = 20
Depth[,6] = 23
Depth[,7] = 25
Depth[,8] = 40
Depth[,9] = 45
Depth[,10] = 50

## Add some pinnacles that would promote aggregation
Depth[3,8] = 33
Depth[4,8] = 33
Depth[3,9] = 33
Depth[4,9] = 33

Depth[1,7] = 44
Depth[1,8] = 44
Depth[1,9] = 44
Depth[1,10] = 50

Depth[2,7] = 39
Depth[2,8] = 39
Depth[2,9] = 39
Depth[2,10] = 39

Depth[4,7] = 39
Depth[4,8] = 39
Depth[4,9] = 39
Depth[4,10] = 50
Depth[2,10] = 50
Depth[3,10] = 50
Depth[2,7] = 39
Depth[3,7] = 39

Depth[7,8] = 33
Depth[8,8] = 33
Depth[7,7] = 33
Depth[8,7] = 33

create_ibm_layer(label = "offshore_Depth", type = "numeric", filename = paste0(LayerDir,"Depth.ibm"), Matrix = Depth, proportion = F)
layers_to_add_to_config = c(layers_to_add_to_config, "Depth.ibm")

## Inshore preference
## ideal spawning areas
Inshore_pref = matrix(0,nrow = nrow, ncol = ncol);
Inshore_pref[,1] = rnorm(nrow, 20, 4)
Inshore_pref[,2] = rnorm(nrow, 15, 4)
Inshore_pref[,3] = rnorm(nrow, 12, 4)
Inshore_pref[,4] = rnorm(nrow, 10, 3)
Inshore_pref[,5] = rnorm(nrow, 3, 2)
Inshore_pref[,5:ncol] = rnorm(nrow * 6, 1, 2)

create_ibm_layer(label = "inshore_pref", type = "numeric", filename = paste0(LayerDir,"Inshore_pref.ibm"), Matrix = Inshore_pref, proportion = F)
layers_to_add_to_config = c(layers_to_add_to_config, "Inshore_pref.ibm")

## SSB is an inshore event, any stragglers on the shelf won't contribute
SSB = matrix(0,nrow = nrow, ncol = ncol);
SSB[,1:5] = 1
create_ibm_layer(label = "SSB_layer", type = "integer", filename = paste0(LayerDir,"SSB_layer.ibm"), Matrix = SSB, proportion = F)
layers_to_add_to_config = c(layers_to_add_to_config, "SSB_layer.ibm")

## Recruits pop up in the inshore as well
recruits = matrix(0,nrow = nrow, ncol = ncol);
recruits[,1:4] = 1 / (4*nrow)
create_ibm_layer(label = "recruit_layer", type = "numeric", filename = paste0(LayerDir,"recruitment_layer.ibm"), Matrix = recruits, proportion = T)
layers_to_add_to_config = c(layers_to_add_to_config, "recruitment_layer.ibm")


## Now add all the !include statements to the config.ibm file
## delete all the other includes before running this loop
#config_file = paste0(BaseDir,"/ibm/config.ibm")
#for (i in 1:length(layers_to_add_to_config)) {
#  line = paste0("!include Layers/",layers_to_add_to_config[i])
#  write(line, file = config_file, append = T)
#}
