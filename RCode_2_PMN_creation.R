#################################################_
### Building a Predictive Multilayer Network 
### Project: Predictive Multilayer Networks
### Author: Tyler N. McFadden
###
### Code accompanying manuscript: 
### McFadden TN, Dirzo R. Harnessing multilayer networks to predict 
### metacommunity responses to global environmental change.
### Submitted to Ecology.
###
###
### About this script:
###   This script contains 8 sections constituting a single pipeline for building
###   a Predictive Multilayer Network from land cover and elevation input data,
###   and for calculating network connectance. Once you have confirmed that all parameters
###   are correct, the entire script can be run at once in the 'Jobs' section of RStudio.
###   The process of running this script will take up to several hours.
###
###   Code for calculating additional network metrics (betweenness centrality and robustness) 
###   is provided in separate R scripts.
###
###   Prior to running this script:
###       - Confirm working directory is correct
###       - Confirm the landscape input file is correct
###       - Confirm all packages have been installed
###       - Confirm parameters at the beginning of each intralayer link section
###             are correct (sections 3, 4, and 5)
###
###   The working directory should contain the following files to begin with:
###       - landscape input file in .tif format
###       - elevation raster (elev_map_sim.tif)
###       - edge list for null network (null_edges_1.csv, null_edges_2.csv, null_edges_3.csv -- file was separated into 3 to comply with GitHub file size limits)
###       - this R script
###
#################################################_

# Set working directory
setwd("/Users/tylermcfadden/Desktop/Stanford/Multilayer_network/Manuscript/Code")

# Load libraries
library(landscapeR)
library(raster)
library(sp)
library(rgdal)
library(reshape2)
library(foreach)
library(doParallel)
library(gdistance)
library(dplyr)
library(tidyverse)
library(igraph)
library(gtools)
library(Matrix)

# Define landscape input file
landscape.file.name <- "starting_landscape.tif"

#################################################_
##### 1. Habitat Metrics #####
###
###
### About this section:
###   This section imports a raster landscape (tiff file) and derives from it a series
###   of habitat metrics based on land cover. These metrics will later be used as inputs for
###   species distribution models.
#################################################_


##########################################################_
##### Calculating habitat metrics (inputs for SDMs) #####_
##########################################################_

# Import landscape file
rs <- raster(landscape.file.name)

# Create blank raster with same dimensions and CRS as landscape
r <- setValues(rs, 0)


##### Land cover calculations #####_

# Calculate land cover within radius
x <- coordinates(rs) # get the coordinates for each pixel in landscape
lc50 <- raster::extract(rs, x, buffer = 50, small = TRUE)
lc100 <- raster::extract(rs, x, buffer = 100, small = TRUE) 
lc250 <- raster::extract(rs, x, buffer = 250, small = TRUE) 
lc500 <- raster::extract(rs, x, buffer = 500, small = TRUE) 


### Forest cover in 100m radius 
# For loop to get percent cover of given class (or multiple classes) for each pixel
lc.prop <- NULL
for(i in 1:length(lc100)){
  lc.prop[i] <- sum(unlist(lc100[i]) %in% c(3,4))/length(unlist(lc100[i])) # change number to change cover class. change lc file name to change radius
}

# Plot percent cover 
r2 <- raster(nrows=100, ncols=100, xmn=656000, xmx=659000, ymn=5588000, 
             ymx=5591000, resolution=30, crs = CRS("+init=epsg:5362"))
lc.matrix <- as.matrix(lc.prop)
r2 <- setValues(r, lc.matrix)
plot(r2)

# Export raster
writeRaster(r2, "forest_cover_100m.tif", format = "GTiff")


### Forest cover in 250m radius 
# For loop to get percent cover of given class (or multiple classes) for each pixel
lc.prop <- NULL
for(i in 1:length(lc250)){
  lc.prop[i] <- sum(unlist(lc250[i]) %in% c(3,4))/length(unlist(lc250[i])) # change number to change cover class. change lc file name to change radius
}

# Plot percent cover 
r2 <- raster(nrows=100, ncols=100, xmn=656000, xmx=659000, ymn=5588000, 
             ymx=5591000, resolution=30, crs = CRS("+init=epsg:5362"))
lc.matrix <- as.matrix(lc.prop)
r2 <- setValues(r, lc.matrix)
plot(r2)

# Export raster
writeRaster(r2, "forest_cover_250m.tif", format = "GTiff")


### Forest cover in 500m radius 
# For loop to get percent cover of given class (or multiple classes) for each pixel
lc.prop <- NULL
for(i in 1:length(lc500)){
  lc.prop[i] <- sum(unlist(lc500[i]) %in% c(3,4))/length(unlist(lc500[i])) # change number to change cover class. change lc file name to change radius
}

# Plot percent cover 
r2 <- raster(nrows=100, ncols=100, xmn=656000, xmx=659000, ymn=5588000, 
             ymx=5591000, resolution=30, crs = CRS("+init=epsg:5362"))
lc.matrix <- as.matrix(lc.prop)
r2 <- setValues(r, lc.matrix)
plot(r2)

# Export raster
writeRaster(r2, "forest_cover_500m.tif", format = "GTiff")


### Primary forest cover in 100m radius 
# For loop to get percent cover of given class (or multiple classes) for each pixel
lc.prop <- NULL
for(i in 1:length(lc100)){
  lc.prop[i] <- sum(unlist(lc100[i]) %in% c(4))/length(unlist(lc100[i])) # change number to change cover class. change lc file name to change radius
}

# Plot percent cover 
r2 <- raster(nrows=100, ncols=100, xmn=656000, xmx=659000, ymn=5588000, 
             ymx=5591000, resolution=30, crs = CRS("+init=epsg:5362"))
lc.matrix <- as.matrix(lc.prop)
r2 <- setValues(r, lc.matrix)
plot(r2)

# Export raster
writeRaster(r2, "primary_forest_100m.tif", format = "GTiff")


### Primary forest cover in 250m radius 
# For loop to get percent cover of given class (or multiple classes) for each pixel
lc.prop <- NULL
for(i in 1:length(lc250)){
  lc.prop[i] <- sum(unlist(lc250[i]) %in% c(4))/length(unlist(lc250[i])) # change number to change cover class. change lc file name to change radius
}

# Plot percent cover 
r2 <- raster(nrows=100, ncols=100, xmn=656000, xmx=659000, ymn=5588000, 
             ymx=5591000, resolution=30, crs = CRS("+init=epsg:5362"))
lc.matrix <- as.matrix(lc.prop)
r2 <- setValues(r, lc.matrix)
plot(r2)

# Export raster
writeRaster(r2, "primary_forest_250m.tif", format = "GTiff")


### Primary forest cover in 500m radius 
# For loop to get percent cover of given class (or multiple classes) for each pixel
lc.prop <- NULL
for(i in 1:length(lc500)){
  lc.prop[i] <- sum(unlist(lc500[i]) %in% c(4))/length(unlist(lc500[i])) # change number to change cover class. change lc file name to change radius
}

# Plot percent cover 
r2 <- raster(nrows=100, ncols=100, xmn=656000, xmx=659000, ymn=5588000, 
             ymx=5591000, resolution=30, crs = CRS("+init=epsg:5362"))
lc.matrix <- as.matrix(lc.prop)
r2 <- setValues(r, lc.matrix)
plot(r2)

# Export raster
writeRaster(r2, "primary_forest_500m.tif", format = "GTiff")


### Identifying forest edges 
# For loop to identify forest pixels immediately adjacent to non-forest
fp.edge <- NULL
for(i in 1:length(lc50)){
  fp.edge[i] <- ifelse((((rs[i] > 2) & (1 %in% unlist(lc50[i]))) | ((rs[i] > 2) & (2 %in% unlist(lc50[i])))), 1, 0)
}

# Plot forest edges
r3 <- raster(nrows=100, ncols=100, xmn=656000, xmx=659000, ymn=5588000, 
             ymx=5591000, resolution=30, crs = CRS("+init=epsg:5362"))
fp.edge.matrix <- as.matrix(fp.edge)
r3 <- setValues(r, fp.edge.matrix)
plot(r3)

# Export raster
writeRaster(r3, "forest_edges.tif", format = "GTiff") # forest cells adjacent to non-forest cells


### Identifying non-forest cells adjacent to forest cells
# For loop to identify non-forest pixels immediately adjacent to non-forest
fp.edge <- NULL
for(i in 1:length(lc50)){
  fp.edge[i] <- ifelse((((rs[i] < 3) & (3 %in% unlist(lc50[i]))) | ((rs[i] < 3) & (4 %in% unlist(lc50[i])))), 1, 0)
}

# Plot forest edges
r3 <- raster(nrows=100, ncols=100, xmn=656000, xmx=659000, ymn=5588000, 
             ymx=5591000, resolution=30, crs = CRS("+init=epsg:5362"))
fp.edge.matrix <- as.matrix(fp.edge)
r3 <- setValues(r, fp.edge.matrix)
plot(r3)

# Export raster
writeRaster(r3, "forest_non_forest_edge.tif", format = "GTiff") # non-forest cells adjacent to forest cells


### Identifying forest interior
# For loop to identify forest interior (>100m from forest edge)
for.int <- NULL
for(i in 1:length(lc100)){
  for.int[i] <- ifelse((((rs[i] > 2) & !(1 %in% unlist(lc100[i])) & !(2 %in% unlist(lc100[i])))), 1, 0)
}

# Plot forest interior
r4 <- raster(nrows=100, ncols=100, xmn=656000, xmx=659000, ymn=5588000, 
             ymx=5591000, resolution=30, crs = CRS("+init=epsg:5362"))
for.int.matrix <- as.matrix(for.int)
r4 <- setValues(r, for.int.matrix)
plot(r4)

# Export raster
writeRaster(r4, "forest_interior.tif", format = "GTiff") # any forest >100m from non-forest


### Identifying primary forest interior
# For loop to identify forest interior (primary forest >100m from any forest edge)
for.int <- NULL
for(i in 1:length(lc100)){
  for.int[i] <- ifelse((((rs[i] == 4) & !(1 %in% unlist(lc100[i])) & !(2 %in% unlist(lc100[i])))), 1, 0)
}

# Plot forest interior
r4 <- raster(nrows=100, ncols=100, xmn=656000, xmx=659000, ymn=5588000, 
             ymx=5591000, resolution=30, crs = CRS("+init=epsg:5362"))
for.int.matrix <- as.matrix(for.int)
r4 <- setValues(r, for.int.matrix)
plot(r4)

# Export raster
writeRaster(r4, "primary_forest_interior.tif", format = "GTiff") # primary forest >100m from non-forest



#################################################_
#####
##### 2. Species Distribution Models #####
###
###
### About this section:
###   This section imports land cover and habitat metric data (produced in the 
###   Landscape simulation script and Habitat metric section of this script), and uses habitat 
###   preferences of seven hypothetical species to predict their probability of occurrence across
###   the landscape
#################################################_


##############################################_
### Habitat associations
##############################################_

# Create dataframe containing the variables of interest for predicting species distributions
x <- coordinates(raster(landscape.file.name)) # cell number and utms (easting then northing)
hab.metrics1 <- as.data.frame(x) # create matrix
colnames(hab.metrics1) <- c("easting", "northing") # change column names

# add rows for each variable listed above
hab.metrics1$cover.type <- melt(t(as.matrix(raster(landscape.file.name))))[,3] # add column for landcover typer
hab.metrics1$elevation <- melt(t(as.matrix(raster("elev_map_sim.tif"))))[,3] 
hab.metrics1$forest.cover100 <- melt(t(as.matrix(raster("forest_cover_100m.tif"))))[,3] 
hab.metrics1$forest.cover250 <- melt(t(as.matrix(raster("forest_cover_250m.tif"))))[,3] 
hab.metrics1$forest.cover500 <- melt(t(as.matrix(raster("forest_cover_500m.tif"))))[,3] 
hab.metrics1$primary.forest500 <- melt(t(as.matrix(raster("primary_forest_500m.tif"))))[,3]
hab.metrics1$forest.interior <- melt(t(as.matrix(raster("forest_interior.tif"))))[,3]
hab.metrics1$forest.edge <- melt(t(as.matrix(raster("forest_edges.tif"))))[,3] # forest cells adjacent to non-forest
hab.metrics1$non.forest.edge <- melt(t(as.matrix(raster("forest_non_forest_edge.tif"))))[,3] # non-forest cells adjacent to forest

head(hab.metrics1)

# Export csv of habitat metrics dataframe
write.csv(hab.metrics1, "habitat_metrics.csv")

#### Assigning habitat preferences

# hummingbird 1 - forest specialist trapliner
hab.metrics1$hummer1 <- ifelse(hab.metrics1$cover.type >2, # prob of occurence is zero in cropland and pasture, multiplied by 0.6 for secondary forest, and multiplied by 0.9 for primary forest (1.5 times more likely to occur in primary vs secondary forest)
                               ifelse(hab.metrics1$cover.type == 3, 0.6 * (1 / (1 + exp(1)^(-20*(hab.metrics1$forest.cover250 - 0.5))))
                                      , 0.9 * (1 / (1 + exp(1)^(-20*(hab.metrics1$forest.cover250 - 0.5))))), 0)

h1 <- raster(nrows=100, ncols=100, xmn=656000, xmx=659000, ymn=5588000, 
             ymx=5591000, resolution=30, crs = CRS("+init=epsg:5362"))
hummer1.matrix <- as.matrix(hab.metrics1$hummer1)
h1 <- setValues(h1, hummer1.matrix)
plot(h1)



# hummingbird 2 - territorial forest generalist
hab.metrics1$hummer2 <- 0.9 / (1 + exp(1)^(-10*(hab.metrics1$forest.cover100 - 0.2))) # occurence prob has logistic relationship with forest cover at 100m
plot(hab.metrics1$hummer2 ~ hab.metrics1$forest.cover100)

h2 <- raster(nrows=100, ncols=100, xmn=656000, xmx=659000, ymn=5588000, 
             ymx=5591000, resolution=30, crs = CRS("+init=epsg:5362"))
hummer2.matrix <- as.matrix(hab.metrics1$hummer2)
h2 <- setValues(h2, hummer2.matrix)
plot(h2)


# Bumblebee - generalist bee, needs forest for nesting
hab.metrics1$bumblebee <- pmax((0.9 / (1 + exp(1)^(-20*(hab.metrics1$primary.forest500 - 0.1)))),
                               (0.9 / (1 + exp(1)^(-10*(hab.metrics1$forest.cover500 - 0.5))))
) # selects the greater value of two logistic relationships
# the first indicates that occurence increases with primary forest cover at 500m. Requires some small amount for nesting habitat
# the second indicates occurence increases with secondary forest cover, though more secondary forest is needed to get same nesting habitat
bb <- raster(nrows=100, ncols=100, xmn=656000, xmx=659000, ymn=5588000, 
             ymx=5591000, resolution=30, crs = CRS("+init=epsg:5362"))
bb.matrix <- as.matrix(hab.metrics1$bumblebee)
bb <- setValues(bb, bb.matrix)
plot(bb)


# Plant 1 - epiphytic bromeliad
hab.metrics1$plant1 <- ifelse(hab.metrics1$cover.type == 4, 0.9, 0) # only occurs in primary forest
p1 <- raster(nrows=100, ncols=100, xmn=656000, xmx=659000, ymn=5588000, 
             ymx=5591000, resolution=30, crs = CRS("+init=epsg:5362"))
p1.matrix <- as.matrix(hab.metrics1$plant1)
p1 <- setValues(p1, p1.matrix)
plot(p1)


# Plant 2 - vine
hab.metrics1$plant2 <- ifelse(hab.metrics1$cover.type >2, # only occurs in forest. 1.5 times more likely at forest edge than interior (edge defined as directly adjacent to non-forest)
                              ifelse(hab.metrics1$forest.edge, 0.9, 0.6), 0)

p2 <- raster(nrows=100, ncols=100, xmn=656000, xmx=659000, ymn=5588000, 
             ymx=5591000, resolution=30, crs = CRS("+init=epsg:5362"))
p2.matrix <- as.matrix(hab.metrics1$plant2)
p2 <- setValues(p2, p2.matrix)
plot(p2)


# Plant 3 - understory herb
test <- 1 / (1 + exp(1)^(-0.01*(hab.metrics1$elevation - 300))) 
plot(test ~ hab.metrics1$elevation)
hab.metrics1$plant3 <- ifelse(hab.metrics1$cover.type >2, # occurs only in forest. increases with elevation. Twice as likely in forest interior vs edge (interior = >100m from edge)
                              ifelse(hab.metrics1$forest.interior, 0.9 * (1 / (1 + exp(1)^(-0.01*(hab.metrics1$elevation- 300))))
                                     , 0.45 * (1 / (1 + exp(1)^(-0.05*(hab.metrics1$elevation - 300))))), 0)
p3 <- raster(nrows=100, ncols=100, xmn=656000, xmx=659000, ymn=5588000, 
             ymx=5591000, resolution=30, crs = CRS("+init=epsg:5362"))
p3.matrix <- as.matrix(hab.metrics1$plant3)
p3 <- setValues(p3, p3.matrix)
plot(p3)


# Plant 4 - shade-intolerant shrub # this one needs some works still
hab.metrics1$plant4 <- ifelse(hab.metrics1$cover.type != 1, # occurs in pasture (more likely with more forest within 250m) or at low levels in forest edges
                              ifelse(hab.metrics1$forest.edge ==1, 0.5, 0)
                              , 0.2 + 0.6 / (1 + exp(1)^(-10*(hab.metrics1$forest.cover250 - 0.1))))
p4 <- raster(nrows=100, ncols=100, xmn=656000, xmx=659000, ymn=5588000, 
             ymx=5591000, resolution=30, crs = CRS("+init=epsg:5362"))
p4.matrix <- as.matrix(hab.metrics1$plant4)
p4 <- setValues(p4, p4.matrix)
plot(p4)


# Map of 100% occurrence probability for connectance 'null model'
r <- raster(nrows=100, ncols=100, xmn=656000, xmx=659000, ymn=5588000, 
            ymx=5591000, resolution=30, crs = CRS("+init=epsg:5362"))
complete.occurrence <- setValues(r, 1)
plot(complete.occurrence)


# Export species distribution rasters
writeRaster(h1, "h1_distribution.tif", format = "GTiff") 
writeRaster(h2, "h2_distribution.tif", format = "GTiff") 
writeRaster(bb, "bb_distribution.tif", format = "GTiff") 
writeRaster(p1, "p1_distribution.tif", format = "GTiff") 
writeRaster(p2, "p2_distribution.tif", format = "GTiff") 
writeRaster(p3, "p3_distribution.tif", format = "GTiff") 
writeRaster(p4, "p4_distribution.tif", format = "GTiff")
writeRaster(complete.occurrence, "complete.occurrence.tif", format = "GTiff")

#################################################_
#####
##### 3. Intralayer links - Hummingbird 1 #####
###
###
### About this section:
###   This section uses the gdistance package to calculate least cost paths for a given
###   species across a hypothetical landscape. The input is a species distribution model
###   in tif format in which each cell contains data on abundance or probability of occurrence.
###   The output is an edgelist reporting the least cost path between each pair of raster cells.
###   Took ~9 minutes to run in the RStudio jobs section using 7 cores.
###   Prior to running code, confirm working directory define parameters, and specify the 
###   import and export file names.
#################################################_

# Define parameters, and specify import/export file names
buff <- 1000   # buffer width based on species movement or dispersal distance (1000 for H1, 500 for H2 and BB)
numcores <- 7  # number of cores to be used for parallel computing
aoi <- 1:10000  # Define range of raster cells to be used in connectivity calculations. Only
# works for a continuous range of numbers. If range does not start at 1, 
# see additional instructions in 'Format data output' section.
import.file.name <- "h1_distribution.tif" # file in tif format
export.file.name <- "h1_leastcost_edges.csv" # file in csv format


# Start parallel computation
registerDoParallel(numcores) 


#################################################_
##### Create and correct transition matrix #####_
#################################################_

# start with raster layer for the species distribution model
sdm <- raster(import.file.name)
plot(sdm)

# add small number to all cells so that there are no cells with prob = 0
nz <- sdm + 0.001 
plot(nz)

# create Transition matrix
trans.nz <- transition(nz, transitionFunction=mean, directions=8) # use mean of adjacent cell values to calculate transition. Calculate for all 8 adjacent cells
trans.nz
image(transitionMatrix(trans.nz))
plot(raster(trans.nz))

# correct the values - diagonal cells should be less connected than adjacent cells
trans.nzc <- geoCorrection(trans.nz, scl = TRUE)
trans.nzc
plot(raster(trans.nzc))


####################################################_
##### Calculate least cost paths between cells #####_
####################################################_

# Parallel loop for computing least cost paths. Set cell range (i) and buffer distance
system.time(
  lcps <- foreach(i = aoi) %dopar% { 
    costDistance(trans.nzc, SpatialPoints(cbind((coordinates(raster(trans.nzc))[i,1]),  
                                                coordinates(raster(trans.nzc))[i,2])), 
                 SpatialPoints(xyFromCell(raster(trans.nzc), raster::extract(raster(trans.nzc), 
                                                                             cbind((coordinates(raster(trans.nzc))[i,1]),  
                                                                                   coordinates(raster(trans.nzc))[i,2]),
                                                                             buffer = buff, cellnumbers = TRUE)[[1]][,1])))
  })


##############################################_
##### Format data output #####_
##############################################_

# For loop for re-naming distance object row and column names to correspond with raster 
#    cell numbers. Set cell range (i) and buffer distance to match definitions for lcps.
#    For rowname/colnames functions, set the n in '[[i - n]]' to the starting i minus 1
#    This function only works for consecutive sequences of i.
system.time(
  for (i in aoi) {
    rownames(lcps[[i - 0]]) <- i
    colnames(lcps[[i - 0]]) <- raster::extract(raster(trans.nzc), 
                                               cbind((coordinates(raster(trans.nzc))[i,1]),  
                                                     coordinates(raster(trans.nzc))[i,2]), buffer = buff, 
                                               cellnumbers = TRUE)[[1]][,1]
  })


# View first part of the distance matrix to visually confirm formatting is correct
as.matrix(lcps[[825]]) 

# Reformat data into an edge list 
#     remember that for all unmentioned edges the distance is infinate, not zero.
#     also need to remember to remove rows that are self connections (e.g. cell 1 to 1, dist = 0)
datalist = list() # create empty list

for (i in aoi) { # forloop to melt data into edge list
  dat <- melt(as.matrix(lcps[[i]]), varnames = c("From", "To"))
  datalist[[i]] <- dat # compile edge lists for each cell into list
}

lcps.edges <- do.call(rbind, datalist) # combine data into a single edge list
colnames(lcps.edges) <- c("From", "To", "LCP")

####################################################_
##### Calculate Euclidean distance between cells #####_
####################################################_

x <- coordinates(raster(import.file.name)) # cell number and utms (easting then northing)

euclid.dist <- as.matrix(dist(x, method = "euclidean", diag = FALSE, upper = FALSE, p = 2))
ed.melt <- melt(euclid.dist)
colnames(ed.melt) <- c("From", "To", "EucDist")
joined.edges <- left_join(lcps.edges, ed.melt, by = c("From","To")) # took approx 20 minutes


# Export edge list to csv file
write.csv(joined.edges, export.file.name, row.names = FALSE)

# End parallel computing
stopImplicitCluster()

#################################################_
#####
##### 4. Intralayer links - Hummingbird 2 #####
###
###
### About this section:
###   This section uses the gdistance package to calculate least cost paths for a given
###   species across a hypothetical landscape. The input is a species distribution model
###   in tif format in which each cell contains data on abundance or probability of occurrence.
###   The output is an edgelist reporting the least cost path between each pair of raster cells.
###   Took ~9 minutes to run in the RStudio jobs section using 7 cores.
###   Prior to running code, confirm working directory define parameters, and specify the 
###   import and export file names.
#################################################_

# Define parameters, and specify import/export file names
buff <- 500   # buffer width based on species movement or dispersal distance (1000 for H1, 500 for H2 and BB)
numcores <- 7  # number of cores to be used for parallel computing
aoi <- 1:10000  # Define range of raster cells to be used in connectivity calculations. Only
# works for a continuous range of numbers. If range does not start at 1, 
# see additional instructions in 'Format data output' section.
import.file.name <- "h2_distribution.tif" # file in tif format
export.file.name <- "h2_leastcost_edges.csv" # file in csv format


# Start parallel computation
registerDoParallel(numcores) 


#################################################_
##### Create and correct transition matrix #####_
#################################################_

# start with raster layer for the species distribution model
sdm <- raster(import.file.name)
plot(sdm)

# add small number to all cells so that there are no cells with prob = 0
nz <- sdm + 0.001 
plot(nz)

# create Transition matrix
trans.nz <- transition(nz, transitionFunction=mean, directions=8) # use mean of adjacent cell values to calculate transition. Calculate for all 8 adjacent cells
trans.nz
image(transitionMatrix(trans.nz))
plot(raster(trans.nz))

# correct the values - diagonal cells should be less connected than adjacent cells
trans.nzc <- geoCorrection(trans.nz, scl = TRUE)
trans.nzc
plot(raster(trans.nzc))


####################################################_
##### Calculate least cost paths between cells #####_
####################################################_

# Parallel loop for computing least cost paths. Set cell range (i) and buffer distance
system.time(
  lcps <- foreach(i = aoi) %dopar% { 
    costDistance(trans.nzc, SpatialPoints(cbind((coordinates(raster(trans.nzc))[i,1]),  
                                                coordinates(raster(trans.nzc))[i,2])), 
                 SpatialPoints(xyFromCell(raster(trans.nzc), raster::extract(raster(trans.nzc), 
                                                                             cbind((coordinates(raster(trans.nzc))[i,1]),  
                                                                                   coordinates(raster(trans.nzc))[i,2]),
                                                                             buffer = buff, cellnumbers = TRUE)[[1]][,1])))
  })


##############################################_
##### Format data output #####_
##############################################_

# For loop for re-naming distance object row and column names to correspond with raster 
#    cell numbers. Set cell range (i) and buffer distance to match definitions for lcps.
#    For rowname/colnames functions, set the n in '[[i - n]]' to the starting i minus 1
#    This function only works for consecutive sequences of i.
system.time(
  for (i in aoi) {
    rownames(lcps[[i - 0]]) <- i
    colnames(lcps[[i - 0]]) <- raster::extract(raster(trans.nzc), 
                                               cbind((coordinates(raster(trans.nzc))[i,1]),  
                                                     coordinates(raster(trans.nzc))[i,2]), buffer = buff, 
                                               cellnumbers = TRUE)[[1]][,1]
  })


# View first part of the distance matrix to visually confirm formatting is correct
as.matrix(lcps[[825]]) 

# Reformat data into an edge list 
#     remember that for all unmentioned edges the distance is infinate, not zero.
#     also need to remember to remove rows that are self connections (e.g. cell 1 to 1, dist = 0)
datalist = list() # create empty list

for (i in aoi) { # forloop to melt data into edge list
  dat <- melt(as.matrix(lcps[[i]]), varnames = c("From", "To"))
  datalist[[i]] <- dat # compile edge lists for each cell into list
}

lcps.edges <- do.call(rbind, datalist) # combine data into a single edge list
colnames(lcps.edges) <- c("From", "To", "LCP")

####################################################_
##### Calculate Euclidean distance between cells #####_
####################################################_

x <- coordinates(raster(import.file.name)) # cell number and utms (easting then northing)

euclid.dist <- as.matrix(dist(x, method = "euclidean", diag = FALSE, upper = FALSE, p = 2))
ed.melt <- melt(euclid.dist)
colnames(ed.melt) <- c("From", "To", "EucDist")
joined.edges <- left_join(lcps.edges, ed.melt, by = c("From","To")) # took approx 20 minutes


# Export edge list to csv file
write.csv(joined.edges, export.file.name, row.names = FALSE)

# End parallel computing
stopImplicitCluster()
#################################################_
#####
##### 5. Intralayer links - Bumblebee #####
###
###
### About this section:
###   This section uses the gdistance package to calculate least cost paths for a given
###   species across a hypothetical landscape. The input is a species distribution model
###   in tif format in which each cell contains data on abundance or probability of occurrence.
###   The output is an edgelist reporting the least cost path between each pair of raster cells.
###   Took ~9 minutes to run in the RStudio jobs section using 7 cores.
###   Prior to running code, confirm working directory define parameters, and specify the 
###   import and export file names.
#################################################_

# Define parameters, and specify import/export file names
buff <- 500   # buffer width based on species movement or dispersal distance (1000 for H1, 500 for H2 and BB)
numcores <- 7  # number of cores to be used for parallel computing
aoi <- 1:10000  # Define range of raster cells to be used in connectivity calculations. Only
# works for a continuous range of numbers. If range does not start at 1, 
# see additional instructions in 'Format data output' section.
import.file.name <- "bb_distribution.tif" # file in tif format
export.file.name <- "bb_leastcost_edges.csv" # file in csv format


# Start parallel computation
registerDoParallel(numcores) 


#################################################_
##### Create and correct transition matrix #####_
#################################################_

# start with raster layer for the species distribution model
sdm <- raster(import.file.name)
plot(sdm)

# add small number to all cells so that there are no cells with prob = 0
nz <- sdm + 0.001 
plot(nz)

# create Transition matrix
trans.nz <- transition(nz, transitionFunction=mean, directions=8) # use mean of adjacent cell values to calculate transition. Calculate for all 8 adjacent cells
trans.nz
image(transitionMatrix(trans.nz))
plot(raster(trans.nz))

# correct the values - diagonal cells should be less connected than adjacent cells
trans.nzc <- geoCorrection(trans.nz, scl = TRUE)
trans.nzc
plot(raster(trans.nzc))


####################################################_
##### Calculate least cost paths between cells #####_
####################################################_

# Parallel loop for computing least cost paths. Set cell range (i) and buffer distance
system.time(
  lcps <- foreach(i = aoi) %dopar% { 
    costDistance(trans.nzc, SpatialPoints(cbind((coordinates(raster(trans.nzc))[i,1]),  
                                                coordinates(raster(trans.nzc))[i,2])), 
                 SpatialPoints(xyFromCell(raster(trans.nzc), raster::extract(raster(trans.nzc), 
                                                                             cbind((coordinates(raster(trans.nzc))[i,1]),  
                                                                                   coordinates(raster(trans.nzc))[i,2]),
                                                                             buffer = buff, cellnumbers = TRUE)[[1]][,1])))
  })


##############################################_
##### Format data output #####_
##############################################_

# For loop for re-naming distance object row and column names to correspond with raster 
#    cell numbers. Set cell range (i) and buffer distance to match definitions for lcps.
#    For rowname/colnames functions, set the n in '[[i - n]]' to the starting i minus 1
#    This function only works for consecutive sequences of i.
system.time(
  for (i in aoi) {
    rownames(lcps[[i - 0]]) <- i
    colnames(lcps[[i - 0]]) <- raster::extract(raster(trans.nzc), 
                                               cbind((coordinates(raster(trans.nzc))[i,1]),  
                                                     coordinates(raster(trans.nzc))[i,2]), buffer = buff, 
                                               cellnumbers = TRUE)[[1]][,1]
  })


# View first part of the distance matrix to visually confirm formatting is correct
as.matrix(lcps[[825]]) 

# Reformat data into an edge list 
#     remember that for all unmentioned edges the distance is infinate, not zero.
#     also need to remember to remove rows that are self connections (e.g. cell 1 to 1, dist = 0)
datalist = list() # create empty list

for (i in aoi) { # forloop to melt data into edge list
  dat <- melt(as.matrix(lcps[[i]]), varnames = c("From", "To"))
  datalist[[i]] <- dat # compile edge lists for each cell into list
}

lcps.edges <- do.call(rbind, datalist) # combine data into a single edge list
colnames(lcps.edges) <- c("From", "To", "LCP")

####################################################_
##### Calculate Euclidean distance between cells #####_
####################################################_

x <- coordinates(raster(import.file.name)) # cell number and utms (easting then northing)

euclid.dist <- as.matrix(dist(x, method = "euclidean", diag = FALSE, upper = FALSE, p = 2))
ed.melt <- melt(euclid.dist)
colnames(ed.melt) <- c("From", "To", "EucDist")
joined.edges <- left_join(lcps.edges, ed.melt, by = c("From","To")) # took approx 20 minutes


# Export edge list to csv file
write.csv(joined.edges, export.file.name, row.names = FALSE)

# End parallel computing
stopImplicitCluster()

#################################################_
#####
##### 6. Interlayer Links #####
###
###
### About this section:
###   This section uses species distributions and simulated interaction preferences to infer
###   interlayer edges (plant-pollinator interactions), and exports the edgelist as a csv file.
#################################################_


# Create global interaction matrix
glob.net <- matrix(c(50,0,100,25,100,60,0,40,0,80,0,80), nrow = 3, byrow = TRUE)
colnames(glob.net) <- c("p1", "p2", "p3", "p4")
rownames(glob.net) <- c("h1", "h2", "bb")

# Define function for importing prob of occurrence matrices and converting prob < 0.1 to prob = 0
import.and.constrain <- function(filename) {
  ifelse(as.matrix(raster(filename)) < 0.1, 0, as.matrix(raster(filename)))
}

# Import and constrain occurrence matrices
h1.10 <- import.and.constrain("h1_distribution.tif")
h2.10 <- import.and.constrain("h2_distribution.tif")
bb.10 <- import.and.constrain("bb_distribution.tif")
p1.10 <- import.and.constrain("p1_distribution.tif")
p2.10 <- import.and.constrain("p2_distribution.tif")
p3.10 <- import.and.constrain("p3_distribution.tif")
p4.10 <- import.and.constrain("p4_distribution.tif")

# Define function for calculating interlayer edges
# pollinator = pollinator occurrence matrix
# plant = plant occurrence matrix
# gnr = pollinator row in glob.net
# gnc = plant column in glob.net
# S.layer = number of source layer
# D.layer = number of destination layer
inter.edge <- function(pollinator, plant, gnr, gnc, S.layer, D.layer) {
  co.occur <- pollinator * plant # probability of co-occurrence
  co.occur <- co.occur * glob.net[gnr,gnc]/sum(glob.net[gnr,]) # scale by foraging preferences in global network, define row and column numbers
  co.occur.melt <- melt(t(co.occur))[ ,3] # melt matrix
  df <- data.frame("S.node" = c(which(co.occur.melt > 0)),"S.layer" = S.layer, # create edgelist
                   "D.node" = c(which(co.occur.melt > 0)),"D.layer" = D.layer)
  df$weight <- subset(co.occur.melt, co.occur.melt > 0) # add weights to edgelist
  df
}

# h1 interlayer
h1.p1 <- inter.edge(h1.10, p1.10, 1, 1, 1, 4)
h1.p3 <- inter.edge(h1.10, p3.10, 1, 3, 1, 6)
h1.p4 <- inter.edge(h1.10, p4.10, 1, 4, 1, 7)

# h2 interlayer
h2.p1 <- inter.edge(h2.10, p1.10, 2, 1, 2, 4)
h2.p2 <- inter.edge(h2.10, p2.10, 2, 2, 2, 5) 
h2.p4 <- inter.edge(h2.10, p4.10, 2, 4, 2, 7)

# bb interlayer
bb.p2 <- inter.edge(bb.10, p2.10, 3, 2, 3, 5) 
bb.p4 <- inter.edge(bb.10, p4.10, 3, 4, 3, 7)

# save dataframe of interlayer links
interlayer3 <- rbind(h1.p1, h1.p3, h1.p4, h2.p1, h2.p2, h2.p4, bb.p2, bb.p4)

# Export interlayer edgelist
write.csv(interlayer3, "interlayer_edges.csv", row.names = FALSE)



#####
##### 7. Edge List Creation #####
###
###
### About this section:
###   This section imports the interlayer and least cost edgelists calculated 
###   previously, defines Kernal Utilization Distributions (homeranges) for each
###   pollinator species, and uses these to scale the least cost links. The final
###   product is a dataframe/csv file of the supragraph edges.
#################################################_



# Import edge list parts
interlayer3<- read.csv("interlayer_edges.csv")
h1.intra <- read.csv("h1_leastcost_edges.csv")
h2.intra <- read.csv("h2_leastcost_edges.csv")
bb.intra <- read.csv("bb_leastcost_edges.csv")


############################################################_
##### Define pollinator Kernal Utilization Distributions #####_
############################################################_

# Fit a gamma distribution to pilot hummingbird movement data (TN McFadden, unpublished data. More info available upon request). 
# Here we provide the relevant parameters derived from the hummingbird movement data.
plot(pgamma(0:500, shape = 1.4106, rate = 0.0304), xlab = "Distance from homerange center", 
     ylab = "Cumulative proportion of time spent")
plot(dgamma(0:500, shape = 1.4106, rate = 0.0304), xlab = "Distance from homerange center", 
     ylab = "Proportion of time spent")

# What proportion of time is spent beyond r radius? For example, 30m
1 - pgamma(0:500, shape = 1.4106, rate = 0.0304)[30] 


### Create hypothetical KUDs for each of the three species

# H1
plot(dgamma(0:1000, shape = 0.9, rate = 0.01), xlab = "Distance from homerange center", 
     ylab = "Proportion of time spent")
plot(pgamma(0:1000, shape = 0.9, rate = 0.01))
h1.pgamma <- pgamma(0:1000, shape = 0.9, rate = 0.01)
1 - h1.pgamma[1000]

# H2
plot(dgamma(0:500, shape = 1.4106, rate = 0.0304), xlab = "Distance from homerange center", 
     ylab = "Proportion of time spent")
plot(pgamma(0:500, shape = 1.4106, rate = 0.0304), xlab = "Distance from homerange center", 
     ylab = "Proportion of time spent")
h2.pgamma <- pgamma(0:500, shape = 1.4106, rate = 0.0304)

# BB
plot(dgamma(0:500, shape = 1.4106, rate = 0.02), xlab = "Distance from homerange center", 
     ylab = "Proportion of time spent")
plot(pgamma(0:500, shape = 1.4106, rate = 0.02), xlab = "Distance from homerange center", 
     ylab = "Proportion of time spent")
bb.pgamma <- pgamma(0:500, shape = 1.4106, rate = 0.02)


############################################################_
##### Convert intralayer distances to edge weights #####_
############################################################_

########  Notes on this section  ################################_
###
### Wij = f(EucDist) * 1/LCP, where f(EucDist) is the proportion of time a pollinator spends beyond
###               the given radius (radius == the Euclidean distance between point i and j). 
###
##########################################################################_

# Write a function for calculating Weight
Wij <- function(pgamma, EucDist, LCP) {
  (1 - pgamma[EucDist]) * (1/LCP)
}

# Format interlayer dataframes
h1.intra <- subset(h1.intra, h1.intra$EucDist > 0)
h2.intra <- subset(h2.intra, h2.intra$EucDist > 0)
bb.intra <- subset(bb.intra, bb.intra$EucDist > 0)

h1.intra$weights <- Wij(h1.pgamma, h1.intra$EucDist, h1.intra$LCP)
h2.intra$weights <- Wij(h2.pgamma, h2.intra$EucDist, h2.intra$LCP)
bb.intra$weights <- Wij(bb.pgamma, bb.intra$EucDist, bb.intra$LCP)

# Scale intralayer link weights so they aren't so small and so there is a lower cuttoff
#         I multiplied all weights by 10 and then set any weight < 0.01 to zero
h1.intra$weights <- h1.intra$weights * 10
h1.intra$weights[h1.intra$weights < 0.01] <- 0 

h2.intra$weights <- h2.intra$weights * 10
h2.intra$weights[h2.intra$weights < 0.01] <- 0 

bb.intra$weights <- bb.intra$weights * 10
bb.intra$weights[bb.intra$weights < 0.01] <- 0 


###############################################_
##### Combine edge lists to a single file #####_
###############################################_

# create edge list for each pollinator intralayer
edge.h1.h1 <- data.frame("S.node" = h1.intra$From,"S.layer" = 1,"D.node" = h1.intra$To,"D.layer" = 1, "weight" = h1.intra$weights)
edge.h2.h2 <- data.frame("S.node" = h2.intra$From,"S.layer" = 2,"D.node" = h2.intra$To,"D.layer" = 2, "weight" = h2.intra$weights)
edge.bb.bb <- data.frame("S.node" = bb.intra$From,"S.layer" = 3,"D.node" = bb.intra$To,"D.layer" = 3, "weight" = bb.intra$weights)

# save dateframe of intralayer links
intralayer <- rbind(edge.h1.h1, edge.h2.h2, edge.bb.bb) 

# combine intra and inter links to a single dataframe
path.multilayer <- rbind(intralayer, interlayer3) 

# remove zeros and NAs
path.multilayer <- subset(path.multilayer, path.multilayer$weight > 0) 

# Export as csv file
write.csv(path.multilayer, "edge_list.csv", row.names = FALSE)


#####
##### 8. Multilayer Connectance #####
###
###
### About this section:
###   First, creates an adjacency matrix from current network. Then creates a 'null'
###   adjacency matrix from the "null_edge_list.csv", which is a null network in which
###   all species occur universally across the landscape. We will use this as our 
###   'null model' for evaluating connectance, since this null network will contain 
###   the maximum possible number of links, given the biological constraints of 
###   species-specific home range KUDs, and species-specific floral preferences.
###
#################################################_


### Import and tidy data
edge.list <- read.csv("edge_list.csv")
# Add column specifiying intra vs inter-layer. Number will need to be changed depending on n pollinators in network
edge.list <- mutate(edge.list, type = ifelse(D.layer < 4, "intralayer", "interlayer")) 
edge.list$S.layer <- chartr("123456789", "ABCDEFGHI", edge.list$S.layer) # use letters instead of numbers to label layers
edge.list$D.layer <- chartr("123456789", "ABCDEFGHI", edge.list$D.layer)

# Combine node and layer labels to a single variable describing the node-layer tuple
edge.list <- edge.list %>% 
  unite(source_tuple, S.layer, S.node)

edge.list <- edge.list %>% 
  unite(dest_tuple, D.layer, D.node)

# Convert edgelist to a adjacency matrix
el.matrix <- as.matrix(edge.list)
g <- graph.edgelist(el.matrix[,1:2], directed = FALSE)
E(g)$weight = as.numeric(el.matrix[,3])
adj <- get.adjacency(g, type = c("both"), attr='weight', sparse=TRUE) 

rowSums(adj) # check that rowSums and colSums are the same
colSums(adj)



### Import and tidy null edge list
null1 <- read.csv("null_edges_1.csv")
null2 <- read.csv("null_edges_2.csv")
null3 <- read.csv("null_edges_3.csv")
edge.list <- rbind(null1, null2, null3)
# Add column specifiying intra vs inter-layer. Number will need to be changed depending on n pollinators in network
edge.list <- mutate(edge.list, type = ifelse(D.layer < 4, "intralayer", "interlayer")) 
edge.list$S.layer <- chartr("123456789", "ABCDEFGHI", edge.list$S.layer) # use letters instead of numbers to label layers
edge.list$D.layer <- chartr("123456789", "ABCDEFGHI", edge.list$D.layer)

# Combine node and layer labels to a single variable describing the node-layer tuple
edge.list <- edge.list %>% 
  unite(source_tuple, S.layer, S.node)

edge.list <- edge.list %>% 
  unite(dest_tuple, D.layer, D.node)

# Convert edgelist to a adjacency matrix
el.matrix <- as.matrix(edge.list)
g <- graph.edgelist(el.matrix[,1:2], directed = FALSE)
E(g)$weight = as.numeric(el.matrix[,3])
null.adj <- get.adjacency(g, type = c("both"), attr='weight', sparse=TRUE) 

rowSums(null.adj) # check that rowSums and colSums are the same
colSums(null.adj)



### Calculate connectance

# Weighted connectance (sum of all link weights in real network divided by sum 
#     of all link weights in null network)
weighted.connectance <- (sum(adj)/2) / (sum(null.adj)/2)   # must divide by two since adj is a symmetrical matrix
sum.weights <- sum(adj)/2

# Binary connectance (number of links in real network divided by number of links
#     in null network)
binary.connectance <- (length(which(adj > 0))/2) / (length(which(null.adj > 0))/2)
n.links <- length(which(adj > 0))/2


### Calculate number of connected plant and pollinator nodes

# Oder rows and columns by names, so that nodes are in order of layer (A,B,C...)
adj <- adj[mixedorder(rownames(adj)), mixedorder(colnames(adj))]

# Create separate matrix for intra-layer interactions
intra.edges <- filter(edge.list, type == "intralayer")
intra.edges <- intra.edges[ , -4]

el.matrix <- as.matrix(intra.edges)
g <- graph.edgelist(el.matrix[,1:2], directed = FALSE)
E(g)$weight = as.numeric(el.matrix[,3])
intra.adj <- get.adjacency(g, type = c("both"), attr='weight', sparse=TRUE) 

head(adj) # full adjacency matrix (supra-graph), with rows and columns sorted alphanumerically
head(intra.adj) # intralayer adjacency matrix, unsorted

starting.pollinators <- length(intersect(names(which(rowSums(adj[1:ncol(intra.adj), 1:ncol(intra.adj)]) > 0)), 
                                         names(which(rowSums(adj[1:ncol(intra.adj), (1+ncol(intra.adj)):ncol(adj)]) > 0))))
starting.plants <- sum(colSums(adj[1:ncol(intra.adj), (1+ncol(intra.adj)):ncol(adj)]) > 0)
starting.alive <- sum(starting.plants, starting.pollinators)

### Create results dataframe
results <- data.frame("Weighted_connectance" = weighted.connectance,
                      "Binary_connectance" = binary.connectance,
                      "n_pollinator_nodes" = starting.pollinators,
                      "n_plant_nodes" = starting.plants,
                      "n_nodes" = starting.alive,
                      "n_links" = n.links,
                      "sum_link_weights" = sum.weights)

write.csv(results, "connectance_results_table.csv")



