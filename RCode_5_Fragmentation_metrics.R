#################################################_
### Forest fragmentation
### Project: Predictive Multilayer Networks
### Author: Tyler N. McFadden
###
### Code accompanying manuscript: 
### McFadden TN, Dirzo R. Harnessing multilayer networks to predict 
### metacommunity responses to global environmental change.
### Submitted to Ecology Letters.
###
###
### About this script:
###   This script uses the landscapemetrics package (R implentation of FRAGSTATS) to calculate 
###   several standard metrics of forest fragmentation. 
###   The input is a map of land cover in tif format in which cell values 3 and 4
###   represent secondary and primary forest cover, respectively.
###   Prior to running code, confirm working directory define parameters, and specify the 
###   import and export file names.
#################################################_


# Define parameters, and specify import/export file names
# Set working directory
setwd("/Users/tylermcfadden/Desktop/Stanford/Multilayer_network/Manuscript/Code")

# Define landscape input file
landscape.file.name <- "starting_landscape.tif"

# Load libraries
library(landscapeR)
library(raster)
library(foreach)
library(doParallel)
library(gdistance)
library(reshape2)
library(dplyr)
library(landscapemetrics)

# Function for calculating mean, SE, n
meanSEn <- function(x) {
  mean <- mean(x, na.rm = TRUE)
  se <- sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))
  n <- length(x[!is.na(x)])
  c(mean,se,n)
}


#################################################_
##### Calculate landscape metics #####_
#################################################_

# start with raster layer for the landscape
landscape <- raster(landscape.file.name)
forest <- subs(landscape, data.frame(id=c(1,2,3,4), v=c(0,0,1,1))) # convert to map of forest cover
forest[is.na(forest[])] <- 0 # set any NA cells to 0
plot(forest)

# Measuring landscape fragmentation metrics via landscapemetrics package
check_landscape(forest) # check that dataset is properly formatted
list_lsm() # list of all possible metrics

forest.cover <- sum(values(forest))/10000

patch.area <- lsm_p_area(forest, directions = 8) # area is reported in ha
patch.area <- subset(patch.area, patch.area$class == 1)
largest.patch <- max(patch.area$value)
mean.patch <- meanSEn(patch.area$value)[1]
SE.patch <- meanSEn(patch.area$value)[2]
n.patches <- meanSEn(patch.area$value)[3]
patch.density <- n.patches/9 # patch density per 100 ha

nneighbor <- lsm_p_enn(forest, directions = 8, verbose = TRUE)
nneighbor <- subset(nneighbor, nneighbor$class == 1)
mean.nnd <- meanSEn(nneighbor$value)[1] # mean + SE nearest neighbor
SE.nnd <- meanSEn(nneighbor$value)[2]

perim <- lsm_p_perim(forest, directions = 8)
perim <- subset(perim, perim$class == 1)
edge.length <- sum(perim$value) # total edge length in m

edge.ratio <- lsm_p_para(forest, directions = 8)
edge.ratio <- subset(edge.ratio, edge.ratio$class == 1)
mean.edge.ratio <- meanSEn(edge.ratio$value)[1] # mean + SE edge ratio (perimeter/area in m/m2)
SE.edge.ratio <- meanSEn(edge.ratio$value)[2]

# create a dataframe to store results
metrics <- data.frame("Scenario" = landscape.file.name, forest.cover,largest.patch, mean.patch,SE.patch,
                                n.patches, patch.density,mean.nnd,SE.nnd,
                               edge.length,mean.edge.ratio,SE.edge.ratio)
 

write.csv(metrics, "fragmentation.metrics.csv")


