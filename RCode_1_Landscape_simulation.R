#################################################_
### Landscape Simulation and Scenario Building
### Project: Predictive Multilayer Networks
### Author: Tyler N. McFadden
###
### Code accompanying manuscript: 
### McFadden TN, Dirzo R. Harnessing multilayer networks to predict 
### metacommunity responses to global environmental change.
### Submitted to Ecological Monographs.
###
### About this script:
###   This script creates a simulated fragmented landscape containing four 
###   land cover types and contains code for how to create a simulated elevation profile, though
###   the actual elevation profile used is provided as a separate file. Code is also
###   provided for simulating the 20 deforestation and reforestation scenarios.
#################################################_

# Set working directory
setwd("/Users/tylermcfadden/Desktop/Stanford/Multilayer_network/Manuscript/Code")

# Load libraries
library(landscapeR)
library(raster)
library(sp)
library(rgdal)
library(geoR)

# Use random number generation methods from old version of R (pre R 3.6.0)
RNGkind(sample.kind = "Rounding")


#################################################_
##### Create a simulated land cover map #####
#################################################_

# Set the basic raster parameters - 30m cells, using UTM coordinate system
r <- raster(nrows=100, ncols=100, xmn=656000, xmx=659000, ymn=5588000, 
            ymx=5591000, resolution=30, crs = CRS("+init=epsg:5362"))
r <- setValues(r, 0)
summary(r)

# Assign land cover (1 = pasture, 2= cropland, 3 = secondary forest, 4 = primary forest)
set.seed(123)
pasture <- makeClass(r, 70, abs(rnorm(70, mean = 100, sd = 100)), val = 1)
set.seed(123)
crops <- makeClass(pasture, 20, abs(rnorm(20, mean = 20, sd = 10)), val = 2, bgr = c(0,1))
set.seed(123)
secondary <- makeClass(crops, 100, abs(rnorm(100, mean = 120, sd = 100)), val = 3)
set.seed(123)
primary <- makeClass(secondary, 5, abs(rnorm(5, mean = 200, sd = 100)), val = 4, bgr = c(0,1,2,3))
rs <- rmSingle(primary)

plot(rs)

# Export raster
writeRaster(rs, "starting_landscape.tif", format = "GTiff")



########################################
### Simulate Elevation across landscape
########################################

#### This code was used to create a realistic map of elevation. However, it is based on random
# number gereration, so the exact map is not reproducible here. 
sim <- grf(10000, grid="reg", cov.pars=c(1, 0.25))
image(sim, col=gray(seq(1, .1, l=30)))

elev <- raster(sim)
elev <- setValues(elev, 100*(values(elev) + 2.5))
plot(elev)

elev.utm <- raster(nrows=100, ncols=100, xmn=656000, xmx=659000, ymn=5588000, 
                   ymx=5591000, resolution=30, crs = CRS("+init=epsg:5362"))
elev.matrix <- as.matrix(values(elev))
elev.utm <- setValues(r, elev.matrix)
plot(elev.utm)

# The actual elevation map used in analysis is available here:
elev.raster <- raster("elev_map_sim.tif") # import elevation raster
plot(elev.raster)


#################################################_
##### Create alternative land cover scenarios #####
#################################################_

# View original landscape
plot(rs)

# Count the number of forest cells on landscape (primary or secondary)
vals <- getValues(rs)
nfcells <- length(subset(vals, vals == 4 | vals == 3))

### Deforestation scenarios (primary and secondary forest converted to pasture)

# Create a function for converting forest to pasture
forest.loss <- function(orig.landscape, prop.change){
  set.seed(123)
  expandClass(orig.landscape, 1, (prop.change * nfcells), bgr = c(3,4))
}

# Create rasters for each deforestation scenario. Defoestation occured in 10% increments
writeRaster(forest.loss(rs, 0.1), "scen.loss.10.tif", format = "GTiff")
writeRaster(forest.loss(rs, 0.2), "scen.loss.20.tif", format = "GTiff")
writeRaster(forest.loss(rs, 0.3), "scen.loss.30.tif", format = "GTiff")
writeRaster(forest.loss(rs, 0.4), "scen.loss.40.tif", format = "GTiff")
writeRaster(forest.loss(rs, 0.5), "scen.loss.50.tif", format = "GTiff")
writeRaster(forest.loss(rs, 0.6), "scen.loss.60.tif", format = "GTiff")
writeRaster(forest.loss(rs, 0.7), "scen.loss.70.tif", format = "GTiff")
writeRaster(forest.loss(rs, 0.8), "scen.loss.80.tif", format = "GTiff")
writeRaster(forest.loss(rs, 0.9), "scen.loss.90.tif", format = "GTiff")
writeRaster(forest.loss(rs, 1), "scen.loss.100.tif", format = "GTiff")


### Forest restoration scenarios (pasture converted to secondary forest)

# Create a function for converting pasture to forest
forest.rest <- function(orig.landscape, prop.change){
  set.seed(123)
  expandClass(orig.landscape, 3, (prop.change * nfcells), bgr = c(1))
}

# Create rasters for each reforestation scenario. Refoestation occured in 10% increments
writeRaster(forest.rest(rs, 0.1), "scen.rest.10.tif", format = "GTiff")
writeRaster(forest.rest(rs, 0.2), "scen.rest.20.tif", format = "GTiff")
writeRaster(forest.rest(rs, 0.3), "scen.rest.30.tif", format = "GTiff")
writeRaster(forest.rest(rs, 0.4), "scen.rest.40.tif", format = "GTiff")
writeRaster(forest.rest(rs, 0.5), "scen.rest.50.tif", format = "GTiff")
writeRaster(forest.rest(rs, 0.6), "scen.rest.60.tif", format = "GTiff")
writeRaster(forest.rest(rs, 0.7), "scen.rest.70.tif", format = "GTiff")
writeRaster(forest.rest(rs, 0.8), "scen.rest.80.tif", format = "GTiff")
writeRaster(forest.rest(rs, 0.9), "scen.rest.90.tif", format = "GTiff")
writeRaster(forest.rest(rs, 1), "scen.rest.100.tif", format = "GTiff")





