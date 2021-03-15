#################################################_
### Betweenness Centrality
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
### This script uses the igraph package to calculate betweenness for all
### cells in the starting landscape. Betweenness can be visualized by plotting
### landscapes with cells colored by betweenness percentile. 
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
library(tidyverse)
library(igraph)
library(gtools)
library(Matrix)
library(viridis)

##################################_
##### Import and format data #####
##################################_

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

# Convert edgelist to an adjacency matrix
el.matrix <- as.matrix(edge.list)
g <- graph.edgelist(el.matrix[,1:2], directed = FALSE)
E(g)$weight = as.numeric(el.matrix[,3])
adj <- get.adjacency(g, type = c("both"), attr='weight', sparse=TRUE) 

# Oder rows and columns by names, so that nodes are in order of layer (A,B,C...)
adj <- adj[mixedorder(rownames(adj)), mixedorder(colnames(adj))]

# check that rowSums and colSums are the same
rowSums(adj)
colSums(adj)



##################################_
##### Centrality via igraph #####
##################################-


graph <- induced_subgraph(g, vids = V(g)[1:length(V(g))])

# Took 147 minutes to run whole network with true betweenness. weights as 1/cost
system.time(
  betweenness <- betweenness(graph, directed = FALSE, weights = 1/E(graph)$weight))

centrality <- as.data.frame(betweenness)
centrality$node <- rownames(centrality)
centrality <- centrality %>% separate(node, into = c("Layer", "Cell"))

ggplot(data = centrality) +
  geom_histogram(aes(betweenness, color = Layer), bins = 300) 


# sum betweenness for nodes occupying same cell
cell.centrality <- centrality %>% group_by(Cell) %>% summarise(sum(betweenness))
cell.centrality$Cell <- as.numeric(cell.centrality$Cell)
cell.centrality <- arrange(cell.centrality, Cell)
summary(cell.centrality$`sum(betweenness)`)


# Import landscape file
rs <- raster("starting_landscape.tif")

# Create raster with same dimensions but with betweenness values
r <- setValues(rs, cell.centrality$`sum(betweenness)`)
plot(r)

# Plot betweenness percentiles with better color scale
cuts <- quantile(cell.centrality$`sum(betweenness)`, c(0,0.5, 0.75, 0.9,.95,0.99, 1))
pal <- scales::viridis_pal(option = "A", direction = 1, begin = 0, end = 1)(6)
par(mar = c(2, 2, 2, 2))
par(mfrow = c(1,1))
par(xpd = NA, tck = NA)
plot(r, breaks=cuts, col = pal, legend = FALSE, box = FALSE, axes = FALSE)
par(xpd = TRUE)
legend(x = 659000, y = 5591000, legend = c("0-50%", "50-75%", "75-90%", "90-95%", "95-99%", "99-100%"), 
       fill = pal,
       cex = 1, bty = "n")

# Plot landscape's land cover and bumblebee distributions for comparison
par(mar = c(0.5, 0.5, 0.5, 0.5))
par(mfrow = c(2,1))
par(xpd = FALSE)
plot(rs,legend = FALSE, col = rev(terrain.colors(6)), axes = FALSE, box = FALSE)
par(xpd = TRUE)
legend(x = 659000, y = 5591000, legend = c("Primary forest", "Secondary forest", "Crops", "Pasture"), 
       fill = c("#00A600FF", "#63C600FF", "#E6E600FF", "#EEB99FFF"),
       cex = 0.9, bty = "n")
plot(raster("bb_distribution.tif"), legend = TRUE, box = FALSE, axes = FALSE)


# Export betweenness data to csv files
write.csv(centrality, "betweenness.csv")
write.csv(cell.centrality, "betweenness_by_cell.csv")




