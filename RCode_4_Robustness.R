#################################################_
### Multilayer Robustness
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
###   This script imports and formats the PMN edgelist and implements an extinction cascade 
###   algorithm to measure network robustness. The algorithm removes one pollinator 
###   node at a time and counts the number of surviving plant and pollinator nodes left
###   in the network. A plant goes secondarily extinct if it loses all of its pollinators 
###   (total interlayer link weights = 0). A pollinator goes secondarily extinct if it 
###   loses all of its plants (total interlayer link weights = 0) or if its spatial 
###   connections are reduced below a certain threshold (total intralayer link weights 
###   < threshold). The intralayer threshold is defined by the user as the intralayer 
###   links weights of the lowest X percentile of pollinator nodes in the ‘starting_landscape’ 
###   network. The order of pollinator node removal is determined as a function of the pollinator’s 
###   probability of occurrence at that location. Conceptually, the entire process takes 
###   place as follows:
###     1.	Pollinator node randomly selected
###          -	Probability of selected node being removed is a function of occurrence probability.
###     2.	Pollinator node is removed 
###     3.	Pollinator nodes with intralayer link weights below threshold go secondarily extinct and are removed
###     4.	Plants and pollinators with interlayer equal to zero go secondarily extinct and are removed
###     5.	Number of surviving plant and pollinator nodes are counted.
###     6.	Repeat steps 1-5 until no surviving nodes remain
###
###   The script is arranged as a series of nested loops, so that the user can run the
###   entire script one time per scenario. The outmost loop (j; line 133, section = "Nested loop for running multiple replicates") 
###   indicates the number of replicates (currently 5). Each replicate uses a different node 
###   removal order. The next loop (k; line 163, section = "Multialyer robustness") indicates which 
###   intralayer extinction threshold values to use. Currently it uses 3 different values 
###   for each replicate (1st, 2.5th, and 5th percentile).
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

# Define lengths of each subnetwork
adj.length <- ncol(adj)

intra.length <- sum(length(which(startsWith(colnames(adj), "A"))),
                    length(which(startsWith(colnames(adj), "B"))),
                    length(which(startsWith(colnames(adj), "C"))))

inter.length <- sum(length(which(startsWith(colnames(adj), "D"))),
                    length(which(startsWith(colnames(adj), "E"))),
                    length(which(startsWith(colnames(adj), "F"))),
                    length(which(startsWith(colnames(adj), "G"))))


####################################_
##### Removal order options  #####
####################################_

# Simple removal orders (not used here)
#removals <- sample(colnames(adj)[1:intra.length]) # random order of pollinator nodes
#removals <- names(sort(colSums(adj)[1:intra.length], decreasing = TRUE)) # ordered removal of pollinator nodes from most to least connected


### Create removal order based on occurrence probability

# create dataframe of tuple name and tuple occurrence probability
h1.values <- data.frame("layer" = "A", "cell" = 1:10000, "p.occur" = values(raster("h1_distribution.tif")))
h2.values <- data.frame("layer" = "B", "cell" = 1:10000, "p.occur" = values(raster("h2_distribution.tif")))
bb.values <- data.frame("layer" = "C", "cell" = 1:10000, "p.occur" = values(raster("bb_distribution.tif")))
p1.values <- data.frame("layer" = "D", "cell" = 1:10000, "p.occur" = values(raster("p1_distribution.tif")))
p2.values <- data.frame("layer" = "E", "cell" = 1:10000, "p.occur" = values(raster("p2_distribution.tif")))
p3.values <- data.frame("layer" = "F", "cell" = 1:10000, "p.occur" = values(raster("p3_distribution.tif")))
p4.values <- data.frame("layer" = "G", "cell" = 1:10000, "p.occur" = values(raster("p4_distribution.tif")))
node.p.occur <- rbind(h1.values, h2.values, bb.values, p1.values, p2.values, p3.values, p4.values)
node.p.occur <- node.p.occur %>% 
  unite(source_tuple, layer, cell)
node.names <- data.frame("source_tuple" = colnames(adj))
tuple.occur <- left_join(node.names, node.p.occur, by = c("source_tuple" = "source_tuple"))


#########################################################_
##### Nested loop for running multiple replicates #####
#########################################################_

# Define results matrix
robust.results <- matrix(NA, nrow = 0, ncol = 4) 
for (j in 1:5) {
  
# This loop to selects extinction order based on occurrence probability. 
#     Takes about a minute to run for whole dataset

selections <- tuple.occur$source_tuple # define 'selections' list
selections <- selections[1:intra.length] # shorten to only pollinator nodes
removals <- matrix(data = NA, nrow = ncol(adj[,1:intra.length]), ncol = 1) # create empty matrix for results (this will be the final removal order)
system.time(
  for (i in 1:intra.length) {
    repeat {
      node <- sample(selections, 1) # select a node
      if(rbinom(n = 1, size = 1, prob =    # will the selected node go primarily extinct?
                (1 - tuple.occur$p.occur[tuple.occur$source_tuple == node])) == 1) {
        break
      } # else print(node) # printing the nodes that don't going extinct is optional, but useful if there is a bug
    }
    removals[i] <- node; # if yes, add node to removals list
    selections <- selections[-which(selections == node)] # and delete from starting selections list
  })

length(removals)
length(unique(removals)) # check that all nodes in removals list are unique



#####################################################_
##### Multilayer robustness #####
#####################################################_

for (k in 1:3) {

# Define removal matrix as upper rectangle of edgelist matrix
removal.matrix <- adj[1:intra.length, ]
tail(rownames(removal.matrix)) # check that matrix is what I expect
tail(colnames(removal.matrix))

# Decide on intralayer weight threshold for pollinator secondary extinctions
        # currently the loop uses 3 values, the 1st, 2.5th, and 5th percentiles of the intralayer weights in original'starting landscape' network
        # Alternative values based on percentiles can be derived using this code (the number 0.01 indicates we are finding the value of the 1st percentile): 
        # quantile(colSums(removal.matrix[,1:intra.length]), 0.01)
threshold <- c(0.1396736, 3.271279, 11.73237)[k]

# Remove rows and columns of pollinator nodes which do not have interlayer edges
removal.matrix <- removal.matrix[-which(rowSums(removal.matrix[, (1+intra.length):adj.length]) == 0),
                                 -which(rowSums(removal.matrix[, (1+intra.length):adj.length]) == 0)]

# Make sure everything looks good
head(colnames(removal.matrix))
length(rownames(removal.matrix))
length(colnames(removal.matrix))

# Remove rows and columns of pollinators nodes with total intralayer edges < extinction threshold (lowest X% of nodes in starting landscape network)
removal.matrix <- removal.matrix[!rownames(removal.matrix) %in% 
                 intersect(rownames(removal.matrix), 
                           names(which(colSums(removal.matrix) < threshold))), 
               !colnames(removal.matrix) %in% 
                 intersect(rownames(removal.matrix), 
                           names(which(colSums(removal.matrix) < threshold)))]

length(rownames(removal.matrix))
length(colnames(removal.matrix))

# Check removal order
removals[1:100]

# Create matrix for storing results from extinction loop
results.mat <- matrix(NA, nrow = length(removals), ncol = 5) 

# Extinction loop
system.time(
for (i in 1:length(removals)) {
    if(length(intersect(removals[i], rownames(removal.matrix))) != 1) {
      next
    } 
  # count starting columns and rows
  start.ncol <- ncol(removal.matrix)
  start.nrow <- nrow(removal.matrix)
  if(nrow(removal.matrix) < 3) {break}  # break this loop if removal matrix has 2 rows, otherwise it will throw an error. 
  # save name of removed node
  node <- removals[i]
  # remove row and column of selected node (primary extinction)
  removal.matrix <- removal.matrix[!rownames(removal.matrix) %in% removals[i], !colnames(removal.matrix) %in% removals[i]]
  if(is.integer(nrow(removal.matrix))) {} else break 
  # remove any pollinator rows and columns with colsums < extinction threshold (secondary extinctions of pollinators due to loss of intralayer links)
  removal.matrix <- removal.matrix[!rownames(removal.matrix) %in% 
                                     intersect(rownames(removal.matrix), 
                                          names(which(colSums(removal.matrix) < threshold))), 
                                   !colnames(removal.matrix) %in% 
                                     intersect(rownames(removal.matrix), 
                                               names(which(colSums(removal.matrix) < threshold)))]
  if(is.integer(nrow(removal.matrix))) {} else break 
  # remove any plant columns with colsums = zero (secondary extinctions of plants due to loss of intralayer links)
  removal.matrix <- removal.matrix[ , !colnames(removal.matrix) %in% 
                                      intersect(colnames(removal.matrix), 
                                          names(which(colSums(removal.matrix) == 0)))]
  if(is.integer(nrow(removal.matrix))) {} else break 
  # count ending columns and rows
  end.ncol <- ncol(removal.matrix)
  end.nrow <- nrow(removal.matrix)
  # save results to matrix
  results.mat[i, 1] <- node
  results.mat[i, 2] <- start.ncol
  results.mat[i, 3] <- start.nrow
  results.mat[i, 4] <- end.ncol
  results.mat[i, 5] <- end.nrow
}) # end of extinction (i) loop

length(unique(results.mat[,1]))



# Tidy results and add columns
results.mat <- results.mat[rowSums(is.na(results.mat)) != ncol(results.mat), ]
results <- as.data.frame(results.mat[,2:5])
results$V1 <- as.numeric(as.character(results$V1))
results$V2 <- as.numeric(as.character(results$V2))
results$V3 <- as.numeric(as.character(results$V3))
results$V4 <- as.numeric(as.character(results$V4))
colnames(results) <- c("starting.ncol", "starting.nrow", "ending.ncol", "ending.nrow")
results <- rbind(results, c(results$ending.ncol[nrow(results.mat)], results$ending.nrow[nrow(results.mat)], 0, 0))
results$starting.poll <- results$starting.nrow
results$starting.plant <- results$starting.ncol - results$starting.nrow
results$ending.poll <- results$ending.nrow
results$ending.plant <- results$ending.ncol - results$ending.nrow
results$n.poll.extinct <- cumsum(results$starting.poll - results$ending.poll)
results$n.plant.extinct <- cumsum(results$starting.plant - results$ending.plant)
results$n.primary.extinct <- as.numeric(rownames(results))
results$sec.poll.extinct <- results$n.poll.extinct - results$n.primary.extinct



# Plot results

# secondary extinctions as a function of primary extinctions
plot1 <- ggplot(data = results) +
  geom_point(aes(x = n.primary.extinct/starting.poll[1], y = sec.poll.extinct/starting.poll[1]), color = "blue", size = 0.03) +
  #geom_point(aes(x = n.primary.extinct/starting.poll[1], y = n.poll.extinct/starting.poll[1]), color = "blue", size = 0.03) +
  geom_point(aes(x = n.primary.extinct/starting.poll[1], y = n.plant.extinct/starting.plant[1]), color = "green", size = 0.03) +
  xlab("Proportion of pollinator nodes removed") +
  ylab("Proportion nodes lost to secondary extinctions") +
  xlim(0,1)
  
ggsave(plot1, file = paste0("plot1_", "thresh_", threshold, "_rep_", j,".png"), width = 14, height = 10, units = "cm")

# surviving nodes as a function of nodes removed
plot2 <- ggplot(data = results) +
  geom_point(aes(x = n.primary.extinct/starting.poll[1], y = ending.poll/starting.poll[1], color = "blue"), size = 1.5) +
  geom_segment(aes(x = 0, xend = 1, y = 1, yend = 0, color = "blue"), linetype = "dashed", size = 1.1) +
  geom_point(aes(x = n.primary.extinct/starting.poll[1], y = ending.plant/starting.plant[1], color = "green"), size = 1.5) +
  xlab("Proportion of pollinator nodes removed") +
  ylab("Proportion of nodes surviving") +
  geom_segment(aes(x = 0, xend = 1, y = 1, yend = 1, color = "green"), linetype = "dashed", size = 1.1) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25)) +
  scale_color_viridis(discrete = TRUE, option = "D", begin = 0, end = 0.4, direction = 1, 
                      name = "", 
                      labels = c("Pollinators", "Plants")) +
  theme(legend.position = c(0.90, 0.80)) +
  theme(legend.background = element_rect(fill="white")) +
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 15)) +
  guides(colour = guide_legend(override.aes = list(size = 2)))

ggsave(plot2, file = paste0("plot2_", "thresh_", threshold, "_rep_", j,".png"), width = 14, height = 10, units = "cm")


# Calculate robustness metrics

# To single guild extictions
robustness.poll2poll <- sum(results$ending.poll) / sum(seq(from = results$starting.poll[1], to = 0, by = -1))
robustness.plant2poll <- sum(results$ending.plant/results$starting.plant[1]) / results$starting.poll[1]

# Save results
robust.results <- rbind(robust.results, 
                        c(threshold, j, robustness.poll2poll, robustness.plant2poll))
#robust.results[j, 1] <- threshold
#robust.results[j, 2] <- j
#robust.results[j, 3] <- robustness.poll2poll
#robust.results[j, 4] <- robustness.plant2poll

} # end of threshold (k) loop
} # end of replicate (j) loop

write.csv(robust.results, "robust.results.csv")



