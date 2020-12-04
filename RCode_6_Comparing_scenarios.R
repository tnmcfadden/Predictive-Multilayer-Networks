############################################_
### Comparing network metrics across scenarios
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
###   Provides plots and linear models comparing PMN connectance and robustness across scenarios. 
###
############################################_

# set working directory
setwd("/Users/tylermcfadden/Desktop/Stanford/Multilayer_network/Manuscript/Code")

# load libraries
library(raster)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
library(viridis)

# load data
all.results <- read.csv("results_all_scenarios.csv")


# Plot weighted connectance vs. change in forest cover
conn <- ggplot(data = all.results) +
  geom_point(aes(x = percent.change, y = Weighted_connectedness, color = ""), size = 6) +
  ylim(0,1) + 
  xlim(-100,100) +
  ylab("Connectance") + 
  xlab("Percent change in forest cover vs. initial scenario") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25)) +
  scale_color_viridis(discrete = TRUE, option = "D", direction = 1) +
  theme(legend.position = "none")
conn 


# Plot pollinator robustness vs. change in forest cover
poll <- ggplot(data = all.results) +
  geom_point(aes(x = percent.change, y = poll.robust, color = as.factor(Threshold)), alpha = 0.3, size = 6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylim(0,1) + 
  xlim(-100,100) +
  ylab("Pollinator robustness") + 
  xlab(element_blank()) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25)) +
  scale_color_viridis(discrete = TRUE, option = "D", direction = -1, 
                      name = "Intralyer extinction\nthreshold", 
                      labels = c("High", "Medium", "Low")) +
  guides(colour = guide_legend(override.aes = list(alpha=1))) +
  theme(legend.position = c(0.90, 0.25)) +
  theme(legend.background = element_rect(fill="white")) +
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 15))
poll  


# Plot plant robustness vs. change in forest cover
plant <- ggplot(data = all.results) +
  geom_point(aes(x = percent.change, y = plant.robust, color = as.factor(Threshold)), alpha = 0.3, size =6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylim(0,1) +
  xlim(-100,100) +
  ylab("Plant robustness") +
  xlab("Percent change in forest cover vs. initial scenario") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25)) +
  scale_color_viridis(discrete = TRUE, option = "D", direction = -1) +
  theme(legend.position = "none")
plant


# Linear models of robustness vs change in forest cover
summary(lm(all.results$poll.robust[all.results$Threshold == "1"] ~ all.results$percent.change[all.results$Threshold == "1"]))
summary(lm(all.results$poll.robust[all.results$Threshold == "2.5"] ~ all.results$percent.change[all.results$Threshold == "2.5"]))
summary(lm(all.results$poll.robust[all.results$Threshold == "5"] ~ all.results$percent.change[all.results$Threshold == "5"]))

summary(lm(all.results$plant.robust[all.results$Threshold == "1"] ~ all.results$percent.change[all.results$Threshold == "1"]))
summary(lm(all.results$plant.robust[all.results$Threshold == "2.5"] ~ all.results$percent.change[all.results$Threshold == "2.5"]))
summary(lm(all.results$plant.robust[all.results$Threshold == "5"] ~ all.results$percent.change[all.results$Threshold == "5"]))

