# Predictive-Multilayer-Networks
This repository contains data and R scripts accompanying McFadden and Dirzo (submitted to Nature Ecology & Evolution) Harnessing multilayer networks to predict metacommunity-wide responses to global environmental change. 

R scripts
- RCode_1_Landscape_simulation.R
- RCode_2_PMN_creation.R
- RCode_3_Betweenness.R
- RCode_4_Robustness.R
- RCode_5_Fragmentation_metrics.R
- RCode_6_Comparing_scenarios.R


Required data files for PMN creation
-	Land cover map: starting_landscape.tif (created in RCode_1_Landscape_simulation.R)
-	Elevation map: elev_map_sim.tif
- Null network for connectance metric: null_edge_list.csv


Example results folder
- 



Instructions

Beginning with the RCode_2_Habitat_metrics script, run the scripts in order to create and analyze the multilayer network. By changing the land cover input file and re-running the entire pipeline in order, one can evaluate the effects of changing land-use scenarios on the robustness of the multilayer network.
