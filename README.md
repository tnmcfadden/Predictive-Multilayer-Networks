# Predictive-Multilayer-Networks
This repository contains data and R scripts accompanying McFadden and Dirzo (in prep.) Harnessing multilayer networks to predict metacommunity-wide responses to global environmental change. 

R scripts
-	RCode_1_Landscape_simulation_2019Nov15.R
-	RCode_2_Habitat_metrics_2019Nov15.R
-	RCode_3_Species_distributions_2019Nov15.R
-	RCode_4_Intralayer_links_2019Nov15.R
-	RCode_5_Interlayer_links_2019Nov15.R
-	RCode_6_Edge_list_creation_2019Nov15.R
-	RCode_7_Robustness_loop_2020April30
-	streamlined_pipeline_2020Jan10updated.R
-	RCode_9_Centrality_2020May03.R

Required data files
-	Land cover map: starting_landscape.tif (created in RCode_1_Landscape_simulation_2019Nov15.R)
-	Elevation map: elev_map_sim.tif

Instructions

Beginning with the RCode_2_Habitat_metrics script, run the scripts in order to create and analyze the multilayer network. By changing the land cover input file and re-running the entire pipeline in order, one can evaluate the effects of changing land-use scenarios on the robustness of the multilayer network.
