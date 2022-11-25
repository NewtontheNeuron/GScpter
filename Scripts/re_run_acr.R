# This is a script to rerun other all_cell_rosters

# The Grin gene analysis
args <- c("grin_only_SDH_DDH_Excit_Inhib", "../../Datasets/neurons_only_2021/clean_neuron_object.RDS")
args <- c("grin_only_human_DH_Excit_Inhib", "../../Datasets/human_spinalcord_2022/top_level_new_annotation.rda")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("loadLibraries.R")
source("JSON_Handler.R")
source("Pre_analysis_functions.R")
source("DotPlot.R")
source("PooledDotPlot.R")
source("ClusterPlot.R")

all_cell_roster <- readRDS("../../Datasets/all_cell_rosters/all_cell_roster_human_grin_only.RDS")
all_cell_roster <- readRDS("../../Datasets/all_cell_rosters/all_cell_roster_grin.RDS")

# Perhaps it is better to modify the load data function
# After this you can do what ever you want

# For grant remove this afterwards
acr <- all_cell_roster
all_cell_roster <- all_cell_roster %>%
  filter(id == "SDH")
mainDP(lbc)
CPR <- createClusterPoolResults(roster = acr, pool.level = "2", scale.method = "z-score")
library(forcats)
CPR <- CPR %>%
  mutate(id = fct_relevel(id, "SDH", "DDH"))
mainPDP(CPR)

lbc <- createListbyCluster("log10")
mainDP(lbc)
CPR <- createClusterPoolResults(pool.level = "2", scale.method = "log10")
mainPDP(CPR)
