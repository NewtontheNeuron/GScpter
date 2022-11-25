# This is a script to act as a main proxy for theresa tam's analysis and to add
# any specialty functions of plots

args <- c("theresa_tam_mouse_SDH_DDH_Excit_Inhib",
          "../../Datasets/neurons_only_2021/clean_neuron_object.RDS")
args <- c("theresa_tam_human_DH_Excit_Inhib",
          "../../Datasets/human_spinalcord_2022/top_level_new_annotation.rda")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("loadLibraries.R")
source("JSON_Handler.R")
source("Pre_analysis_functions.R")
source("DotPlot.R")
source("PooledDotPlot.R")

source("ClusterPlot.R")

RDfile <- load_data(args[2])

# Mouse
extra_pool <- list()
extra_pool[["top"]] <- list("dataset", "age", "final_cluster_assignment", "run", "nCount_RNA")
extra_pool[["1"]] <- list("id", "features.label", "subgr")
extra_pool[["2"]] <- list("id", "features.label")
#extra_pool[["1"]] <- list("class_label", "region_label")

# Human
extra_pool <- list()
extra_pool[["top"]] <- list("nCount_RNA", "nCount_SCT", "batches_combined",
                            "HSC_Subcluster_Annotation", "new_annotation",
                            "top_level_annotation")
extra_pool[["1"]] <- list("subgr", "features.label")
extra_pool[["2"]] <- list("id", "features.label")

# Everything relies on all_cell_roster and cell_roster
all_cell_roster <- returnAllCellRoster(RDfile)
saveRDS(all_cell_roster, "../../Datasets/all_cell_rosters/all_cell_roster_human_theresa_tam.RDS")
all_cell_roster <- readRDS("../../Datasets/all_cell_rosters/all_cell_roster_human_theresa_tam.RDS")
all_cell_roster <- readRDS("../../Datasets/all_cell_rosters/all_cell_roster_mouse_theresa_tam.RDS")

# List by cluster
lbc <- createListbyCluster(scale.method = "log10")
mainDP(lbc)

#pooled dotplot.R
CPR <- createClusterPoolResults(scale.method = "z-score", pool.level = "2")
mainPDP(CPR)
