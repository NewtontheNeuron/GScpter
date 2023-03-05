#TODO: need to take arguments for knowing what data file to use..

#something like this
#args <- c("NMDA_mouse_brain", "C:/Users/no/Documents/Neuroscience - MSc/Winter 2022/BIOL 5502/Rmd_proj/Seurat.ss.rda")
args <- c("human_nmda_ex_inh", "../../Datasets/human_ariel_data/top_level_new_annotation.rda")
args <- c("nmda_sdh_ddh", "NULL")
#args <- commandArgs(trailingOnly = TRUE)
print(paste ("this is the title project name entered: ", args[1]))
print(paste ("this was the given directory:", args[2]))

# The Grin gene analysis
args <- c("grin_only_SDH_DDH_Excit_Inhib", "../../Datasets/neurons_and_glia_2022/final_meta_dataset.rds")
args <- c("grin_only_human_DH_Excit_Inhib", "../../Datasets/human_spinalcord_2022/top_level_new_annotation.rda")
# Unit testing group over
args <- c("grin_only_SDH_DDH_Excit_Inhib_for_unit_testing_groupover", "../../Datasets/neurons_and_glia_2022/final_meta_dataset.rds")

### for kcni
args <- c("combined", "~/Neuroscience - MSc/Summer 2022/kcni_summer_school/Data/mouse_human_practice.rds")

### for marrium
args <- c("marr_1_SDH_DDH_Excit_Inhib", "../../Datasets/neurons_only_2021/clean_neuron_object.RDS")
args <- c("marr_human_DH_Excit_Inhib", "../../Datasets/human_spinalcord_2022/top_level_new_annotation.rda")

### for jessica
args <- c("jess_SDH_DDH_Excit_Inhib", "../../Datasets/neurons_only_2021/clean_neuron_object.RDS")
args <- c("jess_human_DH_Excit_Inhib", "../../Datasets/human_spinalcord_2022/top_level_new_annotation.rda")

### for clare
args <- c("clare_SDH_DDH_Excit_Inhib", "../../Datasets/neuron_and_glia_2022/final_meta_dataset.rds")
# Note switched to final_meta_dataset with Excit-1 not Excit-01
args <- c("clare_human_DH_Excit_Inhib", "../../Datasets/human_spinalcord_2022/top_level_new_annotation.rda")

# Need to set dir
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#source pre_analysis_functions for all the functions
source("loadLibraries.R")
source("JSON_Handler.R")
source("Pre_analysis_functions.R")

#source scripts to get them to run.
# TODO: source a script that sources scripts with all the plots?
source("DotPlot.R")
source("PooledDotPlot.R")
source("ClusterPlot.R")
#source("Quad_BarChart2.r")

#load data
start <- Sys.time()
RDfile <- load_data(args[2])
end <- Sys.time()
(duration <- end - start)
# duration <- c(43.22078, 28.90538, 1126.886, 546.6604)
# start <- c("2022-03-24 08:37:09 EDT", "2022-03-24 10:52:30 EDT", NA, "2022-03-25 08:02:31 EDT")
# condition <- c("Early morning with no browser only MS OneNote", "Nothing running
# but discord and MS OneNote", "Late night lots of things open no Google chrome",
# "Early morning with nothing open")
# data.type <- c("mouse", "mouse", "human", "human")
# duration <- data.frame(dur = duration, start = start, condition = conditions, data = data.type)

# Which extra groupings do you want?
# For extra_pool top remove the cluster area
extra_pool <- list()
extra_pool[["top"]] <- list("dataset", "age", "final_cluster_assignment", "run", "nCount_RNA")
extra_pool[["1"]] <- list("id", "features.label", "subgr")
extra_pool[["2"]] <- list("id", "features.label")
extra_pool[["3"]] <- list("features.label")
#extra_pool[["1"]] <- list("class_label", "region_label")

# Human dataset extrapool
extra_pool <- list()
extra_pool[["top"]] <- list("nCount_RNA", "nCount_SCT", "batches_combined",
                            "HSC_Subcluster_Annotation", "new_annotation",
                            "top_level_annotation")
extra_pool[["1"]] <- list("id", "features.label")
extra_pool[["2"]] <- list("features.label")

# Everything relies on all_cell_roster
all_cell_roster <- returnAllCellRoster(RDfile)
saveRDS(all_cell_roster, "../../Datasets/all_cell_rosters/all_cell_roster_human_jess.RDS")
all_cell_roster <- readRDS("../../Datasets/all_cell_rosters/all_cell_roster_human_jess.RDS")
all_cell_roster <- readRDS("../../Datasets/all_cell_rosters/all_cell_roster_jessica.RDS")

#call each script while passing the data as a value.

#dotplot.R
lbc <- createListbyCluster(scale.method = "z-score")
mainDP(lbc, height = 2000, width = 6700, legend.position = "bottom",
       legend.box = "horizontal", legend.title.angle = 0,
       legend.margin = margin(l = 100))

#pooled dotplot.R
CPR <- createClusterPoolResults(scale.method = "z-score", pool.level = "2")
mainPDP(CPR, height = 1500, width = 1500, factor.order = c("SDH", "DDH"))

# The cluste range plot


##### For human
lbc <- createListbyCluster(scale.method = "log10")
mainDP(lbc, height = 2000, width = 6700, legend.position = "bottom",
       legend.box = "horizontal", legend.title.angle = 0,
       legend.margin = margin(l = 100))

#pooled dotplot.R
CPR <- createClusterPoolResults(scale.method = "log10")
mainPDP(CPR, height = 1500, width = 1500)


