#TODO: need to take arguments for knowing what data file to use..

#something like this
#args <- c("NMDA_mouse_brain", "C:/Users/no/Documents/Neuroscience - MSc/Winter 2022/BIOL 5502/Rmd_proj/Seurat.ss.rda")
#args <- c("human_sex_nmda", "../../Data/human_ariel_data/top_level_new_annotation.rda")
args <- c("nmda_sdh_ddh", "NULL")
#args <- commandArgs(trailingOnly = TRUE)
print(paste ("this is the title project name entered: ", args[1]))
print(paste ("this was the given directory:", args[2]))

# The Grin gene analysis
args <- c("nmda_sdh_ddh", "../../Datasets/neurons_only_2021/clean_neuron_object.RDS")

### for kcni
args <- c("combined", "~/Neuroscience - MSc/Summer 2022/kcni_summer_school/Data/mouse_human_practice.rds")

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
RDSfile <- load_data(args[2])
end <- Sys.time()
duration <- end - start
duration
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
#extra_pool[["1"]] <- list("class_label", "region_label")

# For kcni
colnames(RDSfile@meta.data)
unique(RDSfile@meta.data$subclass_label)
extra_pool <- list()
extra_pool[["top"]] <- list("donor_sex_label")

# Everything relies on all_cell_roster and cell_roster
all_cell_roster <- returnAllCellRoster(RDSfile)
saveRDS(all_cell_roster, "../../Datasets/all_cell_rosters/all_cell_roster_grin.RDS")

#call each script while passing the data as a value.

#dotplot.R
lbc <- createListbyCluster(RDSfile)
mainDP(lbc)

#pooled dotplot.R
CPR <- createClusterPoolResults("1", "z-score")
mainPDP(CPR)

# The cluste range plot


#quadbar_chart2.R
#ListByCluster <- createListbyCluster(RDSfile)
mainQBC(CPR)

