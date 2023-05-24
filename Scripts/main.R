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
args <- c("grin_only_SDH_Excit_Inhib", "../../Datasets/neurons_and_glia_2022/final_meta_dataset.rds")
# Unit testing group over
args <- c("grin_only_SDH_DDH_Excit_Inhib_for_unit_testing_groupover", "../../Datasets/neurons_and_glia_2022/final_meta_dataset.rds")

### for kcni
args <- c("combined", "~/Neuroscience - MSc/Summer 2022/kcni_summer_school/Data/mouse_human_practice.rds")

### for marrium
args <- c("marr_1_SDH_DDH_Excit_Inhib", "../../Datasets/neurons_and_glia_2022/final_meta_dataset.rds")
args <- c("marr_human_DH_Excit_Inhib", "../../Datasets/human_spinalcord_2022/top_level_new_annotation.rda")

### for jessica
args <- c("jess_SDH_DDH_Excit_Inhib", "../../Datasets/neurons_and_glia_2022/final_meta_dataset.rds")
args <- c("jess_human_DH_Excit_Inhib", "../../Datasets/human_spinalcord_2022/top_level_new_annotation.rda")

### for clare
args <- c("clare_SDH_DDH_Excit_Inhib", "../../Datasets/neuron_and_glia_2022/final_meta_dataset.rds")
# Note switched to final_meta_dataset with Excit-1 not Excit-01
args <- c("clare_human_DH_Excit_Inhib", "../../Datasets/human_spinalcord_2022/top_level_new_annotation.rda")

# Laurence
args <- c("grik_SDH_DDH_Excit_Inhib", "../../Datasets/neurons_and_glia_2022/final_meta_dataset.rds")


#### Start up ####
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

RDfile <- load_data(args[2])
assay <- "raw"
assay <- "RNA"
slot <- "data"
#Todo: show user a list of all the possible assays and then the user picks one

# Which groupings do you want?
extra_pool <- list()
extra_pool[["top"]] <- list("dataset", "age", "final_cluster_assignment", "run", "nCount_RNA")
extra_pool[["1"]] <- list("id", "features.label", "subgr")
extra_pool[["2"]] <- list("id", "features.label")
extra_pool[["3"]] <- list("features.label")
extra_pool[["4"]] <- list("features.label", "subgr")
extra_pool[["5"]] <- list("features.label", "subgr", "age")

user_order <- list()
user_order[["1"]] <- list("SDH Excitatory", "SDH Inhibitory",
                       "DDH Excitatory", "DDH Inhibitory")
user_order[["2"]] <- list("SDH", "DDH")
user_order[["4"]] <- list("Excitatory", "Inhibitory")

# Human dataset extrapool
extra_pool <- list()
extra_pool[["top"]] <- list("nCount_RNA", "nCount_SCT", "batches_combined",
                            "HSC_Subcluster_Annotation", "new_annotation",
                            "top_level_annotation")
extra_pool[["1"]] <- list("id", "features.label")
extra_pool[["2"]] <- list("features.label") # TODO: because of pool and paste you either have to throw and error or promt a different selection
extra_pool[["3"]] <- list("features.label", "subgr")
extra_pool[["4"]] <- list("features.label", "batches_combined")

#### All Cell Roster ####
all_cell_roster <- returnAllCellRoster(RDfile)
saveRDS(all_cell_roster, "../../Datasets/all_cell_rosters/all_cell_roster_human_jess.RDS")
all_cell_roster <- readRDS("../../Datasets/all_cell_rosters/all_cell_roster_human_jess.RDS")
all_cell_roster <- readRDS("../../Datasets/all_cell_rosters/all_cell_roster_jessica.RDS")

#call each script while passing the data as a value.

#### Factor recode and relevel ####
# For IASP poster
# Recode the cluster names as ES1 or ID1 etc.
all_cell_roster <- all_cell_roster %>%
  as_tibble() %>%
  mutate(orig.name = cluster,
         cluster = case_when(
           grepl("Excit-", cluster) ~ str_replace(cluster, "Excit-", "E"),
           grepl("Inhib-", cluster) ~ str_replace(cluster, "Inhib-", "I"),
         ))

all_cell_roster$cluster <- factor(all_cell_roster$cluster,
                                  levels = c("E1", "E2", "E3", "E4", "E5",
                                             "E6", "E7", "E8", "E9", "E10", 
                                             "E11", "E12", "E13", "E14", "E15",
                                             "E16", "E17", "E18", "E19", "E20",
                                             "E21", "E22", "E23", "E24", "E25", 
                                             "E26", "E27", "E28", "E29", "E30",
                                             "E31", "E32", "E33", "E34", "E35",
                                             "E36", "I1", "I2", "I3", "I4", "I5",
                                             "I6", "I7", "I8", "I9", "I10",
                                             "I11", "I12", "I13", "I14", "I15",
                                             "I16", "I17", "I18", "I19", "I20"))

#### dotplot.R ####
lbc <- createListbyCluster(scale.method = "zsoflog1pCPM", pre.expm1 = F)
mainDP(lbc, height = 2600, width = 7600, legend.position = "bottom",
       legend.box = "horizontal", legend.title.angle = 0,
       legend.margin = margin(l = 0), base.name = "mouse_FullDotPlot_zsoflog1pCPM")
mainDP(lbc, height = 2100, width = 7400, legend.position = "bottom",
       legend.box = "horizontal", legend.title.angle = 0,
       legend.margin = margin(l = 0), base.name = "mouse_FullDotPlot_zsoflog1pCPM_nolabs",
       rm.labs = T)
mainDP(lbc, height = 2200, width = 7050, legend.position = "right",
       legend.box = "horizontal", legend.title.angle = 90, transp = T,
       legend.margin = margin(l = 0, t = 50), base.name = "canacn2023/mouse_fdp_zsoflog1pCPM",
       max.dot.size = 15, legend.key.height = unit(1.7, "line"),
       legend.text = element_text(angle = 0, vjust = 0.5),
       legend.spacing.x = unit(0.2, "line"),
       legend.box.just = "bottom", saveorret = F) +
  guides(color = guide_colorbar(barwidth = 0.5, barheight = 10,
                                ticks = F, label.position = "left",
                                label.theme = element_text(face = "bold", color = "green4"),
                                title.position = "left", nbin = 3),
         size = guide_legend(label.position = "left", title.position = "right",
                             label.hjust = 1.5, label.vjust = 0,
                             label.theme = element_text(color = "skyblue3", face = "bold",
                                                        size = 20))) +
  scale_color_viridis_c(breaks = c(-1, 0, 1), option = 13) +
  scale_size(limits = c(0, 100), range = c(0, 15))
mainDP(lbc, height = 2200, width = 6400, legend.position = "right",
       legend.box = "horizontal", legend.title.angle = 0, transp = T,
       legend.margin = margin(l = 0, t = - 15), base.name = "canacn2023/mouse_fdp_zsoflog1pCPM",
       max.dot.size = 15, legend.key.width = unit(3.5, "line"),
       legend.text = element_text(angle = 0, hjust = 0.5),
       legend.box.just = "right", legend.spacing.x = unit(3, "line"))

#### pooled dotplot.R ####
CPR <- createClusterPoolResults(scale.method = "log10", pool.level = "1")
mainPDP(CPR, base.name = "PooledDotPlot", height = 1500,
        width = 2000, factor.order = unlist(user_order[["1"]]))

CPR <- createClusterPoolResults(scale.method = "zsoflog1pCPM", pool.level = "2",
                                pre.expm1 = F)
mainPDP(CPR, base.name = "mouse_Pooled_SDHvDDH_zsoflog1pCPM", height = 1500,
        width = 1500, factor.order = unlist(user_order[["2"]]),
        legend.margin = margin(b = -55))
mainPDP(CPR, base.name = "mouse_Pooled_SDHvDDH_zsoflog1pCPM_nolabs", height = 1500,
        width = 1200, factor.order = unlist(user_order[["2"]]),
        legend.margin = margin(), rm.labs = T)
mainPDP(CPR, base.name = "canacn2023/mouse_Pooled_SDHvDDH_zsoflog1pCPM", height = 2200,
        width = 1500, factor.order = unlist(user_order[["2"]]),
        legend.margin = margin(b = -55))


CPR <- createClusterPoolResults(scale.method = "log10", pool.level = "3")
CPR$group.label <- "Dorsal Horn"
mainPDP(CPR, base.name = "Pooled_Dorsal_Horn", height = 1500,
        width = 1500)

CPR <- createClusterPoolResults(scale.method = "zsoflog1pCPM", pool.level = "4",
                                pre.expm1 = F)
mainPDP(CPR, base.name = "Pooled_DH_EvI", height = 1500,
        width = 2000, factor.order = unlist(user_order[["4"]]))
mainPDP(CPR, base.name = "canacn2023/PDP_DH_EvI", height = 2200,
        width = 2750, factor.order = unlist(user_order[["4"]]),
        transp = T, max.dot.size = 15, legend.key.height = unit(2, "line"),
        legend.text.align = 1, legend.spacing.y = unit(0, "line"),
        legend.text = element_text(vjust = 0.5))




# ---- Human SC ----
#


all_cell_roster$cluster <- factor(all_cell_roster$cluster,
                                  levels = c("Ex-Dorsal-1", "Ex-Dorsal-2",
                                             "Ex-Dorsal-3", "Ex-Dorsal-4", 
                                             "Ex-Dorsal-5", "Ex-Dorsal-6",
                                             "Ex-Dorsal-7", "Ex-Dorsal-8",
                                             "Ex-Dorsal-9", "Ex-Dorsal-10",
                                             "Ex-Dorsal-11", "Ex-Dorsal-12",
                                             "Inh-Dorsal-1", "Inh-Dorsal-2",
                                             "Inh-Dorsal-3", "Inh-Dorsal-4",
                                             "Inh-Dorsal-5", "Inh-Dorsal-6",
                                             "Inh-Dorsal-7", "Inh-Dorsal-8",
                                             "Inh-Dorsal-9", "Inh-Dorsal-10"))
# For presentations
# Recode the cluster names as ED1 or ID1 etc.
all_cell_roster <- all_cell_roster %>%
  mutate(orig.name = cluster,
         cluster = case_when(
           grepl("Ex-", cluster) ~ str_replace(cluster, "Ex-Dorsal-", "ED"),
           grepl("Inh-", cluster) ~ str_replace(cluster, "Inh-Dorsal-", "ID")
         ))
all_cell_roster$cluster <- factor(all_cell_roster$cluster,
                                  levels = c("ED1", "ED2", "ED3", "ED4", 
                                             "ED5", "ED6", "ED7", "ED8",
                                             "ED9", "ED10", "ED11", "ED12",
                                             "ID1", "ID2", "ID3", "ID4",
                                             "ID5", "ID6", "ID7", "ID8",
                                             "ID9", "ID10"))

# Graphs for human data
lbc <- createListbyCluster(scale.method = "zsoflog1pCPM", pre.expm1 = F)
mainDP(lbc, height = 2500, width = 6000, legend.position = "bottom",
       legend.box = "horizontal", legend.title.angle = 0, base.name = "human_DotPlot_zsoflog1pCPM",
       legend.margin = margin(l = 0), add.label = F, saveorret = T)
mainDP(lbc, height = 2000, width = 5800, legend.position = "bottom",
       legend.box = "horizontal", legend.title.angle = 0, base.name = "human_DotPlot_zsoflog1pCPM_nolabs",
       legend.margin = margin(l = 0), add.label = F, saveorret = T, rm.labs = T)
mainDP(lbc, height = 2200, width = 4300, legend.position = "right",
       legend.box = "horizontal", legend.title.angle = 90, transp = T,
       legend.margin = margin(l = 0, t = 50), base.name = "canacn2023/human_fdp_zsoflog1pCPM",
       max.dot.size = 15, legend.key.height = unit(1.5, "line"),
       legend.text = element_text(angle = 0),
       legend.spacing.x = unit(0.5, "line"),
       legend.box.just = "top")
mainDP(lbc, height = 2800, width = 3500, legend.position = "right",
       legend.box = "vertical", legend.title.angle = 0, transp = T,
       legend.margin = margin(l = 0, t = - 15), base.name = "canacn2023/human_fdp_zsoflog1pCPM",
       max.dot.size = 15, legend.key.width = unit(1.5, "line"),
       legend.text = element_text(angle = 0, hjust = 0.5),
       legend.box.just = "left", legend.spacing.x = unit(0.5, "line"),
       legend.direction = "horizontal")

# TODO: Fix uncentered legends for jess

#pooled dotplot.R
all_cell_roster$id <- "Dorsal Horn"
CPR <- createClusterPoolResults(scale.method = "zsoflog1pCPM", pool.level = "1")
mainPDP(CPR, "human_PooledDorsalHorn_zsoflog1pCPM", height = 1500, width = 1500,
        legend.margin = margin(t = 55))
mainPDP(CPR, "human_PooledDorsalHorn_zsoflog1pCPM_nolabs", height = 1500, width = 1200,
        legend.margin = margin(), rm.labs = T)

CPR <- createClusterPoolResults(scale.method = "zsoflog1pCPM", pool.level = "3",
                                pre.expm1 = F)
mainPDP(CPR, "canacn2023/PDP_EvI", height = 2200,
        width = 2800, factor.order = unlist(user_order[["4"]]),
        transp = T, max.dot.size = 15, legend.key.height = unit(2, "line"),
        legend.text.align = 1, legend.spacing.y = unit(1, "line"),
        legend.text = element_text(vjust = 0.5))

CPR <- createClusterPoolResults(scale.method = "log10", pool.level = "4")
mainPDP(CPR, "Samples_PDP_DH_count_SCT_lab", height = 2300, width = 2300)

# ClusterRange Plot
mainCluP(pool.level = "3")
