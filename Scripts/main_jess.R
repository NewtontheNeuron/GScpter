# This is a default main script where you run the functoins specific to your
# analysis

# Select your research prompt
args <- c("jess_SDH_DDH_Excit_Inhib", "../../Datasets/neurons_and_glia_2022/final_meta_dataset.rds")

#### Start up ####
# Need to set dir
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#source pre_analysis_functions for all the functions
source("loadLibraries.R")
source("JSON_Handler.R")
source("Pre_analysis_functions.R")
source("DotPlot.R")
source("PooledDotPlot.R")
source("ClusterPlot.R")

loadconfig(args[1])
RDfile <- load_data(args[2])
assay <- "raw"
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

#### All Cell Roster ####
all_cell_roster <- returnAllCellRoster(RDfile)

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
                                  levels = c("E1", "E14", "E16", "E18", "E2", "E3", 
                                             "E8", "E9", "E10", "E12", "E15",
                                             "E4", "E5", "E6", "E13",
                                             "E19", "E20", "E21", "E22", "E23",
                                             "E24", "E25", "E26", "E27", "E28",
                                             "E29", "E30", "E31", "E32", "E33",
                                             "E34", "E35", "E36",
                                             "I5", "I9", "I3", "I4", "I11", "I12",
                                             "I13", "I1", "I2", "I6", "I7", "I8",
                                             "I10", "I14", "I15", "I16", "I17",
                                             "I18", "I19", "I20", "I21"))

#### dotplot.R ####
lbc <- createListbyCluster(scale.method = "zsoflog1pCPM", pre.expm1 = F)
gain <- mainDP(lbc, height = 2200, width = 7050, legend.position = "right",
       legend.box = "horizontal", legend.title.angle = 90, transp = T,
       legend.margin = margin(l = 0, t = 0),
       max.dot.size = 15, legend.box.margin = margin(l = 5, t = 80),
       legend.text = element_text(angle = 0, face = "bold"),
       legend.spacing.x = unit(0.85, "line"), #title = element_text(margin = margin(l=300)),
       legend.box.spacing = unit(1, "line"), global_size = 30,
       legend.box.just = "top", saveorret = F,
       axis.line = element_line(color = "black", linewidth = 1),
       axis.text = element_text(face = "bold"), title = element_text(face = "bold")) +
  guides(color = guide_colorbar(barwidth = 1.2, barheight = 9,
                                ticks = T, label.position = "left",
                                title.position = "top"),
         size = guide_legend(label.position = "left", title.position = "top")) +
  scale_color_viridis_c(breaks = c(-1, 0, 1), option = 5) +
  scale_size(limits = c(0, 100), range = c(0, 15), breaks = c(0, 5, 10, 15, 25, 100))
ggsave("../Output/Jessica/mouse_fdp_zsoflog1pCPM.png",
       plot = gain, height = 3200,
       width = 13600, device = "png", units = "px", dpi = 600, type = "cairo")
# Then without labels
gain <- mainDP(lbc, height = 2200, width = 7050, legend.position = "right",
               legend.box = "horizontal", legend.title.angle = 90, transp = T,
               legend.margin = margin(l = 0, t = 0), rm.labs = T,
               max.dot.size = 15, legend.box.margin = margin(l = 5, t = 0),
               legend.text = element_text(angle = 0, face = "bold"),
               legend.spacing.x = unit(0.85, "line"), #title = element_text(margin = margin(l=300)),
               legend.box.spacing = unit(1, "line"), global_size = 30,
               legend.box.just = "top", saveorret = F,
               axis.line = element_line(color = "black", linewidth = 1),
               axis.text = element_text(face = "bold"), title = element_text(face = "bold")) +
  guides(color = guide_colorbar(barwidth = 1.2, barheight = 9,
                                ticks = T, label.position = "left",
                                title.position = "top"),
         size = guide_legend(label.position = "left", title.position = "top")) +
  scale_color_viridis_c(breaks = c(-1, 0, 1), option = 5) +
  scale_size(limits = c(0, 100), range = c(0, 15), breaks = c(0, 5, 10, 15, 25, 100))
ggsave("../Output/Jessica/mouse_fdp_zsoflog1pCPMnolabs.png",
       plot = gain, height = 3200,
       width = 13600, device = "png", units = "px", dpi = 600, type = "cairo")

# Then without labels on y alone
gain <- mainDP(lbc, height = 2200, width = 7050, legend.position = "right",
               legend.box = "horizontal", legend.title.angle = 90, transp = T,
               legend.margin = margin(l = 0, t = 0), rm.labs = F,
               max.dot.size = 15, legend.box.margin = margin(l = 5, t = 0),
               legend.text = element_text(angle = 0, face = "bold"),
               legend.spacing.x = unit(0.85, "line"), #title = element_text(margin = margin(l=300)),
               legend.box.spacing = unit(1, "line"), global_size = 30,
               legend.box.just = "top", saveorret = F,
               axis.line = element_line(color = "black", linewidth = 1),
               axis.text.y = element_text(face = "bold"), title = element_text(face = "bold")) +
  guides(color = guide_colorbar(barwidth = 1.2, barheight = 9,
                                ticks = T, label.position = "left",
                                title.position = "top"),
         size = guide_legend(label.position = "left", title.position = "top")) +
  scale_color_viridis_c(breaks = c(-1, 0, 1), option = 5) +
  scale_size(limits = c(0, 100), range = c(0, 15), breaks = c(0, 5, 10, 15, 25, 100)) +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank()
  )
ggsave("../Output/Jessica/mouse_fdp_zsoflog1pCPMnolabs.png",
       plot = gain, height = 3800,
       width = 13600, device = "png", units = "px", dpi = 600, type = "cairo")

#### pooled dotplot.R ####
# Mouse and human: GRIN expressoin in SDH vs. DDH neurons
# of the dorsal horn
CPR <- createClusterPoolResults(scale.method = "zsoflog1pCPM", pool.level = "2",
                                pre.expm1 = F)
CPR$group.label <- factor(CPR$group.label, levels = c("SDH", "DDH"))
CPR
mainPDP(CPR, base.name = "PDP_SDHvDDH", height = 2200,
        width = 2750, factor.order = unlist(user_order[["2"]]),
        transp = T, max.dot.size = 20, yieldplot = T,
        legend.key.height = unit(2, "line"),
        legend.text.align = 1, legend.spacing.x = unit(0.85, "line"),
        legend.box.spacing = unit(1, "line"), global_size = 30,
        legend.box.just = "top",
        legend.box.margin = margin(t = 70),
        legend.text = element_text(angle = 0, face = "bold", vjust = 0.5),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.text = element_text(face = "bold"),
        title = element_text(face = "bold")) +
  guides(color = guide_colorbar(barwidth = 1.2, barheight = 9,
                                ticks = T, label.position = "left",
                                title.position = "top"),
         size = guide_legend(label.position = "left", title.position = "top")) +
  scale_color_viridis_c(breaks = c(-1, 0, 1), option = 5) +
  scale_size(limits = c(0, 100), range = c(0, 20)) -> power
ggsave("../Output/Jessica/mouse_PDP_SDHvDDH_zsoflog1pCPM.png",
       plot = power, height = 3200, width = 4400, device = "png", units = "px",
       dpi = 600, type = "cairo")
# Without labels
mainPDP(CPR, base.name = "PDP_SDHvDDH", height = 2200,
        width = 2750, factor.order = unlist(user_order[["2"]]),
        transp = T, max.dot.size = 15, yieldplot = T,
        legend.key.height = unit(2, "line"),
        legend.text.align = 1, legend.spacing.x = unit(0.85, "line"),
        legend.box.spacing = unit(1, "line"), global_size = 30,
        legend.box.just = "top", rm.labs = T,
        legend.box.margin = margin(t = 0),
        legend.text = element_text(angle = 0, face = "bold", vjust = 0.5),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.text = element_text(face = "bold"),
        title = element_text(face = "bold")) +
  guides(color = guide_colorbar(barwidth = 1.2, barheight = 9,
                                ticks = T, label.position = "left",
                                title.position = "top"),
         size = guide_legend(label.position = "left", title.position = "top")) +
  scale_color_viridis_c(breaks = c(-1, -0.5, 0, 0.5, 1), option = 5) +
  scale_size(limits = c(0, 100), range = c(0, 15), breaks = c(0, 5, 10, 50, 100)) -> power
ggsave("../Output/Jessica/mouse_PDP_SDHvDDH_zsoflog1pCPMnolabs.png",
       plot = power, height = 3800, width = 2600, device = "png", units = "px",
       dpi = 600, type = "cairo")



# ---- Human SC ----
#

args <- c("jess_human_DH_Excit_Inhib", "../../Datasets/human_spinalcord_2022/top_level_new_annotation.rda")
source("loadLibraries.R")
source("JSON_Handler.R")
source("Pre_analysis_functions.R")
source("DotPlot.R")
source("PooledDotPlot.R")
source("ClusterPlot.R")

RDfile <- load_data(args[2])
assay <- "RNA"
slot <- "data"

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
gain <- mainDP(lbc, height = 2200, width = 4300, legend.position = "right",
       legend.box = "horizontal", legend.title.angle = 90, transp = T,
       legend.margin = margin(l = 0, t = 0),
       max.dot.size = 15, legend.box.margin = margin(l = 5, t = 70),
       legend.text = element_text(angle = 0, face = "bold"), global_size = 30,
       legend.spacing.x = unit(0.85, "line"),
       legend.box.just = "top", saveorret = F,
       axis.line = element_line(color = "black", linewidth = 1),
       axis.text = element_text(face = "bold"), title = element_text(face = "bold")) +
  guides(color = guide_colorbar(barwidth = 1.2, barheight = 9,
                                ticks = T, label.position = "left",
                                title.position = "top"),
         size = guide_legend(label.position = "left", title.position = "top"),
         label.vjust = -1) +
  scale_color_viridis_c(breaks = c(-1, 0, 1), option = 5) + # 13
  scale_size(limits = c(0, 100), range = c(0, 15))
ggsave("../Output/Jessica_human/human_fdp_zsoflog1pCPM.png",
       plot = gain, height = 2000,
       width = 4100, device = "png", units = "px", dpi = 300, type = "cairo")
# Without labels on y only
gain <- mainDP(lbc, height = 2200, width = 4300, legend.position = "right",
               legend.box = "horizontal", legend.title.angle = 90, transp = T,
               legend.margin = margin(l = 0, t = 0), rm.labs = F,
               max.dot.size = 15, legend.box.margin = margin(l = 5, t = 70),
               legend.text = element_text(angle = 0, face = "bold"), global_size = 30,
               legend.spacing.x = unit(0.85, "line"),
               legend.box.just = "top", saveorret = F,
               axis.line = element_line(color = "black", linewidth = 1),
               axis.text.y = element_text(face = "bold"), title = element_text(face = "bold")) +
  guides(color = guide_colorbar(barwidth = 1.2, barheight = 9,
                                ticks = T, label.position = "left",
                                title.position = "top"),
         size = guide_legend(label.position = "left", title.position = "top"),
         label.vjust = -1) +
  scale_color_viridis_c(breaks = c(-1, 0, 1), option = 5) + # 13
  scale_size(limits = c(0, 100), range = c(0, 15)) +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank()
  )
ggsave("../Output/Jessica_human/human_fdp_zsoflog1pCPMnolabs.png",
       plot = gain, height = 3800,
       width = 8200, device = "png", units = "px", dpi = 600, type = "cairo")


####Human Dorsal Horn Pool####
all_cell_roster$overallgroup <- "Dorsal Horn"
all_cell_roster
extra_pool[["6"]] <- list("features.label", "overallgroup")
CPR <- createClusterPoolResults(scale.method = "zsoflog1pCPM", pool.level = "6",
                                pre.expm1 = F)
CPR
greats <- mainPDP(CPR, height = 2200,
                  width = 1750, transp = T, max.dot.size = 15, yieldplot = T,
                  legend.key.height = unit(2, "line"),
                  legend.text.align = 1, legend.spacing.x = unit(0.85, "line"),
                  legend.box.spacing = unit(1, "line"), global_size = 30,
                  legend.box.just = "top",
                  legend.box.margin = margin(l = 0),
                  legend.text = element_text(angle = 0, face = "bold", vjust = 0.5),
                  axis.line = element_line(color = "black", linewidth = 1),
                  axis.text = element_text(face = "bold"),
                  title = element_text(face = "bold")) +
  guides(color = guide_colorbar(barwidth = 1.2, barheight = 9,
                                ticks = T, label.position = "left",
                                title.position = "top"),
         size = guide_legend(label.position = "left", title.position = "top")) +
  scale_color_viridis_c(breaks = c(-1, -0.5, 0, 0.5, 1), option = 5) +
  scale_size(limits = c(0, 100), range = c(0, 15))
ggsave("../Output/Jessica_human/human_PDP_DH_zsoflog1pCPM.png",
       plot = greats, height = 4400, width = 3900, device = "png", units = "px",
       dpi = 600, type = "cairo")
# Without labs
greats <- mainPDP(CPR, height = 2200,
                  width = 1750, transp = T, max.dot.size = 15, yieldplot = T,
                  legend.key.height = unit(2, "line"),
                  legend.text.align = 1, legend.spacing.x = unit(0.85, "line"),
                  legend.box.spacing = unit(1, "line"), global_size = 30,
                  legend.box.just = "top", rm.labs = T,
                  legend.box.margin = margin(l = 0),
                  legend.text = element_text(angle = 0, face = "bold", vjust = 0.5),
                  axis.line = element_line(color = "black", linewidth = 1),
                  axis.text = element_text(face = "bold"),
                  title = element_text(face = "bold")) +
  guides(color = guide_colorbar(barwidth = 1.2, barheight = 9,
                                ticks = T, label.position = "left",
                                title.position = "top"),
         size = guide_legend(label.position = "left", title.position = "top")) +
  scale_color_viridis_c(breaks = c(-1, -0.5, 0, 0.5, 1), option = 5) +
  scale_size(limits = c(0, 100), range = c(0, 15))
ggsave("../Output/Jessica_human/human_PDP_DH_zsoflog1pCPMnolabs.png",
       plot = greats, height = 3800, width = 2500, device = "png", units = "px",
       dpi = 600, type = "cairo")

