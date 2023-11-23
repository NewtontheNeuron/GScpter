# Quick run script for an analysis using the Yadav et al. (2023) single cell RNA
# sequencing data. Just run everything and you will get everything
runwidth <- c(8200, 8200, 8200, 4400, 4400)
runheight <- 3800

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("loadLibraries.R")
source("JSON_Handler.R")
source("Pre_analysis_functions.R")
source("DotPlot.R")
source("PooledDotPlot.R")
source("ClusterPlot.R")

loadconfig("yadav_2023_human_DH_template")
RDfile <- load_data("../../Datasets/human_spinalcord_2022/top_level_new_annotation.rda")
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
extra_pool[["2"]] <- list("features.label")
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

mainDP(lbc, transp = T, height = runheight, width = runwidth[1], dpi = 600)

# Then without labels
gain <- mainDP(lbc, saveimage = F, rm.labs = "xy", transp = T)
ggsave(outnameform("human_fdp_zsoflog1pCPMnolabs.png"),
       plot = gain, height = runheight, width = runwidth[2], device = "png",
       units = "px", dpi = 600, type = "cairo")

# Without labels on y only
gain <- mainDP(lbc, saveimage = F, rm.labs = "y", transp = T)
ggsave(outnameform("mouse_fdp_zsoflog1pCPMnolabs_yonly.png"),
       plot = gain, height = runheight, width = runwidth[3], device = "png",
       units = "px", dpi = 600, type = "cairo")


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
ggsave(outnameform("human_PDP_DH_zsoflog1pCPM.png"),
       plot = greats, height = runheight, width = runwidth[4], device = "png", units = "px",
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
ggsave(outnameform("human_PDP_DH_zsoflog1pCPMnolabs.png"),
       plot = greats, height = runheight, width = runwidth[5], device = "png", units = "px",
       dpi = 600, type = "cairo")

