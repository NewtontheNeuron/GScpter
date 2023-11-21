# Quick run script for an analysis using the Russ et al. (2021) single cell RNA
# sequencing data. Just run everything and you will get everything
runwidth <- c(13600, 13600, 13600, 4400, 4400)
runheight <- 3800

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("loadLibraries.R")
source("JSON_Handler.R")
source("Pre_analysis_functions.R")
source("DotPlot.R")
source("PooledDotPlot.R")
source("ClusterPlot.R")

loadconfig("russ_2021_mouse_SDH_DDH_template")
RDfile <- load_data("../../Datasets/neurons_and_glia_2022/final_meta_dataset.rds")
assay <- "raw"
slot <- "data"

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

#### Factor recode and relevel ####
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
mainDP(lbc, transp = T, height = runheight, width = runwidth[1], dpi = 600)

# Then without labels
gain <- mainDP(lbc, saveimage = F, rm.labs = "xy", transp = T)
ggsave(outnameform("mouse_fdp_zsoflog1pCPMnolabs.png"),
       plot = gain, height = runheight, width = runwidth[2], device = "png",
       units = "px", dpi = 600, type = "cairo")

# Then without labels on y alone
gain <- mainDP(lbc, saveimage = F, rm.labs = "y", transp = T)
ggsave(outnameform("mouse_fdp_zsoflog1pCPMnolabs_yonly.png"),
       plot = gain, height = runheight, width = runwidth[3], device = "png",
       units = "px", dpi = 600, type = "cairo")



#### pooled dotplot.R ####
# Mouse: GRIN expression in SDH vs. DDH neurons
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
ggsave(outnameform("mouse_PDP_SDHvDDH_zsoflog1pCPM.png"),
       plot = power, height = runheight, width = runwidth[4], device = "png", units = "px",
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
ggsave(outnameform("mouse_PDP_SDHvDDH_zsoflog1pCPMnolabs.png"),
       plot = power, height = runheight, width = runwidth[5], device = "png", units = "px",
       dpi = 600, type = "cairo")

