#### DotPlot
# This script is for creating the unpooled dotplot showing
# the average expression on a color scale and percent expressed
# as the size of the dot for all the clusters gene combinations
# It is similar to the DotPlot done by Seurat
# The DotPlots are split by subgroup

# Function to create the DotPlot
createDotPlot <- function(lbc, transp = F, legend.position, legend.box,
                          legend.title.angle, legend.margin, global_size,
                          add.label, max.dot.size, legend.key.size, ...){
  lbc %>%
    ggplot() +
    geom_point(aes(y = features.label, x = cluster, color = avg.exp.scaled, size = pct.exp)) +
    {if (add.label) geom_label(aes(y = features.label, x = cluster, label = signif(avg.exp, digits = 2)))} +
    labs(x = "Cluster", y = "Gene", color = "Avg exp scaled", size = "% Expressing") +
    scale_size(range = c(0, max.dot.size)) +
    scale_color_viridis_c(option = "plasma") +
    cowplot::theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,# 90 0.5
                                     hjust = 1, size = global_size,
                                     color = "black"),
          axis.text.y = element_text(angle = 0, vjust = 0.5,
                                     hjust = 1, size = global_size,
                                     color = "black"),
          axis.title = element_text(size = global_size, face = "bold"),
          #axis.title.x = element_text(margin = margin(t = -10)),
          legend.key.size = legend.key.size,
          legend.text = element_text(size = global_size/1.75),
          legend.title = element_text(size = global_size, angle = legend.title.angle),
          legend.position = legend.position,
          legend.box = legend.box,
          legend.box.margin = legend.margin,
          legend.spacing.x = unit(1.5, "line"),
          plot.background = element_rect(fill = ifelse(transp == T, "transparent", "white"))) +
    theme(...)
}

mainDP <- function(lbc, transp = F, width = NA, height = NA,
                   legend.position = "right", legend.box = "horizontal",
                   legend.title.angle = 90, legend.margin = margin(),
                   base.name = "DotPlot", global_size = 36, add.label = F,
                   saveorret = T, rm.labs = F, max.dot.size = 20,
                   legend.key.size = unit(1.5, "line"), ...){
  
  ifelse(class(legend.margin)[1] != "margin", print("You nead to use a ggplot margin object"),0)
  Plot <- createDotPlot(lbc, transp, legend.position, legend.box, legend.title.angle,
                        legend.margin, global_size, add.label, max.dot.size,
                        legend.key.size, ...)
  
  # Remove labels
  if(rm.labs) {
    Plot <- Plot +
      theme(
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()
      )
  }
  
  # width = 10000 for 4 subgroups, 4500 for 2 subgroups, width = 2250
  ifelse(saveorret,
         save_image(base.name, Plot,
                    height = ifelse(!is.na(height),height, 2700),
                    width = ifelse(!is.na(width), width, 7200),
                    device = "png"),
         return(Plot))
  
}
# For the iasp poster it is 1800 for both, 4900 on width for human and 8600 for
# the mouse.


