#### DotPlot
# This script is for creating the unpooled dotplot showing
# the average expression on a color scale and percent expressed
# as the size of the dot for all the clusters gene combinations
# It is similar to the DotPlot done by Seurat
# The DotPlots are split by subgroup

mainDP <- function(lbc, transp = F, width = NA, height = NA, dpi = 300,
                   legend.position = "right", legend.box = "horizontal",
                   legend.title.angle = 90, legend.box.margin = margin(l = 5, t = 80),
                   base.name = "DotPlot", global_size = 30, add.label = F,
                   saveimage = T, rm.labs = "", max.dot.size = 15,
                   legend.key.size = unit(1.5, "line"), legendbelow = F, ...){
  
  ifelse(class(legend.box.margin)[1] != "margin", print("You nead to use a ggplot margin object"), 0)
  
  # Create the Dot plot in one function
  Plot <- lbc %>%
    ggplot() +
    geom_point(aes(y = features.label, x = cluster, color = avg.exp.scaled, size = pct.exp)) +
    {if (add.label) geom_label(aes(y = features.label, x = cluster, label = signif(avg.exp, digits = 2)))} +
    labs(x = "Cluster", y = "Gene", color = "Avg exp scaled", size = "% Expressing") +
    scale_size(limits = c(0, 100), range = c(0, 15), breaks = c(0, 5, 10, 15, 25, 100)) +
    scale_color_viridis_c(breaks = c(-1, 0, 1), option = 5) +
    cowplot::theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, face = "bold",
                                     hjust = 1, size = global_size,
                                     color = "black"),
          axis.text.y = element_text(angle = 0, vjust = 0.5, face = "bold",
                                     hjust = 1, size = global_size,
                                     color = "black"),
          axis.title = element_text(size = global_size, face = "bold"),
          axis.line = element_line(color = "black", linewidth = 1),
          title = element_text(face = "bold"),
          legend.key.size = legend.key.size,
          legend.text = element_text(angle = 0, face = "bold", size = global_size/1.75),
          legend.title = element_text(size = global_size, angle = legend.title.angle),
          legend.position = legend.position,
          legend.box = legend.box,
          legend.box.margin = legend.box.margin,
          legend.spacing.x = unit(0.85, "line"),
          legend.box.spacing = unit(1, "line"),
          plot.background = element_rect(fill = ifelse(transp == T, "transparent", "white"))) +
    theme(...)
  
  
  # Below are the fine tuning elements such as removing labels, changing legends
  #
  # Remove labels
  if(rm.labs == "xy") {
    Plot <- Plot +
      theme(
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()
      )
  } else if (rm.labs == "y") {
    Plot <- Plot +
      theme(
        axis.title = element_blank(),
        axis.text.y = element_blank()
      )
  }
  
  # Right legends or bottom legends
  if(legendbelow == F) {
    Plot <- Plot +
      theme(
        legend.position = "right",
        legend.box = "horizontal",
        legend.box.margin = margin(l = 5, t = 80),
        legend.box.just = "top",
      )
  } else {
    Plot <- Plot +
      theme(
        legend.position = "bottom",
        legend.box = "vertical",
        legend.box.margin = margin(r = 5, b = 80),
        legend.box.just = "top",
      )
  }
  
  # Add guides
  Plot <- Plot +
    guides(color = guide_colorbar(barwidth = 1.2, barheight = 9, ticks = T,
                                  label.position = "left", title.position = "top"),
           size = guide_legend(label.position = "left", title.position = "top"))
  
  # External responding statement at the end of the function to either save it or not
  #
  # Save image function
  ifelse(saveimage,
         save_image(base.name, Plot,
                    height = ifelse(!is.na(height),height, 2700),
                    width = ifelse(!is.na(width), width, 7200),
                    device = "png", dpi = dpi),
         return(Plot))
  
}


