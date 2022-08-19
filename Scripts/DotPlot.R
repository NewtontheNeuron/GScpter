#### DotPlot
# This script is for creating the unpooled dotplot showing
# the average expression on a color scale and percent expressed
# as the size of the dot for all the clusters gene combinations
# It is similar to the DotPlot done by Seurat
# The DotPlots are split by subgroup

# Function to create the DotPlot
createDotPlot <- function(lbc_filtered){
  # Setting a gloabal size factor
  global_size <- 40
  # Creating the ggplot DotPlot
  lbc_filtered %>%
    ggplot(aes(y = features.label, x = cluster, color = avg.exp.scaled, size = pct.exp)) +
    geom_point() +
    labs(x = "Cluster", y = "Gene", color = "Avg exp scaled", size = "% Expressing") +
    scale_size(range = c(0, 20)) +
    scale_color_viridis_c(option = "plasma") +
    cowplot::theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,# 45 1
                                     hjust = 1, size = global_size,
                                     color = "black"),
          axis.text.y = element_text(angle = 0, vjust = 0.5,
                                     hjust = 1, size = global_size,
                                     color = "black"),
          axis.title = element_text(size = global_size, face = "bold"),
          legend.key.size = unit(2, "line"),
          legend.text = element_text(size = global_size/1.75),
          legend.title = element_text(size = global_size, angle = 90),
          legend.box.background = element_rect(color = "white"),
          #legend.position = "bottom",
          legend.box = "horizontal",
          legend.box.margin = margin(b = -120), # -190
          legend.spacing.x = unit(0.5, "line"),
          plot.background = element_rect(fill = "white"),
          plot.margin = margin(b = 75),
          panel.background = element_rect(fill = "white"))
}

mainDP <- function(lbc){
  # TODO: Loop through the subgroups and plot and save the DotPlots
  
  Plot <- createDotPlot(lbc)
  
  # Save the image, you can set the prefered width and height
  # width = 10000 for 4 subgroups, 4500 for 2 subgroups, width = 2250
  save_image('DotPlot', Plot, height = 2100, width = 4500)
}


