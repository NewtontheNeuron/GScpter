#### DotPlot
# This script is for creating the unpooled dotplot showing
# the average expression on a color scale and percent expressed
# as the size of the dot for all the clusters gene combinations
# It is similar to the DotPlot done by Seurat
# The DotPlots are split by subgroup

# Function to create the DotPlot
createDotPlot <- function(lbc, transp = F, legend.position, legend.box,
                          legend.title.angle, legend.margin){
  # Setting a gloabal size factor
  global_size <- 36
  # Creating the ggplot DotPlot
  lbc %>%
    ggplot(aes(y = features.label, x = cluster, color = avg.exp.scaled, size = pct.exp)) +
    geom_point() +
    labs(x = "Cluster", y = "Gene", color = "Avg exp scaled", size = "% Expressing") +
    scale_size(range = c(0, 20)) +
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
          legend.key.size = unit(1.5, "line"),
          legend.text = element_text(size = global_size/1.75),
          legend.title = element_text(size = global_size, angle = legend.title.angle),
          #legend.box.background = element_rect(color = "white"),
          legend.position = legend.position,
          legend.box = legend.box,
          #legend.box.margin = margin(l = -80),
          legend.box.margin = legend.margin, # -190 -125
          legend.spacing.x = unit(1.5, "line"),
          plot.background = element_rect(fill = ifelse(transp == T, "transparent", "white")))
}

mainDP <- function(lbc, transp = F, width = NA, height = NA,
                   legend.position = "right", legend.box = "horizontal",
                   legend.title.angle = 90, legend.margin = margin()){
  # TODO: Loop through the subgroups and plot and save the DotPlots
  ifelse(class(legend.margin)[1] != "margin", print("You nead to use a ggplot margin object"),0)
  Plot <- createDotPlot(lbc, transp, legend.position, legend.box, legend.title.angle,
                        legend.margin)
  
  # Save the image, you can set the prefered width and height
  # width = 10000 for 4 subgroups, 4500 for 2 subgroups, width = 2250
  save_image('DotPlot', Plot, height = ifelse(!is.na(height), height, 2700),
             width = ifelse(!is.na(width), width, 7200), device = "png")
}
# For the iasp poster it is 1800 for both, 4900 on width for human and 8600 for
# the mouse.


