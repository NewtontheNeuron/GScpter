library(rstudioapi)

#set working directory to the one this file is currently in
setwd(dirname(getActiveDocumentContext()$path))

source("Pre_analysis_functions.R")

#function not returning anything.
ListbyClusterAll <- returnListbyClusterAll()

# Plotting all the relevant clusters form the data by
# using ClusterPoolAll
Gene <- ListbyClusterAll$features.label
Cluster <- ListbyClusterAll$id
AvgExpScaled <- ListbyClusterAll$avg.exp.scaled
markers <- Gene %>% unique()

Plot <- ListbyClusterAll %>% 
  filter(features.label %in% markers) %>% 
    mutate(`% Expressing` = ListbyClusterAll$pct.exp) %>% 
    ggplot(aes(y=Gene, x = Cluster, color = AvgExpScaled, size = `% Expressing`)) +
    geom_point() + 
    scale_size(range = c(0, 25)) +
    scale_color_viridis_c(option = "plasma") + 
    cowplot::theme_cowplot() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 25), 
            axis.title = element_text(size = 25, face="bold"), legend.key.size = unit(2.1, "line"), 
            legend.text = element_text(size = 17), legend.title = element_text(size = 20)) +
    theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size = 25))
Plot

# Save the image
# You will have to resize the Rstudio box
# or set the prefered width and height
# >>>> input required >>>>

save_image('DotPlot',Plot)