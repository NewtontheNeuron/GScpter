library(rstudioapi)

#set working directory to the one this file is currently in
setwd(dirname(getActiveDocumentContext()$path))

source("Pre_analysis_functions.R")

#function not returning anything.
ClusterPoolResults <- returnClusterpoolResult()

# Plot the pooled dotplot
Gene <- ClusterPoolResults$features.label
Cluster <- ClusterPoolResults$id
AvgExpScaled <- ClusterPoolResults$avg.exp.re.scaled
markers <- Gene %>% unique()

Plot <- ClusterPoolResults %>% 
  filter(features.label %in% markers) %>% 
    mutate(`% Expressing` = ClusterPoolResults$pct.exp) %>% 
      ggplot(aes(y=Gene, x = Cluster, color = AvgExpScaled, size = `% Expressing`)) + 
        geom_point() + 
        scale_size(range = c(0, 20)) +
        scale_color_viridis_c(option = "plasma") + 
        cowplot::theme_cowplot() + 
        theme(axis.title = element_text(size=20,face="bold"), legend.key.size=unit(1, "line")) +
        theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5, size=15)) +
        theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5, size=15))
Plot

# Save the image
# You will have to resize the Rstudio box
# or set the preferred width and height
# >>>> input required >>>>

save_image('PooledDotPlot',Plot, width = 5000, height = 2500)
