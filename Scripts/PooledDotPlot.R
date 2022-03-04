mainPDP <- function(ClusterPoolResults){
  
  # Plot the pooled dotplot
  ClusterPoolResults$ClusterAndSubgroup <- paste(ClusterPoolResults$id, ClusterPoolResults$SubGroup) 
  Gene <- ClusterPoolResults$features.label
  Subgroup <- ClusterPoolResults$ClusterAndSubgroup
  AvgExpScaled <- ClusterPoolResults$avg.exp.z.scaled
  markers <- Gene %>% unique()
  
  Plot <- ClusterPoolResults %>% 
    filter(features.label %in% markers) %>% 
    mutate(`% Expressing` = ClusterPoolResults$pct.exp) %>% 
    ggplot(aes(y=Gene, x = Subgroup, color = AvgExpScaled, size = `% Expressing`)) + 
    geom_point() + 
    scale_size(range = c(0, 20)) +
    scale_color_viridis_c(option = "plasma") + 
    cowplot::theme_cowplot() + 
    theme(axis.title = element_text(size=20,face="bold"), legend.key.size=unit(1, "line")) +
    theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0.5, size=15)) +
    theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5, size=15)) +
    theme(panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"))
  Plot
  
  # Save the image
  save_image('PooledDotPlot', Plot) 
}

#uncomment this function if you want to load the data
#RDSfile <- load_data()

#ClusterPoolResults <- createClusterPoolResults(RDSfile)

#mainPDP(ClusterPoolResults)
