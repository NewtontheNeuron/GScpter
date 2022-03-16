mainPDP <- function(CPR){
  
  # Plot the pooled dotplot
  CPR$ClusterAndSubgroup <- paste(CPR$id, CPR$SubGroup) 
  Gene <- CPR$features.label
  Subgroup <- CPR$ClusterAndSubgroup
  AvgExpScaled <- CPR$avg.exp.z.scaled
  markers <- Gene %>% unique()
  
  Plot <- CPR %>% 
    filter(features.label %in% markers) %>% 
    mutate(`% Expressing` = CPR$pct.exp) %>% 
    ggplot(aes(y=Gene, x = Subgroup, color = AvgExpScaled, size = `% Expressing`)) + 
    geom_point() + 
    scale_size(range = c(0, 20)) +
    scale_color_viridis_c(option = "plasma") + 
    cowplot::theme_cowplot() + 
    theme(axis.title = element_text(size = 20, face = "bold"),
          legend.key.size=unit(1, "line"),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5, size=15)) + # changed -45 angle to 0
    theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5, size=15)) +
    theme(panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"))
  Plot
  
  # Save the image
  save_image('PooledDotPlot', Plot) 
}
