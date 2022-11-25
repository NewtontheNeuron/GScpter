mainPDP <- function(CPR, transp = F, height = NA, width = NA,
                    factor.order = c()){
  global_size <- 20

  CPR$ClusterAndSubgroup <- paste(CPR$id, CPR$subgr)
  # TODO: this should be done at the Cluster pool results level

  
  Plot <- CPR %>%
    mutate(ClusterAndSubgroup = ifelse(length(factor.order) > 0,
                       fct_relevel(ClusterAndSubgroup, factor.order),
                       ClusterAndSubgroup)) %>%
    ggplot(aes(y = features.label, x = ClusterAndSubgroup, color = avg.exp.scaled, size = pct.exp)) + 
    geom_point() +
    labs(x = "Group", y = "Gene", color = "Avg exp scaled", size = "% Expressing") +
    scale_size(range = c(0, 20)) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
    scale_color_viridis_c(option = "plasma") + 
    cowplot::theme_cowplot() + 
    theme(axis.title = element_text(size = global_size, face = "bold"),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5,
                                     size=global_size, color = "black"),
          axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5,
                                     size = global_size, color = "black"),
          legend.key.size = unit(1, "line"),
          legend.text = element_text(size = global_size/1.75),
          legend.title = element_text(size = global_size, angle = 90),
          legend.box = "horizontal",
          legend.spacing.x = unit(0.5, "line"),
          plot.background = element_rect(fill = ifelse(transp == T, "transparent", "white")))

  Plot
  
  # Save the image
  save_image('PooledDotPlot', Plot, height = ifelse(!is.na(height), height, 2400),
             width = ifelse(!is.na(width), width, 3000)) 
}
