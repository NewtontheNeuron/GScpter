#### Cluster Plot
# This is a plot that shows the raw unlog scaled expression levels
# of the genes for each cluster for each cell with in each clusterpool
# The black bar indicates the median

# plot of expresson of every cell with color scale for the cluster
createClusterPlot <- function(cpdata){
  cpdata[[1]] %>%
    ggplot(mapping = aes(features.label, raw_counts, group = cluster)) + # Human had to use raw_counts instead of
    geom_jitter(position = position_jitterdodge(dodge.width = 1, # expm1(raw counts), which also changed the label
                                                jitter.width = 0), # from unlog-scaled to log-scaled
                aes(color = cluster)) + # Using raw_counts will work with mouse as well and is
    stat_summary(fun = median, geom = "crossbar", # Probably a safe bet, or you could set a parameter in the function
                 position = position_dodge(width = 1)) +
    #scale_color_viridis_d(option = "plasma") +
    scale_fill_hue(h = c(180, 270)) +
    scale_y_continuous(expand = c(0,0), limits = c(-1, 12)) + # human: had to change the limit
    labs(x = "Gene", y = "Raw log-scaled expression", # solution rounded 5% above the highest value make a
         color = "Cluster", title = names(cpdata)) + # setLimit function
    cowplot::theme_cowplot() + 
    theme(legend.position = "bottom",
          axis.title = element_text(size = 15, face = "bold"))
}

mainCluP <- function(cell_roster){
  # use plot arkv logic
  plot_arkv <- list()
  average <- list()
  # Loop through the cluster pools and generate the graph
  for (index in 1:length(cell_roster)) {
    plot_arkv[[names(cell_roster[index])]] <- createClusterPlot(cell_roster[index])
  }
  
  # Now add them together
  masterplot <- plot_arkv[[1]] + plot_arkv[[2]] #+ plot_arkv[[3]] + plot_arkv[[4]] #+
    #plot_arkv[[5]]
  # human: had to reduce this to 2 not 4 plots
  
  # Now save image
  save_image("ClusterPlot", masterplot, height = 2500, width = 4100)
}

# Run main to get the plots
mainCluP(cell_roster = cell_roster)
