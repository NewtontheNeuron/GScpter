#### Cluster Plot
# This is a plot that shows the raw unlog scaled expression levels
# of the genes for each cluster for each cell with in each clusterpool
# The black bar indicates the median

# plot of expresson of every cell with color scale for the cluster
mainCluP <- function(roster = all_cell_roster, pool.level = "1", transp = F,
                     dot.global.size = 15){
  any_overgrouped <- list.any(roster[unlist(extra_pool[[pool.level]])], is.list(.))
  if (any_overgrouped) {
    roster <- group_expand(roster, pool.level)
  }
  roster$group.label <- pool_level_paste(roster, extra_pool[[pool.level]])
  
  clusterplot <- roster %>%
    ggplot(mapping = aes(features.label, raw_counts, group = cluster)) + # Human had to use raw_counts instead of
    geom_jitter(position = position_jitterdodge(dodge.width = 1, # expm1(raw counts), which also changed the label
                                                jitter.width = 0), # from unlog-scaled to log-scaled
                aes(color = cluster)) + # Using raw_counts will work with mouse as well and is
    stat_summary(fun = mean, geom = "crossbar", # Probably a safe bet, or you could set a parameter in the function
                 position = position_dodge(width = 1)) +
    scale_color_viridis_d() +
    scale_y_continuous(expand = c(0,0),
                       limits = c(0, p5max(roster$raw_counts))) + # human: had to change the limit
    labs(x = "Gene", y = "Log-scaled expression", color = "Cluster") +
    facet_wrap(~group.label) +
    cowplot::theme_cowplot() + 
    theme(axis.title = element_text(size = dot.global.size, face = "bold"),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5,
                                     size=dot.global.size, color = "black"),
          axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5,
                                     size = dot.global.size, color = "black"),
          legend.key.size = unit(1, "line"),
          legend.text = element_text(size = dot.global.size/1.75),
          legend.title = element_text(size = dot.global.size, angle = 90),
          legend.box = "horizontal",
          legend.position = "bottom",
          legend.spacing.x = unit(0.5, "line"),
          plot.background = element_rect(fill = ifelse(transp == T, "transparent", "white")))

  save_image("ClusterPlot", clusterplot, height = 2500, width = 5000)
}
