# This is a script that creates a dot plot from CPR or clusterpoolresults data

mainPDP <- function(CPR, base.name = "PooledDotPlot", transp = F, height = NA, width = NA,
                    factor.order = c(), global_size = 30, add.label = F,
                    legend.margin = margin(), rm.labs = "", max.dot.size = 15, yieldplot = F, ...){
  
  Plot <- CPR %>%
    ggplot() + 
    geom_point(aes(y = features.label, x = group.label, color = avg.exp.scaled, size = pct.exp)) +
    {if (add.label) geom_label(aes(y = features.label, x = group.label, label = signif(avg.exp, digits = 2)))} +
    labs(x = "Group", y = "Gene", color = "Avg exp scaled", size = "% Expressing") +
    scale_size(range = c(0, max.dot.size)) +
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
          legend.margin = legend.margin,
          plot.background = element_rect(fill = ifelse(transp == T, "transparent", "white"))) +
    theme(...)

  Plot
  
  
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
  
  # Add guides
  Plot <- Plot +
    guides(color = guide_colorbar(barwidth = 1.2, barheight = 9, ticks = T,
                                  label.position = "left", title.position = "top"),
           size = guide_legend(label.position = "left", title.position = "top"))
  
  
  # External responding statement at the end of the function to either save it or not
  #
  # Save image function
  ifelse(yieldplot, return(Plot),
         save_image(base.name, Plot, height = ifelse(!is.na(height), height, 2400),
                    width = ifelse(!is.na(width), width, 3000)))
}
