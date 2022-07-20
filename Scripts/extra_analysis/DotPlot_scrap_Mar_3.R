library(rstudioapi)

#set working directory to the one this file is currently in
setwd(dirname(getActiveDocumentContext()$path))
source("Pre_analysis_functions.R")

main <- function(){
  
  # Start a plot arkv
  plot_arkv <- list()
  
  # Let us loop through list by cluster
  for (wrk_cp in 1:length(ListByCluster)) {
    Gene <- ListByCluster[[wrk_cp]]$features.label
    Cluster <- ListByCluster[[wrk_cp]]$id
    AvgExpScaled <- ListByCluster[[wrk_cp]]$avg.exp.scaled
    markers <- Gene %>% unique()
    
    plt <- ListByCluster[[wrk_cp]] %>% 
      filter(features.label %in% markers) %>% 
      mutate(`% Expressing` = ListByCluster[[wrk_cp]]$pct.exp) %>% 
      ggplot(aes(y=Gene, x = Cluster, color = AvgExpScaled, size = `% Expressing`)) +
      geom_point() + 
      scale_size(range = c(0, 25)) +
      scale_color_viridis_c(option = "plasma") + 
      cowplot::theme_cowplot() + 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 25), 
            axis.title = element_text(size = 25, face="bold"), legend.key.size = unit(2.1, "line"), 
            legend.text = element_text(size = 17), legend.title = element_text(size = 20)) +
      theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size = 25)) + 
      theme(panel.background = element_rect(fill = "white"),
            plot.background = element_rect(fill = "white"))
    
    # Save it to the plot arkv
    plot_arkv[[wrk_cp]] <- plt
  }
  
  # Now create the plot and stack them from top to bottom
  Plot <- ggpubr::ggarrange(Plot1, Plot2, ncol = 1, nrow = 2, legend = "bottom") +
    theme(legend.key = element_rect(fill = "white"),
          legend.text = element_text(color = "white"),
          legend.title = element_text(color = "white")) +
    guides(color = guide_legend(override.aes = list(color = NA)))
  
  # Use patchwork
  Plot <- Plot1 + Plot2 + plot_layout(nrow = 2, guides = "collect")
  
  # Save the image
  # You will have to resize the Rstudio box
  # or set the prefered width and height
  #width = 10000 for 4 subgroups, 4500 for 2 subgroups, width = 2250
  save_image('DotPlot', Plot, height = 6000, width = 7000)
}

#run this function if you want to load the data
#load_data()

main()

