library(rstudioapi)

#set working directory to the one this file is currently in
#setwd(dirname(getActiveDocumentContext()$path))

source("Pre_analysis_functions.R")

#run this function if you want to load the data
#load_data()

ClusterPoolResults <- returnClusterpoolResult()

# Plot the pooled dotplot
ClusterPoolResults$ClusterAndSubgroup <- paste(ClusterPoolResults$id, ClusterPoolResults$SubGroup) 
Gene <- ClusterPoolResults$features.label
Subgroup <- ClusterPoolResults$ClusterAndSubgroup
AvgExpScaled <- ClusterPoolResults$avg.exp.re.scaled
markers <- Gene %>% unique()

Plot <- ClusterPoolResults %>% 
  filter(features.label %in% markers) %>% 
    mutate(`% Expressing` = ClusterPoolResults$pct.exp) %>% 
      ggplot(aes(y=Gene, x = Subgroup, color = AvgExpScaled, size = `% Expressing`)) + 
        geom_point() + 
        #scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
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
# You will have to resize the Rstudio box
# or set the preferred width and height
save_image('PooledDotPlot', Plot)
