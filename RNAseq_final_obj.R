# Specify the libraries that are used.

library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(ggdendro)
library(cowplot)
library(patchwork)
library(dplyr)
library(stringr)
library(data.table)
library(tibble)
library(viridisLite)
library(Cairo)

# Function for adding the file to the workspace and saving the workspace.
# This is only necessary outside Rstudio. Rstudio
# does the same thing. The first function adds new
# data into the working environment (workspace).
# The next function saves the data to the .RData file.
# The third function loads the .RData file.
# You do not need to add new data and save it twice.
# Most of the time you will only need load_data(), the third function.
# >>>> input required >>>>
using_RStudio <- toupper(readline(prompt = "Have you loaded the data before? Y/N "))

while (using_RStudio != "Y" & using_RStudio != "N"){
  print("Please enter either 'Y' or 'N'")
  using_RStudio <- readline(prompt = "HAve you loaded the data before? Y/N ")
}

newdata <- function(){
  filename <- file.choose()
  clean_neuron_object <<- readRDS(filename)
}
savedata <- function(){
  save(clean_neuron_object, file = '.RData')
}
load_data <- function(){
  load(file = '.RData')
}


if (using_RStudio == "N"){
  newdata() # Do this once
  #savedata() # Do this once
}

load_data()

# Setting the genes to be investigated
# >>>> input required >>>>
features <- c("rna_Grin1", "rna_Grin2a", "rna_Grin2b", "rna_Grin2c", "rna_Grin2d", "rna_Grin3a", "rna_Grin3b")

# Extract the dotplot information from Seurat into an object.
# The dotplot information involves the average expression and
# the percent expressed.
b <- DotPlot(clean_neuron_object, features = features)

# Start and define the data frame that will contain the results.
# This step ensures that the values from the dot plot can be 
# Stored in the data frame in the correct place with the correct
# character type.
##ClusterPoolResults <- data.frame(avg.exp=numeric(), pct.exp=numeric(), features.plot=character(), id=character(), avg.exp.scaled=numeric(), features.label=character())
##ClusterPoolResults$features.plot <- as.character(ClusterPoolResults$features.plot)
##ClusterPoolResults$id <- as.character(ClusterPoolResults$id)
##ClusterPoolResults$features.label <- as.character(ClusterPoolResults$features.label)

# The list of clusters that are included in each group
# I could add the functionality in the future to compare an 'n'
# amount of clusterpools.
# >>>> input required >>>>
ClusterPool1 <- c('Excit-5', 'Excit-6', 'Excit-20', 'Excit-21', 'Excit-22', 'Excit-23', 'Excit-24', 'Excit-25', 'Excit-26', 'Excit-27', 'Excit-29', 'Excit-30', 'Excit-31', 'Excit-32', 'Excit-34', 'Excit-35', 'Excit-36')
ClusterPool2 <- c('Inhib-3', 'Inhib-6', 'Inhib-8', 'Inhib-12', 'Inhib-14', 'Inhib-15', 'Inhib-16', 'Inhib-18', 'Inhib-19', 'Inhib-20', 'Inhib-21')
ClusterPool3 <- c("Excit-01", "Excit-02", "Excit-03", "Excit-08", "Excit-09", "Excit-10", "Excit-12", "Excit-14", "Excit-15", "Excit-16", "Excit-18", "Excit-04", "Excit-05", "Excit-13", "Excit-19")
ClusterPool4 <- c("Inhib-01", "Inhib-02", "Inhib-03", "Inhib-04", "Inhib-05", "Inhib-06", "Inhib-07", "Inhib-09", "Inhib-10","Inhib-11","Inhib-12", "Inhib-13")
# First I make the rest of the code take other thngs
ClusterPoolAll <- c(ClusterPool1, ClusterPool2, ClusterPool3, ClusterPool4)

# Include a way to label clusterpool1 and 2 for example
# Excitatory vs. Inhibitory. Place them in the order they
# appear above. For instance the cluster pool 1 name should
# come first.
# >>>> input required >>>>
clusterpool_names <- c('Excitatory', 'Inhibitory')
clusterpool_subgroup <- c('SDH', 'DDH', 'SDH', 'DDH')

# Set project name
# >>>> input required >>>>
project_name <- readline(prompt = "Set your project name (use DDH for now): ")

# Filtering the data set based on the clusters into
# 2 clusterpools and a clusterpool containing all of
# the clusters and their average expression and percent
# expressed data.
ListByCluster1 <- b$data[b$data$id %in% ClusterPool1,]
ListByCluster2 <- b$data[b$data$id %in% ClusterPool2,]
ListByCluster3 <- b$data[b$data$id %in% ClusterPool3,]
ListByCluster4 <- b$data[b$data$id %in% ClusterPool4,]
ListbyClusterAll <- b$data[b$data$id %in% ClusterPoolAll,]

# Remove the 'rna_' label on most genes in the Levine Dataset
# Clean the features into presentable labels without the 'rna_' tag.
clean_label_list <- str_remove(features, 'rna_')

# -------------------------------------------------------------------
# Function for pooling the average expression and percent 
# expressed from the dotplot data. The function also re-scales avg.exp
# using a z-score system similar to Seurat's z-score system.
PoolnShare <- function(id1, id2, subgr1, subgr2, subgr3, subgr4){
  
  #Defining the full Data Table, it will clear every time
  ClusterPoolResults <- data.frame(avg.exp=numeric(), pct.exp=numeric(), features.plot=character(), id=character(), avg.exp.scaled=numeric(), features.label=character(), subgr=character())
  ClusterPoolResults$features.plot <- as.character(ClusterPoolResults$features.plot)
  ClusterPoolResults$id <- as.character(ClusterPoolResults$id)
  ClusterPoolResults$features.label <- as.character(ClusterPoolResults$features.label)
  ClusterPoolResults$subgr <- as.character(ClusterPoolResults$subgr)
  
  #Separates the data based on the cluster pools
  #ListByCluster1 <- b$data[b$data$id %in% cp1,]
  #ListByCluster2 <- b$data[b$data$id %in% cp2,]
  
  #For repeating function
  PoolAllRepeat <- 1

  #Average, Scale and save
  PoolAll <- function (cluster, identity, subgr){
    for(i in features){
      ListByGene <- cluster[str_detect(row.names(cluster), i), ]
      NewListItem <- data.frame(colMeans(ListByGene[1]), colMeans(ListByGene[2]), Col3=i, 
        Col4=identity, colMeans(ListByGene[5]), Col6=clean_label_list[match(str_remove(i, 'rna_'), clean_label_list, SubGroup=subgr)])
      NewListItem$Col3 <- as.character(NewListItem$Col3)
      NewListItem$Col4 <- as.character(NewListItem$Col4)
      NewListItem$Col6 <- as.character(NewListItem$Col6)
      NewListItem$SubGroup <- as.character(NewListItem$SubGroup)
      ClusterPoolResults[nrow(ClusterPoolResults) + 1, ] <<- NewListItem
      # The purpose of this 
      if(PoolAllRepeat == 1){
        row.names(ClusterPoolResults)[nrow(ClusterPoolResults)] <<- i
      } else {
        row.names(ClusterPoolResults)[nrow(ClusterPoolResults)] <<- paste(i, PoolAllRepeat)
      }
    }
    PoolAllRepeat <<- PoolAllRepeat + 1
  }
  
  # Running the function PoolAll, (If there were extra pools put them here before the 'Rescale' bellow)
  # Turn this into a for loop
  PoolAll(ListByCluster1, id1, subgr1)
  PoolAll(ListByCluster2, id2, subgr2)
  PoolAll(ListByCluster3, id1, subgr3)
  PoolAll(ListByCluster4, id2, subgr4)
  
  #Rescaling average expression based on z scores of the full new dataset
  ClusterPoolResults[, ncol(ClusterPoolResults) + 1] <- 
    data.frame(avg.exp.re.scaled = (ClusterPoolResults[, 1] - colMeans(ClusterPoolResults[1]))/sd(ClusterPoolResults$avg.exp))
  return(ClusterPoolResults)
}

# Running the pool and share function
ClusterPoolResults <- PoolnShare(clusterpool_names[1], clusterpool_names[2], clusterpool_subgroup[1], clusterpool_subgroup[2], clusterpool_subgroup[3], clusterpool_subgroup[4])

# Function for saving images with specific folder,
# filename, and date. If the folder for the project name
# does not exist you will have to make it.
save_image <- function(base_filename, set_height = 2500, set_width = 1500){
  dir.create(project_name)
  curr_date <- format(Sys.Date(),"%b_%d%_%Y" )
  filename <- sprintf("%s/%s%s%s.png", project_name, base_filename, project_name, curr_date)
  ggsave(filename, plot = Plot, device = "png", height = set_height, width = set_width, units = "px", type = "cairo")
}

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
# or set the prefered width and height
# >>>> input required >>>>
save_image('PooledDotPlot')

# --------------------------------------------------
# Ploting all the relevant clusters form the data by
# using ClusterPoolAll
ListbyClusterAll[, ncol(ListbyClusterAll) + 1] <- data.frame(features.label = clean_label_list)
Gene <- ListbyClusterAll$features.label
Cluster <- ListbyClusterAll$id
AvgExpScaled <- ListbyClusterAll$avg.exp.scaled
markers <- Gene %>% unique()

Plot <- ListbyClusterAll %>% filter(features.label %in% markers) %>% 
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
save_image('DotPlot')

# -----------------------------------------------
# Bargraphs for the average expression and percent
# expressed data, and also for both clusterpools
SE <- function (x) {
  sd(x)/sqrt(length(x))
}
AvgExpPar <- function (x) {
  x %>% group_by(features.plot) %>% summarise(std.err = SE(avg.exp), avg.exp = mean(avg.exp), lower = avg.exp - std.err, upper = avg.exp + std.err)
}
PctExpPar <- function (x) {
  x %>% group_by(features.plot) %>% summarise(std.err = SE(pct.exp), pct.exp = mean(pct.exp), lower = pct.exp - std.err, upper = pct.exp + std.err)
}

Plot_4_Bar <- function (avg1, pct1, avg2, pct2, c1, c2){
  n_feat <- "rna_Grin2c" #length(features) + 0.1
  
  Bween_pool <- function(method){
    ListbyClusterAll <- b$data[b$data$id %in% ClusterPoolAll,]
    GeneStatResults <- data.frame(t=numeric(), df=numeric(), p.value=numeric(), name=character())
    GeneStatResults$name <- as.character(GeneStatResults$name)
    for (i in features){
      ListbyGene <- ListbyClusterAll[str_detect(row.names(ListbyClusterAll), i), ]
      ListByCluster1 <- ListbyGene[ListbyGene$id %in% ClusterPool1,]
      ListByCluster2 <- ListbyGene[ListbyGene$id %in% ClusterPool2,]
      if (method == "avg"){
        ClusterStatResults <- t.test(ListByCluster1$avg.exp, ListByCluster2$avg.exp)
      } else if (method == "pct"){
        ClusterStatResults <- t.test(ListByCluster1$pct.exp, ListByCluster2$pct.exp)
      }
      ClusterStatResults <- append(ClusterStatResults, i, after = length(ClusterStatResults))
      x <- ClusterStatResults[[11]]
      NewStatResult <- data.frame(t=ClusterStatResults$statistic, df=ClusterStatResults$parameter, p.value=ClusterStatResults$p.value, name=x)
      NewStatResult$name <- as.character(NewStatResult$name)
      GeneStatResults[nrow(GeneStatResults) + 1, ] <- NewStatResult
    }
    GeneStatResults[ , "p.value.adj"] <- as.numeric()
    for(i in 1:length(GeneStatResults[,ncol(GeneStatResults)])) {
      GeneStatResults[i ,ncol(GeneStatResults)] <- (GeneStatResults[i,]$p.value)*nrow(GeneStatResults)
    }
    return(GeneStatResults)
  }
  #rsource is residing source, i.e. where the significance star resides
  Plot_Bween <- function(wplot, pvalues, method, rsource){
    for(i in 1:nrow(pvalues)){
      Ext_gene <- rsource[rsource$features.plot %in% pvalues$name[i],]
      if(method == "avg"){
        value <- apply(Ext_gene[1], 2, max, na.rm=TRUE) + (apply(Ext_gene[1], 2, max, na.rm=TRUE)*0.10)
      } else if(method == "pct"){
        value <- apply(Ext_gene[2], 2, max, na.rm=TRUE) + (apply(Ext_gene[2], 2, max, na.rm=TRUE)*0.10)
      }
      if(pvalues$p.value[i] < 0.05){
        stars <- "*"
        if(pvalues$p.value[i] < 0.01){
          stars <- "**"
          if(pvalues$p.value[i] < 0.001){
            stars <- "***"
          }
        }
      }
      if(pvalues$p.value.adj[i] < 0.05){
        wplot <- wplot + annotate("text", label = stars, x = pvalues$name[i], y = value, size = 6, colour = "black")
      }
      else {
        #wplot <- wplot + annotate("text", label = "ns", x = pvalues$name[i], y = value, size = 3, colour = "black")
      }
    }
    return(wplot)
  }
  
  Run_ANOVA <- function (cluster, method, anova.p.val=NULL) {
    if(method == "avg"){
      ANOVA <- aov(avg.exp ~ features.plot, data = cluster)
    } else if (method == "pct"){
      ANOVA <- aov(pct.exp ~ features.plot, data = cluster)
    }
    if (isTRUE(anova.p.val)){
      anova <- summary(ANOVA)[[1]][["Pr(>F)"]][1]
      return(signif(anova, digits = 2))
    } else if(is.null(anova.p.val)){
      
    }
    Pthoc <- dunnettT3Test(ANOVA)$p.value
    return(Pthoc)
  }
  
  Loop_Cell <- function (x){
    pvalues <- data.frame(p.value=numeric(), x=character(), xend=character())
    pvalues$x <- as.character(pvalues$x)
    pvalues$xend <- as.character(pvalues$xend)
    for(i in 1:ncol(x)){
      for (j in 1:nrow(x)){
        Working_Cell <- x[j,i]
        if (is.na(Working_Cell) == FALSE && Working_Cell > 0.05){
          wpvalues <- data.frame(p.value=Working_Cell, x=colnames(x)[i], xend=rownames(x)[j])
          wpvalues$x <- as.character(wpvalues$x)
          wpvalues$xend <- as.character(wpvalues$xend)
          pvalues[nrow(pvalues) + 1, ] <- wpvalues
        }
      }
    }
    return(pvalues)
  }
  
  Position_ANOVA <- function(source, method){
    if(method == "avg"){
      pvpos <- apply(source[3], 2, max, na.rm=TRUE) - (apply(source[3], 2, max, na.rm=TRUE)*0.05)
    } else if (method == "pct"){
      increment <- 5
      pvpos <- apply(source[3], 2, max, na.rm=TRUE) - (apply(source[3], 2, max, na.rm=TRUE)*0.05)
    }
    return(pvpos)
  }
  
  highest_mean <- function(source, method, wplot, pvalues){
    if(method == "avg"){
      increment <- 1
      start_point <- apply(source[1], 2, max, na.rm=TRUE) + (apply(source[1], 2, max, na.rm=TRUE)*0.05)
    } else if (method == "pct"){
      increment <- 5
      start_point <- apply(source[2], 2, max, na.rm=TRUE) + (apply(source[2], 2, max, na.rm=TRUE)*0.05)
    }
    for(i in 1:nrow(pvalues)){
      value <- start_point + (i - 1)*increment
      wplot <- wplot + annotate("segment", x = pvalues$x[i], xend = pvalues$xend[i], y = value, yend = value, colour = "black")
    }
    return(wplot)
  }
  
  exp.def <- function (source, method) {
    if(method == "avg"){
      source$avg.exp
    } else if (method == "pct"){
      source$pct.exp
    }
  }
  avg.lab <- "Average Expression"
  pct.lab <- "%Expressed"
  Plot_details <- function (a,b,c,d,e,method){
    ggplot(a, aes(features.plot, d, fill = factor(features.plot))) + 
      #theme(legend.position = "none") + 
      geom_col(color="black", show.legend = FALSE) + 
      scale_fill_viridis_d(option = "plasma") + 
      #scale_fill_discrete(fill = factor(a$features.plot), guide=FALSE) +
      geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width=0.9), width = 0.50) +
      geom_point(data = b, y = c, fill = "white", color="black") +
      scale_y_continuous(expand = c(0,0), limits = c(0, NA)) + 
      labs(x="Gene", y=e) +
      cowplot::theme_cowplot() + 
      scale_x_discrete(labels=clean_label_list) +
      theme(axis.title = element_text(size=20,face="bold")) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=20)) +
      theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=20)) +
      annotate("text", label=paste("p < ", Run_ANOVA(b, method, anova.p.val = TRUE)), x = 4.5, y = Position_ANOVA(a, method), size = 8.5, color = "black")
  }
  #Use 1 function to run,correct and draw t tests. Use 5% less of the larget value to plcae the anova pvalue. Export the post hocs as a table.
  Plot_1 <- Plot_details(avg1, c1, exp.def(c1, "avg"), exp.def(avg1, "avg"), avg.lab, "avg") 
  #Plot_1 <- highest_mean(c1, "avg", Plot_1, Loop_Cell(Run_ANOVA(c1, "avg")))
  Plot_2 <- Plot_details(avg2, c2, exp.def(c2, "avg"), exp.def(avg2, "avg"), avg.lab, "avg")
  Plot_2 <- Plot_Bween(Plot_2, Bween_pool("avg"), "avg", c2)
  #Plot_2 <- highest_mean(c2, "avg", Plot_2, Loop_Cell(Run_ANOVA(c2, "avg")))
  Plot_3 <- Plot_details(pct1, c1, exp.def(c1, "pct"), exp.def(pct1, "pct"), pct.lab, "pct")
  #Plot_3 <- highest_mean(c1, "pct", Plot_3, Loop_Cell(Run_ANOVA(c1, "pct")))
  Plot_4 <- Plot_details(pct2, c2, exp.def(c2, "pct"), exp.def(pct2, "pct"), pct.lab, "pct")
  Plot_4 <- Plot_Bween(Plot_4, Bween_pool("pct"), "pct", c2)
  #Plot_4 <- highest_mean(c2, "pct", Plot_4, Loop_Cell(Run_ANOVA(c2, "pct")))
  Plot <- Plot_1 + Plot_2 + Plot_3 + Plot_4
  return(Plot)
}
Plot <- Plot_4_Bar(AvgExpPar(ListByCluster1), PctExpPar(ListByCluster1), AvgExpPar(ListByCluster2), PctExpPar(ListByCluster2), ListByCluster1, ListByCluster2)
Plot

# Save the image
# You will have to resize the Rstudio box
# or set the prefered width and height
# >>>> input required >>>>
save_image('BarPlot')
