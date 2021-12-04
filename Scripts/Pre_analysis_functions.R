# Function for adding the file to the workspace and saving the workspace.
# This is only necessary outside Rstudio or if loading the data for this first time.
# >>>> input required >>>>

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
library(rstudioapi)
library(datasets)
#set working directory to the one this file is currently in
setwd(dirname(getActiveDocumentContext()$path))

#get JSON_Handler functions.
source("JSON_Handler.R")

load_data <- function(){
  RStudio <- toupper(readline(prompt = "Are you using RStudio to run this? "))
  if (RStudio != "Y"){
    loadedBefore <- toupper(readline(prompt = "Have you loaded the RDS data before? Y/N "))
    
    while (loadedBefore != "Y" & loadedBefore != "N"){
      print("Please enter either 'Y' or 'N'")
      loadedBefore <- toupper(readline(prompt = "Have you loaded the data before? Y/N "))
    }

    if (loadedBefore == "N"){
      #was newdata and savedata functions below
      print("Loading data into R... this might take a while.")
      filename <- file.choose()
      clean_neuron_object <<- readRDS(filename)
      save(clean_neuron_object, file = '.RData')
    } 
    
    #was loadfile function
    load(file = '.RData')
    print("RDS file loaded!")
  }
  
  print("Good to go!")
}

load_data()

#get the list of features from JSON
features <- returnFeatures()
project_name <- returnProjectName()

# Extract the dotplot information from Seurat into an object.
# The dotplot information involves the average expression and
# the percent expressed.
b <- DotPlot(clean_neuron_object, features = features)

# Start and define the data frame that will contain the results.
# This step ensures that the values from the dot plot can be 
# Stored in the data frame in the correct place with the correct
# character type.
ClusterPoolResults <<- data.frame(avg.exp=numeric(), pct.exp=numeric(), features.plot=character(), id=character(), avg.exp.scaled=numeric(), features.label=character())
ClusterPoolResults$features.plot <- as.character(ClusterPoolResults$features.plot)
ClusterPoolResults$id <- as.character(ClusterPoolResults$id)
ClusterPoolResults$features.label <- as.character(ClusterPoolResults$features.label)

#Get named clusterpools
#TODO: remove hard coded inhibitory/excitatory and DDH.
clusterpool_names <- c("inhibitory", "excitatory")
ClusterPoolAll <- c()

for (name in clusterpool_names){
  Clusterpool <- returnClusterpoolGenes(name, "DDH")
  ClusterPoolAll <- c(ClusterPoolAll, Clusterpool)
}

# Filtering the data set based on the clusters into
# 2 clusterpools and a clusterpool containing all of
# the clusters and their average expression and percent
# expressed data.

ListbyClusterAll <- b$data[b$data$id %in% ClusterPoolAll,]

#ListByCluster1 <- ListbyClusterAll %>%
#    filter(grepl("Excit-", id)) %>%
#    as_tibble()

ListByCluster1 <- b$data[b$data$id %in% grep("Excit-", ClusterPoolAll, value = TRUE),]
ListByCluster2 <- b$data[b$data$id %in% grep("Inhib-", ClusterPoolAll, value = TRUE),]

# Remove the 'rna_' label on most genes in the Levine Dataset
# Clean the features into presentable labels without the 'rna_' tag.
clean_label_list <- str_remove(features, 'rna_')
#------

PoolnShare <- function(cp1, cp2, id1, id2){
  
  #Defining the full Data Table, it will clear every time
  ClusterPoolResults <<- data.frame(avg.exp=numeric(), pct.exp=numeric(), features.plot=character(), id=character(), avg.exp.scaled=numeric(), features.label=character())
  ClusterPoolResults$features.plot <- as.character(ClusterPoolResults$features.plot)
  ClusterPoolResults$id <- as.character(ClusterPoolResults$id)
  ClusterPoolResults$features.label <- as.character(ClusterPoolResults$features.label)
  
  #For repeating function
  PoolAllRepeat <- 1
  
  #Average, Scale and save
  PoolAll <- function (cluster, identity){
    for(i in features){
      ListByGene <- cluster[str_detect(row.names(cluster), i), ]
      NewListItem <- data.frame(colMeans(ListByGene[1]), colMeans(ListByGene[2]), Col3=i, 
        Col4=identity, colMeans(ListByGene[5]), Col6=clean_label_list[match(str_remove(i, 'rna_'), clean_label_list)])
      NewListItem$Col3 <- as.character(NewListItem$Col3)
      NewListItem$Col4 <- as.character(NewListItem$Col4)
      NewListItem$Col6 <- as.character(NewListItem$Col6)
      ClusterPoolResults[nrow(ClusterPoolResults) + 1, ] <<- NewListItem
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
  PoolAll(ListByCluster1, id1)
  PoolAll(ListByCluster2, id2)
  
  #Rescaling average expression based on z scores of the full new dataset
  ClusterPoolResults[, ncol(ClusterPoolResults) + 1] <- 
    data.frame(avg.exp.re.scaled = (ClusterPoolResults[, 1] - colMeans(ClusterPoolResults[1]))/sd(ClusterPoolResults$avg.exp))
  return(ClusterPoolResults)
}


# Function for pooling the average expression and percent 
# expressed from the dotplot data. The function also re-scales avg.exp
# using a z-score system similar to Seurat's z-score system.
#
# TODO: Make this function take a vector of clusterpool all instead of indivdual 2 scores + 2 clusterpool names.

# Running the pool and share function
#poolnShare not working.
ClusterPool1 <- returnClusterpoolGenes(clusterpool_names[1], "DDH")
ClusterPool2 <- returnClusterpoolGenes(clusterpool_names[2], "DDH")

ClusterPoolResults <- PoolnShare(ClusterPool1, ClusterPool2 , clusterpool_names[1], clusterpool_names[2])

#str(ClusterPoolResults)

#function to set the width and height of plot saved as an image into global values
#Parameters: 
#   listOfClusters : vector of names of clusters
#   listOfGenes: a vector of names of genes used.
setHeightWidthImage <- function(listOfClusters, listOfGenes){
  numClusters <- length(listOfClusters) / length(listOfGenes)
  numGenes <- length(listOfGenes)
  
  if (numClusters > 0 && numClusters <= 5){
    plot_width <<- (numClusters * 210) + 1000
  } else if (numClusters > 5){
    plot_width <<- (numClusters * 160) + 1000
  }
  plot_height <<- (numGenes * 250) + 150
  
  
  message(sprintf("number of Clusters: %2.f", numClusters))
  message(sprintf("plot height (in px): %2.f", plot_width))
  
  message(sprintf("number of genes formatted: %2.f", numGenes))
  message(sprintf("plot width (in px) %2.f", plot_height))
}

# Function for saving images with specific folder,
# filename, and date. If the folder for the project name
# does not exist you will have to make it.
save_image <- function(base_filename, Plot,  width = 2000, height = 2000){
  dir.create(sprintf("../Output/%s", project_name), showWarnings = FALSE)
  curr_date <- format(Sys.Date(),"%b_%d%_%Y" )
  filename <- sprintf("../Output/%s/%s%s.png", project_name, base_filename, curr_date)
  ggsave(filename, plot = Plot, device = "png", height = height, width = width, units = "px", type = "cairo")
}

returnClusterpoolResult <- function(){
  return(ClusterPoolResults)
}

returnListbyClusterAll <- function(){
    ListbyClusterAll[, ncol(ListbyClusterAll) + 1] <- data.frame(features.label = clean_label_list)
    return(ListbyClusterAll)
}

SE <- function (x) {
  sd(x)/sqrt(length(x))
}

AvgExpPar <- function (x) {
  x %>% 
    group_by(features.plot) %>% 
    summarise(std.err = SE(avg.exp), avg.exp = mean(avg.exp), lower = avg.exp - std.err, upper = avg.exp + std.err)
}

PctExpPar <- function (x) {
  x %>% 
    group_by(features.plot) %>% 
    summarise(std.err = SE(pct.exp), pct.exp = mean(pct.exp), lower = pct.exp - std.err, upper = pct.exp + std.err)
}