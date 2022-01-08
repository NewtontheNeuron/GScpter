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
  RStudio <- toupper(readline(prompt = "Are you using RStudio to run this? Y/N "))
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
      save(clean_neuron_object, file = '../.RData')
    } 
    
    #was loadfile function
    # This is good but when newdata and savedata were run then
    # we do not need to run load because it would already have clean_neuron_object from
    # above.
    if (loadedBefore == "Y"){
      load(file = '../.RData') # We could also place .RData in the Data folder
      # We also need some form of error handling when the file does not exist so that it prompts
      # user the .RData file does not exist please ensure that it is in the correct folder or
      # do you wish to exit
      print("RDS file loaded!") 
    }
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

clusterpool_names <- returnClusterpool_names()
clusterpool_subgroup <- returnClusterpool_subgroups()

ClusterPoolAll <- c()
ListByCluster <- list()

index <- 1
for (name in clusterpool_names){
  for (subgroup in clusterpool_subgroup){
    Clusterpool <- returnClusterpoolGenes(name, subgroup)
    #create master list of all gene names
    ClusterPoolAll <- c(ClusterPoolAll, Clusterpool)
    #create key name
    key <- paste("clusterpool", index, sep = "")
    #key value pair into ListByCluster
    ListByCluster[[key]] <- b$data[b$data$id %in% Clusterpool,]
    index <- index + 1
  }
}

ListbyClusterAll <- b$data[b$data$id %in% ClusterPoolAll,]

# Remove the 'rna_' label on most genes in the Levine Dataset
# Clean the features into presentable labels without the 'rna_' tag.
clean_label_list <- str_remove(features, 'rna_')
#------
PoolnShare <- function(id, subgr){
  # Start and define the data frame that will contain the results.
  # This step ensures that the values from the dot plot can be 
  # Stored in the data frame in the correct place with the correct
  # character type.
  ClusterPoolResults <<- data.frame(avg.exp=numeric(), pct.exp=numeric(), features.plot=character(), id=character(), avg.exp.scaled=numeric(), features.label=character())
  ClusterPoolResults$features.plot <- as.character(ClusterPoolResults$features.plot)
  ClusterPoolResults$id <- as.character(ClusterPoolResults$id)
  ClusterPoolResults$features.label <- as.character(ClusterPoolResults$features.label)
  ClusterPoolResults$SubGroup <- as.character(ClusterPoolResults$SubGroup)
  
  #For repeating function
  PoolAllRepeat <- 1
  
  #Average, Scale and save
  PoolAll <- function (cluster, identity, subgr){
    for(i in features){
      ListByGene <- cluster[str_detect(row.names(cluster), i), ]
      NewListItem <- data.frame(colMeans(ListByGene[1]), colMeans(ListByGene[2]), Col3=i, 
        Col4=identity, colMeans(ListByGene[5]), Col6=clean_label_list[match(str_remove(i, 'rna_'), clean_label_list)], SubGroup= subgr)
      NewListItem$Col3 <- as.character(NewListItem$Col3)
      NewListItem$Col4 <- as.character(NewListItem$Col4)
      NewListItem$Col6 <- as.character(NewListItem$Col6)
      NewListItem$SubGroup <- as.character(NewListItem$SubGroup)
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
  for (subgr_index in 1:length(subgr)){
    for (id_index in 1:length(id)){
      PoolAll(ListByCluster[[subgr_index]], id[[id_index]], subgr[[subgr_index]])
    }
  }

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

#create id vector and subgroups vector for poolnshare function.
id <- clusterpool_names
subgr <- clusterpool_subgroup
numberOfGroups <- length(id) * length(subgr)

ClusterPoolResults <- PoolnShare(id, subgr)

#function to set the width and height of plot saved as an image into global values
resizeImage <- function(base_filename, numOfClusters, numOfGroups, height, width){

  if (height == 1 && width == 1){
    if (base_filename == "Quad_BarPlot"){
      height <- 5000
      width <- 5000
    } else if (base_filename == "PooledDotPlot"){
      print(numOfClusters)
      print(numOfGroups)
      height <- 500 + (200 * ( numOfClusters / numOfGroups ))
      width <- 800 + (300 * numOfGroups)
    } else if (base_filename == "DotPlot"){
      height <- 5000
      width <- 5000
    } else {
      print("Not a recognized plot name. Setting Height and Width of Plot to 5000px * 5000px")
      print("If you want to add the plot to resizing, change resizeImage() function.")
      height <- 5000
      width <- 5000
    }
  } else if (height == 1){
    if (base_filename == "Quad_BarPlot"){
      height <- 5000
    } else if (base_filename == "PooledDotPlot"){
      height <- 500 + (200 * ( numOfClusters / numOfGroups ))
    } else if (base_filename == ""){
      height <- 5000
    } else {
      print("Not a recognized plot name. Setting Height of Plot to 5000px")
      print("If you want to add the plot to resizing, change resizeImage() function to include your plot name.")
      height <- 5000
    }
  } else if (width == 1){
    if (base_filename == "Quad_BarPlot"){
      width <- 5000
    } else if (base_filename == 'PooledDotPlot'){
      print("hi")
      width <- 800 + (300 * numOfGroups)
    } else if (base_filename == ""){
      width <- 5000
    } else {
      print("Not a recognized plot name. Setting Width of Plot to 5000px")
      print("If you want to add the plot to resizing, change resizeImage() function to include your plot name.")
      width <- 5000
    }
  } 
  print(sprintf("%s plot height (in px): %2.f, and width (in px): %2.f", base_filename, height, width))
  return(c(height, width)) 
}

# Function for saving images with specific folder,
# filename, and date. If the folder for the project name
# does not exist you will have to make it.
save_image <- function(base_filename, Plot, height = 1, width = 1){
  dimensions = resizeImage(base_filename, length(ClusterPoolResults$features.label) - 1, numberOfGroups, height, width)
  dir.create(sprintf("../Output/%s", project_name), showWarnings = FALSE)
  curr_date <- format(Sys.Date(),"%b_%d%_%Y" )
  filename <- sprintf("../Output/%s/%s%s.jpg", project_name, base_filename, curr_date)
  ggsave(filename, plot = Plot, device = "jpg", height = dimensions[1], width = dimensions[2], units = "px", type = "cairo")
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
