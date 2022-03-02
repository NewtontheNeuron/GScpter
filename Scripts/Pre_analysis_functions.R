# Function for adding the file to the workspace and saving the workspace.
# This is only necessary outside Rstudio or if loading the data for
# this first time.
# >>>> input required >>>>

library(Seurat) #* **
library(ggplot2) #~ **
#library(patchwork) #~
#library(tidyverse) #*
#library(ggdendro) #~
#library(cowplot) # I do not think this is needed here
library(dplyr) #~ **
library(stringr) #~ **
#library(data.table) #~ **
#library(tibble) #~ **
#library(viridisLite)
library(Cairo) #* **
library(rstudioapi) #* **
#library(datasets) #*
#set working directory to the one this file is currently in
setwd(dirname(getActiveDocumentContext()$path))

#get JSON_Handler functions.
source("JSON_Handler.R")

load_data <- function() {

  #was newdata and savedata functions below
  print("Loading data into R... this might take a while.")
  filename <- file.choose()
  clean_neuron_object <<- readRDS(filename)
  save(clean_neuron_object, file = '../.RData')
  #was loadfile function
  # This is good but when newdata and savedata were run then
  # we do not need to run load because it would already have clean_neuron_object from
  # above.
  load(file = '../.RData') # We could also place .RData in the Data folder
  # We also need some form of error handling when the file does not exist so
  # that it prompts user the .RData file does not exist please ensure that it
  # is in the correct folder or do you wish to exit.
  print("RDS file loaded!") 
  print("Good to go!")
}

#get the list of features from JSON
features <- returnFeatures()
project_name <- returnProjectName()

# Remove the 'rna_' label on most genes in the Levine Dataset
# Clean the features into presentable labels without the 'rna_' tag.
features_no_key <- str_remove(features, 'rna_')

# Extract the dotplot information from Seurat into an object.
# The dotplot information involves the average expression and
# the percent expressed.
b <- DotPlot(clean_neuron_object, features = features)

clusterpool_names <- returnClusterpool_names()
clusterpool_subgroup <- returnClusterpool_subgroups()

ClusterPoolAll <- c()
ListByCluster <- list()
cell_roster <- list()

index <- 1
for (name in clusterpool_names){
  for (subgroup in clusterpool_subgroup){
    Clusterpool <- returnClusterpoolGenes(name, subgroup)
    
    #create master list of all gene names
    ClusterPoolAll <- c(ClusterPoolAll, Clusterpool)
    
    #create key name
    key <- paste(subgroup, name, sep = " ")
    
    #key value pair into ListByCluster
    ListByCluster[[key]] <- b$data[b$data$id %in% Clusterpool,]
    
    # First grab the indices of the features
    GeneIndicies <- which(unlist(clean_neuron_object@assays$RNA@data@Dimnames[1]) %in% features_no_key)
    
    # Then grab the indicies of the cells for the spcefic Clusterpool
    CellIndicies <- which(unlist(clean_neuron_object@active.ident) %in% Clusterpool)
    
    # Get the relevant data using the the indicies otained above
    raw_data <- GetAssayData(clean_neuron_object, assay = 'RNA', slot = 'data')[GeneIndicies, CellIndicies]
    
    # Dataframe with the cells and the associated clusters
    clusters_n_cells <- as.data.frame(unlist(clean_neuron_object@active.ident)[CellIndicies])
    
    # Now we turn raw_data into a data frame
    cell_roster[[key]] <- raw_data %>% # Using dplyr
      # Change the object to a data frame
      as.data.frame() %>% 
      # Make the row names an actual column, necessary for an acurate pivot
      mutate('features.label' = row.names(raw_data)) %>%
      # Pivot the data frame with the cells and raw counts in separte columns
      pivot_longer(cols = -features.label, names_to = 'cell.barcode', values_to = 'raw_counts') %>%
      # Add a column with the clusters
      mutate('cluster' = clusters_n_cells[match(cell.barcode, row.names(clusters_n_cells)), 1], .keep = "unused") %>%
      # Add a column with the id and subgroup
      mutate('id' = name, 'subgr' = subgroup)
    index <- index + 1
  }
}

ListbyClusterAll <- b$data[b$data$id %in% ClusterPoolAll,]

# Now to make a master roster of the cells
all_cell_roster <- data.table::rbindlist(cell_roster)

#------
PoolnShare <- function(id, subgr){
  # Start and define the data frame that will contain the results.
  # This step ensures that the values from the dot plot can be
  # Stored in the data frame in the correct place with the correct
  # class of variables.
  ClusterPoolResults <<- data.frame(avg.exp = numeric(),
                                    pct.exp = numeric(),
                                    features.plot = character(),
                                    id = character(),
                                    features.label = character(),
                                    SubGroup = character(),
                                    avg.std.err = numeric(),
                                    avg.lower = numeric(),
                                    avg.upper = numeric())
  ClusterPoolResults$features.plot <- as.character(ClusterPoolResults$features.plot)
  ClusterPoolResults$id <- as.character(ClusterPoolResults$id)
  ClusterPoolResults$features.label <- as.character(ClusterPoolResults$features.label)
  ClusterPoolResults$SubGroup <- as.character(ClusterPoolResults$SubGroup)

  # Adding the standard error function
  SE <- function(x) {
  sd(x) / sqrt(length(x))
  }

  # Average, percent, scale and save
  PoolAll <- function (cluster, identity, subgr){
    for(i in features){
      ListByGene <- cluster[str_detect(row.names(cluster), i), ] # Change this to filter by gene
      wfeature <- features_no_key[match(str_remove(i, 'rna_'), features_no_key)] # i.e. the working feature 
      # or the current feature in the loop

      # Now we compute the average expression and percent expressed
      # First grab the indicies of the working feature
      # The average expression and percent expressed is calculated
      # just like it is done in Seurat/utilities.R
      GeneIndicies <- match(wfeature, unlist(clean_neuron_object@assays$RNA@data@Dimnames[1]))
      # Then grab the indicies of the cells in the specific list of clusters
      CellIndicies <- which(unlist(clean_neuron_object@active.ident) %in% unique(ListByCluster[[subgr_index]]$id))
      # Get the relevant data using the the indicies otained above
      Transcript_exp <- GetAssayData(clean_neuron_object, assay = 'RNA', slot = 'data')[GeneIndicies, CellIndicies]
      # Exponentiate the data and then compute the mean
      raw.avg <- mean(expm1(Transcript_exp))
      # Now lets look at save some descriptive statistics for average expression
      avg.std.err <- SE(expm1(Transcript_exp))
      avg.lower <- raw.avg - avg.std.err
      avg.upper<- raw.avg + avg.std.err
      # Now the percent expressed does not require exponentiation
      raw.pct <- length(Transcript_exp[Transcript_exp > 0]) / length(Transcript_exp) * 100
      
      # Add a new row to ClusterPoolResults and place the values in the relevant collumn
      ClusterPoolResults <<- ClusterPoolResults %>% # This method works but is a bit slow
        add_row(avg.exp = raw.avg,                  # In the futrue I could just try adding
                pct.exp = raw.pct,                  # in the variables one data type at a time.
                features.plot = i,
                id = identity,
                features.label = wfeature,
                SubGroup = subgr,
                avg.std.err = avg.std.err,
                avg.lower = avg.lower,
                avg.upper = avg.upper)
      # This might no longer be needed because it is not even used
    }
    # This is because it is just saving the feature as the row name
    # which is like the dot plot output but not necessary,
    # we may even find a more efficent way of accessing the info
    # besides from a dot plot print out.
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
    c((ClusterPoolResults$avg.exp - colMeans(ClusterPoolResults[1]))/sd(ClusterPoolResults$avg.exp))
  # Then set the name
  names(ClusterPoolResults)[ncol(ClusterPoolResults)] <- 'avg.exp.z.scaled'
  return(ClusterPoolResults)
}



# Function for pooling the average expression and percent 
# expressed from the dotplot data. The function also re-scales avg.exp
# using a z-score system similar to Seurat's z-score system.
#
# TODO: Make this function take a vector of clusterpool
# all instead of indivdual 2 scores + 2 clusterpool names.

# Running the pool and share function

#create id vector and subgroups vector for poolnshare function.
id <- clusterpool_names
subgr <- clusterpool_subgroup
numberOfGroups <- length(id) * length(subgr)
numberOfGenes <- length(ClusterPoolAll)

ClusterPoolResults <- PoolnShare(id, subgr) # Error I get nothing in the ClusterPoolResults

#function to set the width and height of plot saved as an image into global values
resizeImage <- function(base_filename, numOfClusters, numOfGroups, numOfGenes, height, width){

  if (height == 1 && width == 1){
    if (base_filename == "Quad_BarPlot"){
      height <- 5000
      width <- 5000
    } else if (base_filename == "PooledDotPlot"){
      height <- 500 + (200 * ( numOfClusters / numOfGroups ))
      width <- 800 + (300 * numOfGroups)
    } else if (base_filename == "DotPlot"){
      print(numOfClusters)
      print(numOfGenes)
      height <- 1000 + (200 * numOfClusters / numOfGroups)
      width <- 1000 + (150 * numOfGenes)
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
      height <- 900 + (200 * ( numOfClusters / numOfGroups ))
    } else if (base_filename == "DotPlot"){
      height <- 500 + (200 * ( numOfClusters / numOfGroups ))
    } else {
      print("Not a recognized plot name. Setting Height of Plot to 5000px")
      print("If you want to add the plot to resizing, change resizeImage() function to include your plot name.")
      height <- 5000
    }
  } else if (width == 1){
    if (base_filename == "Quad_BarPlot"){
      width <- 5000
    } else if (base_filename == 'PooledDotPlot'){
      width <- 800 + (300 * numOfGroups)
    } else if (base_filename == "DotPlot"){
      width <- 1000 + (150 * numOfGenes)
    } else {
      print("Not a recognized plot name. Setting Width of Plot to 5000px")
      print("If you want to add the plot to resizing, change resizeImage() function to include your plot name.")
      width <- 5000
    }
  } 
  print(sprintf("%s plot height (in px): %2.f, and width (in px): %2.f", base_filename, height, width))
  return(c(height, width)) 
}

# Function for saving images with specific folder, filename, and date. 
save_image <- function(base_filename, Plot, height = 1, width = 1){
  #set working directory to output directory
  # If Output does not exist, create it.
  if (!dir.exists("../Output")){
    dir.create("../Output", showWarnings = FALSE)
  }
  # In any cass set the working directory
  setwd("../Output")

  #get proper calculated dimensions
  dimensions = resizeImage(base_filename, length(ClusterPoolResults$features.label) - 1, numberOfGroups, numberOfGenes, height, width)
  
  #get current date to put on graph title
  curr_date <- format(Sys.Date(),"%b_%d%_%Y" )

  dir.create(sprintf("%s", project_name), showWarnings = FALSE)
  #put all quad bar plots in another folder inside the project.
  if (substr(base_filename, 1,3) == "QB_"){
    dir.create(sprintf("%s/Quad_BarPlot", project_name), showWarnings = FALSE)
    filename <- sprintf("%s/Quad_BarPlot/%s_%s.jpg", project_name, base_filename, curr_date)
  } else {
    filename <- sprintf("%s/%s_%s.jpg", project_name, base_filename, curr_date)
  }
  
  #save plot to proper location as a jpg.
  ggsave(filename, plot = Plot, device = "jpg", height = dimensions[1], width = dimensions[2], units = "px", type = "cairo")
}

returnListByCluster <- function(){
  return(ListByCluster)
}

returnClusterpoolResult <- function(){
  return(ClusterPoolResults)
}

returnListbyClusterAll <- function(){
    ListbyClusterAll[, ncol(ListbyClusterAll) + 1] <- data.frame(features.label = features_no_key)
    return(ListbyClusterAll)
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

SE <- function(x) {
  sd(x) / sqrt(length(x))
}

