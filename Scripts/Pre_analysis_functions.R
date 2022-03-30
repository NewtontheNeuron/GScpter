#### Pre-analysis functions
# This script contains functions that are used by multiple graphing scripts
# to extract information from seurat objects and process the information
# tidy the information, and perform basic calculations.

# library(Seurat) #* **
# library(ggplot2) #~ **
# library(dplyr) #~ **
# library(Cairo) #* **
# library(rstudioapi) #* **

# #get the list of features from JSON
# features <- returnFeatures()
# project_name <- returnProjectName()

createCellRoster <- function(clean_neuron_object){
  # Get the presentable genelabels by removing the 'rna_' tag
  features_no_key <- returnCleanLabelList()
  
  # Get the groups and subgroups from the json file
  groups <- returnClusterpool_names()
  subGroups <- returnClusterpool_subgroups()
  
  # TODO: need to find a way to get the correct formating for the genes
  # for example GRIN1 vs. Grin1
  # Now merging meta data and active ident and finding where the clusters are stored
  # First merge @meta.data and @active.ident
  meta_ident <<- as.data.frame(clean_neuron_object@meta.data) %>%
    mutate(cells = row.names(as.data.frame(clean_neuron_object@active.ident)),
           id_plus = as.data.frame(clean_neuron_object@active.ident)[[1]])
  # Next look for the first cluster pool within each collumn
  all_clusters <- returnAllClusters()
  # Then get the names of the collumns that have at least one of the
  # clusters in the clusterpool
  clusters_location <- names(which.max(apply(
    meta_ident, 
    FUN = function (x) sum(x %in% all_clusters), 
    MARGIN = 2))) # Just take any one of the culumns or names (there may be multiple)
  # Now it can be used further down to get the index for the cells i.e. CellIndicies
  # TODO: add the ability to introduce secondary identifiers
  # TODO: turn this into a function

  # Create empty list and vector
  cell_roster <- list()

  for (group in groups){
    for (subGroup in subGroups){

      Clusterpool <- returnClusters(group, subGroup)

      # First grab the indices of the features
      GeneIndicies <- which(unlist(clean_neuron_object@assays$RNA@data@Dimnames[1]) %in% features_no_key)
      # Then grab the indicies of the cells for the spcefic Clusterpool
      CellIndicies <- which(meta_ident[[clusters_location]] %in% Clusterpool)
      # Get the relevant data using the the indicies otained above
      raw_data <- GetAssayData(clean_neuron_object, assay = 'RNA', slot = 'data')[GeneIndicies, CellIndicies]

      # Dataframe with the cells and the associated clusters
      clusters_n_cells <- meta_ident %>%
        select(all_of(clusters_location), id_plus) %>%
        slice(CellIndicies)
      
      key <- paste(subGroup, group, sep = " ")

      # Turn raw_data into a data frame
      cell_roster[[key]] <- raw_data %>% # Using dplyr
      # Change the object to a data frame
      as.data.frame() %>% 
      # Make the row names an actual column, necessary for an acurate pivot
      mutate('features.label' = row.names(raw_data)) %>%
      # Pivot the data frame with the cells and raw counts in separte columns
      pivot_longer(cols = -features.label, names_to = 'cell.barcode', values_to = 'raw_counts') %>%
      # Add a column with the clusters
      mutate('cluster' = clusters_n_cells[match(cell.barcode, row.names(clusters_n_cells)), 1],
             'sample' = clusters_n_cells[match(cell.barcode, row.names(clusters_n_cells)), 2], .keep = "unused") %>%
      # Add a column with the id and subGroup
      mutate('id' = group, 'subgr' = subGroup)
    }
  }
  return(cell_roster)
}

# TODO: These can be separated as utility functions.R
# Create a function for calculating percent expressed
pct_calc <- function (x) {
  length(x[x > 0]) / length(x) * 100
}

# Create a function to calculate the z score
zs_calc <- function (x) {
  (x - mean(x)) / sd(x)
}

# Adding the standard error function
SE <- function(x) {
  sd(x) / sqrt(length(x))
}

returnAllCellRoster <- function(clean_neuron_object){
  cell_roster <<- createCellRoster(clean_neuron_object)

  # Now to make a master roster of the cells
  all_cell_roster <- data.table::rbindlist(cell_roster)

  return(all_cell_roster)
}

createListbyCluster <- function(clean_neuron_object){
  
  # Calculate average expression, percent expressed, and average expression scaled
  # for each cluster and gene combination
  lbc <- all_cell_roster %>%
    # group by cluster and gene combinations
    group_by(cluster, features.label) %>%
    # calculate the average expression and percent expressed
    summarise(avg.exp = mean(expm1(raw_counts)),
              pct.exp = pct_calc(raw_counts)) %>%
    # Ungroup and perform the appropriate scaling
    ungroup(cluster, features.label) %>%
    # This scaling turns to NA any thing 4 standard deviations above the mean
    mutate(avg.exp.scaled = mad_scaled(avg.exp)) # human: median vs zs_scaled

  # lbc_new %>% ggplot(aes(x=cluster, y = log(avg.exp, base = 100), color = features.label)) + geom_point()
  # I had to exclude a very high outlier for Grin2a and shift to a log 10 scale.
  # TODO: look into the scaling.
  # TODO: protocol for creating new graphs
  # TODO: make the standard deviation calculation exclude the extremely large values

  # TODO: Make PoolnShare (deprecated now CPR_new) function take a vector of clusterpool
  # all instead of indivdual 2 scores + 2 clusterpool names. 

  return(lbc)
}

createClusterPoolResults <- function(clean_neuron_object){
  # Using all_cell_roster to calcualte the average expression and percent expressed
  # for each cluster pool (id and subgr) and gene combination
  CPR <- all_cell_roster %>%
    # group by clusterpool (id and subgroup) gene combinations
    group_by(id, subgr, features.label) %>%
    # perform the following calculations within the groups and store the
    # neccessary information
    summarise(avg.exp = mean(expm1(raw_counts)),
              pct.exp = pct_calc(raw_counts),
              avg.std.err = SE(expm1(raw_counts)),
              avg.lower = avg.exp - avg.std.err,
              avg.upper = avg.exp + avg.std.err) %>%
    # Ungroup the data table and caluclate the appropriate scaling
    ungroup(id, subgr, features.label) %>%
    mutate(avg.exp.z.scaled = mad_scaled(avg.exp)) # human: median vs zs_scaled

  return(CPR)
}

#function to set the width and height of plot saved as an image into global values
getImageDimensions <- function(base_filename, height, width){

  #get amount of data to fill dimensions.
  numClusterpools <- length(returnClusterpool_names())
  numSubgroups <- length(returnClusterpool_subgroups())

  numOfGroups <- numClusterpools * numSubgroups
  numOfGenes <- length(returnAllClusters())
  numOfClusters <- length(returnFeatures())

  if (height == 1 && width == 1){
    if (base_filename == "Quad_BarPlot"){
      height <- 5000
      width <- 5000
    } else if (base_filename == "PooledDotPlot"){
      height <- 1500 + (200 * ( numOfClusters / numOfGroups ))
      width <- 800 + (300 * numOfGroups)
    } else if (base_filename == "DotPlot"){
      height <- 2100 + (200 * numOfClusters / numOfGroups)
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
      height <- 1500 + (200 * ( numOfClusters / numOfGroups ))
    } else if (base_filename == "DotPlot"){
      height <- 2100 + (200 * numOfClusters / numOfGroups)
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
save_image <- function(base_filename, Plot, height = 1, width = 1, dpi = 300){
  
  #set working directory to output directory
  if (!(grepl( "Output", getwd(), fixed = TRUE))){
    setwd("../Output" )
  }
  # In any case set the working directory
  #setwd("../Output")

  #get proper calculated dimensions
  dimensions = getImageDimensions(base_filename, height, width)
  
  #get current date to put on graph title
  curr_date <- format(Sys.Date(),"%b_%d%_%Y" )

  #create folder in output directory with project name
  project_name <- returnProjectName()
  dir.create(project_name, showWarnings = FALSE)

  #put all quad bar plots in another folder inside the project.
  if (substr(base_filename, 1,3) == "QB_"){
    dir.create(sprintf("%s/Quad_BarPlot", project_name), showWarnings = FALSE)
    filename <- sprintf("%s/Quad_BarPlot/%s_%s.jpg", project_name, base_filename, curr_date)
  } else {
    filename <- sprintf("%s/%s_%s.jpg", project_name, base_filename, curr_date)
  }
  
  #save plot to proper location as a jpg.

  ggsave(filename, plot = Plot, device = "jpg",
  height = dimensions[1],
  width = dimensions[2], units = "px",
  type = "cairo", limitsize = FALSE, dpi = dpi) #, dpi = 400)

  #set working directory to what it was before
  setwd("../Scripts" )
}

#call this function to load an RDS data file into an object
load_data <- function(fileLocation){

  print("Loading data into R...")

  if (fileLocation == "NULL"){
    filename <- file.choose()
    clean_neuron_object <- readRDS(filename)
  } else if (str_detect(fileLocation,
                        regex("^.*.(rda|rdata)$", ignore_case = T))) {
    # Start a new environment
    myenv <- new.env()
    # Load the .RData or RDA into the environment
    load(fileLocation, envir = myenv)
    # set the first object in the environment to the proper name
    clean_neuron_object <- get0(ls(myenv)[1], envir = myenv)
    # remove the environment
    rm(myenv)
  }
  else {
    clean_neuron_object <- readRDS(paste(fileLocation, "/clean_neuron_object.RDS", sep = ""))
  }

  print("RDS file loaded!")

  return(clean_neuron_object)
}
