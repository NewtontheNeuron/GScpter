#### Pre-analysis functions
# This script contains functions that are used by multiple graphing scripts
# to extract information from seurat objects and process the information
# tidy the information, and perform basic calculations.


createCellRoster <- function(clean_neuron_object){
  # Get the presentable genelabels by removing the 'rna_' tag
  features_no_key <- returnCleanLabelList()
  
  # Get the groups and subgroups from the json file
  groups <- returnClusterpool_names()
  subGroups <- returnClusterpool_subgroups()
  
  # Next look for the first cluster pool within each collumn
  all_clusters <- returnAllClusters()
  
  # TODO: need to find a way to get the correct formating for the genes
  # for example GRIN1 vs. Grin1
  # Now merging meta data and active ident and finding where the clusters are stored
  # First merge @meta.data and @active.ident
  meta_ident <- as.data.frame(clean_neuron_object@meta.data) %>%
    mutate(cells = row.names(as.data.frame(clean_neuron_object@active.ident)),
           id_plus = as.data.frame(clean_neuron_object@active.ident)[[1]],
           cell.barcode = row.names(.))
  # Then get the names of the collumns that have at least one of the
  # clusters in the clusterpool
  clusters_location <- names(which.max(apply(
    meta_ident, 
    FUN = function (x) sum(x %in% all_clusters), 
    MARGIN = 2)))
  # TODO: add the ability to introduce secondary identifiers
  # TODO: turn this into a function

  # Create empty list and vector
  cell_roster <- list()

  for (group in groups){
    for (subGroup in subGroups){

      Clusterpool <- returnClusters(group, subGroup)

      # First grab the indices of the features
      GeneIndicies <- which(unlist(clean_neuron_object@assays$RNA@data@Dimnames[1]) %in% features_no_key)
      # Then grab the indices of the cells for the spcefic Clusterpool
      CellIndicies <- which(meta_ident[[clusters_location]] %in% Clusterpool)
      # Get the relevant data using the the indicies otained above
      raw_data <- GetAssayData(clean_neuron_object, assay = 'RNA', slot = 'data')[GeneIndicies, CellIndicies]

      # Dataframe with the cells and the associated clusters
      # add all the extra variables that are selected in extra_pool
      clusters_n_cells <- meta_ident %>%
        select(cell.barcode = cells,
               all_of(clusters_location), all_of(unlist(extra_pool[["top"]]))) %>%
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
      full_join(clusters_n_cells, by = "cell.barcode") %>%
      # TODO: what if there is something else named cluster?
      rename("cluster" = all_of(clusters_location)) %>%
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

createListbyCluster <- function(scale.method = "z-score"){
  
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
    mutate(avg.exp.scaled = case_when(
      scale.method == "z-score" ~ zs_calc(avg.exp),
      scale.method == "log10" ~ log10(avg.exp),
      scale.method == "log1p" ~ log1p(avg.exp)
    )) # human: median vs zs_scaled

  # lbc_new %>% ggplot(aes(x=cluster, y = log(avg.exp, base = 100), color = features.label)) + geom_point()
  # I had to exclude a very high outlier for Grin2a and shift to a log 10 scale.
  # TODO: look into the scaling.
  # TODO: protocol for creating new graphs
  # TODO: make the standard deviation calculation exclude the extremely large values

  # TODO: Make PoolnShare (deprecated now CPR_new) function take a vector of clusterpool
  # all instead of indivdual 2 scores + 2 clusterpool names. 

  return(lbc)
}

createClusterPoolResults <- function(roster = all_cell_roster,
                                     pool.level = "1",
                                     scale.method = "z-score", ...){
  # Using all_cell_roster to calcualte the average expression and percent expressed
  # for each cluster pool (id and subgr) and gene combination
  CPR <- roster %>%
    # group by clusterpool (id and subgroup) gene combinations
    group_by(across(unlist(extra_pool[[pool.level]]))) %>%
    # perform the following calculations within the groups and store the
    # neccessary information
    summarise(avg.exp = mean(expm1(raw_counts)),
              pct.exp = pct_calc(raw_counts),
              avg.std.err = SE(expm1(raw_counts)),
              avg.lower = avg.exp - avg.std.err,
              avg.upper = avg.exp + avg.std.err,
              ...) %>%
    # Ungroup the data table and caluclate the appropriate scaling
    ungroup(everything()) %>%
    mutate(avg.exp.scaled = case_when(
      scale.method == "z-score" ~ zs_calc(avg.exp),
      scale.method == "log10" ~ log10(avg.exp),
      scale.method == "log1p" ~ log1p(avg.exp)
    )) # human: median vs zs_scaled

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
  else if (str_detect(fileLocation,
                      regex("^.*.(rds)$", ignore_case = T))){
    clean_neuron_object <- readRDS(fileLocation)
  }

  print("RDS file loaded!")

  return(clean_neuron_object)
}

# The functions bellow will be used to add/remove/configure grouping layers
# This function adds a layer to the provided list of lists
add_layer <- function(layer_list = list(top = list()), newLayerItems = NULL) {
  # the new layer items should be a list and not nothing
  if(is.null(newLayerItems) | !is.list(newLayerItems)){
    return("You need to input a list of items using newLayerItems = ...")
  }
  layer_name <- paste("layer", length(layer_list), sep = "")
  new_list <- setNames(list(newLayerItems), layer_name)
  layer_list <- append(layer_list, new_list)
  return(layer_list)
}
# Utilized by the other functions
# This function fixes the number order of the layers
fix_layer_number <- function(layer_list = list(top = list())){
  names(layer_list) <- lapply(seq_along(layer_list),
                              FUN = function(x) ifelse(names(layer_list)[[x]] != "top",
                                                       paste("layer", x - 1, sep = ""), "top"))
  return(layer_list)
}
# This function removes a layer from the provided list of lists
rm_layer <- function(layer_list = list(top = list()), layer_number = NULL){
  if(is.null(layer_number)){
    return(NULL)
  }
  layer_list <- within(layer_list, rm(paste("layer", layer_number, sep = "")))
  layer_list <- fix_layer_number(layer_list)
  return(layer_list)
}
# This function changes a specified layer
change_layer <- function(layer_list = list(top = list()),
                         layer_number = NULL,
                         newLayerItems = NULL){
  if(is.null(layer_number) | is.list(newLayerItems)){
    return("You need to specify the layer number that you want to change and
           the list of items that will replace it")
  }
  layer_list[[paste("layer", layer_number, sep = "")]] <- newLayerItems
  return(layer_list)
}
# This function reorders the layer based on the type of reordering
reorder_layer <- function(layer_list = list(top = list()),
                          layer_number = NULL,
                          type = c("top", "bottom", "up", "down")) {
  if(is.null(layer_number)){
    return(NULL)
  }
  layer_name <- paste("layer", layer_number, sep = "")
  
  # Get the above and bellow layers
  split_layers <- function(){
    layers_switched <- which(names(layer_list) %in% c(layer_name, layer_m))
    layers_above <- names(layer_list)[-c(1, layers_switched,
                                         (layer_number + 1):length(layer_list))]
    layers_below <- names(layer_list)[-c(1, layers_switched,
                                          which(names(layer_list) %in% layers_above))]
    return(list(above = layers_above, below = layers_below))
  }
  if(type == "top"){
    layer_list <- layer_list[c("top",
                               layer_name,
                               names(layer_list)[-c(1, which(names(layer_list) %in% layer_name))]
                               )]
  } else if (type == "bottom") {
    layer_list <- layer_list[c("top",
                               names(layer_list)[-c(1, which(names(layer_list) %in% layer_name))],
                               layer_name)]
  } else if (type == "up" & layer_number > 1) {
    layer_m <- paste("layer", layer_number - 1, sep = "")
    unused <- split_layers()
    layer_list <- layer_list[c("top",
                               unused$above,
                               layer_name,
                               layer_m,
                               unused$below)]
  } else if (type == "down" & layer_number < length(layer_list)) {
    layer_m <- paste("layer", layer_number + 1, sep = "")
    unused <- split_layers()
    layer_list <- layer_list[c("top",
                               unused$above,
                               layer_m,
                               layer_name,
                               unused$below)]
  }
  layer_list <- fix_layer_number(layer_list)
  return(layer_list)
}
# function for ignoring a level
ignore_level <- function(listoflevels = c(), nameoflevel = "") {
  listoflevels <- listoflevels[-which(listoflevels %in% nameoflevel)]
  if(!exists(ignored_levels)){
    ignored_levels <<- c()
  }
  ignored_levels <<- c(ignored_levels, nameoflevel)
  return(listoflevels)
}
# function to reinstate levels
reinstate_level <- function(listoflevels = c(), nameoflevel = "") {
  listoflevels <- c(listoflevels, nameoflevel)
  if(!exists(ignored_levels)){
    ignored_levels <<- c()
  }
  ignored_levels <<- ignored_levels[-which(ignored_levels %in% nameoflevel)]
}
# TODO: function to reorder levels

