createCellRoster <- function(clean_neuron_object){

  features_no_key <- returnCleanLabelList()

  groups <- returnClusterpool_names()
  subGroups <- returnClusterpool_subgroups()

  for (group in groups){
    for (subgroup in subGroups){

      Clusterpool <- returnClusters(group, subGroup)

      GeneIndices <- which(unlist(clean_neuron_object@assays$RNA@data@Dimnames[1]) %in% features_no_key)
      CellIndicies <- which(unlist(clean_neuron_object@active.ident) %in% Clusterpool)

      raw_data <- GetAssayData(clean_neuron_object, assay = 'RNA', slot = 'data')[GeneIndicies, CellIndicies]

      clusters_n_cells <- as.data.frame(unlist(clean_neuron_object@active.ident)[CellIndicies])

      key <- paste(subgroup, name, sep = " ")

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
    }
  }  
}

returnAllCellRoster <- function(){
  cell_roster <- createCellRoster()

  all_cell_roster <- data.table::rbindlist(cell_roster)

  return(all_cell_roster)
}

createListbyCluster <- function(clean_neuron_object){
  
  #get gene names as features and project name
  features <- returnFeatures()

  #get dotplot information, information about average expression and percent expressed
  b <- getDotPlot(clean_neuron_object, features)


  groups <- returnClusterpool_names()
  subGroups <- returnClusterpool_subgroups()

  ListByCluster <- list()

  for (group in groups){
    for (subgroup in subGroups){

      #create key name from subgroup and group so it looks good on the graph.
      Clusterpool <- returnClusters(group, subgroup)
      key <- paste(subgroup, group, sep = " ")
      ListByCluster[[key]] <- b$data[b$data$id %in% Clusterpool,]
    }
  }  

  return(ListByCluster)
}

createListByClusterAll <- function(clean_neuron_object){
  
  #get gene names as features and project name
  features <- returnFeatures()

  #get dotplot information, information about average expression and percent expressed.
  b <- getDotPlot(clean_neuron_object, features)

  ClusterPoolAll <- returnAllClusters()

  ListbyClusterAll <- b$data[b$data$id %in% ClusterPoolAll,]

  clean_label_list <- returnCleanLabelList()
  ListbyClusterAll[, ncol(ListbyClusterAll) + 1] <- data.frame(features.label = clean_label_list)

  return(ListbyClusterAll)

}

createClusterPoolResults <- function(clean_neuron_object){

  #get group, subgroup and features from JSON as a list
  groups <- returnClusterpool_names()
  subGroups <- returnClusterpool_subgroups()
  features <- returnFeatures()

  #get dotplot information, information about average expression and percent expressed.

  #get ListByCluster, will be necessary to extract information later on.
  ListByCluster <- createListbyCluster(clean_neuron_object)

  #create base dataframe that we will be adding data into.
  ClusterPoolResults <- data.frame(avg.exp = numeric(),
                                  pct.exp = numeric(),
                                  features.plot = character(),
                                  id = character(),
                                  features.label = character(),
                                  SubGroup = character(),
                                  avg.std.err = numeric(),
                                  avg.lower = numeric(),
                                  avg.upper = numeric())
  
  #for all possible combinations of group and subgroup, add a row with relevant information
  for (group in groups){
    for (subGroup in subGroups){
      for (feature in features){

        feature_noRNA_ <- str_remove(feature, 'rna_')
        
        currentCluster <- ListByCluster[[paste(subGroup, group, sep = " ")]]

        #get wanted index of wanted data
        GeneIndices <- match(feature_noRNA_, unlist(clean_neuron_object@assays$RNA@data@Dimnames[1]))
        CellIndicies <- which(unlist(clean_neuron_object@active.ident) %in% unique(currentCluster$id))

        Transcript_exp <- GetAssayData(clean_neuron_object, assay = 'RNA', slot = 'data')[GeneIndices, CellIndicies]

        raw.avg <- mean(expm1(Transcript_exp))
        raw.pct <- length(Transcript_exp[Transcript_exp > 0]) / length(Transcript_exp) * 100

        avg.std.err <- SE(expm1(Transcript_exp))
        avg.lower <- raw.avg - avg.std.err
        avg.upper<- raw.avg + avg.std.err


        ClusterPoolResults <- ClusterPoolResults %>%
        add_row(avg.exp = raw.avg,
                pct.exp = raw.pct, 
                features.plot = feature, 
                id = group, 
                features.label = feature_noRNA_, 
                SubGroup = subGroup,  
                avg.std.err = avg.std.err, 
                avg.lower = avg.lower,
                avg.upper = avg.upper)
      }
    }
  }


  #I just shortened the name so it could fit on one line.
  CPR <- ClusterPoolResults

  CPR[, ncol(CPR) + 1] <- data.frame(avg.exp.re.scaled = (CPR[, 1] - colMeans(CPR[1]))/sd(CPR$avg.exp))
  
  names(CPR)[ncol(CPR)] <- 'avg.exp.z.scaled'

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
save_image <- function(base_filename, Plot, height = 1, width = 1){
  
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
    dir.create(sprintf("%s/QBC", project_name), showWarnings = FALSE)
    filename <- sprintf("%s/QBC/%s_%s.jpg", project_name, base_filename, curr_date)
  } else {
    filename <- sprintf("%s/%s_%s.jpg", project_name, base_filename, curr_date)
  }
  
  #save plot to proper location as a jpg.
  ggsave(filename, plot = Plot, device = "jpg", height = dimensions[1], width = dimensions[2], units = "px", type = "cairo", limitsize = FALSE)

  #set working directory to what it was before
  setwd("../Scripts" )
}

#run DotPlot seurat function to get wanted dataframe using genes as features.
getDotPlot <- function(clean_neuron_object, features){
  dotPlot <- DotPlot(clean_neuron_object, features = features)
  return(dotPlot)
}

#calculate the standard error
SE <- function (x) {
  sd(x)/sqrt(length(x))
}

#return average expression
AvgExpPar <- function (x) {
  x %>% 
    group_by(features.plot) %>% 
    summarise(std.err = SE(avg.exp), avg.exp = mean(avg.exp), lower = avg.exp - std.err, upper = avg.exp + std.err)
}

#return percent expressed
PctExpPar <- function (x) {
  x %>% 
    group_by(features.plot) %>% 
    summarise(std.err = SE(pct.exp), pct.exp = mean(pct.exp), lower = pct.exp - std.err, upper = pct.exp + std.err)
}

#call this function to load an RDS data file into an object
load_data <- function(fileLocation){
  
  if (fileLocation == "NULL"){
    filename <- file.choose()
    object <- readRDS(filename)
  } else if (str_detect(fileLocation,
                        regex("^.*.(rda|rdata)$", ignore_case = T))) {
    myenv <- new.env()
    load(fileLocation, envir = myenv)
    object <- get0(ls(myenv)[1], envir = myenv)
    rm(myenv)
  } else if (str_detect(fileLocation, regex("^.*.(rds)$", ignore_case = T))){
    object <- readRDS(fileLocation)
  } else {
    return(NULL)
  }
  
  return(object)
}

# The functions bellow will be used to add/remove/configure grouping layers
# This function adds a layer to the provided list of lists
create_listoflayers <- function() {
  list_of_layers <<- list(top = list())
}
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
  layer_name <- paste("layer", layer_number, sep = "")
  print(layer_name)
  layer_list <- layer_list[-which(names(layer_list) == layer_name)]
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



