#### Pre-analysis functions
# This script contains functions that are used by multiple graphing scripts
# to extract information from seurat objects and process the information
# tidy the information, and perform basic calculations.


returnAllCellRoster <- function(RDfile){
  # Get the presentable genelabels by removing the 'rna_' tag
  features_no_key <- returnCleanLabelList()
  
  # Get the groups and subgroups from the json file
  groups <- returnClusterpool_names()
  subGroups <- returnClusterpool_subgroups()
  
  # Next look for the first cluster pool within each column
  all_clusters <- returnAllClusters()
  
  # TODO: need to find a way to get the correct formatting for the genes
  # for example GRIN1 vs. Grin1
  # Now merging meta data and active ident and finding where the clusters are stored
  # First merge @meta.data and @active.ident
  meta_ident <- as.data.frame(RDfile@meta.data) %>%
    mutate(cells = row.names(as.data.frame(RDfile@active.ident)),
           id_plus = as.data.frame(RDfile@active.ident)[[1]],
           cell.barcode = row.names(.))
  # Then get the names of the collumns that have at least one of the
  # clusters in the clusterpool
  clusters_location <- names(which.max(apply(
    meta_ident, 
    FUN = function (x) sum(x %in% all_clusters), 
    MARGIN = 2)))
  # TODO: add the ability to introduce secondary identifiers
  # TODO: turn this into a function

  # Stores information on the clusters and the associated groups
  # TODO: find a json package friendly way of doing this
  Clusterpool <- list()
  for (group in groups){
    for (subGroup in subGroups){
      key <- paste(subGroup, group, sep = " ")
      # TODO: Space delineation not so good
      Clusterpool[[key]] <- returnClusters(group, subGroup)
    }
  }
  # First grab the indices of the features
  GeneIndicies <- which(unlist(RDfile@assays[[which(names(RDfile@assays) == assay)]]@data@Dimnames[1]) %in% features_no_key)
  # Then grab the indices of the cells for the spcefic Clusterpool
  CellIndicies <- which(meta_ident[[clusters_location]] %in% unique(all_clusters))
  # Get the relevant data using the the indicies otained above
  raw_data <- GetAssayData(RDfile, assay = assay, slot = slot)[GeneIndicies, CellIndicies]

  # Dataframe with the cells and the associated clusters
  # add all the extra variables that are selected in extra_pool
  clusters_n_cells <- meta_ident %>%
    select(cell.barcode = cells,
           all_of(clusters_location), all_of(unlist(extra_pool[["top"]]))) %>%
    slice(CellIndicies) %>%
    mutate(cats = lapply(.[[clusters_location]],
                         function (x) list.which(Clusterpool, x %in% .)) %>>%
             lapply(function (x) names(Clusterpool)[x]),
           id = lapply(cats, function (x) unique(str_split_i(x, " ", i = 2))),
           subgr = lapply(cats, function (x) unique(str_split_i(x, " ", i = 1)))) %>%
    select(-cats)
  
  # Turn raw_data into a data frame
  all_cell_roster <- raw_data %>% # Using dplyr
  # Change the object to a data frame
  as.data.frame() %>% 
  # Make the row names an actual column, necessary for an accurate pivot
  mutate('features.label' = row.names(raw_data)) %>%
  # Pivot the data frame with the cells and raw counts in separate columns
  pivot_longer(cols = -features.label, names_to = 'cell.barcode', values_to = 'raw_counts') %>%
  # Add a column with the clusters
  full_join(clusters_n_cells, by = "cell.barcode") %>%
  # TODO: what if there is something else named cluster?
  rename("cluster" = all_of(clusters_location))
  
  return(all_cell_roster)
}
# TODO: The next steps are to figure out how to retrieve the information for grouping
# Perhaps create a centralized function for this.

# TODO: These can be separated as utility functions.R
# Create a function for calculating percent expressed
pct_calc <- function (x) {
  length(x[x > 0]) / length(x) * 100
}

# Create a function to calculate the z score
zs_calc <- function (x) {
  (x - mean(x)) / sd(x)
}

# Counts per million scale
# Change RNA to raw
run_CPM <-  function (roster, scale.method, pre.exp = F) {
  if (scale.method %in% c("CPM", "log1mCPM", "log1pCPM",
                          "zsoflog1pCPM")) {
    roster %>%
      mutate(raw_counts = case_when(
        pre.exp == F ~ (raw_counts / nCount_RNA) * 10^6,
        pre.exp == T ~ expm1(raw_counts / nCount_RNA) * 10^6))
    
  } else if (scale.method == "countspertotal") {
    
    roster %>%
      mutate(raw_counts = raw_counts * nCount_RNA / 10^6)
    
  } else {
    roster
  }
}

# Adding the standard error function
SE <- function(x) {
  sd(x) / sqrt(length(x))
}

# Function for getting 5% over the max of a set of values
p5max <- function(x) {
  round(max(x) * 1.05)
}

# Function to wrap text
wrapper <- function(x, ...) 
{
  paste(strwrap(x, ...), collapse = "\n")
}

# Calculate average expression, percent expressed, and average expression scaled
# for each cluster and gene combination
createListbyCluster <- function(scale.method = "z-score", pre.expm1 = T){
  lbc <- all_cell_roster %>%
    run_CPM(scale.method) %>%
    # group by cluster and gene combinations
    group_by(cluster, features.label) %>%
    # calculate the metrics
    summarise(avg.exp = ifelse(pre.expm1, mean(expm1(raw_counts)), mean(raw_counts)),
              pct.exp = pct_calc(raw_counts),
              med.exp = median(raw_counts),
              .groups = "keep") %>%
    ungroup(cluster, features.label) %>%
    mutate(avg.exp.scaled = case_when(
      scale.method == "z-score" ~ zs_calc(avg.exp),
      scale.method == "log10" ~ log10(avg.exp),
      scale.method == "log1p" ~ log1p(avg.exp),
      scale.method %in% c("CPM", "noscale") ~ avg.exp,
      scale.method == "log1mCPM" ~ log10(avg.exp - 1),
      scale.method == "log1pCPM" ~ log10(avg.exp + 1),
      scale.method == "zsoflog1pCPM" ~ zs_calc(log10(avg.exp + 1))
    ),
    avg.exp.scaled = ifelse(is.infinite(avg.exp.scaled), NA, avg.exp.scaled))

  # TODO: protocol for creating new graphs

  return(lbc)
}

# Function for strategically expanding overlapping groups
# to make the data more tidy.
group_expand <- function (roster, pool.level = "1", overgrouped = NULL) {
  if (is.null(overgrouped)) {
    overgrouped <- list.which(roster, is.list(.) & .name %in%
                                unlist(extra_pool[[pool.level]]))
  }
  
  for (over in overgrouped) {
    roster <- roster %>%
      unnest_longer(colnames(roster[over]))
  }
  return(roster)
}

# Function that pastes from specific columns in a list of column names.
pool_level_paste <- function(pool.result, pool.list) {
  pool.list <- pool.list[-which(pool.list %in% "features.label")]
  pool.result <- pool.result[, unlist(pool.list)]
  group.label <- apply(pool.result, 1, paste, collapse = " ")
  return(group.label)
}

# Function that groups the information for all_cell_roster at any level
# and provides the averaged expression and percent expressed for that group
createClusterPoolResults <- function(roster = all_cell_roster,
                pool.level = "1",
                scale.method = "z-score",
                pre.expm1 = T, ...) {
  
  any_overgrouped <- list.any(roster[unlist(extra_pool[[pool.level]])], is.list(.))
  
  if (any_overgrouped) {
    roster <- group_expand(roster, pool.level)
  }
  
  CPR <- roster %>%
    run_CPM(scale.method) %>%
    group_by(across(unlist(extra_pool[[pool.level]]))) %>%
    summarise(avg.exp = ifelse(pre.expm1, mean(expm1(raw_counts)), mean(raw_counts)),
              pct.exp = pct_calc(raw_counts),
              avg.std.err = SE(expm1(raw_counts)),
              avg.lower = avg.exp - avg.std.err,
              avg.upper = avg.exp + avg.std.err,
              .groups = "keep",
              ...) %>%
    ungroup(everything()) %>%
    mutate(avg.exp.scaled = case_when(
      scale.method == "z-score" ~ zs_calc(avg.exp),
      scale.method == "log10" ~ log10(avg.exp),
      scale.method == "log1p" ~ log1p(avg.exp),
      scale.method %in% c("CPM", "noscale") ~ avg.exp,
      scale.method == "log1mCPM" ~ log10(avg.exp - 1),
      scale.method == "log1pCPM" ~ log10(avg.exp + 1),
      scale.method == "zsoflog1pCPM" ~ zs_calc(log10(avg.exp + 1))),
      group.label = pool_level_paste(., extra_pool[[pool.level]]))
  
  return(CPR)
  }
  # TODO: Table is missing other descriptors

# Function to set the width and height of plot saved as an image into global values
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
save_image <- function(base_filename, Plot, height = 1, width = 1,
                       dpi = 300, bg, device = "png"){
  
  #set working directory to output directory
  if (!(grepl( "Output", getwd(), fixed = TRUE))){
    setwd("../Output" )
  }

  #get proper calculated dimensions
  dimensions <- getImageDimensions(base_filename, height, width)
  
  #get current date to put on graph title
  curr_date <- format(Sys.Date(),"%b_%d%_%Y" )

  #create folder in output directory with project name
  project_name <- returnProjectName()
  dir.create(project_name, showWarnings = FALSE)
  
  #Set the file name
  filename <- sprintf("%s/%s_%s.%s", project_name, base_filename, curr_date, device)
  
  #save plot to proper location as a jpg.
  type <- ifelse(device == "emf", NA, "cairo")

  ggsave(filename, plot = Plot, device = device,
         height = dimensions[1],
         width = dimensions[2], units = "px",
         type = type, limitsize = FALSE, dpi = dpi)

  #set working directory to what it was before
  setwd("../Scripts" )
}

#call this function to load an RDS data file into an object
load_data <- function(fileLocation){

  print("Loading data into R...")

  if (fileLocation == "NULL"){
    filename <- file.choose()
    rdfile <- readRDS(filename)
  } else if (str_detect(fileLocation,
                        regex("^.*.(rda|rdata)$", ignore_case = T))) {
    # Start a new environment
    myenv <- new.env()
    # Load the .RData or RDA into the environment
    load(fileLocation, envir = myenv)
    # set the first object in the environment to the proper name
    rdfile <- get0(ls(myenv)[1], envir = myenv)
    # remove the environment
    rm(myenv)
  }
  else if (str_detect(fileLocation,
                      regex("^.*.(rds)$", ignore_case = T))){
    rdfile <- readRDS(fileLocation)
  }

  print("RDS file loaded!")

  return(rdfile)
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

