library(Seurat)

# Set wd to the dir of the current script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load('../../Data/human_ariel_data/top_level_new_annotation.rda')
cell_class <- read.csv('../../Data/human_ariel_data/top_level_annotations_metadata.csv')
str(cell_class)
levels(cell_class$x)

# only run the human stuff
human <- integrated_top_level_obj
clean_neuron_object <- integrated_top_level_obj
rm(integrated_top_level_obj)

# Run both mouse and human
human <- clean_neuron_object
filename <- file.choose()
mouse <- readRDS(filename)
clean_neuron_object <- mouse
gc()

clean_neuron_object <- human
rm(human)

# Find where the clusters are (done outside the loop)
relevant_info <- as.data.frame(clean_neuron_object@meta.data) %>%
  mutate(cells = row.names(as.data.frame(clean_neuron_object@active.ident)),
         id_plus = as.data.frame(clean_neuron_object@active.ident)[[1]])
  
  
clusters_found <- names(which(apply(relevant_info, FUN = function (x) "Excit-01" %in% x, MARGIN = 2) == TRUE))

# this is done inside the loop

CellIndicies <- which(relevant_info[[clusters_found[1]]] %in% c("Excit-01"))
