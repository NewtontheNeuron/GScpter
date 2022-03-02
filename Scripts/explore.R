library(Seurat)

# Set wd to the dir of the current script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load('../../Data/human_ariel_data/top_level_new_annotation.rda')
cell_class <- read.csv('../../Data/human_ariel_data/top_level_annotations_metadata.csv')
str(cell_class)
levels(cell_class$x)
