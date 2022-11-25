# This Script is for the integration of the mouse and human spinal cord
# single cell RNA datasets
# I will be following quite synonymously with the Seurat package and the scripts
# That I have gotten from the kcni summer school

library(Seurat)
library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Import both data sets
source("../Pre_analysis_functions.R")
Human <- load_data("../../../Datasets/human_spinalcord_2022/top_level_new_annotation.rda")
Mouse <- load_data("../../../Datasets/neuron_and_glia_2022/final_meta_dataset.rds")

# Create the homologous gene names
# Take all the mouse genes and turn them into their human homolog names
library(homologene)
mouse_genes <- Mouse@assays$RNA@counts@Dimnames[[1]]
human_genes <- Human@assays$RNA@counts@Dimnames[[1]]
mouse_bc <- Mouse@assays$RNA@counts@Dimnames[[2]]
human_bc <- Human@assays$RNA@counts@Dimnames[[2]]
homolo <- mouse2human(mouse_genes, db = homologeneData2)
rm(Human, Mouse)
gc()
# Remove the homologues that are not present in the human dataset
homoloB <- homolo[-which(!(homolo$humanGene %in% human_genes)), ]
# Merge mouse genes with homoloB
mgdf <- data.frame(mouseGene = mouse_genes, order = seq_along(1:length(mouse_genes)))
homoloC <- full_join(homolo[,1:2], mgdf, by = "mouseGene")
# Remove the extra human gene homologues on a first come first serve basis
duplicates <- homoloC %>% group_by(order) %>% count() %>% filter(n > 1, !is.na(order))
duplicates
torm <- c()
for (dup in seq_along(duplicates[[1]])){
  cat("\n\n-----start-----\n")
  print(dup)
  ndup <- duplicates[dup, ]
  print(ndup)
  torm <- c(torm, which(homoloC$order == ndup$order)[2:ndup$n])
  print(torm)
  cat("-----end-----")
}
homoloD <- homoloC[-torm, ]
homoloD %>% group_by(order) %>% count() %>% filter(n > 1)
# Order homoloD based on how it was in mouse genes
homoloD <- homoloD[order(homoloD$order), ]
# There should not be any NAs in the mouseGene column
any(is.na(homoloD$mouseGene))
# mosueGene should be equal to mouse_genes
all(homoloD$mouseGene == mouse_genes)
# Create a new column with the mouse genes but the homologous version for those
# that have a human homologue
homoloE <- homoloD %>%
  mutate(newMGenes = case_when(
    is.na(humanGene) ~ mouseGene,
    !is.na(humanGene) ~ humanGene
  ))
homoloE
homoloE[which(homoloE$newMGenes == "TUBA3E"),]
data.frame(table(homoloE$newMGenes))

# I need to remove the duplicates formed when a human gene maps to more than one
# mouse gene
# I can do so on a first come first searved basis
duplicates2 <- data.frame(table(homoloE$newMGenes))
duplicates2 <- duplicates2[which(duplicates2$Freq > 1),]
duplicates2
tocg <- c()
for(dup2 in 1:nrow(duplicates2)){
  ndup2 <- duplicates2[dup2, ]
  eachdup <- homoloE[which(homoloE$newMGenes %in% ndup2$Var1),]
  if (any(tolower(eachdup$mouseGene) == tolower(eachdup$newMGene))) {
    tocg <- c(tocg,
              row.names(eachdup[which(
                tolower(eachdup$mouseGene) != tolower(eachdup$newMGene)),
                ]))
  } else if (all(tolower(eachdup$mouseGene) != tolower(eachdup$newMGene))) {
    tocg <- c(tocg, row.names(eachdup[2:nrow(eachdup),]))
  }
}
tocg
homoloF <- homoloE %>%
  mutate(newMGenes = case_when(
    row.names(homoloE) %in% tocg ~ mouseGene,
    !(row.names(homoloE) %in% tocg) ~ newMGenes
  ))
homoloF
homoloF[which(homoloF$newMGenes == "TUBA3E"),]
which(data.frame(table(homoloF$newMGenes))$Freq > 1)
data.frame(table(homoloF$newMGenes))[20405,]
homoloF[which(homoloF$newMGenes == "PISD"),]
# One of them did not change because PISD and Pisd are apparently different genes
# in the mouse. So we can change it directly.
homoloF[which(homoloF$mouseGene == "Pisd"),]$newMGenes <- "Pisd"
homoloF[which(homoloF$mouseGene == "Pisd"),]
which(data.frame(table(homoloF$newMGenes))$Freq > 1)
# Excelent now there are no duplicates

# Grab the human and mouse matricies and change the mouse genes then attempte the merge
Human <- load_data("../../../Datasets/human_spinalcord_2022/top_level_new_annotation.rda")
Mouse <- load_data("../../../Datasets/neuron_and_glia_2022/final_meta_dataset.rds")
cMouse <- Mouse@assays$RNA@data
mMouse <- Mouse@meta.data %>% mutate(species = "mouse")
rm(Mouse)
cHuman <- Human@assays$RNA@data
mHuman <- Human@meta.data %>% mutate(species = "human")
rm(Human)
rm(duplicates, homolo, homoloB, homoloC, homoloD, mgdf, ndup, dup, torm,
   duplicates2, dup2, ndup2, homoloE, tocg, eachdup)
gc()
cMouse@Dimnames[[1]] <- homoloF$newMGenes

# B12
# Create the mouse and human seurat objects
sMouse <- CreateSeuratObject(cMouse, meta.data = mMouse)
sHuman <- CreateSeuratObject(cHuman, meta.data = mHuman)
rm(cMouse, cHuman, mMouse, mHuman)
gc()
# Create the merged seurat object
# Normalize and find variable features for each of the seurat objects separate
# Use only the genes that are in common
homoloG <- filter(.data = homoloF, !is.na(humanGene))
homoloG
# Next remove the extra mouse genes that the human genome does not have
duplicates3 <- homoloG %>% group_by(humanGene) %>% count() %>% filter(n > 1, !is.na(order))
duplicates3
torm2 <- c()
for (dup3 in seq_along(duplicates3[[1]])){
  cat("\n\n-----start-----\n")
  print(dup3)
  ndup3 <- duplicates3[dup3, ]
  print(ndup3)
  torm2 <- c(torm2, which(homoloG$humanGene == ndup3$humanGene)[2:ndup3$n])
  print(torm2)
  cat("-----end-----")
}
homoloH <- homoloG[-torm2, ]
homoloH
# Create the list of seurat objects
data.list <- list(sMouse, sHuman)
rm(sMouse, sHuman)
gc()
data.list <- lapply(data.list, function(x) {
  #x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10^6)
  x <- x[homoloH$newMGenes,]
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# Select features that are variable across both seurat objects
features_b <- SelectIntegrationFeatures(object.list = data.list)
rm(human_bc, mouse_bc, human_genes, mouse_genes, dup3, ndup3, torm2, duplicates3,
   homoloF, homoloG, homoloH)
gc
# The following step takes the largest amount of time It elapsed 6hs, 8m, and 10s
anchors <- FindIntegrationAnchors(object.list = data.list, anchor.features = features_b)
str(anchors)
rm(features_b, data.list)
gc()
# Save point: anchors
saveRDS(anchors, "../../../Datasets/integrated_mouse_human/mouse_human_int_anchors.RDS")
anchors <- readRDS("../../../Datasets/integrated_mouse_human/mouse_human_int_anchors.RDS")

# this command creates an 'integrated' data assay
# Make sure that only the anchorset is in the the Global Evirnoment for memory reasons.
mouse_human_int <- IntegrateData(anchorset = anchors)
rm(anchors)
gc()
# Save point: integrated data
saveRDS(mouse_human_int, "../../../Datasets/integrated_mouse_human/mouse_human_int_preprocess.RDS")
mouse_human_int <- readRDS("../../../Datasets/integrated_mouse_human/mouse_human_int_preprocess.RDS")
# scale the data
mouse_human_int <- ScaleData(mouse_human_int)
# Run the PCAs with a standard number dimensions
mouse_human_int <- RunPCA(mouse_human_int, pcs = 50)
# Show the Elbow plot
ElbowPlot(mouse_human_int, ndims = 50)
# Clustering
mouse_human_int <- FindNeighbors(mouse_human_int, dims = 1:30)
mouse_human_int <- FindClusters(mouse_human_int, resolution = 0.5)
table(mouse_human_int@meta.data$seurat_clusters)
# UMAP
mouse_human_int <- RunUMAP(mouse_human_int, reduction = "pca", dims = 1:30)
DimPlot(mouse_human_int, reduction = "umap")
DimPlot(mouse_human_int, reduction = "umap", group.by = "species")
# It looks like there is batch effect
View(mouse_human_int@meta.data)
mouse_human_int@meta.data[ncol(mouse_human_int@meta.data) + 1] <-
  ifelse(!is.na(mouse_human_int@meta.data$final_coarse_types),
         mouse_human_int@meta.data$final_coarse_types,
         mouse_human_int@meta.data$top_level_annotation)
any(is.na(mouse_human_int@meta.data[ncol(mouse_human_int@meta.data)]))
colnames(mouse_human_int@meta.data)[ncol(mouse_human_int@meta.data)] <- "overall_annotations"
# Mutate similar cells
mouse_human_int@meta.data <- mouse_human_int@meta.data %>%
  mutate(overall_annotations = case_when(
    overall_annotations == "Astrcytes" ~ "Astrocytes",
    overall_annotations == "OPCs" ~ "OPC",
    overall_annotations == "Ependymal Cells" ~ "Ependymal",
    overall_annotations == "Schwann Cells" ~ "Schwann",
    overall_annotations == overall_annotations ~ overall_annotations
  ))
# Redo the Dimplot
DimPlot(mouse_human_int, reduction = "umap", group.by = "overall_annotations") +
  scale_color_viridis_d()
DimPlot(mouse_human_int, reduction = "umap", group.by = "orig.ident")




# It turns out that you do not need to merge at this step proceed to step B12
merged_mat <- merge.Matrix(cMouse, cHuman, by.x = cMouse@Dimnames[[1]],
                           by.y = cHuman@Dimnames[[1]], all.x = T, all.y = T)
saveRDS(merged_mat, "../../../Datasets/integrated_mouse_human/merged_mtx.RDS")
rm(cMouse, cHuman, homoloF)
gc()
merged_md <- bind_rows(mMouse, mHuman)
merged_md
saveRDS(merged_md, "../../../Datasets/integrated_mouse_human/merged_md.RDS")


# Start from here next time
merged_mat <- readRDS("../../../Datasets/integrated_mouse_human/merged_mtx.RDS")
merged_md <- readRDS("../../../Datasets/integrated_mouse_human/merged_md.RDS")
# Create the Seurat object
mouse_human_int <- CreateSeuratObject(counts = merged_mat, meta.data = merged_md)
# This did not work
homologeneDataVeryNew <- updateHomologene()
new <- mouse2human(mouse_genes, db = homologeneDataVeryNew)
homolo$humanGene == new$humanGene

# Acquire the data slot matrix, which contains log-normalized values of the counts slot
# , for each dataset. Also save the metadata. Then exchange the mouse genes for its
# human homologous names

mHuman <- as.data.frame(Human@meta.data) %>% mutate(species = "human")
cMouse <- as.data.frame(Mouse@assays$RNA@data)
# To account for the genes that are not homologous to the human genome for instance
# mouse gm genes
mgdf <- data.frame(mouseGene = mouse_genes, order = seq_along(1:length(mouse_genes)))
hmdf <- data.frame(humanGene = human_genes)
new_homolo <- full_join(homolo[,1:2], mgdf, by = "mouseGene")
# I dont want the human genes that I dont have
# new_homolo[which(new_homolo$mouseGene == "Sgk3"),]
new_homolo <- right_join(new_homolo, hmdf, by = "humanGene")
new_homolo <- new_homolo[order(new_homolo$order), ]
# We will have to remove the the extra human gene
# Algorithm: remove them on a first come first serve basis
duplicates <- new_homolo %>% group_by(order) %>% count() %>% filter(n > 1, !is.na(order))
a <- new_homolo %>% filter(!is.na(order))
torm <- c()
for (dup in seq_along(duplicates[[1]])){
  cat("\n\n-----start-----\n")
  print(dup)
  ndup <- duplicates[dup, ]
  print(ndup)
  print(a[which(a$order == ndup$order)[2:ndup$n],])
  torm <- c(torm, which(a$order == ndup$order)[2:ndup$n])
  print(torm)
  cat("----end-----")
}
# Now remove the extra genes
homolo2 <- a[-torm, ]
homolo2 %>% group_by(order) %>% count() %>% filter(n > 1)
rm(a, duplicates, ndup, mgdf, hmdf, homolo, new_homolo, human_genes, mouse_genes, dup, torm)
# homolo2 should have the master list of genes in common between the mouse and human
# Attempt to merge the mouse count matrix to only the genes in the data frame
# Remove the pre analysis functions
mMouse <- Mouse@assays$RNA@data
mdMouse <- Mouse@meta.data %>% mutate(species = "mouse")
rm(Mouse)
gc()

head(mMouse)[, 1]
head(homolo2)
homolo2 <- homolo2[order(homolo2$order), ]
head(homolo2)

row.names(new_homolo) <- new_homolo$order
new_homolo$mouseGene == mouse_genes

# Order homolo based on how it is in the dimnames of the mouse dataset
homolo[head(match(mouse_genes, homolo$mouseGene)), 2]

Mouse@assays$RNA@data@Dimnames[[1]]

mMouse <- as.data.frame(Mouse@meta.data) %>% mutate(species = "mouse")

rm(Human, Mouse)

cMouse <- merge(cMouse, new[,1:2], by.x = "row.names", by.y = "mouseGene", all.x = T, all.y = F)
row.names(cMouse) <- humanGene
cMouse <- cMouse[, 1:(ncol(cMouse) - 1)]

# Test example matrix
library(Matrix)
i <- c(1,5,2,4,2,2,8)
j <- c(2,5,3,2,4,2,4)
x <- rpois(7,2)
y <- rpois(7,3)
dn1 <- list(mouse_genes[1:8], mouse_bc[1:5])
dn2 <- list(human_genes[1:5], human_bc[1:8])

M1 <- sparseMatrix(i, j, x = x, dimnames = dn1)
M2 <- sparseMatrix(j, i, x = y, dimnames = dn2)

M1
M2

# Lets try a merge even if things dont line up well
M3 <- merge.Matrix(M1, M2, by.x = M1@Dimnames[[1]], by.y = M2@Dimnames[[1]])
M3

# Lets try manipulating the actual code from:
# https://github.com/cvarrichio/Matrix.utils/blob/master/R/Matrix.utils.R#L348

merge.Matrix <- function(x, y, by.x, by.y, all.x = TRUE, all.y = TRUE,out.class=class(x),
                       fill.x=ifelse(is(x,'sparseMatrix'),FALSE,NA),fill.y=fill.x,...)
{
  requireNamespace('grr')
  if(is.null(dim(x)))
    return(grr::matches(by.x,by.y,all.x,all.y,indexes=FALSE))
  indices <- grr::matches(by.x, by.y, all.x, all.y, nomatch = NULL)
  x <- rbind(x, fill.x)
  x <- as(grr::extract(x,indices$x),out.class)
  
  y <- rbind(y, fill.y)
  if(!is.null(colnames(x)) & !is.null(colnames(y)))
    colnames(y)[colnames(y) %in% colnames(x)] <- paste('y',colnames(y)[colnames(y) %in% colnames(x)],sep='.')
  y <- as(grr::extract(y,indices$y),out.class)
  
  result <- cbind2(x, y)
  row.names(result) <- ifelse(row.names(result) == "fill.x", row.names(y), row.names(x))
  return(result)
}

# Lets try a merging with the new function
M3 <- merge.Matrix(M1, M2, by.x = M1@Dimnames[[1]], by.y = M2@Dimnames[[1]])
M3
