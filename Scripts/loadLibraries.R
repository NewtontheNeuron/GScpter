# Function for adding the file to the workspace and saving the workspace.
# This is only necessary outside Rstudio or if loading the data for this first time.
options(repos = list(CRAN="http://cran.rstudio.com/"))

list.of.packages <- c("ggplot2", "Seurat", "patchwork", "tidyverse", "ggdendro", "dplyr", "tidyr", "stringr", "data.table", "tibble", "viridisLite", "Cairo", "datasets", "rjson", "rstudioapi")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(ggdendro)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(tibble)
library(viridisLite)
library(Cairo)
library(datasets)
library(rjson)
library(rstudioapi)