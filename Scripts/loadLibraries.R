# Function for adding the file to the workspace and saving the workspace.
# This is only necessary outside Rstudio or if loading the data for this first time.
options(repos = list(CRAN="http://cran.rstudio.com/"))

list.of.packages <- c("ggplot2", "Seurat", "dplyr", "tidyr", "stringr",
                      "data.table", "cowplot", "Cairo", "rjson", "rstudioapi",
                      "forcats", "rlist", "pipeR")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(Cairo)
library(rjson)
library(rstudioapi)
library(forcats)
library(rlist)
library(pipeR)
