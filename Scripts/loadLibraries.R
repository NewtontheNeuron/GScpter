options(repos = list(CRAN="http://cran.rstudio.com/"))

list.of.packages <- c("ggplot2", "Seurat", "dplyr", "tidyr", "stringr", "Cairo", "rjson", "rstudioapi")
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