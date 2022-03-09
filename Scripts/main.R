#TODO: need to take arguments for knowing what data file to use..

#something like this
args <- commandArgs(trailingOnly = TRUE)
print(paste ("this is the title project name entered: ", args[1]))
print(paste ("this was the given directory:", args[2]))


#source pre_analysis_functions for all the functions
source("loadLibraries.R")
source("JSON_Handler.R")
source("Pre_analysis_functions.R")

#source scripts to get them to run.
source("DotPlot.R")
source("PooledDotPlot.R")
source("Quad_BarChart2.r")

#load data
RDSfile <- load_data(args[2])

#call each script while passing the data as a value.

#dotplot.R
ListbyClusterAll <- createListByClusterAll(RDSfile)
mainDP(ListByClusterAll)

#pooled dotplot.R
ClusterPoolResults <- createClusterPoolResults(RDSfile)
mainPDP(ClusterPoolResults)

#quadbar_chart2.R
ListByCluster <- createListbyCluster(RDSfile)
mainQBC(ListByCluster)
