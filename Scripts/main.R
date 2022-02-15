#TODO: need to take arguments for knowing what data file to use..

#something like this
args <- commandArgs(trailingOnly = TRUE)
print(paste ("this is the title project name entered: ", args[1]))

#source pre_analysis_functions that will call load data
source("Pre_analysis_functions.R")

#source scripts to get them to run.
source("DotPlot.R")
source("PooledDotPlot.R")
source("Quad_BarChart2.r")


#call each script while passing a value.
ListbyClusterAll <- returnListbyClusterAll()
mainDP(ListByClusterAll)

ClusterPoolResults <- returnClusterpoolResult()
mainPDP(ClusterPoolResults)

ListByCluster <- returnListByCluster()
mainQBC(ListByCluster)