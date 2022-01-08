library(rstudioapi)
library(rjson)

#set working directory to the one this file is currently in
setwd(dirname(getActiveDocumentContext()$path))

#get yaml/json data. Only works if current directory is changed.
Data <- fromJSON(file = "../Data/JSON/Example_2_Subgroups.json")

#return list of features from JSON object.
returnFeatures <- function(){
    return(Data$features)
}

returnClusterpool_names <- function(){
    return(Data$clusterpool_names)
}

returnClusterpool_subgroups <- function(){
    return(Data$subgroup_names)
}

#return list of genes to look at according to specification from JSON object.
#TODO: Add error handling (i.e. what happens when the cell that was called doesn't exist.)
returnClusterpoolGenes <- function(Type, Name){
    listOfGenes <- Data[["clusterpools"]][[Type]][[Name]]
    return(listOfGenes)
}

#return name of the projection from JSON object
returnProjectName <- function(){
    return(Data$project_name)
}

