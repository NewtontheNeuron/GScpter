#set working directory to the one this file is currently in
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#library(rjson)

#get yaml/json data. Only works if current directory is changed.
Data <- fromJSON(file = paste("../Data/JSON/", args[1], ".json", sep = ""))

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
returnClusters <- function(Type, Name){
    listOfClusters <- Data[["clusterpools"]][[Type]][[Name]]
    return(listOfClusters)
}

returnAllClusters <- function(){
  
    #return all clusters in every clusterpool.
    groups <- returnClusterpool_names()
    subgroups <- returnClusterpool_subgroups()
    
    allClusters <- c()
    
    for (group in groups){
      for (subgroup in subgroups){
        Clusterpool <- returnClusters(group, subgroup)
        allClusters <- c(allClusters, Clusterpool )
      }
    }
    return(allClusters)
  
}

#return name of the projection from JSON object
returnProjectName <- function(){
    return(Data$project_name)
}

returnCleanLabelList <- function(){
  features <- returnFeatures()
  clean_label_list <- str_remove(features, 'rna_')
  return(clean_label_list)
}

