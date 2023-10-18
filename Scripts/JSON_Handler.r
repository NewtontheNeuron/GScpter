# JSON_Handler.r contains functions pertained to gathering information from
# config files


loadconfig <- function(fileaddres = ".") {
  
  if (file.exists(paste(fileaddres, ".json", sep = ""))) {
    
    Data <<- fromJSON(file = paste(fileaddres, ".json", sep = ""))
    
  } else if (!file.exists(fileaddres)) {
    
    projectaddress <- paste("../Data/JSON/", fileaddres, ".json", sep = "")
    
    if (!file.exists(projectaddress)) {
      
      return("Please double check that the file or file address is spelt properly.")
      
    } else {
      
      Data <<- fromJSON(file = projectaddress)
      
    }
    
  } else {
    
    Data <<- fromJSON(file = fileaddres)
    
  }
}

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

