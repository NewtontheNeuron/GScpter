#get needed libraries
source("loadLibraries.R")

#hard coded for clean_neuron_object and top_level_new_annotation.rda
main <- function(){
    args <- commandArgs(trailingOnly = TRUE)
    if (str_detect(args[1], ".RDS")) {
        print(args[1])

        return
        clean_neuron_object <- readRDS(args[1])
        createJSONRodent(clean_neuron_object)
    } else {
        print(args[1])

        return
        load(file = args[1])

        createJSONHuman(integrated_top_level_obj)
    }
}

createJSONRodent <- function(RObject){
    #for rodent dataset:

    #all clusters in the data
    list_of_clusters <- levels(RObject@active.ident)
    unique_list_of_clusters <- sort(unique(list_of_clusters))

    #all genes in the data
    list_of_genes <- RObject@assays$RNA@data@Dimnames[1]
    unique_list_of_genes <- unique(list_of_genes)

    data <- vector(mode="list", length = 2)
    data[[1]] <- unique_list_of_clusters
    data[[2]] <- unique_list_of_genes

    json_data <- toJSON(data)

    write(json_data, "../gui_data.json")

}

createJSONHuman <- function(RObject){
    #for human dataset:

    #all clusters in the data
    list_of_clusters <- RObject@meta.data$new_annotation
    unique_list_of_clusters <- sort(unique(list_of_clusters))

    #all genes in the data
    list_of_genes <- RObject@assays$RNA@data@Dimnames[1]
    unique_list_of_genes <- unique(list_of_genes)

    data <- vector(mode="list", length = 2)
    data[[1]] <- unique_list_of_clusters
    data[[2]] <- unique_list_of_genes

    json_data <- toJSON(data)

    write(json_data, "../gui_data.json")

}

main()


