#shinyapp/Scripts/app.R

library(shiny)
library(data.tree)
#(jsonlite)

#source pre_analysis_functions for all the functions
source("loadLibraries.R")
source("JSON_Handler.r")
source("Pre_analysis_functions.R")

#source scripts to get them to run.
source("DotPlot.R")
source("PooledDotPlot.R")

#### Define elements ####
#(inputId = "num", "Enter number of values:", min = 1, max = 100, value = 50)
#select type of data
selectIn <- selectInput(inputId = "select", label = "Mouse or Human Data?", choices = c("Mouse", "Human"))

#text in
textInGene <- selectizeInput(inputId = "geneText", label = "Gene:", choices = NULL)
textInClusterpool <- textInput(inputId = "clusterpoolText", label = "Clusterpool:")
textInSubgroup <- textInput(inputId = "subgroupText", label = "Subgroup:")
textInCluster <- textInput(inputId = "clusterText", label = "Cluster: ")

# tree display and set
listClusterpool <- selectInput(inputId = "clusterpoolList", label = "Clusterpool", choices = NULL, multiple = TRUE, selectize = FALSE)
listSubgroup <- selectInput(inputId = "subgroupList", label = "Subgroup", choices = NULL, multiple = TRUE, selectize = FALSE)
clusterSelect <- checkboxGroupInput(inputId = "selectCluster", label = "Clusters:", choices = 1:12, inline = TRUE, width = "400px")
myTree <- Node$new("pools")

# Adding variables of grouping
listofvariables <- selectInput(inputId = "listofvariables", label = "Variables:",
                               choices = NULL, multiple = TRUE, selectize = FALSE)
listoflevels <- selectInput(inputId = "listoflevels", label = "Levels:", choices = NULL,
                            multiple = TRUE, selectize = FALSE)
# Add or Remove grouping layers
addLayer <- actionButton(inputId = "addVariable", "Add layer(s)")
changeLayer <- actionButton(inputId = "cngVariable", "Change layer with selected variables")
removeLayer <- actionButton(inputId = "rmVariable", "Remove layer(s)")

# recode variables
recodeInput <- textInput(inputId = "recodeInput", label = NULL)
recodeLevel <- actionButton(inputId = "recodeLevel", "Recode level")

# list of grouping layers
listoflayers <- selectInput(inputId = "listoflayers", label = "Grouping layers:",
                            choices = NULL, multiple = TRUE, selectize = FALSE)

# layer movement buttons
layer_to_top <- actionButton(inputId = "layertotop", label = "top")
layer_to_bottom <- actionButton(inputId = "layertobottom", label = "bottom")
layer_up <- actionButton(inputId = "layertoup", label = "up")
layer_down <- actionButton(inputId = "layertodown", label = "down")

# level movement buttons
level_to_top <- actionButton(inputId = "leveltotop", label = "top")
level_to_bottom <- actionButton(inputId = "leveltobottom", label = "bottom")
level_up <- actionButton(inputId = "leveltoup", label = "up")
level_down <- actionButton(inputId = "leveltodown", label = "down")

#data file in
options(shiny.maxRequestSize = 10 * 10^9)
fileIn <- fileInput(inputId = "file", label = "Input Data File")

#action buttons
geneButton <- actionButton(inputId = "geneButton", "Add")
clusterpoolButton <- actionButton(inputId = "clusterpoolButton", "Add")
subgroupButton <- actionButton(inputId = "subgroupButton", "Add")
clusterButton <- actionButton(inputId = "clusterButton", "Add")
button <- actionButton(inputId = "button", "Generate Plots")

#data output
geneOut <- htmlOutput(outputId = "genePrint")
clusterpoolOut <- textOutput("clusterpoolPrint")
subgroupOut <- textOutput("subgroupPrint")
clusterOut <- textOutput("clusterPrint")

#plots to output
textOut <- htmlOutput("fileprogress")
dotplotOut <- plotOutput(outputId = "dotplot")
pooleddotplotOut <- plotOutput(outputId = "pooleddotplot")

ui <- fluidPage(
  
  titlePanel("Clusterpool Comparison Tool"),
  # Link to the css style sheet
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "stylesheet.css")
  ),
  # The row divider for step 1 of the process
  h2("Step 1"),
  fluidRow(
    # The upload RDS file button
    # TODO: configure clusters/main groups
    column(6, fileIn,
           textOut),
    # The use existing dataset button
    column(6, selectIn)
  ),
  # The row divider for step 2 of the process
  fluidRow(
    column(12,
      # The selection and input
      h2("Step 2"),
      flowLayout(
        textInGene,
        div(style = "margin-top: 25px", geneButton),
      ),
      div(style = "background-color: #d1d1d1; width: 100%; height: 75px", geneOut),
      # The output
      flowLayout(
        textInClusterpool,
        div(style = "margin-top: 25px", clusterpoolButton),
        listClusterpool),
      flowLayout(
        textInSubgroup,
        div(style = "margin-top: 25px", subgroupButton),
        listSubgroup),
      flowLayout(clusterSelect)
    )
  ),
  # The row divider for step 3 of the process
  verticalLayout(
    h2("Step 3"),
    p("Configure more layers of grouping"),
    flowLayout(
      verticalLayout(
        listofvariables,
        addLayer
      ),
      verticalLayout(
        fluidRow(
          column(7, listoflevels),
          column(5,
                 level_to_top,
                 level_up,
                 level_down,
                 level_to_bottom)
        ),
        fluidRow(
          column(6, recodeInput),
          column(6, recodeLevel)
        )
      ),
      verticalLayout(
        fluidRow(
          column(7, listoflayers),
          column(5,
                 layer_to_top,
                 layer_up,
                 layer_down,
                 layer_to_bottom)
        ),
        changeLayer,
        removeLayer
      )
    )
  ),
  # The row divider for step 4 of the process
  h2("Finalize"),
  button
)

#### Functions ####
# Function to update the layer box when things change
updatelayerbox <- function(session){
  display_l <<- lapply(seq_along(list_of_layers),
                      function(x) paste(names(list_of_layers[x]), ":",
                                        paste(unlist(list_of_layers[[x]]),
                                              collapse = " ")))
  display_l <<- display_l[-1]
  updateSelectInput(session, "listoflayers", choices = unlist(display_l))
}

JSON_Object = {}

#### Server ####
#inputs are things that the user can change / interact with
#output are things that a user can see like a plot 
server <- function(input, output, session) {
  
  #lists of data
  myValues <- reactiveValues()
  
  # Load the data for everyone to use
  observeEvent(input$file, {
    output$fileprogress <- renderUI({
      div(class = "loading-file", "The file is loading. Please wait.")
    })
    # Running the function to load the file
    RDfile <<- load_data(input$file$datapath)
    genes <<- unlist(RDfile@assays$RNA@counts@Dimnames)
    clusternames <<- unique(unlist(as.data.frame(RDfile@active.ident)[1]))
    variables <<- names(RDfile@meta.data)
    updateSelectizeInput(session, 'geneText', choices = genes,
                         server = TRUE, # multiple = TRUE,
                         options = list(maxOptions = 7, placeholder = "Type in your gene..."))
    updateCheckboxGroupInput(session, "selectCluster", choices = clusternames,
                             inline = TRUE)
    updateSelectInput(session, "listofvariables", choices = variables)
    # Creating list of layers
    create_listoflayers()
    # Creating modified levels list
    create_mod_lev_lyr()
    # Change the message
    #if(exists(x = RDfile)) {
    output$fileprogress <- renderUI({
      div(class = "loaded-file", sprintf("File loaded: %s", input$file$datapath))
    })
    #}
  })
  
  #add gene to list
  observeEvent(input$geneButton, {
    myValues$genes <- c(myValues$genes, input$geneText)
    geneOut <<- tagAppendChild(geneOut,
                               div(class = "gene_boxes",
                                            myValues$genes[length(myValues$genes)]))
    output$genePrint <- renderUI(geneOut)
  })
  
  #add clusterpool to list
  observeEvent(input$clusterpoolButton, {
    myValues$clusterpools <- c( myValues$clusterpools, input$clusterpoolText)
    updateSelectInput(session, "clusterpoolList", choices = myValues$clusterpools)
    assign(input$clusterpoolText, myTree$AddChild(input$clusterpoolText))
    print(myTree, "clus")
  })

  #add subgroup to list
  observeEvent(input$subgroupButton, {
    myValues$subgroups <- c( myValues$subgroups, input$subgroupText)
    updateSelectInput(session, "subgroupList", choices = myValues$subgroups)
    myTree$Do(function(node) assign(input$subgroupText,
                                    node$AddChild(input$subgroupText, clus = NA)),
              filterFun = function(x) x$level == 2)
    print(myTree, "clus")
  })
  
  # fetch clusters
  observeEvent(input$clusterpoolList, {
    observeEvent(input$subgroupList, {
      wcp <<- input$clusterpoolList
      wsg <<- input$subgroupList
      print(wcp)
      print(wsg)
      if(length(wcp) != 1 | length(wsg) != 1){
        print("I cannot just get one thing. Do you want to remove them?")
      }
      clusterstobe <- myTree$Get(filterFun = function(self) all(self$path == c("pools", wcp, wsg)),
                                 attribute = "clus")
      print(as.vector(clusterstobe))
      updateCheckboxGroupInput(session, "selectCluster",
                               selected = as.vector(clusterstobe))
    })
  })
  # Save the input
  #TODO: Error when nothing is selected
  observeEvent(input$selectCluster, {
    print(c("pools", wcp, wsg))
    clusters <- input$selectCluster
    #print(myTree)
    myTree$Set(clus = list(clusters),
               filterFun = function(self) all(self$path == c("pools", wcp, wsg)))
    print(myTree, "clus", "path")
  })

  #when button is pressed, generate all plots
  observeEvent(input$button, {
    #dotplot.R
    output$dotplot <- renderPlot({
      title <- "Dot Plot"
      #RDSfile <- load_data(input$file$datapath)
      ListByClusterAll <- createListByClusterAll(RDfile)
      mainDP(ListByClusterAll)
    })
    
    #pooled dotplot.R
    output$pooleddotplot <- renderPlot({
      
      JSON_Object[["features"]] = myValues$genes
      #write to json so data can properly be loaded by json_handler.R
      #write_json(JSON_Object, file = ShinyData.json)
      
      
      title <- "Pooled Dot Plot"
      #RDSfile <- load_data(input$file$datapath)
      ClusterPoolResults <- createClusterPoolResults(RDfile)
      mainPDP(ClusterPoolResults)
    }) 
  })
  
  # Step 3 add, remove variables to grouping layers and recode levels
  # The variables section should be populated when the file is loaded
  observeEvent(input$listofvariables, {
    if(length(input$listofvariables) == 1) {
      levels <- unique(RDfile@meta.data[[input$listofvariables]])
      if(length(levels) < 100) {
        updateSelectInput(session, "listoflevels", choices = levels)
      }
    } else {
      updateSelectInput(session, "listoflevels", choices = NULL)
    }
  })
  
  # Add a layer based on the selected variables
  observeEvent(input$addVariable, {
    if(length(input$listofvariables) > 0) {
      #str(input$listofvariables)
      list_of_layers <<- add_layer(layer_list = list_of_layers,
                                   newLayerItems = list(input$listofvariables))
      updatelayerbox(session)
    }
  })
  
  # Remove a layer based on the selected layer
  observeEvent(input$rmVariable, {
    if(length(input$listoflayers) == 1) {
      index <- which(display_l %in% input$listoflayers)
      list_of_layers <<- rm_layer(layer_list = list_of_layers,
                                  layer_number = index)

      updatelayerbox(session)
    }
  })
  # change a layer based on the selected variables
  observeEvent(input$cngVariable, {
    print("Change function - list of layers then variables")
    print(input$listoflayers)
    print(input$listofvariables)
    if(length(input$listoflayers) == 1 & length(input$listofvariables) > 0){
      index <- which(display_l %in% input$listoflayers)
      print("Change function - index")
      print(index)
      print("Change function - display l")
      print(display_l)
      print(list(input$listofvariables))
      cng_lr <- change_layer(list_of_layers,
                             layer_number = index,
                             newLayerItems = list(input$listofvariables))
      print("Change function - index from changelayer paf.r")
      print(cng_lr)
      print("Change function - the list of layers")
      print(list_of_layers)
      list_of_layers <<- cng_lr$main
      print(cng_lr$log)
      updatelayerbox(session)
    }
  })
  
  # Movement of a layer to the top, bottom, 1 increment up or down of the stack
  observeEvent(input$layertotop, {
    index <- which(display_l %in% input$listoflayers)
    list_of_layers <<- reorder_layer(list_of_layers, index, "top")
    updatelayerbox(session)
  })
  observeEvent(input$layertobottom, {
    index <- which(display_l %in% input$listoflayers)
    list_of_layers <<- reorder_layer(list_of_layers, index, "bottom")
    updatelayerbox(session)
  })
  observeEvent(input$layertoup, {
    index <- which(display_l %in% input$listoflayers)
    list_of_layers <<- reorder_layer(list_of_layers, index, "up")
    updatelayerbox(session)
  })
  observeEvent(input$layertodown, {
    index <- which(display_l %in% input$listoflayers)
    list_of_layers <<- reorder_layer(list_of_layers, index, "down")
    updatelayerbox(session)
  })
  
  # Movement of a level to the top, bottom, 1 increment up or down of the stack
  observeEvent(input$leveltotop, {
    levels <- as.list(levels(RDfile@meta.data[[input$listofvariables]]))
    level <- input$listoflevels
    variable <- input$listofvariables
    modified_levels[variable] <<- reorder_level(levels, level, "top")
    updateSelectInput(session, "listoflevels", choices = levels)
  })
  observeEvent(input$layertobottom, {
    levels <- as.list(unique(RDfile@meta.data[[input$listofvariables]]))
    level <- input$listoflevels
    variable <- input$listofvariables
    modified_levels[[variable]] <<- reorder_level(levels, level, "bottom")
    updatelayerbox(session)
  })
  observeEvent(input$layertoup, {
    levels <- as.list(unique(RDfile@meta.data[[input$listofvariables]]))
    level <- input$listoflevels
    variable <- input$listofvariables
    modified_levels[[variable]] <<- reorder_level(levels, level, "up")
    updatelayerbox(session)
  })
  observeEvent(input$layertodown, {
    levels <- as.list(unique(RDfile@meta.data[[input$listofvariables]]))
    level <- input$listoflevels
    variable <- input$listofvariables
    modified_levels[[variable]] <<- reorder_level(levels, level, "down")
    updatelayerbox(session)
  })
  
  # TODO: create the back-end for this part.
  # TODO: ? help hover box
  # TODO: Need to know that something happens after you add clup or subgr
  # TODO: Example of layering
  # TODO: remove level
  # TODO: the layers start over after hitting change layer
  # TODO: Error when there is no cluster pools selected or present wcp not found
  # TODO: Crash on layer movement when there is no layer selected.

  
}

shinyApp(ui = ui, server = server)
