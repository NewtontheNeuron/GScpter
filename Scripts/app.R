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
                            multiple = FALSE)
# Add or Remove grouping layers
addLayer <- actionButton(inputId = "addVariable", "Add layer(s)")
removeLayer <- actionButton(inputId = "rmVariable", "Remove layer(s)")

# recode variables
recodeInput <- textInput(inputId = "recodeInput", label = NULL)
recodeLevel <- actionButton(inputId = "recodeLevel", "Recode level")

# list of grouping layers
listoflayers <- selectInput(inputId = "listoflayers", label = "Grouping layers:",
                            choices = NULL, multiple = TRUE, selectize = FALSE)


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
        listoflevels,
        fluidRow(
          column(6, recodeInput),
          column(6, recodeLevel)
        )
      ),
      verticalLayout(
        listoflayers,
        removeLayer
      )
    )
  ),
  # The row divider for step 4 of the process
  h2("Finalize"),
  button
)

JSON_Object = {}
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
    updateCheckboxGroupInput(session, "selectCluster", choices = clusternames)
    updateSelectInput(session, "listofvariables", choices = variables)
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
    geneOut <<- tagAppendChild(geneOut, div(class = "gene_boxes", myValues$genes[length(myValues$genes)]))
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
    myTree$Do(function(node) assign(input$subgroupText, node$AddChild(input$subgroupText, clus = NA)), filterFun = function(x) x$level == 2)
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
      clusterstobe <- myTree$Get(filterFun = function(self) all(self$path == c("pools", wcp, wsg)), attribute = "clus")
      print(as.vector(clusterstobe))
      updateCheckboxGroupInput(session, "selectCluster", selected = as.vector(clusterstobe))
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
  # TODO: create the backend for this part.

  
}

shinyApp(ui = ui, server = server)
