#shinyapp/Scripts/app.R

library(shiny)
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
textInGene <- textInput(inputId = "geneText", label = "Gene:")
textInClusterpool <- textInput(inputId = "clusterpoolText", label = "Clusterpool:")
textInSubgroup <- textInput(inputId = "subgroupText", label = "Subgroup:" )
textInCluster <- textInput(inputId = "clusterText", label = "Cluster: ")

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
geneOut <- textOutput("genePrint")
clusterpoolOut <- textOutput("clusterpoolPrint")
subgroupOut <- textOutput("subgroupPrint")
clusterOut <- textOutput("clusterPrint")

#plots to output
textOut <- textOutput("text")
dotplotOut <- plotOutput(outputId = "dotplot")
pooleddotplotOut <- plotOutput(outputId = "pooleddotplot")

ui <- fluidPage(
  
  titlePanel("Clusterpool Comparison Tool"),
  
  
  #fluid row columns need to add up to 12, integers only.
  fluidRow(
    
    column(3,
        selectIn,
        textInGene,
        geneOut,
        textInClusterpool,
        clusterpoolOut,
        textInSubgroup,
        subgroupOut,
        textInCluster,
        clusterOut,
        fileIn,
    ),
    column(1,
        br(),
        br(),
        br(),
        br(),
        br(),
        geneButton,
        br(),
        br(),
        br(),
        clusterpoolButton,
        br(),
        br(),
        br(),
        subgroupButton,
        br(),
        br(),
        br(),
        clusterButton
    ),
    column(8,
        button,
        dotplotOut, pooleddotplotOut
           
    )
  ),
  
  textOut, 
)

JSON_Object = {}
#inputs are things that the user can change / interact with
#output are things that a user can see like a plot 
server <- function(input, output) {
  
  #lists of data
  myValues <- reactiveValues()
  
  #add gene to list
  observeEvent(input$geneButton, {
    myValues$genes <- c(myValues$genes, input$geneText)
    output$genePrint <- renderText(myValues$genes)
  })
  
  #add clusterpool to list
  observeEvent(input$clusterpoolButton, {
    myValues$clusterpools <- c( myValues$clusterpools, input$clusterpoolText)
    output$clusterpoolPrint <- renderText(myValues$clusterpools)
  })

  #add subgroup to list
  observeEvent(input$subgroupButton, {
    myValues$subgroups <- c( myValues$subgroups, input$subgroupText)
    output$subgroupPrint <- renderText(myValues$subgroups)
  })
  
  output$text <- renderText({ sprintf("File loaded: %s", input$file$datapath) })

  #when button is pressed, generate all plots
  observeEvent(input$button, {
    #dotplot.R
    output$dotplot <- renderPlot({
      title <- "Dot Plot"
      RDSfile <- load_data(input$file$datapath)
      ListByClusterAll <- createListByClusterAll(RDSfile)
      mainDP(ListByClusterAll)
    })
    
    #pooled dotplot.R
    output$pooleddotplot <- renderPlot({
      
      JSON_Object[["features"]] = myValues$genes
      #write to json so data can properly be loaded by json_handler.R
      #write_json(JSON_Object, file = ShinyData.json)
      
      
      title <- "Pooled Dot Plot"
      RDSfile <- load_data(input$file$datapath)
      ClusterPoolResults <- createClusterPoolResults(RDSfile)
      mainPDP(ClusterPoolResults)
    }) 
  })

  
}

shinyApp(ui = ui, server = server)
