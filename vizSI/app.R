#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyjs)
library(shinycssloaders)
library(shinyBS)
library(shinydashboard)

library(stringr)
library(data.table)
library(tidyverse)
library(GenomicRanges)
library(ggforce)

library(DBI)
library(ggplot2)
library(gridExtra)

# library(jtoolsscree)

# setwd("./vizSI")
source("./get_missplicing.R")

ui <- fluidPage(
  
  useShinyjs(),
  
  titlePanel("vizSI - visualization of Spliced Introns"),
  
    tabsetPanel(
      
      id = "maintabset", 
      type = "tabs",
      tabPanel(title = "Gene Search",
               
               sidebarLayout(
                 sidebarPanel = sidebarPanel(
                   div(
                     id = "inputPanel",
                     h3("Gene of interest:"),
                     hr(),
                     selectizeInput(inputId = "gene", 
                                    label = "Gene:", 
                                    choices = NULL, 
                                    multiple = TRUE,
                                    options = list(
                                      placeholder = "Search by gene",
                                      maxItems = 1,
                                      options = list(create = FALSE)),
                                    selected = NULL),
                     hr(),
                     selectizeInput(inputId = "geneTissue", "Tissue:", 
                                    choices = tissue_GTEx_choices_alphabetical,
                                    multiple = F, 
                                    options = list(maxItems = 1), 
                                    selected = "Brain-FrontalCortex_BA9"),
                     hr(),
                     actionButton(inputId = "geneButton", label = "Accept"),
                     hidden(
                       p(id = "intronID", ""),
                       p(id = "novelID", "")
                       )
                     ),
                   div(
                     id = "intronPanel",
                     uiOutput("intronPanelOutput")
                     ),
                   width = 3
                ),
                mainPanel = mainPanel(
                  uiOutput("geneOutput") %>% withSpinner(color="#0dc5c1"),
                  uiOutput("intronGeneDetail"),
                  
                  bsModal(id = "modalDetailsIntron", 
                          title = NULL, 
                          trigger = NULL, 
                          size = "large", 
                          uiOutput("modalIntronDetail")),
                  bsModal(id = "modalDetailsNovel", 
                          title = NULL, 
                          trigger = NULL, 
                          size = "large", 
                          uiOutput("modalNovelDetail")),
                  width = 9
               )
      )),
      tabPanel(title = "Plots",
               sidebarLayout(
                sidebarPanel(
                  h3("Plot generation:"),
                  hr(),
                  selectizeInput(inputId = "tissueDistances", 
                                 "Tissue:",
                                 choices = tissue_GTEx_choices_alphabetical,
                                 multiple = F, 
                                 options = list(maxItems = 1), 
                                 selected = "Brain-FrontalCortex_BA9"),
                  #numericInput(inputId = "bp", label = "Distance (in bp):", value = 30, min = 1, max = 100000),
                  hr(),
                  actionButton(inputId = "acceptBtnPlots", label = "Accept"),
                  width = 3
                ),
      
                # Show a plot of the generated distribution
                mainPanel(
                  plotOutput("distancesOutput", height = 600) %>% withSpinner(color="#0dc5c1"),
                  plotOutput("moduloOutput", height = 600) %>% withSpinner(color="#0dc5c1"),
                  plotOutput("missplicingOutput", height = 600) %>% withSpinner(color="#0dc5c1"),
                  width = 9
                )
              )
      )
    
  )
)



server <- function(input, output, session) {
   

  ## Fill the dropdown with the list of genes --------------------------------------------------------------------------------------
  
  updateSelectizeInput(session, 'gene', choices = genes, server = TRUE, selected = "GBA")
  
  
  
  ## Get the list of novel junctions attached to an annotated intron - shown in a different tab ------------------------------------
  
  observe({
    
    cdata <- parseQueryString(session$clientData$url_search)
    
    # print(paste0("Intron '", input$intronID,"' details:"))
    # print(cdata[['intron']])
    
    if (!is.null(cdata[['intron']])) {

      hideElement(id = "inputPanel")
      showElement(id = "intronPanel")
      
      
      ## Text output with the info from the selected intron 
      
      output$intronPanelOutput <- renderUI({
        
        tagList(
          h3("Data selection:"),
          hr(),
          p(strong("Gene: "), paste0(cdata[['gene']],".")),
          p(strong("Intron: "),paste0("'", cdata[['coordinates']],"'.")),
          p(strong("Length: "),paste0(cdata[['length']]," bp.")),
          p(strong("Mis-splicing type: "),paste0("'", cdata[['type']],"' splice site.")),
          hr(),
          p(strong("Project: "),paste0("'", cdata[['tissue']],"'."))
        )

        
      })
      
      
      output$intronGeneDetail <- renderUI({
      
        tagList(
        
          h2(paste0("Novel events from intron '", cdata[['coordinates']],"':")),
          
          DT::datatable(get_novel_data(intron_ID = cdata[['intron']],
                                       tissue = cdata[['tissue']]), 
                        options = list(pageLength = 20,
                                       order = list(8, 'desc'),
                                       columnDefs = list(list(visible=FALSE, targets=c(2))),
                                       autoWidth = F,
                                       rowGroup = list(dataSrc = 1),
                                       rowCallback = DT::JS("function(row, data) {
                                                              var onclick_f = 'Shiny.setInputValue(\"novelID\",\"' + data[2] + '\");$(\"#modalDetailsNovel\").modal(\"show\");';
                                                              //var num = '<a id=\"goA\" role=\"button\" onclick = ' + onclick_f + ' >' + data[8] + '</a>';
                                                              var num = data[8];
                                                              $('td:eq(7)', row).html(num);
                                                            }"
                                       )),
                        #extensions = 'RowGroup', 
                        width = "100%",
                        rownames = FALSE)
        )
        
      })
    }
    
  })
  
  ## Popup modal window with the detail of the individuals presenting an annotated intron ------------------------------------------
  
  output$modalIntronDetail <- renderUI({
    
    #cdata <- parseQueryString(session$clientData$url_search)
    
    tagList(
      
      DT::datatable(get_intron_details(intron_id = input$intronID,
                                       tissue = input$geneTissue), 
                    options = list(pageLength = 20,
                                   autoWidth = T,
                                   order = list(2, 'desc')),
                    width = "100%",
                    rownames = FALSE)
    )
    
  })
  
  ## Popup modal window with the detail of the individuals presenting a novel junction ---------------------------------------------
  
  output$modalNovelDetail <- renderUI({
    
    cdata <- parseQueryString(session$clientData$url_search)
    
    tagList(
      
      DT::datatable(get_novel_details(novel_id = input$novelID,
                                      tissue = cdata[['tissue']]), 
                    options = list(pageLength = 20,
                                   autoWidth = T,
                                   order = list(2, 'desc')),
                    width = "100%",
                    rownames = FALSE)
    )
    
  })
  
  
  
  
  ## Get all annotated introns from the selected gene -----------------------------------------------------------------------------

  geneSearchUI <- eventReactive(input$geneButton, {
    
    table3_caption = paste0("Mean values have been calculated across all samples from the selected tissue. MSR_D = Mean MisSplicing Ratio at Donor (5'ss) positions. 
                            MSR_A = Mean MisSplicing Ratio at Acceptor (3'ss) positions.")
    
    tagList(

      
      h1(input$gene),
      h2("Annotated introns (ensembl v104 - March 2021):"),
      br(),
      
      DT::datatable(get_gene_intron_data(input$gene, input$geneTissue), 
                    options = list(pageLength = 20,
                                   columnDefs = list(list(visible=FALSE, targets=c(2))),
                                   order = list(1, 'asc'),
                                   rowGroup = list(dataSrc = 1),
                                   autoWidth = F,
                                   rowCallback = DT::JS(
                                     "function(row, data) {
                                        
                                        var href = encodeURI('https://soniagarciaruiz.shinyapps.io/vizsi/?intron=' + data[2] + '&coordinates=' + data[0] + '&gene=' + data[10] + '&type=' + data[1] + '&length=' + data[3] + '&tissue=' + $(\"#geneTissue\").val());
                                        var num = '<a id=\"goA\" role=\"button\" target=\"_blank\" href=' + href + '>' + data[0] + '</a>';
                                        $('td:eq(0)', row).html(num);


                                        var onclick_f = 'Shiny.setInputValue(\"intronID\",\"' + data[2] + '\");$(\"#modalDetailsIntron\").modal(\"show\");';
                                        //var num = '<a id=\"goA\" role=\"button\" onclick = ' + onclick_f + ' >' + data[7] + '</a>';
                                        var num = data[7];
                                        $('td:eq(6)', row).html(num);

                                       
                                     }"
                                   )
                                   ),
                    extensions = 'RowGroup', 
                    width = "100%",
                    #selection = 'none',
                    rownames = F,
                    caption = htmltools::tags$caption(
                      style = 'caption-side: bottom;',
                      table3_caption
                    )
                    ),
      br(), 
      br(),
      br(),
      br()
    )
  })
  
  
  ## What happens after the user clicks the 'Accept' button -----------------------------------------------------------------------
  
  output$geneOutput = renderUI({
    # showElement(id = "inputPanel")
    # hideElement(id = "intronPanel")
    geneSearchUI()
  })
  
  
  
  ####################################################
  ################# PLOTS TABS #######################
  ####################################################
  
  ## DISTANCES
  plotDistances <- eventReactive(input$acceptBtnPlots, {
    plot_distances(input$tissueDistances)
  })
  plotModulo <- eventReactive(input$acceptBtnPlots, {
    plot_modulo(input$tissueDistances)
  })
  plotMissplicing <- eventReactive(input$acceptBtnPlots, {
                   plot_missplicing(input$tissueDistances)
  })


  output$distancesOutput = renderPlot(expr = {
    plotDistances()
  }, width = 600, height = 600)
  output$moduloOutput = renderPlot(expr = {
    plotModulo()
  }, width = 600, height = 600)
  output$missplicingOutput = renderPlot(expr = {
    plotMissplicing()
  }, width = 600, height = 600)
  
}

# Run the application 
shinyApp(ui = ui, server = server)

