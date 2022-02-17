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

library(DBI)
library(ggplot2)
library(gridExtra)

# setwd("./vizSI")
source("./get_missplicing.R")

ui <- function(){
  
  fluidPage(
    
    fluidRow(
      tagList(
        #hidden(
          div(id = "searchDiv",
            
            br(),br(),br(),br(),br(),
            h2(id = "title", "intronDB"),
            br(),
            wellPanel(id = "searchRow",
                      selectizeInput(inputId = "gene", 
                                     label = "Gene:", 
                                     choices = NULL, 
                                     multiple = TRUE,
                                     selected = NULL,
                                     options = list(
                                       placeholder = "Search by gene",
                                       maxItems = 1,
                                       options = list(create = FALSE))),
                      hr(),
                      selectizeInput(inputId = "geneTissue", "Tissue:", choices = tissue_GTEx_choices_alphabetical,
                                     multiple = F, options = list(maxItems = 1), selected = "Brain-FrontalCortex_BA9"),
                      hr(),
                      actionButton(inputId = "geneButton", label = "Accept"),
                      hidden(
                        p(id = "intronID", "")
                        )
            ),
            tags$style(type="text/css", "#searchDiv {text-align:left; margin:auto; width:30%;} #title{text-align:center; margin:auto; width:30%;}")
            ),
        #),
        br(),
        br()
    )),
    fluidRow(
      column(1),
      column(10,
             tagList(
               div(id = "resultsDiv",
                   uiOutput("geneOutput") %>% withSpinner(color="#0dc5c1"),
                   tags$style(type="text/css", "#resultsDiv {text-align:left; margin:auto; width:100%;}")),
               hidden(
                 div(id = "results2Div",
                     uiOutput("intronGeneDetail"),
                     bsModal(id = "modalIntronGene", 
                             title = "Intron detail:", 
                             trigger = NULL, 
                             size = "large", 
                             uiOutput("modalNovelDetail")),
                     tags$style(type="text/css", "#results2Div {text-align:left; margin:auto; width:100%;}"))
               )
      )),
      column(1)
      
    )
  )
}
#ui = (htmlOutput("page"))

# ui <- fluidPage(
#   useShinyjs(),
#   
#   # Application title
#   titlePanel("vizSI - visualization of Spliced Introns"),
#    
#   # Sidebar with a slider input for number of bins 
#   tabsetPanel(
#     
#     id = "maintabset", type = "tabs",
# 
#     tabPanel("Intron Search",
#       
#       sidebarLayout(
#        sidebarPanel(
#          h3("Intron of interest:"),
#          selectizeInput(inputId = "chr", "Chr:", choices = chr_choices,  multiple = F, options = list(maxItems = 1), selected = "10"),
#          numericInput(inputId = "start", label = "Start:", value = 133366994, min = 1, max = 1000000000),
#          numericInput(inputId = "end", label = "End:", value = 133368922, min = 1, max = 1000000000),
#          selectizeInput(inputId = "strand", "Strand:", choices = strand_choices,  multiple = F, options = list(maxItems = 1), selected = "-"),
#          hr(),
#          selectizeInput(inputId = "tissue", "Tissue:", choices = tissue_GTEx_choices_alphabetical,  multiple = F, 
#                         options = list(maxItems = 1), selected = "Brain-FrontalCortex_BA9"),
#          checkboxInput(inputId = "all_tissues", label = "All Tissues", value = FALSE),
#          hr(),
#          actionButton(inputId = "update", label = "Accept"), 
#          hidden(
#            p(id = "intronID", ""),
#            p(id = "novelID", "")
#          ),
#          width = 3
#        ),
#         
#         # Show a plot of the generated distribution
#         mainPanel(
#       
#           uiOutput("intronOutput") %>% withSpinner(color="#0dc5c1"), 
#           bsModal(id = "modalIntron", title = "Intron detail:", trigger = "goA", 
#                   size = "large", uiOutput("intronDetail")),
#           bsModal(id = "modalNovel", title = "Novel juntion detail:", trigger = "goB", 
#                   size = "large", uiOutput("novelDetail")),
#           width = 9
#           
#         )
#      )
#    ),
#    tabPanel("Gene Search",
#             
#             sidebarLayout(
#               sidebarPanel(
#                 h3("Gene of interest:"),
#                 selectizeInput(inputId = "gene", "Gene:", choices = NULL,  multiple = F, 
#                                options = list(maxItems = 1), selected = "ENSG00000145335"),
#                 hr(),
#                 selectizeInput(inputId = "geneTissue", "Tissue:", choices = tissue_GTEx_choices_alphabetical,  
#                                multiple = F, options = list(maxItems = 1), selected = "Brain-FrontalCortex_BA9"),
#                 hr(),
#                 actionButton(inputId = "geneButton", label = "Accept"), 
#                 width = 3
#               ),
#               
#               # Show a plot of the generated distribution
#               mainPanel(
#                 uiOutput("geneOutput") %>% withSpinner(color="#0dc5c1"), 
#                 width = 9
#               )
#             )
#    )
#    # ,
#    # tabPanel("Plots",
#    # 
#    #          sidebarLayout(
#    #            sidebarPanel(
#    #              selectizeInput(inputId = "tissueDistances", "Tissue:", 
#    #                             choices = tissue_GTEx_choices_alphabetical,  
#    #                             multiple = F, options = list(maxItems = 1), selected = "Brain-FrontalCortex_BA9"),
#    #              numericInput(inputId = "bp", label = "Distance (in bp):", value = 30, min = 1, max = 100000),
#    #              hr(),
#    #              actionButton(inputId = "acceptBtnDist", label = "Accept"),
#    #              width = 3
#    #            ),
#    # 
#    #            # Show a plot of the generated distribution
#    #            mainPanel(
#    #              plotOutput("plotsOutput") %>% withSpinner(color="#0dc5c1"),
#    #              width = 9
#    #            )
#    #          )
#    # )
#   )
# )



server <- function(input, output, session) {
   
  #shinyjs::disable(id = "all_tissues")
  
  # observe({
  #   cdata <- parseQueryString(session$clientData$url_search)
  #   print(cdata[['intron']])
  #   
  #   # if (!is.null(cdata[['intron']])) {
  #   #   show(input$searchDiv)
  #   # } else {
  #   #   hide(input$searchDiv)
  #   # }
  # })
  #shinyjs::hide(id = "searchRow")
  
  #clickedNovel <- NULL
  
  ############################################################
  ##################### INTRON SEARCH ########################
  ############################################################
  
  # getPageTitle <- eventReactive(input$update, {
  #   
  #   if (input$all_tissues) {
  #     return(paste0("All GTEx tissues"))
  #   } else {
  #     return(paste0("", names(tissue_GTEx_choices_alphabetical)[tissue_GTEx_choices_alphabetical == input$tissue], ""))
  #   }
  # 
  # })
  # 
  # getIntronData <- eventReactive(input$update, {
  #   
  #   withProgress(message = 'Obtaining intron data ...', value = 0.5, {
  #     
  #     # validate(
  #     #   need(input$start < input$end, "Error: the starting position must be smaller than the ending position.")
  #     # )
  #     
  #     tissue <- input$tissue
  #     if (input$all_tissues) {
  #       tissue <- NULL
  #       }
  #    
  #     get_missplicing_ratio(intron_coordinates = GRanges(seqnames = input$chr,
  #                                                        ranges = IRanges(start = input$start,
  #                                                                         end = input$end),
  #                                                        strand = input$strand),
  #                           tissue = tissue,
  #                           all_tissues = input$all_tissues)
  #   })
  #   
  # })
  # 
  # getNovelData <- eventReactive(input$update, {
  #   
  #   withProgress(message = 'Obtaining novel junction data ...', value = 0.5, {
  #     tissue <- input$tissue
  #     
  #     if (input$all_tissues) {
  #       tissue <- NULL
  #     }
  #     
  #     get_novel_data(intron_coordinates = GRanges(seqnames = input$chr,
  #                                                 ranges = IRanges(start = input$start, end = input$end),
  #                                                 strand = input$strand),
  #                    tissue = tissue,
  #                    all_tissues = input$all_tissues)
  #   })
  #   
  # })
  # 
  # intronSearchUI <- eventReactive(input$update, {
  #   
  #   table1_caption = NULL
  #   table2_caption = NULL
  #   
  #   if (input$all_tissues) {
  #     table1_caption = paste0("Mean values have been calculated across all samples from all tissues. MSR_D = Mean Mis-splicing Ratio at Donor (5'ss) positions. 
  #                             MSR_A = Mean Mis-splicing Ratio at Acceptor (3'ss) positions.")
  #     table2_caption = paste0("Mean values have been calculated across all samples from all tissues. Mean_MSR = Mean Mis-splicing Ratio.") 
  #   } else {
  #     table1_caption = paste0("Mean values have been calculated across all samples from the current tissue. MSR_D = Mean Mis-splicing Ratio at Donor (5'ss) positions. 
  #                             MSR_A = Mean Mis-splicing Ratio at Acceptor (3'ss) positions.")
  #     table2_caption = paste0("Mean values have been calculated across all samples from the current tissue. Mean_MSR = Mean Mis-splicing Ratio.")
  #   }
  #   
  # 
  #   tagList(
  #     
  #     ## Intron section
  #     h1(getPageTitle()),
  #     hr(),
  #     h3("Intron of interest:"),
  #     br(),
  #     DT::datatable(getIntronData(),
  #                   options = list(pageLength = 20,
  #                                  autoWidth = F,
  #                                  rowGroup = list(dataSrc = 1),
  #                                  rowCallback = DT::JS(
  #                                    "function(row, data) {
  #                                       var onclick_f = 'Shiny.setInputValue(\"intronID\",\"' + data[0] + '\");$(\"#modalIntron\").modal(\"show\");';
  # 
  #                                       var num = '<a id=\"goA\" role=\"button\" onclick = ' + onclick_f + ' >' + data[10] + '</a>';
  #                                       $('td:eq(10)', row).html(num);
  #                                    }"
  #                                  )
  #                                  ),
  #                   #callback = DT::JS("table.order([0, 'asc']).draw();"),
  #                   extensions = 'RowGroup',
  #                   selection = 'none',
  #                   elementId = "table1",
  #                   rownames = FALSE,
  #                   height = 200,
  #                   width = 1250,
  #                   caption = htmltools::tags$caption(
  #                     style = 'caption-side: bottom;',
  #                     table1_caption
  #                   )
  #                   ),
  #     br(),
  #     br(),
  #     hr(),
  #     ## Novel junction section
  #     h3("Novel junctions attached to the intron of interest:"),
  #     br(),
  #     DT::datatable(getNovelData(), 
  #                   options = list(pageLength = 20,
  #                                  autoWidth = F,
  #                                  order = list(1, 'desc'),
  #                                  rowGroup = list(dataSrc = 1),
  #                                  rowCallback = DT::JS(
  #                                    "function(row, data) {
  #                                       var onclick_f = 'Shiny.setInputValue(\"novelID\",\"' + data[0] + '\");$(\"#modalNovel\").modal(\"show\");';
  # 
  #                                       var num = '<a id=\"goB\" role=\"button\" onclick = ' + onclick_f + ' >' + data[11] + '</a>';
  #                                       $('td:eq(11)', row).html(num);
  #                                    }"
  #                                  )),
  #                   extensions = 'RowGroup',
  #                   selection = 'none',
  #                   elementId = "table2",
  #                   rownames = FALSE,
  #                   width = 1250,
  #                   
  #                   caption = htmltools::tags$caption(
  #                     style = 'caption-side: bottom;',
  #                     table2_caption
  #                   ))
  #     
  #     
  #   )
  #   
  # })
  # 
  # output$intronOutput = renderUI({
  #   intronSearchUI()
  # })
  # 
  # isolate(observeEvent(input$all_tissues, {
  # 
  #   if (input$all_tissues) {
  #     shinyjs::disable(id = "tissue")
  #   } else {
  #     shinyjs::enable(id = "tissue")
  #   }
  # 
  # }))
  # 
  # 
  # output$intronDetail <- renderUI({
  #   
  #   tagList(
  #     h2(paste0("Intron '", input$intronID,"' details:")),
  #     DT::datatable(get_intron_details(intron_id = input$intronID,
  #                                      tissue = input$tissue), 
  #                   options = list(pageLength = 20,
  #                                  autoWidth = T),
  #                   width = 800,
  #                   rownames = FALSE)
  #   )
  #   
  # })
  # 
  # output$novelDetail <- renderUI({
  # 
  #   tagList(
  #     h2(paste0("Novel junction '", input$novelID,"' details:")),
  #     DT::datatable(get_novel_details(novel_id = input$novelID,
  #                                     tissue = input$tissue), 
  #                   options = list(pageLength = 20,
  #                                  autoWidth = T),
  #                   width = 800,
  #                   rownames = FALSE)
  #   )
  #   
  # })
  
  
  ############################################################
  ###################### GENE SEARCH #########################
  ############################################################
  
  
  observe({
    
    cdata <- parseQueryString(session$clientData$url_search)
    
    #print(cdata[['intron']])
    
    if (!is.null(cdata[['intron']])) {
      
      #show(input$results2Div)
      #
      #hide(input$searchRow)
      hideElement(id = "searchRow")
      output$intronGeneDetail <- renderUI({
      hideElement(id = "searchRow")
        tagList(
          hideElement(id = "searchRow",asis = T),
          h2(paste0("Novel events from intron '", cdata[['coordinates']],"':")),
          DT::datatable(get_novel_data(intron_ID = cdata[['intron']],
                                       tissue = "Brain-FrontalCortex_BA9"), 
                        options = list(pageLength = 20,
                                       autoWidth = F,
                                       rowCallback = DT::JS("function(row, data) {
                                                              var onclick_f = 'Shiny.setInputValue(\"intronID\",\"' + data[0] + '\");$(\"#modalIntronGene\").modal(\"show\");';
                                                              var num = '<a id=\"goA\" role=\"button\" onclick = ' + onclick_f + ' >' + data[11] + '</a>';
                                                              $('td:eq(11)', row).html(num);
                                                            }"
                                       )),
                        width = "100%",
                        rownames = FALSE)
        )
        
      })
    }
    
  })
  
  output$modalNovelDetail <- renderUI({
    
    tagList(
      #h2(paste0("Intron '", input$intronID,"' details:")),
      DT::datatable(get_intron_details(intron_id = input$intronID,
                                       tissue = "Brain-FrontalCortex_BA9"), 
                    options = list(pageLength = 20,
                                   autoWidth = T),
                    width = 800,
                    rownames = FALSE)
    )
    
  })
  
  updateSelectizeInput(session, 'gene', choices = genes, server = TRUE)
  
  

  geneSearchUI <- eventReactive(input$geneButton, {
    
    table3_caption = paste0("Mean values have been calculated across all samples from the selected tissue. MSR_D = Mean MisSplicing Ratio at Donor (5'ss) positions. 
                            MSR_A = Mean MisSplicing Ratio at Acceptor (3'ss) positions.")
    
    tagList(
      ## Intron section
      h1(input$gene),
      h2("Annotated introns:"),
      br(),
      DT::datatable(get_gene_intron_data(input$gene, input$geneTissue), 
                    options = list(pageLength = 20,
                                   columnDefs = list(list(visible=FALSE, targets=c(2))),
                                   order = list(1, 'asc'),
                                   rowGroup = list(dataSrc = 1),
                                   autoWidth = F,
                                   rowCallback = DT::JS(
                                     "function(row, data) {
                                        
                                        var href = encodeURI('http://localhost:8787/p/7463/?intron=' + data[2] + '&coordinates=' + data[0]);
                                        var num = '<a id=\"goA\" role=\"button\" target=\"_blank\" href=' + href + '>' + data[0] + '</a>';
                                        $('td:eq(0)', row).html(num);

                                     }"
                                   )
                                   ),
                    extensions = 'RowGroup', 
                    #style = "bootstrap4",
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
  
  output$geneOutput = renderUI({
    geneSearchUI()
  })
  
  
  
  ############################################################
  ################# DISTANCES & MODULO #######################
  ############################################################
  
  # ## DISTANCES
  # plotsUI <- eventReactive(input$acceptBtnDist, {
  #   ptlist <- list(plot_distances(input$tissueDistances, input$bp),
  #                  plot_modulo(input$tissueDistances, input$bp),
  #                  plot_missplicing(input$tissueDistances))
  #   
  #   grid.arrange(grobs = ptlist,
  #                ncol = 1,
  #                nrow = length(ptlist))
  #   
  #   })
  # 
  # 
  # output$plotsOutput = renderPlot(expr = {
  #   plotsUI()
  # }, width = 600, height = 1800)
 
  
}

# Run the application 
shinyApp(ui = ui, server = server)

