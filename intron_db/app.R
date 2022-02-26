#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

### start, end only stored but calculates the third

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

library(sandwich)
library(ggstance)

# library(jtoolsscree)

# setwd("./intron_db")
source("./get_missplicing.R")
con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")


ui <- navbarPage(
  
  useShinyjs(),
  #shinyFeedback::useShinyFeedback(),
  
  title = "IDB - Intron DataBase",
  id = "intron_db",
  selected = "one",
  theme = bslib::bs_theme(bootswatch = "cosmo",
                          version = 4),
  
  ################################################
  ## PANEL 'ONE'
  ################################################
  
  tabPanel(title = "Intron Search",
           value = "one",
           
           tags$head(tags$style(HTML(".shiny-split-layout > div { overflow: visible;}"))),
           
           sidebarLayout(
             
             sidebarPanel = sidebarPanel(
               id = "sidebar_panel_tab1",
               
               div(
                 id = "geneInputPanel_tab1",
                 
                 # ## Option 1
                 # hr(),
                 # span(strong("Search:")),
                 # shiny::radioButtons(inputId = "radiobutton_introntype_tab1",
                 #                     label = "",
                 #                     choiceNames = c("Introns",
                 #                                     "Novel Junctions"),
                 #                     choiceValues = c("introns",
                 #                                      "novel"),
                 #                     selected = "introns"),
                 
                 strong(span("Using Ensembl Release 105 (Dec 2021)")),
                 
                 hr(),
                 
                 ## Option 1
                 span(strong("Search by:")),
                 shiny::radioButtons(inputId = "radiobutton_searchtype_tab1",
                                     label = "",
                                     choiceNames = c("Coordinates (hg38)",
                                                     "Gene"),
                                     choiceValues = c("radio_bycoordinates_tab1",
                                                      "radio_bygene_tab1"),
                                     selected = "radio_bygene_tab1"),
                 
                 
                 splitLayout(id="chr_strand_tab1",
                             
                   shiny::selectInput(inputId = "chr_tab1",
                                      label = "Chr",
                                      choices = NULL,
                                      multiple = F),
                   shiny::selectInput(inputId = "strand_tab1",
                                      label = "Strand",
                                      choices = c("+", "-"),
                                      selected = "+",
                                      multiple = F)
                   ),
                 
                 splitLayout(id="start_end_tab1",
                   numericInput(inputId = "start_tab1",
                                label = "Start",
                                value = 44905842),
                   numericInput(inputId = "end_tab1",
                                label = "End",
                                value = 44906601)
                   ),
                 
                
                 selectizeInput(inputId = "gene_tab1",
                                label = "Gene:",
                                choices = NULL,
                                multiple = F,
                                options = list(
                                  placeholder = "Search by gene",
                                  options = list(create = FALSE)),
                                selected = NULL),
                 
                 hr(),
                 
                 
                 
                 ## Option 4
                 
                 p(strong("Sample selection:")),
                 shiny::selectizeInput(inputId = "data_bases_tab1",
                                       label = "Project:",
                                       choices = NULL,
                                       multiple = TRUE,
                                       options = list(
                                         placeholder = "",
                                         maxItems = 15,
                                         options = list(create = FALSE)),
                                       selected = NULL),
                 shiny::selectizeInput(inputId = "clusters_tab1",
                                       label = "Samples:",
                                       choices = NULL,
                                       multiple = TRUE,
                                       options = list(
                                         placeholder = "",
                                         maxItems = 15,
                                         options = list(create = FALSE)),
                                       selected = NULL),
                 
                 hr(),
                 
                 
                 ## Option 6
                 p(strong("Support for novel annotation:")),
                 sliderInput(inputId = "threshold_tab1", 
                             label = "% individuals:",
                             min = 0, 
                             max = 95,
                             value = 1,
                             step = 10,
                             round = T),
                 
                 shiny::span("NOTE: this is the minimum % individuals in which any given novel junction attached to the intron of interest is required to be found.", style = "font-size:90%;"),
                 hr(),
                 
                 
                 
                 
                 ## Option 5
                 p(strong("Additional filters:")),
                 splitLayout(
                   checkboxInput(inputId = "clinvar_tab1", 
                                 label = "ClinVar", value = FALSE),
                   checkboxInput(inputId = "mane_tab1", 
                                 label = "MANE Select", value = FALSE)
                 ),
                 
                
                 
                 
                 ## BUTTON
                 hr(),
                 actionButton(inputId = "geneButton_tab1", label = "Accept"),
                 hidden(
                   p(id = "intronID_tab1", ""),
                   p(id = "db_tab1", ""),
                   p(id = "cluster_tab1", ""),
                   p(id = "novelID_tab1", "")
                   )
                 ),
             div(
               id = "intronPanel_tab1",
               uiOutput("intronPanelOutput_tab1")
             ),
               
               width = 3
              ),
             
              mainPanel = mainPanel(
                id = "main_tab1",
                uiOutput("geneOutput_tab1"), # %>% withSpinner(color="#0dc5c1"),
                uiOutput("intronGeneDetail_tab1"),
  
                bsModal(id = "modalIntronPlot_tab1",
                        title = "MANE transcript visualization",
                        trigger = NULL,
                        size = "large",
                        plotOutput("modalIntronPlot_tab1")),
                # bsModal(id = "modalDetailsIntron",
                #         title = NULL,
                #         trigger = NULL,
                #         size = "large",
                #         uiOutput("modalIntronDetail")),
                # bsModal(id = "modalDetailsNovel",
                #         title = NULL,
                #         trigger = NULL,
                #         size = "large",
                #         uiOutput("modalNovelDetail")),
                width = 9
             )
             )),
  
  
  ################################################
  ## PANEL 'TWO'
  ################################################
  
  tabPanel(title = "Novel Junction Search",
       
           value = "two",
           
           tags$head(tags$style(HTML(".shiny-split-layout > div { overflow: visible;}"))),
           sidebarLayout(
             
             sidebarPanel = sidebarPanel(
               id = "sidebar_panel_tab2",
               div(
                 id = "geneInputPanel_tab2",
                 
                 ## Option 1
                 # hr(),
                 # span(strong("Search:")),
                 # shiny::radioButtons(inputId = "radiobutton_introntype_tab2",
                 #                     label = "",
                 #                     choiceNames = c("Introns",
                 #                                     "Novel Junctions"),
                 #                     choiceValues = c("introns",
                 #                                      "novel"),
                 #                     selected = "introns"),
                 
                 
                 strong(span("Using Ensembl Release 105 (Dec 2021)")),
                 hr(),
                 
                 ## Option 1
                 span(strong("Search by:")),
                 shiny::radioButtons(inputId = "radiobutton_searchtype_tab2",
                                     label = "",
                                     choiceNames = c("Coordinates (hg38)",
                                                     "Gene"),
                                     choiceValues = c("radio_bycoordinates_tab2",
                                                      "radio_bygene_tab2"),
                                     selected = "radio_bygene_tab2"),
                 
                 
                 splitLayout(id="chr_strand_tab2",
                             
                             shiny::selectInput(inputId = "chr_tab2",
                                                label = "Chr",
                                                choices = NULL,
                                                multiple = F),
                             shiny::selectInput(inputId = "strand_tab2",
                                                label = "Strand",
                                                choices = c("+", "-"),
                                                selected = "+",
                                                multiple = F)
                 ),
                 
                 splitLayout(id="start_end_tab2",
                             numericInput(inputId = "start_tab2",
                                          label = "Start",
                                          value = 44906263),
                             numericInput(inputId = "end_tab2",
                                          label = "End",
                                          value = 44906601)
                 ),
                 
                 
                 selectizeInput(inputId = "gene_tab2",
                                label = "Gene:",
                                choices = NULL,
                                multiple = F,
                                options = list(
                                  placeholder = "Search by gene",
                                  maxItems = 1,
                                  options = list(create = FALSE)),
                                selected = NULL),
                 
                 hr(),
                 
                 
                 
                 ## Option 4
                 
                 p(strong("Sample selection:")),
                 shiny::selectizeInput(inputId = "data_bases_tab2",
                                       label = "Project:",
                                       choices = NULL,
                                       multiple = TRUE,
                                       options = list(
                                         placeholder = "",
                                         maxItems = 15,
                                         options = list(create = FALSE)),
                                       selected = NULL),
                 shiny::selectizeInput(inputId = "clusters_tab2",
                                       label = "Samples:",
                                       choices = NULL,
                                       multiple = TRUE,
                                       options = list(
                                         placeholder = "",
                                         maxItems = 15,
                                         options = list(create = FALSE)),
                                       selected = NULL),
                 
                 hr(),
                 
                 
                 ## Option 6
                 p(strong("Support for novel annotation:")),
                 sliderInput(inputId = "threshold_tab2", 
                             label = "% individuals:",
                             min = 1, 
                             max = 90,
                             value = 1,
                             step = 10,
                             round = T),
                 
                 shiny::span("NOTE: this is the minimum % individuals in which the novel junction is required to be found.", style = "font-size:90%;"),

                 
                 ## BUTTON
                 hr(),
                 actionButton(inputId = "geneButton_tab2", label = "Accept"),
                 hidden(
                   p(id = "intronID", ""),
                   p(id = "novelID", "")
                 )
               ),
               div(
                 id = "intronPanel_tab2",
                 uiOutput("intronPanelOutput_tab2")
               ),
               
               width = 3
             ),
             
             mainPanel = mainPanel(
               id = "main_tab2",
               uiOutput("geneOutput_tab2"),
               uiOutput("intronGeneDetail_tab2"),
               
               # bsModal(id = "modalDetailsIntron",
               #         title = NULL,
               #         trigger = NULL,
               #         size = "large",
               #         uiOutput("modalIntronDetail")),
               # bsModal(id = "modalDetailsNovel",
               #         title = NULL,
               #         trigger = NULL,
               #         size = "large",
               #         uiOutput("modalNovelDetail")),
               width = 9
             )
           )),
  
  # tabPanel(title = "Novel Annotation Search",
  #          value = "two",
  #          sidebarLayout(
  #            sidebarPanel = sidebarPanel(
  #              div(
  #                h3("Search by:"),
  #                ##
  #                shiny::radioButtons(inputId = "radiobutton_searchtype",
  #                                    label = "",
  #                                    choiceNames = c("Novel junction coordinates",
  #                                                    "Gene",
  #                                                    "Novel discovery across all databases"),
  #                                    choiceValues = c("radio_bycoordinates",
  #                                                     "radio_bygene",
  #                                                     "radio_discovery"),
  #                                    selected = character(0)),
  #                
  #                ## Option 1
  #                shiny::selectInput(inputId = "chr",
  #                                   label = "Chr",
  #                                   choices = NULL,
  #                                   multiple = F),
  #                numericInput(inputId = "start",
  #                             label = "Start",
  #                             value = 87894110),
  #                numericInput(inputId = "end",
  #                             label = "End",
  #                             value = 87925512),
  #                shiny::selectInput(inputId = "strand",
  #                                   label = "Strand",
  #                                   choices = c("+", "-"),
  #                                   selected = "+",
  #                                   multiple = F),
  #       
  #                ## Option 2
  #                selectizeInput(inputId = "gene_tab2",
  #                               label = "Gene:",
  #                               choices = NULL,
  #                               multiple = F,
  #                               options = list(
  #                                 placeholder = "Search by gene",
  #                                 maxItems = 1,
  #                                 options = list(create = FALSE)),
  #                               selected = NULL),
  #                
  #                
  #                hr(),
  #                
  #                ## slider number of individuals
  #                h4("% individuals:"),
  #                sliderInput(inputId = "threshold_tab2",
  #                            label = "",
  #                            min = 10, max = 90,
  #                            value = 10, step = 10),
  #                
  #                hr(),
  #                
  #                ## slider number of individuals
  #                h4("Data bases:"),
  #                shiny::selectizeInput(inputId = "data_bases",
  #                                       label = "",
  #                                       choices = NULL,
  #                                       multiple = TRUE,
  #                                       options = list(
  #                                         placeholder = "Select the databases to look into",
  #                                         maxItems = 15,
  #                                         options = list(create = FALSE)),
  #                                       selected = NULL),
  #                
  #                hr(),
  #                actionButton(inputId = "novelAnnotationButton", label = "Accept")
  #                
  #                ),
  #              width = 3
  #         
  #          ),
  #          mainPanel = mainPanel(
  #                uiOutput("novelAnnotationOutput") %>% withSpinner(color="#0dc5c1"),
  #                uiOutput("novel_annotation_Detail"),
  #                width = 9
  #              )
  #         )
  #     ),
  
  
  ################################################
  ## PANEL 'THREE'
  ################################################
  
  tabPanel(title = "Datasets",
           value = "three",
           fluidPage(
             fluidRow(
               # column(4, 
               #        div(
               #          h3(strong("GTEx")),
               #          shiny::tags$ul(
               #            shiny::tags$li("mRNA-Seq expression profiling"),
               #            shiny::tags$li("11 brain tissues"),
               #            shiny::tags$li("Accession number: SRP012682")
               #        )
               #        )
               # ),
               column(4, 
                      div(
                        h3(strong("GTExV8")),
                        shiny::tags$ul(
                          shiny::tags$li("2931 BRAIN samples from 13 brain regions"),
                          shiny::tags$li("Samples were collected primarily for molecular assays including WGS, WES, and RNA-Seq."))
                      )
               ),
               column(8, 
                      plotOutput("GTEx_output")
               )),
             fluidRow(
               # column(4, 
               #        div(
               #          h3(strong("GTEx")),
               #          shiny::tags$ul(
               #            shiny::tags$li("mRNA-Seq expression profiling"),
               #            shiny::tags$li("11 brain tissues"),
               #            shiny::tags$li("Accession number: SRP012682")
               #        )
               #        )
               # ),
               column(4, 
                      div(
                        h3(strong("PD/Control")),
                        shiny::tags$ul(
                          shiny::tags$li("mRNA-Seq expression profiling"),
                          shiny::tags$li("29 Parkinson's Disease and 44 neurologically normal controls from post-mortem human subjects"),
                          shiny::tags$li("All samples are from male donors"),
                          shiny::tags$li("Frontal cortex B9 samples."),
                          shiny::tags$li("Accession number:",  shiny::a("SRP058181",
                                                                        href= URLencode(URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68719"),
                                                                        target="_blank")),
                          shiny::tags$li("Illumina HiSeq 2000 (Homo sapiens)"))
                      )
               ),
               column(8, 
                      plotOutput("PD_output")
               )),
             fluidRow(
               column(4, 
                      div(
                        h3(strong("HD/Control")),
                        shiny::tags$ul(
                          shiny::tags$li("mRNA-Seq Expression profiling"),
                          shiny::tags$li("20 Huntington's Disease and 49 neurologically normal control samples from post-mortem human subjects"),
                          shiny::tags$li("Frontal cortex B9 samples."),
                          shiny::tags$li("Accession number:",
                                         shiny::a("SRP051844",
                                                  href= URLencode(URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810"),
                                                  target="_blank")),
                          shiny::tags$li("Illumina HiSeq 2000 (Homo sapiens)")
                        )
                      )
               ),
               column(8, 
                      plotOutput("HD_output")
               )
             )
             
           )
           )
           
  
      # tabPanel(title = "Plots",
      #          sidebarLayout(
      #           sidebarPanel(
      #             h3("Plot generation:"),
      #             hr(),
      #             selectizeInput(inputId = "tissueDistances", 
      #                            "Tissue:",
      #                            choices = tissue_GTEx_choices_alphabetical,
      #                            multiple = F, 
      #                            options = list(maxItems = 1), 
      #                            selected = "Brain-FrontalCortex_BA9"),
      #             #numericInput(inputId = "bp", label = "Distance (in bp):", value = 30, min = 1, max = 100000),
      #             hr(),
      #             actionButton(inputId = "acceptBtnPlots", label = "Accept"),
      #             width = 3
      #           ),
      # 
      #           # Show a plot of the generated distribution
      #           mainPanel(
      #             plotOutput("distancesOutput", height = 600) %>% withSpinner(color="#0dc5c1"),
      #             #plotOutput("moduloOutput", height = 600) %>% withSpinner(color="#0dc5c1"),
      #             plotOutput("missplicingOutput", height = 600) %>% withSpinner(color="#0dc5c1"),
      #             plotOutput("lmOutput", height = 600) %>% withSpinner(color="#0dc5c1"),
      #             width = 9
      #           )
      #         )
      # )
    
  
)



server <- function(input, output, session) {
   
  
  ##################################################
  ## TAB 'ONE'
  ##################################################
  

  ## Fill the dropdowns and hide/show inputs --------------------------------------------------------------------------------------
  
  updateSelectizeInput(session, 'chr_tab1', choices = chr_choices, server = TRUE, selected = "19")
  updateSelectizeInput(session, 'gene_tab1', choices = genes_choices, server = TRUE, selected = 5547)
  updateSelectizeInput(session, 'data_bases_tab1', choices = db_choices, server = TRUE, selected = "BRAIN")
 
  
  
  shinyjs::hideElement(id = "chr_strand_tab1")
  shinyjs::hideElement(id = "start_end_tab1")
  #shinyjs::hideElement(id = "gene_tab1")
  # shinyjs::disable(id = "geneButton")
  
  
  #shinyjs::hideElement(id = "geneInputPanel_tab1")
  #shinyjs::hideElement(id = "genePanel_tab1")
  #shinyjs::disable(id = "mane_tab1")
  
  
  ## Observers ----------------------------------------------------------------------------------------------------------------------
  
  
  # ## Intron type radio-button (i.e. by intron or by novel junction)
  # 
  # observeEvent(input$radiobutton_introntype_tab1,{
  #   freezeReactiveValue(x = input, "geneButton_tab1")
  #   
  #   if (input$radiobutton_introntype_tab1 == "novel") {
  #     
  #     updateNumericInput(inputId = "start_tab1", value = 0)
  #     updateNumericInput(inputId = "end_tab1", value = 0)
  #     
  #     shinyjs::disable(id = "clinvar_tab1")
  #     updateCheckboxInput(inputId = "clinvar_tab1", value = F)
  #     
  #   } else if (input$radiobutton_introntype_tab1 == "introns") {
  #     
  #     shinyjs::enable(id = "clinvar_tab1")
  #     updateNumericInput(inputId = "start_tab1", value = 0)
  #     updateNumericInput(inputId = "end_tab1", value = 0)
  #     
  #   }
  # })
  
  ## Search type radio-button (i.e. by coordinates or by gene name)
  observeEvent(input$radiobutton_searchtype_tab1,{
    
    #freezeReactiveValue(x = input, "geneButton_tab1")
    
    if (input$radiobutton_searchtype_tab1 == "radio_bygene_tab1") {
      
      shinyjs::showElement(id = "gene_tab1")
      
      shinyjs::hideElement(id = "chr_strand_tab1")
      shinyjs::hideElement(id = "start_end_tab1")
      
      
      shiny::updateSliderInput(session = session, inputId = "threshold_tab1", value = 1)
      shinyjs::enable(id = "data_bases_tab1")
      
      
    } else if (input$radiobutton_searchtype_tab1 == "radio_bycoordinates_tab1") {
      
      shinyjs::hideElement(id = "gene_tab1")
      
      shinyjs::showElement(id = "chr_strand_tab1")
      shinyjs::showElement(id = "start_end_tab1")
      
      
      shiny::updateSliderInput(session = session, inputId = "threshold_tab1", value = 1)
      shinyjs::enable(id = "data_bases_tab1")
      
    }
    
    shinyjs::enable(id = "geneButton_tab1")
    
  })
  
  ## Hierarchical select boxes
  observeEvent(input$data_bases_tab1, {
    
    ## In any case, we stop the reaction of the button
    #freezeReactiveValue(x = input, "geneButton_tab1")

    if (any(input$data_bases_tab1 == "all")) {
      ## Empty the select and add only the 'all' option.
      
      updateSelectizeInput(session = session, inputId = 'clusters_tab1', choices = c("All" = "all"), server = TRUE, selected = "all")
      
    } else {
      
      
      
      query = paste0("SELECT * FROM 'master'")
      con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
      df_all_projects_metadata <- dbGetQuery(conn = con, statement = query) 
      DBI::dbDisconnect(conn = con)
      
      choices <- NULL
      
      # projects <- c("GTEx", "PD", "HD")
      projects <- input$data_bases_tab1
      
      for (db in projects) { # db <- data_bases_tab1[1]
      
        data_bases <- df_all_projects_metadata %>%
          filter(SRA_project == db)
        
        clusters <- data_bases$cluster %>% unique() %>% as.list()
        names(clusters) <- data_bases$diagnosis %>% unique() %>% as.list()
        choices <- c(choices, clusters)
        
      }
      options_selected <- input$clusters_tab1
      
      if (is.null(options_selected)) {
        options_selected <- choices[1]
      }
      updateSelectizeInput(session = session, inputId = 'clusters_tab1', choices = choices, 
                           server = TRUE, selected = options_selected)
      
      
    }
  })
  
  ## Other controls
  # observeEvent(c(input$clusters_tab1, input$mane_tab1, input$clinvar_tab1, input$gene_tab1),{
  #   freezeReactiveValue(x = input, name = "geneButton_tab1")
  # })
  
  
  
  ## Get the list of novel junctions attached to an annotated intron - shown in a different tab ------------------------------------
  
  observe({
    
    cdata <- parseQueryString(session$clientData$url_search)
    
    
    
    if (!is.null(cdata[['intron']])) {

      print(cdata[['intron']])
      removeCssClass(id = "main_tab1", "col-sm-9")
      addCssClass(id = "main_tab1", "col-sm-12")
      
      shinyjs::hide(id = "sidebar_panel_tab1")
      shinyjs::hide(id = "intronPanelOutput_tab1")
      
      updateNavbarPage(session = session, inputId = "intron_db", selected = "one")

      
      output$intronGeneDetail_tab1 <- renderUI({
      
        tagList(
        
          h2(paste0("Novel events attached to intron 'ID=", cdata[['intron']],"':")),
          
          DT::renderDataTable(get_novel_data(intron = URLdecode(URL = cdata[['intron']]),
                                             db = URLdecode(URL = cdata[['db']]),
                                             sample_group = URLdecode(URL = cdata[['cluster']])), 
                              options = list(pageLength = 20,
                                             order = list(9, 'desc'),
                                             columnDefs = list(list(visible=FALSE, targets=c(2))),
                                             autoWidth = F),
                              width = "100%",
                              rownames = FALSE)
        )
        
      })
      
    } else if (!is.null(cdata[['id']])) {
      
      ## Hide elements from tab1
      
      
      removeCssClass(id = "main_tab2", "col-sm-9")
      addCssClass(id = "main_tab2", "col-sm-12")
      
      shinyjs::hide(id = "sidebar_panel_tab2")
      shinyjs::hide(id = "intronPanelOutput_tab2")
      
      updateNavbarPage(session = session, inputId = "intron_db", selected = "two")
      
      
      output$intronGeneDetail_tab2 <- renderUI({
        
        tagList(
          
          h2(paste0("Details of the novel junction 'ID=", cdata[['id']],"' across all projects from the IDB:")),
          
          DT::renderDataTable(search_novel_junction(id = cdata[['id']]), 
                              options = list(pageLength = 20,
                                             order = list(8, 'desc')
                              ),
                              #extensions = 'RowGroup', 
                              width = "100%",
                              rownames = FALSE)
        )
        
      })
    } 
  })
  

  ## Modal Popups --------------------------------------------------------------------
  output$modalIntronPlot_tab1 <- renderPlot({

    plot_transcript_from_intron(intron_id = str_replace_all(string = input$intronID_tab1, pattern = "%20", replacement = " "),
                                db = str_replace_all(string = input$db_tab1, pattern = "%20", replacement = " "),
                                clust = str_replace_all(string = input$cluster_tab1, pattern = "%20", replacement = " ") )


  }, width = "auto", height = "auto")
  # ## Popup modal window with the detail of the individuals presenting an annotated intron
  # 
  # output$modalIntronDetail <- renderUI({
  #   
  #   #cdata <- parseQueryString(session$clientData$url_search)
  #   
  #   tagList(
  #     
  #     DT::renderDataTable(get_intron_details(intron_id = input$intronID,
  #                                            tissue = input$geneTissue), 
  #                         options = list(pageLength = 20,
  #                                        autoWidth = T,
  #                                        order = list(2, 'desc')),
  #                         width = "100%",
  #                         rownames = FALSE)
  #   )
  #   
  # })
  # 
  # ## Popup modal window with the detail of the individuals presenting a novel junction
  # 
  # output$modalNovelDetail <- renderUI({
  #   
  #   cdata <- parseQueryString(session$clientData$url_search)
  #   
  #   tagList(
  #     
  #     DT::renderDataTable(get_novel_details(novel_id = input$novelID,
  #                                     tissue = cdata[['tissue']]), 
  #                         options = list(pageLength = 20,
  #                                        autoWidth = T,
  #                                        order = list(2, 'desc')),
  #                         width = "100%",
  #                         rownames = FALSE)
  #   )
  #   
  # })
  
  
  
  
  ## Get all annotated introns from the selected gene -----------------------------------------------------------------------------

  # observeEvent(input$button1, {
  #   output$button2 <- renderUI({
  #     actionButton("button2", label = "Press Button 2")
  #   })
  # })
  observeEvent(input$geneButton_tab1,  {
    
    output$geneOutput_tab1 = renderUI({
      
    # req(!is.null(input$clusters_tab1),
    #     !is.null(input$data_bases_tab1), 
    #     cancelOutput = T)

    title <- "Annotated introns - Details"

    IDB_data <- search_intron(type = "introns",
                              chr = input$chr_tab1,
                              start = input$start_tab1,
                              end = input$end_tab1,
                              strand = input$strand_tab1,
                              gene = input$gene_tab1,
                              threshold = input$threshold_tab1,
                              search_type = input$radiobutton_searchtype_tab1,
                              data_bases = input$data_bases_tab1,
                              clusters = input$clusters_tab1,
                              mane = input$mane_tab1,
                              clinvar = input$clinvar_tab1)
    
   
    
    if (any(names(IDB_data) == "Message")) {
      tagList(
        h1(title),
        DT::renderDataTable(IDB_data)
      )
      
    } else {
      
      table3_caption = paste0("MSR_D = Mis-splicing Ratio at Donor (5'ss) position. MSR_A = Mis-splicing Ratio at Acceptor (3'ss) position.\nReference annotation used: Ensembl v104 (March 2021)")
      
      if (input$radiobutton_searchtype_tab1 == "radio_bycoordinates_tab1") {
        
        info <- paste0("Reference intron '", IDB_data$ID %>% unique(), "' from '", IDB_data$Gene %>% unique, "' gene.")
        
        # info <- p(strong("Coordinates: "),paste0("'", IDB_data$Coordinates %>% unique(), "'."), 
        #           br(),
        #           strong("Width: "),paste0("'", IDB_data$Width %>% unique(), "'."),
        #           br(),
        #           strong("Ss5score: "),paste0("'", IDB_data$Ss5score %>% unique(), "'."),
        #           br(),
        #           strong("Ss3score: "),paste0("'", IDB_data$Ss3score %>% unique(), "'."),
        #           br(),
        #           strong("ClinVar: "),paste0("'", IDB_data$ClinVar %>% unique(), "'."),
        #           br(),
        #           strong("Gene: "),paste0("'", IDB_data$Gene %>% unique(), "'."))
        # 
        # IDB_data <- IDB_data %>%
        #   select(-Width, -Ss5score, -Ss3score, -ClinVar, -Gene)
        
      } else {
        info <- paste0("All ", str_to_lower(title), " from ", input$gene_tab1, " gene")
      }
      
     
      IDB_data <- IDB_data %>%
        mutate("More" = "See") %>%
        relocate("More", .before = "sample")
      
      
      # print(IDB_data)
      
      tagList(
        
        h1(title),
        br(),
        
        div(info),
        br(),
        
        DT::renderDataTable(IDB_data,
                            options = list(pageLength = 20,
                                           columnDefs = list(list(visible=FALSE, targets=c(14))),
                                           order = list(0, 'asc'),
                                           rowGroup = list(dataSrc = 0),
                                           autoWidth = F,
                                           rowCallback = DT::JS("function(row, data) {
                                           if (data[1] != 'never') {
                                                // It's the intron view
                                                var href = encodeURI('https://soniagarciaruiz.shinyapps.io/intron_db/?intron=' + data[0] + '&db=' + data[12] + '&cluster=' + data[14]);
                                                var num = '<a id=\"goA\" role=\"button\" target=\"_blank\" href=' + href + '> Open missplicing events </a>';
                                                //var num = 'Check missplicing events';
                                                


                                            //var href = encodeURI('https://soniagarciaruiz.shinyapps.io/intron_db/?intron=' + data[2] + '&coordinates=' + data[0] + '&gene=' + data[11] + '&type=' + data[1] + '&clinvar=' + data[10] + '&length=' + data[3] + '&db=' + data[12] + '&cluster=' + data[13]);
                                            //var num = '<a id=\"goA\" role=\"button\" target=\"_blank\" onclick=' + href + ' class = 'button'>' + data[0] + '</a>';
                                            //$('td:eq(0)', row).html(num);


                                            var onclick_f = 'Shiny.setInputValue(\"intronID_tab1\",\"' + encodeURI(data[0]) + '\");Shiny.setInputValue(\"db_tab1\",\"' + encodeURI(data[12]) + '\");Shiny.setInputValue(\"cluster_tab1\",\"' + encodeURI(data[14]) + '\");$(\"#modalIntronPlot_tab1\").modal(\"show\");';
                                            //var onclick_f = 'Shiny.setInputValue(\"intronID_tab1\",\"ppp\");Shiny.setInputValue(\"db_tab1\",\"aaa\");Shiny.setInputValue(\"cluster_tab1\",\"aaa\");$(\"#modalIntronPlot_tab1\").modal(\"show\");';
                                            console.log(onclick_f)
                                            num = num + '<br/><a id=\"goA\" role=\"button\" onclick = ' + onclick_f + ' ><button>Visualize transcript</button></a>';
                                            $('td:eq(13)', row).html(num);
                                            
                                          } else {
                                            var num = 'N/A';
                                            $('td:eq(13)', row).html(num);
                                          }
                                       }"
                                     )
                                     ),
                            extensions = 'RowGroup',
                            width = "100%",
                            #selection = 'none',
                            rownames = F,
                            caption = htmltools::tags$caption(
                              style = 'caption-side: bottom;',
                              table3_caption)
                      ),
        br(), 
        br(),
        br(),
        br()
      
      ) 
    }
    })
  })
  
  
  ## What happens after the user clicks the 'Accept' button -----------------------------------------------------------------------
  
  # output$geneOutput_tab1 = renderUI({
  #   #showElement(id = "geneInputPanel_tab1")
  #   #showElement(id = "genePanel_tab1")
  #   
  #   geneSearchUI()
  #   
  # })
  
  
  # output$db_warning <- renderText({
  #   db_validation()
  # })
  # 
  # output$cluster_warning <- renderText({
  #   cluster_validation()
  # })
  
  ##################################################
  ## TAB 'THREE' - PLOTS
  ##################################################
  
  output$GTEx_output = renderPlot({
    plot_metadata("BRAIN")
    
  })
  output$PD_output = renderPlot({
    plot_metadata("SRP058181")

  })
  output$HD_output = renderPlot({
    plot_metadata("SRP051844")
    
  })
  
  
  ##################################################
  ## TAB 'TWO'
  ##################################################
  
  
  ## Fill the dropdowns and hide/show inputs --------------------------------------------------------------------------------------
  
  updateSelectizeInput(session, 'chr_tab2', choices = chr_choices, server = TRUE, selected = "19")
  updateSelectizeInput(session, 'gene_tab2', choices = genes_choices, server = TRUE, selected = "5547")
  updateSelectizeInput(session, 'data_bases_tab2', choices = db_choices, server = TRUE, selected = "BRAIN")
  
  shinyjs::hideElement(id = "chr_strand_tab1")
  shinyjs::hideElement(id = "start_end_tab1")

  
  

  
  ## Observers ----------------------------------------------------------------------------------------------------------------------
  
  
  # ## Intron type radio-button (i.e. by intron or by novel junction)
  # 
  # observeEvent(input$radiobutton_introntype_tab1,{
  #   freezeReactiveValue(x = input, "geneButton")
  #   
  #   if (input$radiobutton_introntype_tab1 == "novel") {
  #     
  #     updateNumericInput(inputId = "start_tab1", value = 0)
  #     updateNumericInput(inputId = "end_tab1", value = 0)
  #     
  #     shinyjs::disable(id = "clinvar")
  #     updateCheckboxInput(inputId = "clinvar", value = F)
  #     
  #   } else if (input$radiobutton_introntype_tab1 == "introns") {
  #     
  #     shinyjs::enable(id = "clinvar")
  #     updateNumericInput(inputId = "start_tab1", value = 0)
  #     updateNumericInput(inputId = "end_tab1", value = 0)
  #     
  #   }
  # })
  
  ## Search type radio-button (i.e. by coordinates or by gene name)
  
  observeEvent(input$radiobutton_searchtype_tab2,{
    
    #freezeReactiveValue(x = input, "geneButton_tab2")
    
    if (input$radiobutton_searchtype_tab2 == "radio_bygene_tab2") {
      
      shinyjs::showElement(id = "gene_tab2")
      
      shinyjs::hideElement(id = "chr_strand_tab2")
      shinyjs::hideElement(id = "start_end_tab2")
      
      
      shiny::updateSliderInput(session = session, inputId = "threshold_tab2", value = 1)
      shinyjs::enable(id = "data_bases_tab2")
      
      
    } else if (input$radiobutton_searchtype_tab2 == "radio_bycoordinates_tab2") {
      
      shinyjs::hideElement(id = "gene_tab2")
      
      shinyjs::showElement(id = "chr_strand_tab2")
      shinyjs::showElement(id = "start_end_tab2")
      
      
      shiny::updateSliderInput(session = session, inputId = "threshold_tab2", value = 1)
      shinyjs::enable(id = "data_bases_tab2")
      
    }
    
    shinyjs::enable(id = "geneButton_tab2")
    
  })
  
  
  
  ## Hierarchical select boxes
  
  observeEvent(input$data_bases_tab2, {
    
    if (any(input$data_bases_tab2 == "all")) {
      ## Empty the select and add only the 'all' option.
      
      updateSelectizeInput(session = session, inputId = 'clusters_tab2', choices = c("All" = "all"), 
                           server = TRUE, selected = "all")
      
    } else {
      
      # print(input$data_bases_tab2)
      
      query = paste0("SELECT * FROM 'master'")
      con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
      df_all_projects_metadata <- dbGetQuery(conn = con, statement = query) 
      DBI::dbDisconnect(conn = con)
      
      choices <- NULL
      
      # projects <- c("GTEx", "PD", "HD")
      projects <- input$data_bases_tab2
      
      for (db in projects) { # db <- data_bases_tab1[1]
        
        data_bases <- df_all_projects_metadata %>%
          filter(SRA_project == db)
        
        clusters <- data_bases$cluster %>% unique() %>% as.list()
        names(clusters) <- data_bases$diagnosis %>% unique() %>% as.list()
        choices <- c(choices, clusters)
        
      }
      options_selected <- input$clusters_tab2
      
      if (is.null(options_selected)) {
        options_selected <- choices[1]
      }
      updateSelectizeInput(session = session, inputId = 'clusters_tab2', choices = choices, server = TRUE, selected = options_selected)
      
  
      
      
    }
    
    ## In any case, we stop the reaction of the button
    #freezeReactiveValue(x = input, "geneButton_tab2")
    
  })
  
  # observeEvent(input$clusters_tab2,{
  #   freezeReactiveValue(x = input, "geneButton_tab2")
  # })
  
  
  
  ## Get the list of novel junctions attached to an annotated intron - shown in a different tab ------------------------------------
  
  # observe({
  #   
  #   cdata <- parseQueryString(session$clientData$url_search)
  #   
  #   # print(paste0("Intron '", input$intronID,"' details:"))
  #   # print(cdata[['intron']])
  #   
  #   
  #   
  # })
  # 
  ## Required values ---------------------------------------------------------------------------------------------------------------
  # cluster_validation <- eventReactive(input$geneButton, {
  # 
  #   cluster_null <- is.null(input$clusters_tab1)
  #   shinyFeedback::feedbackWarning(inputId = "clusters_tab1", 
  #                                  show = cluster_null, 
  #                                  text = "Please select a cluster",
  #                                  color = "red")
  # }, ignoreInit = T)
  # db_validation <- eventReactive(input$geneButton, {
  #   
  #   db_null <- is.null(input$data_bases_tab1)
  #   shinyFeedback::feedbackWarning(inputId = "data_bases_tab1", 
  #                                  show = db_null, 
  #                                  text = "Please select a DB",
  #                                  color = "red")
  # 
  # }, ignoreInit = T)
  
  
  
  
  # ## Popup modal window with the detail of the individuals presenting an annotated intron ------------------------------------------
  # 
  # output$modalIntronDetail <- renderUI({
  #   
  #   #cdata <- parseQueryString(session$clientData$url_search)
  #   
  #   tagList(
  #     
  #     DT::renderDataTable(get_intron_details(intron_id = input$intronID,
  #                                            tissue = input$geneTissue), 
  #                         options = list(pageLength = 20,
  #                                        autoWidth = T,
  #                                        order = list(2, 'desc')),
  #                         width = "100%",
  #                         rownames = FALSE)
  #   )
  #   
  # })
  # 
  # ## Popup modal window with the detail of the individuals presenting a novel junction ---------------------------------------------
  # 
  # output$modalNovelDetail <- renderUI({
  #   
  #   cdata <- parseQueryString(session$clientData$url_search)
  #   
  #   tagList(
  #     
  #     DT::renderDataTable(get_novel_details(novel_id = input$novelID,
  #                                           tissue = cdata[['tissue']]), 
  #                         options = list(pageLength = 20,
  #                                        autoWidth = T,
  #                                        order = list(2, 'desc')),
  #                         width = "100%",
  #                         rownames = FALSE)
  #   )
  #   
  # })
  
  
  
  
  ## Get all annotated introns from the selected gene -----------------------------------------------------------------------------
  
  observeEvent(input$geneButton_tab2,  {
    
    output$geneOutput_tab2 = renderUI({
    
    req(!is.null(input$clusters_tab2),
        !is.null(input$data_bases_tab2), 
        cancelOutput = T)
    
    
    title <- "Novel junctions - Details"
    
    IDB_data <- search_intron(type = "novel",
                              chr = input$chr_tab2,
                              start = input$start_tab2,
                              end = input$end_tab2,
                              strand = input$strand_tab2,
                              gene = input$gene_tab2,
                              threshold = input$threshold_tab2,
                              search_type = input$radiobutton_searchtype_tab2,
                              data_bases = input$data_bases_tab2,
                              clusters = input$clusters_tab2,
                              mane = F, 
                              clinvar = F)
    
    
    if (any(names(IDB_data) == "Message")) {
      
      tagList(
        h1(title),
        DT::renderDataTable(IDB_data)
      )
      
    } else {
      
      table3_caption = paste0("MSR_D = Mis-splicing Ratio at Donor (5'ss) position. MSR_A = Mis-splicing Ratio at Acceptor (3'ss) position.\nReference annotation used: Ensembl v104 (March 2021)")
      
      # if (input$radiobutton_searchtype_tab2 == "radio_bycoordinates_tab2") {
      #   
      #   info <- p(strong("Coordinates: "),paste0("'", IDB_data$Coordinates %>% unique(), "'."), 
      #             br(),
      #             strong("Width: "),paste0("'", IDB_data$Width %>% unique(), "'."),
      #             br(),
      #             strong("Ss5score: "),paste0("'", IDB_data$Ss5score %>% unique(), "'."),
      #             br(),
      #             strong("Ss3score: "),paste0("'", IDB_data$Ss3score %>% unique(), "'."),
      #             br(),
      #             strong("ClinVar: "),paste0("'", IDB_data$ClinVar %>% unique(), "'."),
      #             br(),
      #             strong("Gene: "),paste0("'", IDB_data$Gene %>% unique(), "'."))
      #   
      #   IDB_data <- IDB_data %>%
      #     select(-Width, -Ss5score, -Ss3score, -ClinVar, -Gene)
      #   
      # } else {
        info <- paste0("All ", str_to_lower(title), " from ", input$gene_tab2, " gene")
        
        
      #}
      
      IDB_data <- IDB_data %>%
        select(-ClinVar) %>%
        mutate("More" = "See") %>%
        relocate("More", .before = "sample")
      
      # print(IDB_data)
      
      tagList(
        
        h1(title),
        br(),
        
        div(info),
        br(),
        
        DT::renderDataTable(IDB_data,
                            options = list(pageLength = 20,
                                           columnDefs = list(list(visible=FALSE, 
                                                                  targets=c(14))),
                                           
                                           order = list(0, 'asc'),
                                           rowGroup = list(dataSrc = 0),
                                           autoWidth = F,
                                           rowCallback = DT::JS(
                                             "function(row, data) {
                                             if (data[1].includes('novel')) {
                                             // It's the novel junction view
                                             
                                             var href = encodeURI('https://soniagarciaruiz.shinyapps.io/intron_db/?id=' + data[0]);
                                             var num = '<a id=\"goA\" role=\"button\" target=\"_blank\" href=' + href + '> Check across the IDB </a>';
                                             //var num = 'Check across the IDB';
                                             $('td:eq(13)', row).html(num);
                                             
                                             
                                             }}"
                                             )),
                            extensions = 'RowGroup', 
                            width = "100%",
                            #selection = 'none',
                            rownames = F,
                            caption = htmltools::tags$caption(
                              style = 'caption-side: bottom;',
                              table3_caption)
                                           ),
        br(), 
        br(),
        br(),
        br()
      )
      }
    
    })
  })
  
  
  ## What happens after the user clicks the 'Accept' button -----------------------------------------------------------------------
  
  # output$geneOutput_tab2 = renderUI({
  #   # showElement(id = "geneInputPanel")
  #   # hideElement(id = "genePanel")
  #   if(!is.null(input$clusters_tab2) &&
  #      !is.null(input$data_bases_tab2) &&
  #      !is.null(input$radiobutton_introntype_tab2)) 
  #     geneSearchUI()
  #   
  # })
  
  
  # output$db_warning <- renderText({
  #   db_validation()
  # })
  # 
  # output$cluster_warning <- renderText({
  #   cluster_validation()
  # })
  
  # 
  # ## Fill the dropdown with the list of genes -------------------------------------------
  # updateSelectizeInput(session, 'gene_tab2', choices = genes, server = TRUE, selected = "PTEN")
  # 
  # ## Fill the dropdown with the list of chr -------------------------------------------
  # updateSelectizeInput(session, 'chr', choices = chr_choices, server = TRUE, selected = "10")
  # 
  # ## Fill the dropdown with the list of data bases -------------------------------------------
  # updateSelectizeInput(session, 'data_bases', choices = db_choices, server = TRUE, selected = "GTEx")
  # 
  # shinyjs::hideElement(id = "chr")
  # shinyjs::hideElement(id = "start")
  # shinyjs::hideElement(id = "end")
  # shinyjs::hideElement(id = "strand")
  # shinyjs::hideElement(id = "gene_tab2")
  # 
  # shinyjs::disable(id = "novelAnnotationButton")
  # 
  # 
  # observeEvent(input$radiobutton_searchtype,{
  #   
  #   freezeReactiveValue(x = input, "novelAnnotationButton")
  #   
  #   if (input$radiobutton_searchtype == "radio_bygene") {
  #     
  #     shinyjs::showElement(id = "gene_tab2")
  #     
  #     shinyjs::hideElement(id = "chr")
  #     shinyjs::hideElement(id = "start")
  #     shinyjs::hideElement(id = "end")
  #     shinyjs::hideElement(id = "strand")
  #     
  #     updateSelectizeInput(session, 'data_bases', choices = db_choices, server = TRUE, selected = "GTEx")
  #     shiny::updateSliderInput(session = session, inputId = "threshold_tab2", value = 10)
  #     shinyjs::enable(id = "data_bases")
  # 
  #     
  #   } else if (input$radiobutton_searchtype == "radio_bycoordinates") {
  #     
  #     shinyjs::hideElement(id = "gene_tab2")
  #     
  #     shinyjs::showElement(id = "chr")
  #     shinyjs::showElement(id = "start")
  #     shinyjs::showElement(id = "end")
  #     shinyjs::showElement(id = "strand")
  #     
  #     updateSelectizeInput(session, 'data_bases', choices = db_choices, server = TRUE, selected = "GTEx")
  #     shiny::updateSliderInput(session = session, inputId = "threshold_tab2", value = 10)
  #     shinyjs::enable(id = "data_bases")
  # 
  #   } else {
  #     
  #     shinyjs::hideElement(id = "chr")
  #     shinyjs::hideElement(id = "start")
  #     shinyjs::hideElement(id = "end")
  #     shinyjs::hideElement(id = "strand")
  #     shinyjs::hideElement(id = "gene_tab2")
  #     
  #     shiny::updateSelectizeInput(session = session, inputId = "data_bases", 
  #                                 choices = c("All"), server = TRUE, selected = "All")
  #     shiny::updateSliderInput(session = session, inputId = "threshold_tab2", value = 90)
  #     shinyjs::disable(id = "data_bases")
  #     
  #   }
  #   
  #   shinyjs::enable(id = "novelAnnotationButton")
  #   
  # }, ignoreInit = TRUE)
  # 
  # 
  # 
  # 
  # 
  # ## Main function ---------------------------------------------------------------------
  # 
  # novelAnnotationSearchUI <- eventReactive(input$novelAnnotationButton, {
  #   
  #   
  #   
  #  
  #   tagList(
  #     
  #     h1("Potential novel annotation:"),
  #     br(),
  #     DT::renderDataTable(get_novel_annotation_data(chr = input$chr,
  #                                                   start = input$start,
  #                                                   end = input$end,
  #                                                   strand = input$strand,
  #                                                   gene = input$gene_tab2,
  #                                                   threshold = input$threshold_tab2,
  #                                                   search_type = input$radiobutton_searchtype,
  #                                                   data_bases = input$data_bases),
  #                         rownames = F),
  #     
  #     br(), 
  #     br(),
  #     br(),
  #     br()
  #     
  #   )
  #   
  # })
  # 
  # ## Accept button
  # 
  # output$novelAnnotationOutput = renderUI({
  #   novelAnnotationSearchUI()
  # })
  
  

  ##################################################
  ## TAB 'THREE' - PLOTS
  ##################################################
  
  
  
  # ## DISTANCES
  # plotDistances <- eventReactive(input$acceptBtnPlots, {
  #   plot_distances(input$tissueDistances)
  # })
  # plotModulo <- eventReactive(input$acceptBtnPlots, {
  #   plot_modulo(input$tissueDistances)
  # })
  # plotMissplicing <- eventReactive(input$acceptBtnPlots, {
  #                  plot_missplicing(input$tissueDistances)
  # })
  # plotLm <- eventReactive(input$acceptBtnPlots, {
  #   plot_lm(input$tissueDistances)
  # })
  # 
  # 
  # output$distancesOutput = renderPlot(expr = {
  #   plotDistances()
  # }, width = 600, height = 600)
  # # output$moduloOutput = renderPlot(expr = {
  # #   plotModulo()
  # # }, width = 600, height = 600)
  # output$missplicingOutput = renderPlot(expr = {
  #   plotMissplicing()
  # }, width = 600, height = 600)
  # 
  # output$lmOutput = renderPlot(expr = {
  #   plotLm()
  # }, width = 600, height = 600)
  
}

# Run the application 
shinyApp(ui = ui, server = server)

