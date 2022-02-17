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
  selected = "one",
  theme = bslib::bs_theme(bootswatch = "cosmo"),
  
  tabPanel(title = "Intron Search",
           value = "one",
           tags$head(tags$style(HTML(".shiny-split-layout > div { overflow: visible;}"))),
           sidebarLayout(
             
             sidebarPanel = sidebarPanel(
               id = "sidebar_panel",
               div(
                 id = "geneInputPanel",
                 
                 ## Option 1
                 hr(),
                 span(strong("Search:")),
                 shiny::radioButtons(inputId = "radiobutton_introntype_tab1",
                                     label = "",
                                     choiceNames = c("Introns",
                                                     "Novel Junctions"),
                                     choiceValues = c("introns",
                                                      "novel"),
                                     selected = "introns"),
                 
                 
                 
                 hr(),
                 
                 ## Option 1
                 span(strong("Search by:")),
                 shiny::radioButtons(inputId = "radiobutton_searchtype_tab1",
                                     label = "",
                                     choiceNames = c("Coordinates",
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
                                value = 87894110),
                   numericInput(inputId = "end_tab1",
                                label = "End",
                                value = 87925512)
                   ),
                 
                
                 selectizeInput(inputId = "gene_tab1",
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
                 sliderInput(inputId = "threshold", 
                             label = "% individuals:",
                             min = 1, 
                             max = 90,
                             value = 1,
                             step = 10,
                             round = T),
                 
                 shiny::span("NOTE: this is the minimum % individuals in which the novel junction has to be found.", style = "font-size:85%;"),
                 hr(),
                 
                 ## Option 5
                 p(strong("Additional filters:")),
                 splitLayout(
                   checkboxInput(inputId = "clinvar", 
                                 label = "ClinVar", value = FALSE),
                   checkboxInput(inputId = "mane", 
                                 label = "MANE Select", value = T)
                 ),
                 
                
                 
                 
                 ## BUTTON
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
                id = "main",
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
  tabPanel(title = "Datasets",
           value = "three",
           fluidPage(
             fluidRow(
               column(4, 
                      div(
                        h3(strong("GTEx")),
                        shiny::tags$ul(
                          shiny::tags$li("mRNA-Seq expression profiling"),
                          shiny::tags$li("11 brain tissues"),
                          shiny::tags$li("Accession number: SRP012682")
                      )
                      )
               ),
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
  
  updateSelectizeInput(session, 'gene_tab1', choices = genes, server = TRUE, selected = "APOE")
  updateSelectizeInput(session, 'data_bases_tab1', choices = db_choices, server = TRUE, selected = "GTEx")
  updateSelectizeInput(session, 'chr_tab1', choices = chr_choices, server = TRUE, selected = "10")
  
  
  shinyjs::hideElement(id = "chr_strand_tab1")
  shinyjs::hideElement(id = "start_end_tab1")
  #shinyjs::hideElement(id = "gene_tab1")
  # shinyjs::disable(id = "geneButton")
  
  
  shinyjs::hideElement(id = "geneInputPanel")
  shinyjs::hideElement(id = "genePanel")
  shinyjs::disable(id = "mane")
  
  
  ## Hierarchical select boxes ----------------------------------------------------------------------------------------------------
  
  observeEvent(input$data_bases_tab1, {

    db_choices_full <- readRDS(file = "./dependencies/db_choices_simplified.rds")
    # db_choices_full <- db_choices_full %>%
    #   mutate(tidy = ifelse(db_choices_full$data_base == "GTEx", paste0(tidy, " (GTEx)"), tidy)) %>%
    #   mutate(tidy = ifelse(db_choices_full$data_base == "PD", paste0(tidy, " (PD/Control)"), tidy)) %>%
    #   mutate(tidy = ifelse(db_choices_full$data_base == "HD", paste0(tidy, " (HD/Control)"), tidy))
      
    # gtex <- readRDS(file = "./dependencies/gtex_tissues_tidy.rds")
    # HD <- readRDS(file = "./dependencies/HDControl_clusters_tidy.rds")
    # PD <- readRDS(file = "./dependencies/PDControl_clusters_tidy.rds")
    # 
    # db_choices_full <- db_choices_full %>%
    #   mutate(cluster = c(gtex$sample[7:17],
    #                      PD$sample,
    #                      HD$sample)) %>%
    #   dplyr::rename(tidy = clusters)
    # saveRDS(object = db_choices_full, file = "./dependencies/db_choices_simplified.rds")
    
    choices <- NULL
    
    # projects <- c("GTEx", "PD", "HD")
    projects <- input$data_bases_tab1
    for (db in projects) { # db <- data_bases_tab1[1]
      data_bases <- db_choices_full %>%
        filter(data_base == db)
      
      clusters <- data_bases$cluster %>% as.list()
      names(clusters) <- data_bases$tidy %>% as.character()
      choices <- c(choices, clusters)
    }
    options_selected <- input$clusters_tab1
    
    if (is.null(options_selected)) {
      options_selected <- choices[1]
    }
    updateSelectizeInput(session = session, inputId = 'clusters_tab1', choices = choices, 
                         server = TRUE, selected = options_selected)
    
    freezeReactiveValue(x = input, "geneButton")
  })
  
  observeEvent(input$clusters_tab1,{
    freezeReactiveValue(x = input, "geneButton")
  })
  
  
  
  ## Observers ----------------------------------------------------------------------------------------------------------------------
  
  
  observeEvent(input$radiobutton_introntype_tab1,{
    freezeReactiveValue(x = input, "geneButton")
    
    if (input$radiobutton_introntype_tab1 == "novel") {
      
      updateNumericInput(inputId = "start_tab1", value = 0)
      updateNumericInput(inputId = "end_tab1", value = 0)
      
      shinyjs::disable(id = "clinvar")
      updateCheckboxInput(inputId = "clinvar", value = F)
      
    } else if (input$radiobutton_introntype_tab1 == "introns") {
      
      shinyjs::enable(id = "clinvar")
      updateNumericInput(inputId = "start_tab1", value = 0)
      updateNumericInput(inputId = "end_tab1", value = 0)

    }
  })
  
  
  
  observeEvent(input$radiobutton_searchtype_tab1,{
    
    freezeReactiveValue(x = input, "geneButton")
    
    if (input$radiobutton_searchtype_tab1 == "radio_bygene_tab1") {
      
      shinyjs::showElement(id = "gene_tab1")
      
      shinyjs::hideElement(id = "chr_strand_tab1")
      shinyjs::hideElement(id = "start_end_tab1")
      
      
      shiny::updateSliderInput(session = session, inputId = "threshold", value = 1)
      shinyjs::enable(id = "data_bases_tab1")
   
      
    } else if (input$radiobutton_searchtype_tab1 == "radio_bycoordinates_tab1") {
      
      shinyjs::hideElement(id = "gene_tab1")
      
      shinyjs::showElement(id = "chr_strand_tab1")
      shinyjs::showElement(id = "start_end_tab1")
      
      
      shiny::updateSliderInput(session = session, inputId = "threshold", value = 1)
      shinyjs::enable(id = "data_bases_tab1")
      
    }
    
    shinyjs::enable(id = "geneButton")
    
  }, ignoreInit = TRUE)
  
  
  # observeEvent(input$enovel,{
  #   if (input$enovel == T) {
  #     shinyjs::enable(id = "threshold")
  #   } else {
  #     shinyjs::disable(id = "threshold")
  #   }
  #   
  # }, ignoreInit = TRUE)
  
  
  ## Get the list of novel junctions attached to an annotated intron - shown in a different tab ------------------------------------
  
  observe({
    
    cdata <- parseQueryString(session$clientData$url_search)
    
    # print(paste0("Intron '", input$intronID,"' details:"))
    # print(cdata[['intron']])
    
    if (!is.null(cdata[['id']])) {
      
      hideElement(id = "geneInputPanel")
      showElement(id = "genePanel")
      
      removeCssClass(id = "main", "col-sm-9")
      addCssClass(id = "main", "col-sm-12")
      
      shinyjs::hide(id = "sidebar_panel")
      shinyjs::hide(id = "intronPanelOutput")
      
      
      output$intronGeneDetail <- renderUI({
        
        tagList(
          
          h2(paste0("Details of the novel junction 'ID=", cdata[['id']],"' across all projects from the IDB:")),
          
          DT::renderDataTable(get_novel_annotation_data(id = cdata[['id']]), 
                              options = list(pageLength = 20,
                                             order = list(8, 'desc')
                                             # columnDefs = list(list(visible=FALSE, targets=c(2))),
                                             # autoWidth = F,
                                             # rowGroup = list(dataSrc = 1), rowCallback = DT::JS("function(row, data) {
                                             #                  var onclick_f = 'Shiny.setInputValue(\"novelID\",\"' + data[2] + '\");$(\"#modalDetailsNovel\").modal(\"show\");';
                                             #                  //var num = '<a id=\"goA\" role=\"button\" onclick = ' + onclick_f + ' >' + data[8] + '</a>';
                                             #                  var num = data[8];
                                             #                  $('td:eq(7)', row).html(num);
                                             #                }"
                                             # )
                                             ),
                              #extensions = 'RowGroup', 
                              width = "100%",
                              rownames = FALSE)
        )
        
      })
    } else if (!is.null(cdata[['intron']])) {

      hideElement(id = "geneInputPanel")
      showElement(id = "genePanel")
      
      removeCssClass(id = "main", "col-sm-9")
      addCssClass(id = "main", "col-sm-12")
      
      shinyjs::hide(id = "sidebar_panel")
      shinyjs::hide(id = "intronPanelOutput")
      
      ## Text output with the info from the selected intron 
      
      # output$intronPanelOutput <- renderUI({
      #   
      #   tagList(
      #     h3("Intron selected:"),
      #     hr(),
      #     
      #     p(strong("Intron: "),paste0("'", cdata[['coordinates']],"'.")),
      #     p(strong("Intron length: "),paste0(cdata[['length']]," bp.")),
      #     p(strong("Intron mis-spliced at: "),paste0("'", cdata[['type']],"' splice site.")),
      #     p(strong("Intron with ClinVar mutations: "),paste0(cdata[['clinvar']],".")),
      #     hr(),
      #     p(strong("Gene: "), paste0(cdata[['gene_tab1']],".")),
      #     
      #     hr(),
      #     p(strong("Sample from: "),paste0("'", cdata[['tissue']],"'."))
      #   )
      # 
      #   
      # })
      
      
      output$intronGeneDetail <- renderUI({
      
        tagList(
        
          h2(paste0("Novel events attached to intron 'ID=", cdata[['intron']],"':")),
          
          DT::renderDataTable(get_novel_data(intron = cdata[['intron']],
                                             db = cdata[['db']],
                                             cluster = cdata[['cluster']]), 
                              options = list(pageLength = 20,
                                             order = list(9, 'desc'),
                                             columnDefs = list(list(visible=FALSE, targets=c(2))),
                                             autoWidth = F#,
                                       #       rowGroup = list(dataSrc = 1), rowCallback = DT::JS("function(row, data) {
                                       #                        var onclick_f = 'Shiny.setInputValue(\"novelID\",\"' + data[2] + '\");$(\"#modalDetailsNovel\").modal(\"show\");';
                                       #                        //var num = '<a id=\"goA\" role=\"button\" onclick = ' + onclick_f + ' >' + data[8] + '</a>';
                                       #                        var num = data[8];
                                       #                        $('td:eq(7)', row).html(num);
                                       #                      }"
                                       # )
                                       ),
                        #extensions = 'RowGroup', 
                        width = "100%",
                        rownames = FALSE)
        )
        
      })
    } else {
      showElement(id = "geneInputPanel")
      hideElement(id = "genePanel")
    }
    
  })
  
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
  
  
  
  
  ## Popup modal window with the detail of the individuals presenting an annotated intron ------------------------------------------
  
  output$modalIntronDetail <- renderUI({
    
    #cdata <- parseQueryString(session$clientData$url_search)
    
    tagList(
      
      DT::renderDataTable(get_intron_details(intron_id = input$intronID,
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
      
      DT::renderDataTable(get_novel_details(novel_id = input$novelID,
                                      tissue = cdata[['tissue']]), 
                          options = list(pageLength = 20,
                                         autoWidth = T,
                                         order = list(2, 'desc')),
                          width = "100%",
                          rownames = FALSE)
    )
    
  })
  
  
  
  
  ## Get all annotated introns from the selected gene -----------------------------------------------------------------------------

  geneSearchUI <- eventReactive(c(input$geneButton,input$threshold), {
    
    req(!is.null(input$clusters_tab1),
        !is.null(input$data_bases_tab1), cancelOutput = T)

    
    
    if (input$radiobutton_introntype_tab1 == "introns") {
      if (input$radiobutton_searchtype_tab1 == "radio_bycoordinates_tab1") 
        title <- "Annotated intron - Details"
      else
        title <- "Annotated introns - Details"
    }
      
    else {
      if (input$radiobutton_searchtype_tab1 == "radio_bycoordinates_tab1") 
        title <- "Novel junction - Details"
      else
        title <- "Novel junctions - Details"
    }
      
    
    
    IDB_data <- get_gene_intron_data(type = input$radiobutton_introntype_tab1,
                                     chr = input$chr_tab1,
                                     start = input$start_tab1,
                                     end = input$end_tab1,
                                     strand = input$strand_tab1,
                                     gene = input$gene_tab1,
                                     threshold = input$threshold,
                                     search_type = input$radiobutton_searchtype_tab1,
                                     data_bases = input$data_bases_tab1,
                                     clusters = input$clusters_tab1,
                                     mane = F, 
                                     clinvar = input$clinvar)
    
    
    if (any(names(IDB_data) == "Message")) {
      tagList(
        h1(title),
        DT::renderDataTable(IDB_data)
      )
      
    } else {
      
      table3_caption = paste0("MSR_D = Mis-splicing Ratio at Donor (5'ss) position. MSR_A = Mis-splicing Ratio at Acceptor (3'ss) position.\nReference annotation used: Ensembl v104 (March 2021)")
      
      if (input$radiobutton_searchtype_tab1 == "radio_bycoordinates_tab1") {
        
        info <- p(strong("Coordinates: "),paste0("'", IDB_data$Coordinates %>% unique(), "'."), 
                  br(),
                  strong("Width: "),paste0("'", IDB_data$Width %>% unique(), "'."),
                  br(),
                  strong("Ss5score: "),paste0("'", IDB_data$Ss5score %>% unique(), "'."),
                  br(),
                  strong("Ss3score: "),paste0("'", IDB_data$Ss3score %>% unique(), "'."),
                  br(),
                  strong("ClinVar: "),paste0("'", IDB_data$ClinVar %>% unique(), "'."),
                  br(),
                  strong("Gene: "),paste0("'", IDB_data$Gene %>% unique(), "'."))
       
        IDB_data <- IDB_data %>%
          select(-Width, -Ss5score, -Ss3score, -ClinVar, -Gene)
        
      } else {
        info <- paste0("All ", str_to_lower(title), " from ", input$gene_tab1, " gene")
        
        if (input$radiobutton_introntype_tab1 == "novel") {
          IDB_data <- IDB_data %>%
            select(-ClinVar)
        }
        
      }
     
      IDB_data <- IDB_data %>%
        mutate("More" = "See") %>%
        relocate("More", .before = "sample")
      
      
      print(IDB_data)
      tagList(
        
        h1(title),
        br(),
        
        div(info),
        br(),
        
        DT::renderDataTable(IDB_data,
                            options = list(pageLength = 20,
                                           columnDefs = list(list(visible=FALSE, 
                                                                  targets=c(15))),
                                           
                                           order = list(0, 'asc'),
                                           rowGroup = list(dataSrc = 0),
                                           autoWidth = F,
                                           rowCallback = DT::JS(
                                       "function(row, data) {
                                          if (data[1].includes('novel')) {
                                            // It's the novel junction view

                                            var href = encodeURI('http://localhost:8787/p/3772/?id=' + data[0]);
                                            var num = '<a id=\"goA\" role=\"button\" target=\"_blank\" href=' + href + '> Check across the IDB </a>';
                                            $('td:eq(14)', row).html(num);


                                          } else if (data[1] != 'never') {
                                            // It's the intron view
              
                                            var href = encodeURI('http://localhost:8787/p/3772/?intron=' + data[0] + '&db=' + data[13] + '&cluster=' + data[15]);
                                            var num = '<a id=\"goA\" role=\"button\" target=\"_blank\" href=' + href + '> Check missplicing events </a>';
$('td:eq(14)', row).html(num);
  
  
                                            //var onclick_f = 'Shiny.setInputValue(\"intronID\",\"' + data[2] + '\");$(\"#modalDetailsIntron\").modal(\"show\");';
                                            //var num = '<a id=\"goA\" role=\"button\" onclick = ' + onclick_f + ' >' + data[7] + '</a>';
                                            //var num = data[7];
                                            //$('td:eq(6)', row).html(num);

                                            //var href = encodeURI('http://localhost:8787/p/3220/?intron=' + data[2] + '&coordinates=' + data[0] + '&gene=' + data[11] + '&type=' + data[1] + '&clinvar=' + data[10] + '&length=' + data[3] + '&db=' + data[12] + '&cluster=' + data[13]);
                                            //var num = '<a id=\"goA\" role=\"button\" target=\"_blank\" href=' + href + '>' + data[0] + '</a>';
                                            //$('td:eq(0)', row).html(num);
  
  
                                            //var onclick_f = 'Shiny.setInputValue(\"intronID\",\"' + data[2] + '\");$(\"#modalDetailsIntron\").modal(\"show\");';
                                            //var num = '<a id=\"goA\" role=\"button\" onclick = ' + onclick_f + ' >' + data[7] + '</a>';
                                            //var num = data[7];
                                            //$('td:eq(6)', row).html(num);
                                          } else {
                                            //var num = data[0];
                                            //$('td:eq(0)', row).html(num);
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
  
  
  ## What happens after the user clicks the 'Accept' button -----------------------------------------------------------------------
  
  output$geneOutput = renderUI({
    # showElement(id = "geneInputPanel")
    # hideElement(id = "genePanel")
    if(!is.null(input$clusters_tab1) &&
       !is.null(input$data_bases_tab1) &&
       !is.null(input$radiobutton_introntype_tab1)) 
    geneSearchUI()
    
  })
  
  
  # output$db_warning <- renderText({
  #   db_validation()
  # })
  # 
  # output$cluster_warning <- renderText({
  #   cluster_validation()
  # })
  
  
  
  # ##################################################
  # ## TAB 'TWO'
  # ##################################################
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

