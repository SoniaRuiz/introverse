#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

### start, end only stored but calculates the third
library(DBI)
library(tidyverse)
library(shiny)
library(shinyjs)
library(shinycssloaders)
library(shinyBS)
library(shinydashboard)

library(stringr)
library(data.table)

library(GenomicRanges)
library(ggforce)


library(ggplot2)
library(gridExtra)

library(sandwich)
library(ggstance)
library(aws.s3)

# library(jtoolsscree)

# con <- DBI::dbConnect(drv = RSQLite::SQLite(),"./dependencies/splicing.sqlite")
# setwd("./intron_db")
# con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")


# s3con <- aws.s3::s3s3connection(object = "s3://intron-db/splicing.sqlite")
# aws.s3::object_exists(object = "s3://intron-db/splicing.sqlite")

source("./get_missplicing.R")
con <- DBI::dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
DBI::dbListTables(con)

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
                 shiny::checkboxInput(inputId = "all_tissues_tab1",
                                      label = "All tissues",
                                      value = F),
                 shiny::selectizeInput(inputId = "data_bases_tab1",
                                       label = "Body region:",
                                       choices = NULL,
                                       multiple = TRUE,
                                       options = list(
                                         placeholder = "",
                                         maxItems = 3,
                                         options = list(create = FALSE)),
                                       selected = NULL),
                 shiny::selectizeInput(inputId = "clusters_tab1",
                                       label = "Samples:",
                                       choices = NULL,
                                       multiple = TRUE,
                                       options = list(
                                         placeholder = "",
                                         #maxItems = 15,
                                         options = list(create = FALSE)),
                                       selected = NULL),
                 
                 hr(),
                 
                 
                 ## Option 6
                 p(strong("Support for novel annotation:")),
                 shiny::checkboxInput(inputId = "novel_annotation_tab1",
                                      label = "Filter introns with potential novel annotation",
                                      value = F),
                 sliderInput(inputId = "threshold_tab1", 
                             label = "% individuals:",
                             min = 10, 
                             max = 90,
                             value = 10,
                             step = 10,
                             round = T),
                 
                 shiny::span(id = "span_threshold_tab1", 
                             "NOTE: this is the minimum % individuals in which any novel junction attached to the intron of interest is required to be found.", style = "font-size:90%;"),
                 hr(),
                 
                 
                 
                 
                 ## Option 5
                 p(strong("Additional filters:")),
                 splitLayout(
                   checkboxInput(inputId = "clinvar_tab1", 
                                 label = "ClinVar", value = FALSE),
                   checkboxInput(inputId = "mane_tab1", 
                                 label = "MANE Select", value = T)
                 ),
                 
                
                 
                 
                 ## BUTTON
                 hr(),
                 actionButton(inputId = "geneButton_tab1", label = "Accept"),
                 hidden(
                   p(id = "novelID_tab1", ""),
                   p(id = "intronID_tab1", ""),
                   p(id = "db_tab1", ""),
                   p(id = "cluster_tab1", "")
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
                uiOutput("geneOutput_tab1"),# %>% withSpinner(color="#0dc5c1"),
                uiOutput("intronGeneDetail_tab1"),
  
                bsModal(id = "modalVisualiseTranscript_tab1",
                        title = NULL,
                        trigger = NULL,
                        size = "large",
                        plotOutput("modalVisualiseTranscript_tab1")),
                bsModal(id = "modalVisualiseTranscriptNovel_tab1",
                        title = NULL,
                        trigger = NULL,
                        size = "large",
                        plotOutput("modalVisualiseTranscriptNovel_tab1")),
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
             ))
)

################################################
## SERVER SIDE
################################################

server <- function(input, output, session) {
   
  
  ##################################################
  ## TAB 'ONE'
  ##################################################
  

  ## Fill the dropdowns and hide/show inputs --------------------------------------------------------------------------------------
  
  updateSelectizeInput(session, 'chr_tab1', choices = chr_choices, server = TRUE, selected = "19")
  updateSelectizeInput(session, 'gene_tab1', choices = genes_choices, server = TRUE, selected = "ENSG00000171862")
  updateSelectizeInput(session, 'data_bases_tab1', choices = db_choices, server = TRUE, selected = "BLOOD")
 
  
  
  shinyjs::hideElement(id = "chr_strand_tab1")
  shinyjs::hideElement(id = "start_end_tab1")
  #shinyjs::hideElement(id = "gene_tab1")
  # shinyjs::disable(id = "geneButton")
  
  
  #shinyjs::hideElement(id = "geneInputPanel_tab1")
  #shinyjs::hideElement(id = "genePanel_tab1")
  shinyjs::disable(id = "threshold_tab1")
  shinyjs::disable(id = "radiobutton_searchtype_tab1")
  
  
  ## Observers ----------------------------------------------------------------------------------------------------------------------
  
  # observeEvent(input$intron_db, {
  #   if (intron$intron_db == "one")
  #     refresh()
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
  
  ## Hierarchical BODY SAMPLES select boxes
  observeEvent(input$data_bases_tab1, {
    
    ## In any case, we stop the reaction of the button
    #freezeReactiveValue(x = input, "geneButton_tab1")

    if (any(input$data_bases_tab1 == "all")) {
      
      ## Empty the select and add only the 'all' option.
      updateSelectizeInput(session = session, inputId = 'clusters_tab1', choices = c("All" = "all"), server = TRUE, selected = "all")
      
      
    } else {
      
      query <- paste0("SELECT * FROM 'master'")
      con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
      df_all_projects_metadata <- dbGetQuery(conn = con, statement = query) 
      DBI::dbDisconnect(conn = con)
      
      choices <- NULL
      
      # projects <- c("GTEx", "PD", "HD")
      projects <- input$data_bases_tab1
      
      for (db in projects) { # db <- data_bases_tab1[1]
      
        data_bases <- df_all_projects_metadata %>%
          filter(SRA_project == db) %>%
          dplyr::arrange(cluster_tidy)  
        
        clusters <- data_bases$cluster %>% unique() %>% as.list()
        names(clusters) <- data_bases$cluster_tidy %>% unique() %>% as.list()
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
  
  ## Checkbox 'All Tissues' 
  observeEvent(input$all_tissues_tab1, {
    if (input$all_tissues_tab1) {
      shinyjs::disable(id = "data_bases_tab1")
      shinyjs::disable(id = "clusters_tab1")
    } else {
      shinyjs::enable(id = "data_bases_tab1")
      shinyjs::enable(id = "clusters_tab1")
    }
    
  })
  
  ## Checkbox 'novel annotation' only
  observeEvent(input$novel_annotation_tab1, {
    if (input$novel_annotation_tab1) {
      shinyjs::enable(id = "threshold_tab1")
      shinyjs::show(id = "span_threshold_tab1")
      shiny::updateSliderInput(session = session, inputId = "threshold_tab1", value = 10)
    } else {
      shinyjs::disable(id = "threshold_tab1")
      shinyjs::hide(id = "span_threshold_tab1")
    }
  })
  
  
  
  ## Get the list of novel junctions attached to an annotated intron - shown in a different tab ------------------------------------
  
  observe({
    
    cdata <- parseQueryString(session$clientData$url_search)
    
    
    
    if (!is.null(cdata[['intron']])) {

      print(cdata[['intron']])
      
      # intron=chr10%253A87894110-87925512%253A%252B&db=BRAIN&cluster=Brain%20-%20Amygdala
      removeCssClass(id = "main_tab1", "col-sm-9")
      addCssClass(id = "main_tab1", "col-sm-12")
      
      shinyjs::hide(id = "sidebar_panel_tab1")
      shinyjs::hide(id = "intronPanelOutput_tab1")
      
      updateNavbarPage(session = session, inputId = "intron_db", selected = "one")

      novel_junctions_from_intron <- get_novel_data_from_intron(intron_id = URLdecode(URL = cdata[['intron']]),
                                                                db = URLdecode(URL = cdata[['db']]),
                                                                sample_group = URLdecode(URL = cdata[['cluster']])) %>% dplyr::mutate(More = "")
      
      output$intronGeneDetail_tab1 <- renderUI({
      
        tagList(
        
          h2(paste0("Novel events attached to intron 'ID=", URLdecode(URL = cdata[['intron']]),"':")),
          
          DT::renderDataTable(novel_junctions_from_intron, 
                              options = list(pageLength = 20,
                                             order = list(9, 'desc'),
                                             columnDefs = list(list(visible=FALSE, targets=c(1))),
                                             autoWidth = F,
                                             rowCallback = DT::JS("function(row, data) {
         
                                                var onclick_f = 'Shiny.setInputValue(\"novelID_tab1\",\"' + encodeURI(data[1]) + '\");Shiny.setInputValue(\"cluster_tab1\",\"' + encodeURI(data[10]) + '\");Shiny.setInputValue(\"db_tab1\",\"' + encodeURI(data[11]) + '\");$(\"#modalVisualiseTranscriptNovel_tab1\").modal(\"show\");';
                                                //var onclick_f = '$(\"#modalVisualiseTranscriptNovel_tab1\").modal(\"show\");';
                                          
                                                console.log(onclick_f)
                                                var num = '<a id=\"goA\" role=\"button\" onclick = ' + onclick_f + ' ><button>Visualize MANE transcript</button></a>';
                                                $('td:eq(11)', row).html(num);
                                                  
                                             }"
                                             )),
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
          
          h2(paste0("Details of the novel junction 'ID=", URLdecode(URL = cdata[['id']]),"' across all projects from the IDB:")),
          
          DT::renderDataTable(get_novel_data_across_idb(novel_id = URLdecode(URL = cdata[['id']])), 
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
  output$modalVisualiseTranscript_tab1 <- renderPlot({

    visualise_transcript(intron_id = str_replace_all(string = input$intronID_tab1, pattern = "%20", replacement = " "),
                                     db = str_replace_all(string = input$db_tab1, pattern = "%20", replacement = " "),
                                     clust = str_replace_all(string = input$cluster_tab1, pattern = "%20", replacement = " ") )


  }, width = "auto", height = "auto")
  
  output$modalVisualiseTranscriptNovel_tab1 <- renderPlot({
    
    visualise_transcript(novel_id = str_replace_all(string = input$novelID_tab1, pattern = "%20", replacement = " "),
                         db = str_replace_all(string = input$db_tab1, pattern = "%20", replacement = " "),
                         clust = str_replace_all(string = input$cluster_tab1, pattern = "%20", replacement = " ") )
    
    
    # ggplot() +
    #   theme_void() +
    #   geom_text(aes(0,0,label='Coming soon...')) + 
    #   theme(text = element_text(element_text(size = "14"))) 
    
    
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


  observeEvent(input$geneButton_tab1,  {
    
    output$geneOutput_tab1 = renderUI({

    title <- "Annotated introns - Details"
    
    threshold <- input$threshold_tab1
    if (!input$novel_annotation_tab1) {
      threshold <- -1
    }

    IDB_data <- main_IDB_search(type = "introns",
                              chr = input$chr_tab1,
                              start = input$start_tab1,
                              end = input$end_tab1,
                              strand = input$strand_tab1,
                              gene = input$gene_tab1,
                              threshold = threshold,
                              search_type = input$radiobutton_searchtype_tab1,
                              all_data_bases = input$all_tissues_tab1,
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
      
      table3_caption = paste0("MSR_D = normalised mis-splicing ratio at the donor (5'ss) position.
                              MSR_A = normalised mis-splicing ratio at the acceptor (3'ss) position.\n
                              Reference annotation used: Ensembl v105 (Dec 2021)")
      
      if (input$radiobutton_searchtype_tab1 == "radio_bycoordinates_tab1") {
        
        info <- paste0("Reference intron '", IDB_data$ID %>% unique(), "' from ", IDB_data$Gene %>% unique, " gene.")
        
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
        print(IDB_data)
        info <- paste0("All ", str_to_lower(title), " from ", IDB_data$Gene %>% unique, " gene")
      }
      
     
      IDB_data <- IDB_data %>%
        mutate("More" = "See")
      
      
      # print(IDB_data)
      
      tagList(
        
        h1(title),
        br(),
        
        div(info),
        br(),
        
        DT::renderDataTable(IDB_data,
                            options = list(pageLength = 20,
                                           columnDefs = list(list(visible=FALSE, targets=c(13))),
                                           order = list(0, 'asc'),
                                           rowGroup = list(dataSrc = 0),
                                           autoWidth = F,
                                           rowCallback = DT::JS("function(row, data) {
                                           if (data[1] != 'never') {
                                           
                                                // It's the intron view
                                                var intron_coord = encodeURIComponent(data[0])
                                                var href = encodeURI('https://soniagarciaruiz.shinyapps.io/intron_db/?intron=' + data[13] + '&db=' + data[12] + '&cluster=' + data[11]);
                                                var num = '<a id=\"goA\" role=\"button\" target=\"_blank\" href=' + href + '> Open missplicing events </a>';
                                                //var num = 'Check missplicing events';
                                                
    
                                                var onclick_f = 'Shiny.setInputValue(\"intronID_tab1\",\"' + encodeURI(data[13]) + '\");Shiny.setInputValue(\"db_tab1\",\"' + encodeURI(data[12]) + '\");Shiny.setInputValue(\"cluster_tab1\",\"' + encodeURI(data[11]) + '\");$(\"#modalVisualiseTranscript_tab1\").modal(\"show\");';
                                                
                                                console.log(onclick_f)
                                                num = num + '<br/><a id=\"goA\" role=\"button\" onclick = ' + onclick_f + ' ><button>Visualize MANE</button></a>';
                                                $('td:eq(13)', row).html(num);
                                            
                                          } else {
                                            var num = 'Intron with correct splicing';
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
  
 
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)

