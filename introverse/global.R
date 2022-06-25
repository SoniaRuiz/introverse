library(DBI)
library(tidyverse)

# Sys.setenv("AWS_ACCESS_KEY_ID"=Sys.getenv("AWS_ACCESS_KEY_ID"),
#            "AWS_SECRET_ACCESS_KEY"=Sys.getenv("AWS_SECRET_ACCESS_KEY"),
#            "AWS_DEFAULT_REGION"=Sys.getenv("AWS_DEFAULT_REGION"))


con <- DBI::dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
DBI::dbListTables(con)

# setwd(dir = "intron_db/")



############################################
## STATS FROM THE iSAdb
############################################

## METADATA --------------------------------

p <- FALSE

if (p) {
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
  
  ## Number of samples used
  df_metadata %>% nrow()
  
  ## Number of tissues used
  df_metadata$tissue %>% unique() %>% length()
  
  ## Number of samples per tissue
  df_metadata %>%
    count(tissue) %>%
    arrange(n)
  
  
  
  ## INTRONS ----------------------------------
  
  
  query <- paste0("SELECT * FROM 'intron'")
  df_intron <- dbGetQuery(con, query) %>% as_tibble()
  
  df_intron$ref_junID %>% unique() %>% length()
  
  
  
  
  
  ## NOVEL ------------------------------------
  
  
  query <- paste0("SELECT * FROM 'novel'")
  df_novel <- dbGetQuery(con, query) %>% as_tibble()
  df_novel
  
  df_novel %>% filter(novel_type == "novel_acceptor") %>% distinct(novel_coordinates) %>% nrow()
  df_novel %>% filter(novel_type == "novel_donor") %>% distinct(novel_coordinates) %>% nrow()
  
  
  
  ## GENES ------------------------------------
  
  ## GET FROM GENE TABLE
  query <- paste0("SELECT * FROM 'gene'")
  df_gene <- dbGetQuery(con, query)
  df_gene$gene_id %>% unique() %>% length()
  
  
  ## MIS-SPLICING ACROSS TISSUES --------------------
  
  SRA_projects <- (df_metadata$SRA_project %>% unique())
  
  
  missplicing_stats <- map_df(SRA_projects, function(db) {
    # db <- SRA_projects[1]
  
    clusters <- df_metadata %>%
      filter(SRA_project == db) %>%
      distinct(cluster) %>%
      pull()
    
    map_df(clusters, function(cluster) { 
    
      # cluster <- clusters[1]   
      
      print(cluster)
      
      query <- paste0("SELECT * FROM '", paste0(cluster, "_", db, "_misspliced"), "'") 
      df_cluster <- dbGetQuery(con, query) %>% as_tibble()
      
      df_cluster$MSR_D %>% median()
      df_cluster$MSR_A %>% median()
      
      missplicing <- df_cluster %>% count(ref_junID) %>% pull(n)
      missplicing %>% min()
      missplicing %>% max
      missplicing %>% mean
      
      return(data.frame(cluster = cluster,
                        MSR_D_mean = df_cluster$MSR_D %>% mean(),
                        MSR_A_mean = df_cluster$MSR_A %>% mean(),
                        MSR_D_median = df_cluster$MSR_D %>% median(),
                        MSR_A_median = df_cluster$MSR_A %>% median(),
                        mean_unique_junc = missplicing %>% mean ))
    })
  })
  
  missplicing_stats %>% 
    as_tibble() %>%
    arrange(MSR_D_mean) %>%
    mutate(across(where(is.numeric), round, 2)) %>%
    dplyr::select(Tissue = cluster,
                  "Mean MSR_D across mis-spliced introns" = MSR_D_mean,
                  "Mean MSR_A across mis-spliced introns" = MSR_A_mean,
                  "Mean number of unique novel junctions per mis-spliced intron" = mean_unique_junc) %>%
    left_join(df_metadata %>%
                count(tissue) %>%
                arrange(n),
              by = c("Tissue" = "tissue")) %>%
    rename("N samples" = n) %>%
    
    write.csv(file = "/home/sruiz/PROJECTS/splicing-project-app/paper/stats.csv")

}