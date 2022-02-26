###################################
## CONNECTION TO THE DB  
###################################

library(DBI)
library(tidyverse)

# Create an ephemeral in-memory RSQLite database
setwd("./intron_db/")
dir.create(file.path("./dependencies"), showWarnings = F, recursive = T)
con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")

dbListTables(con)
# dbDisconnect(conn = con)
## Check the schema


query <- paste0("SELECT * from sqlite_schema")
dbGetQuery(con, query)

###################################
## REMOVE ALL TABLES
###################################
con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=0")
tables <- dbListTables(con)
for (table in tables) {
  if (table != "master" && table != "gene_name" && table != "mane") {
    dbRemoveTable(conn = con, table)
    print(paste0("Table: '", table, "' has been removed!"))
  }
}
dbListTables(con)


###################################
## CREATE MANE TABLE
###################################

# hg_MANE <- rtracklayer::import(con = "/data/references/MANE/MANE.GRCh38.v1.0.ensembl_genomic.gtf")
hg_MANE_tidy <- hg_MANE %>%
  as_tibble() %>%
  select(-source, -score, -phase, -gene_id, -gene_type, -tag, -protein_id, -db_xref,-transcript_type,-exon_id,-exon_number ) %>%
  mutate(transcript_id = transcript_id %>% str_sub(start = 1, end = 15)) %>%
  drop_na()


DBI::dbWriteTable(conn = con,
                  name = "mane",
                  value = hg_MANE_tidy,
                  overwrite = T)

# dbRemoveTable(conn = con, "mane")

###################################
## CREATE MASTER TABLE
###################################
SRA_projects <- c("BRAIN","SRP058181","SRP051844")

df_all_projects_metadata <- map_df(SRA_projects, function(project) { 
  
  print(paste0(Sys.time(), " - getting info from ", project))
  # project <- SRA_projects[1]
  # project <- SRA_projects[2]
  # project <- SRA_projects[3]
  
  
  metadata <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", project, "/raw_data/samples_metadata.rds"))
  
  if (project == "BRAIN") {
    
    
    df_project_metadata <- data.frame(age = metadata$gtex.age,
               rin = metadata$gtex.smrin %>% as.character(),
               gender = metadata$gtex.sex %>% as.character(),
               tissue = metadata$gtex.smtsd,
               cluster = metadata$gtex.smtsd,
               diagnosis = paste0("GTExV8 - ", metadata$gtex.smtsd),
               avg_read_length = metadata$recount_seq_qc.avg_len,
               mapped_read_count = metadata$recount_qc.star.all_mapped_reads,
               SRA_project_tidy = "GTEXv8 - BRAIN",
               SRA_project = metadata$recount_project.project)
  } else {
      
      ## EXTRACT METADATA
      
      df_project_metadata <- map_df(metadata$sra.sample_attributes %>% as.vector, function(characteristic) {
        
        # characteristic <- (metadata$sra.sample_attributes %>% as.vector)[[1]]
        #str_replace(string = characteristic, replacement = "##", pattern = '\\|')
        characteristic <- str_split(string = characteristic, pattern = "\\|", simplify = T)
        
        ## TISSUE
        ind <- str_detect(characteristic, pattern = "tissue;;")
        if (any(ind))
          donor_tissue <- str_sub(string = characteristic[ind], 
                                  start = str_locate(characteristic[ind], pattern = ";;")[[2]] + 1,
                                  end = characteristic[ind] %>% str_count())
        else
          donor_tissue <- NA
        
        ## GENDER
        ind <- str_detect(characteristic, pattern = "gender")
        if (any(ind))
          donor_gender <- str_sub(string = characteristic[ind], 
                                  start = str_locate(characteristic[ind], pattern = ";;")[[2]] + 1,
                                  end = characteristic[ind] %>% str_count()) %>% as.character()
        else
          donor_gender <- NA
        
        ## RIN
        ind <- str_detect(characteristic, pattern = "rin")
        if (any(ind))
          donor_rin <- str_sub(string = characteristic[ind], 
                               start = str_locate(characteristic[ind], pattern = ";;")[[2]] + 1,
                               end = characteristic[ind] %>% str_count())
        else
          donor_rin <- NA
        
        ## AGE
        ind <- str_detect(characteristic, pattern = "death")
        if (any(ind))
          donor_age <- str_sub(string = characteristic[ind], 
                               start = str_locate(characteristic[ind], pattern = ";;")[[2]] + 1,
                               end = characteristic[ind] %>% str_count())
        else
          donor_age <- NA
        
        ind <- str_detect(characteristic, pattern = "diagnosis")
        if (any(ind))
          donor_diagnosis <- str_sub(string = characteristic[ind], 
                                     start = str_locate(characteristic[ind], pattern = ";;")[[2]] + 1,
                                     end = characteristic[ind] %>% str_count())
        else
          donor_diagnosis <- NA
        
        
        
        
        data.frame(age = donor_age,
                   rin = donor_rin,
                   gender = donor_gender %>% as.character(),
                   tissue = donor_tissue,
                   diagnosis = donor_diagnosis,
                   stringsAsFactors = F) %>% 
          return()
      
        
      
    })
    ## CLUSTER
    group <- str_sub(string = metadata$sra.sample_title, 
                     start = 1,
                     end = str_locate(metadata$sra.sample_title, pattern = "_")[[1]]) 
      
      
    df_project_metadata <- df_project_metadata %>%
      mutate(cluster = group,
             avg_read_length = metadata$recount_seq_qc.avg_len,
             mapped_read_count = metadata$recount_qc.star.all_mapped_reads,
             SRA_project = metadata$recount_project.project)
    
    
    if (project == "SRP058181") {
      if (any(df_project_metadata %>%
              filter(SRA_project == project) %>%
              pull(diagnosis) %>% is.na())) {
        df_project_metadata[df_project_metadata$cluster == "C_", "diagnosis"] <- "PD/Control - Neurologically normal"
        df_project_metadata[df_project_metadata$cluster == "P_", "diagnosis"] <- "PD/Control - Parkinson's Disease"
      }
      
      df_project_metadata[, "SRA_project_tidy"] <- "PD/Control"
      
    } else if (project == "SRP051844") {
      if (any(df_project_metadata %>%
              filter(SRA_project == project) %>%
              pull(diagnosis) %>% is.na())) {
        df_project_metadata[df_project_metadata$cluster == "C_", "diagnosis"] <- "Control"
        df_project_metadata[df_project_metadata$cluster == "H_", "diagnosis"] <- "HD"
      }
      df_project_metadata[, "SRA_project_tidy"] <- "HD/Control"
      df_project_metadata <- df_project_metadata %>%
        mutate(diagnosis = paste0("HD/Control - ", diagnosis))
    }
  }
  
  
  df_project_metadata %>% head %>% as_tibble()%>% print
  df_project_metadata %>% 
    return()
  
})

# saveRDS(object = df_all_projects_metadata,
#         file = "./dependencies/df_all_projects_metadata.rds")


DBI::dbWriteTable(conn = con,
                  name = "master",
                  value = df_all_projects_metadata,
                  overwrite = T)


###################################
## CREATE GENES TABLE
###################################

# Query to the DB
gtf_version <- "105"

# Query to the DB
query = paste0("SELECT * FROM 'master'")
con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
df_all_projects_metadata <- dbGetQuery(con, query) 


SRA_projects <- df_all_projects_metadata$SRA_project %>% unique()

gene_names <- NULL

for (db in SRA_projects) {
  
  # db <- SRA_projects[1]
  
  print(paste0(Sys.time(), " --> Working with '", db, "' DataBase..."))
  
  clusters <- df_all_projects_metadata %>%
    filter(SRA_project == db) %>%
    distinct(cluster) %>%
    pull()
  
  base_folder <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", db, "/")
  
  
  for (cluster in clusters) { 
    
    # cluster <- clusters[1]
    
    print(paste0(Sys.time(), " --> ", cluster))
    
    if (file.exists(paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
                      cluster, "/v", gtf_version, "/", cluster, "_db_introns.rds"))) {
      df_introns <- readRDS(file = paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
                                          cluster, "/v", gtf_version, "/", cluster, "_db_introns.rds"))
      
      gene_names <- c(df_introns %>%
                        select(gene_name) %>%
                        unnest(gene_name) %>% 
                        unlist() %>%
                        unname() %>%
                        unique(), gene_names) %>% unique()
    }
    
    if (file.exists( paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
                       cluster, "/v", gtf_version, "/", cluster, "_db_novel.rds"))){
      df_novel <- readRDS(file = paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
                                        cluster, "/v", gtf_version, "/", cluster, "_db_novel.rds"))
      
      gene_names <- c(df_novel %>%
                        select(gene_name) %>%
                        unnest(gene_name) %>% 
                        unlist() %>%
                        unname() %>%
                        unique(), gene_names) %>% unique()
      
    }
  }
}

## CREATE GENE_NAME TABLE 

con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")

# DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=0")
# dbRemoveTable(conn = con, "gene_name")


DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=1")
gene_table_statement <- paste0("CREATE TABLE IF NOT EXISTS 'gene_name'", 
                               "(gene_id INTEGER PRIMARY KEY NOT NULL,
                               gene_name TEXT NOT NULL)")
res <- DBI::dbSendQuery(conn = con, statement = gene_table_statement)
DBI::dbClearResult(res)

## POPULATE GENE_NAME TABLE

gene_names <- gene_names %>% 
  as_tibble() %>%
  drop_na() %>%
  dplyr::rename(gene_name = value) %>%
  tibble::rowid_to_column("gene_id")
  

any(str_detect(gene_names$gene_name, pattern = "c\\("))

DBI::dbAppendTable(conn = con,
                   name = "gene_name", 
                   value = gene_names)


DBI::dbDisconnect(conn = con)


###################################
## CREATE AND POPULATE TABLES
###################################

con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=1")

gtf_version <- "105"

# Query to the DB
query = paste0("SELECT * FROM 'master'")
df_all_projects_metadata <- dbGetQuery(con, query) 

SRA_projects <- (df_all_projects_metadata$SRA_project %>% unique())[1]

## SEARCH gene_id in gene_name table
con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
query = paste0("SELECT * FROM 'gene_name'")
all_genes <- dbGetQuery(con, query)


for (db in SRA_projects) {
  
  # db <- SRA_projects[1]
  
  print(paste0(Sys.time(), " --> Working with '", db, "' DataBase..."))
  
  clusters <- df_all_projects_metadata %>%
    filter(SRA_project == db) %>%
    distinct(cluster) %>%
    pull()

  
  base_folder <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", db, "/")
  
  for (cluster in clusters) { 
    
    print(paste0(Sys.time(), " --> ", cluster))
    # cluster <- clusters[11]
    # cluster <- clusters[1]
    # cluster <- clusters[2]
    
    ###############################
    ## CREATE INTRON TABLE
    ###############################
    #dbRemoveTable(conn = con, paste0(cluster, "_", db, "_db_introns"))
    intron_table_statement <- paste0("CREATE TABLE IF NOT EXISTS '", paste0(cluster, "_", db, "_db_introns"), "'", 
                                     "(ref_junID TEXT PRIMARY KEY NOT NULL, 
                                     ref_ss5score DOUBLE NOT NULL, ref_ss3score DOUBLE NOT NULL, 
                                     ref_type VARCHAR(10) NOT NULL, ref_n_individuals INTEGER NOT NULL, ref_mean_counts DOUBLE NOT NULL, 
                                     ref_missplicing_ratio_tissue_ND DOUBLE NOT NULL, ref_missplicing_ratio_tissue_NA DOUBLE NOT NULL, 
                                     clinvar_type TEXT NOT NULL, 
                                     MANE BOOL NOT NULL,
                                     gene_name INTEGER,
                                     FOREIGN KEY (gene_name) REFERENCES 'gene_name' (gene_id))")
    
    res <- DBI::dbSendQuery(conn = con, statement = intron_table_statement)
    DBI::dbClearResult(res)
    
    
    ###############################
    ## CREATE NOVEL TABLE
    ###############################
    
    # Create the novel junctions table
    novel_table_statement <- paste0("CREATE TABLE IF NOT EXISTS '", paste0(cluster, "_", db, "_db_novel"), "'",
                                    "(ref_junID TEXT NOT NULL, 
                                    novel_junID TEXT NOT NULL, 
                                    novel_ss5score DOUBLE NOT NULL, 
                                    novel_ss3score DOUBLE NOT NULL, 
                                    novel_type TEXT NOT NULL, 
                                    novel_n_individuals INTEGER NOT NULL, 
                                    novel_mean_counts DOUBLE NOT NULL, 
                                    distance INTEGER NOT NULL,
                                    gene_name INTEGER,
                                    PRIMARY KEY (ref_junID, novel_junID),
                                    FOREIGN KEY (ref_junID) REFERENCES '", paste0(cluster, "_", db, "_db_introns"), "' (ref_junID),
                                    FOREIGN KEY (gene_name) REFERENCES 'gene_name' (gene_id))")
    
    res <- DBI::dbSendQuery(conn = con, statement = novel_table_statement)
    DBI::dbClearResult(res)
    
    
    print(paste0(Sys.time(), ". Tables created!"))
    
    ###############################
    ## INSERT DATA
    ###############################
    
    ## INTRON TABLE
    
    if (file.exists(paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
                      cluster, "/v", gtf_version, "/", cluster, "_db_introns.rds"))) {
      
      df_introns <- readRDS(file = paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
                                          cluster, "/v", gtf_version, "/", cluster, "_db_introns.rds"))
      
      df_introns_tidy <- df_introns %>% 
        dplyr::select(-tx_id_junction) %>%
        dplyr::mutate(#strand = strand %>% as.character(),
                      #seqnames = seqnames %>% as.character(),
                      #start = start %>% as.integer(),
                      #end = end %>% as.integer(),
                      #ref_junID = ref_junID %>% as.integer(),
                      ref_mean_counts = round((ref_sum_counts / ref_n_individuals), digits = 2),
                      gene_name = gene_name %>% as.character()) %>%
        filter(u2_intron == T) %>%
        distinct(ref_junID, .keep_all = T) %>%
        select(ref_junID,
               #seqnames, start, end, strand,
               ref_ss5score, ref_ss3score,
               ref_type, ref_n_individuals, ref_mean_counts, 
               ref_missplicing_ratio_tissue_ND, ref_missplicing_ratio_tissue_NA,
               clinvar_type, 
               MANE,
               #protein_coding, 
               gene_name)
      
      ## Merge data
      df_introns_tidy <- merge(x = df_introns_tidy,
                               y = all_genes,
                               by = "gene_name") %>%
        select(-gene_name) %>%
        dplyr::rename(gene_name = gene_id)
      
      
      DBI::dbAppendTable(conn = con,
                         name = paste0(cluster, "_", db, "_db_introns"), 
                         value = df_introns_tidy)
      
      
      print(paste0(Sys.time(), ". Intron table populated!"))
      
      DBI::dbReadTable(conn = con, name = paste0(cluster, "_", db, "_db_introns"))
    }
    
    ## QC
    # df_introns_tidy_test <- df_introns_tidy[1,]
    # df_introns_tidy_test[, "ref_junID"] <- NA
    # DBI::dbAppendTable(conn = con,
    #                    name = paste0(cluster, "_db_introns"), 
    #                    value = df_introns_tidy_test,
    #                    row.names = NULL)
    #  
    #  DBI::dbReadTable(conn = con, name = paste0(cluster, "_db_introns")) %>%
    #    filter(ref_junID %>% is.na())
    
    
    
    ## NOVEL TABLE
    
    
    if (file.exists( paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
                       cluster, "/v", gtf_version, "/", cluster, "_db_novel.rds"))) {
      
      df_novel_gr <- readRDS(file = paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
                                           cluster, "/v", gtf_version, "/", cluster, "_db_novel.rds"))
      
      df_novel_tidy <- df_novel_gr %>% 
        dplyr::mutate(#ref_junID = ref_junID %>% as.integer(),
                      #novel_junID = novel_junID %>% as.integer(),
                      novel_type = novel_type %>% as.character(),
                      #strand = strand %>% as.character(),
                      gene_name = gene_name %>% as.character(),
                      novel_mean_counts = round((novel_sum_counts / novel_n_individuals), digits = 2)#,
                      #seqnames = seqnames %>% as.character(),
                      #start = start %>% as.integer(), end = end %>% as.integer()
                      ) %>%
        select(ref_junID,
               novel_junID,
               #seqnames,start, end, strand,
               novel_ss5score, novel_ss3score, novel_type,
               novel_n_individuals, novel_mean_counts,
               distance, 
               #protein_coding, 
               #gene_id, 
               gene_name)
      
      
      
      ## Merge data
      df_novel_tidy <- merge(x = df_novel_tidy,
                               y = all_genes,
                               by = "gene_name") %>%
        select(-gene_name) %>%
        dplyr::rename(gene_name = gene_id)
      
      
      DBI::dbAppendTable(conn = con,
                         name = paste0(cluster, "_", db, "_db_novel"), 
                         value = df_novel_tidy)
      
      print(paste0(Sys.time(), ". Novel table populated!"))
    }
    # ## QC
    # df_novel_tidy_test <- df_novel_tidy[1,]
    # df_novel_tidy_test[, "ref_junID"] <- NULL
    # df_novel_tidy_test[, "novel_junID"] <- NA
    # DBI::dbAppendTable(conn = con,
    #                    name = paste0(cluster, "_db_novel"), 
    #                    value = df_novel_tidy_test)
    # 
    # DBI::dbReadTable(conn = con, name = paste0(cluster, "_", db, "_db_novel")) %>%
    #   filter(is.na(ref_junID))
    
  }
}
dbListTables(con)
DBI::dbDisconnect(conn = con)



###################################
## CREATE 'GENE' TABLE
###################################

tables <- dbListTables(con)
all_genes <- NULL
for (table in tables) {
  if (table != "master") {
    #table <- tables[1]
    # Query to the DB
    query = paste0("SELECT gene_name FROM '", table, "'")
    con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
    all_genes <- c(all_genes, dbGetQuery(con, query)) %>% unlist %>% unique()
    DBI::dbDisconnect(conn = con)
    
  }
}
saveRDS(object = all_genes %>% as.character(), file = "./dependencies/all_genes_names.rds")
dbListTables(con)
