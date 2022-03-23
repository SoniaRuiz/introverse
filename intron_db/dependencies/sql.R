###################################
## CONNECTION TO THE DB  
###################################

library(DBI)
library(tidyverse)
library(data.table)

# Create an ephemeral in-memory RSQLite database
# setwd("./intron_db/")
#dir.create(file.path("./dependencies"), showWarnings = F, recursive = T)
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
  #if (str_detect(table %>% tolower(), pattern = "brain")) {
    dbRemoveTable(conn = con, table)
    print(paste0("Table: '", table, "' has been removed!"))
  #}
}
dbListTables(con)

SRA_projects <- c("BRAIN", "ADRENAL_GLAND", "KIDNEY", "SMALL_INTESTINE", "SALIVARY_GLAND",
                  "SPLEEN", "LIVER", "BONE_MARROW", "OVARY", "VAGINA", "UTERUS", "BLADDER", "BLOOD",
                  "CERVIX_UTERI", "FALLOPIAN_TUBE", "STOMACH", "ESOPHAGUS",
                  "SKIN", "PANCREAS", "BREAST", "TESTIS", "PITUITARY", "PROSTATE", "BLOOD_VESSEL",
                  "ADIPOSE_TISSUE", "HEART", "MUSCLE", "COLON", "THYROID", "NERVE", "LUNG")

###################################
## CREATE MASTER TABLE
###################################

create_master_table <- function(GTEx = T) {
  
  con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
  
  
  
  df_metadata <- map_df(SRA_projects, function(project) { 
    
    print(paste0(Sys.time(), " - getting info from ", project))
    # project <- SRA_projects[1]
    # project <- SRA_projects[2]
    # project <- SRA_projects[3]
    
    
    metadata <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", 
                                      project, "/raw_data/samples_metadata.rds"))
    
    if (GTEx) {
      
      
      df_project_metadata <- data.frame(age = metadata$gtex.age,
                                        rin = metadata$gtex.smrin %>% as.character(),
                                        gender = metadata$gtex.sex %>% as.character(),
                                        tissue = metadata$gtex.smtsd,
                                        cluster = metadata$gtex.smtsd, #str_remove_all(metadata$gtex.smtsd, pattern = " ") %>% tolower(),
                                        cluster_tidy = paste0("GTExV8 - ", metadata$gtex.smtsd),
                                        avg_read_length = metadata$recount_seq_qc.avg_len,
                                        mapped_read_count = metadata$recount_qc.star.all_mapped_reads,
                                        SRA_project_tidy = paste0("GTEXv8 - ", metadata$recount_project.project),
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
    
    
    #df_project_metadata %>% head %>% as_tibble()%>% print
    df_project_metadata %>% 
      return()
    
  })
  
  # saveRDS(object = df_all_projects_metadata,
  #         file = "./dependencies/df_all_projects_metadata.rds")
  
  
  DBI::dbWriteTable(conn = con,
                    name = "master",
                    value = df_metadata,
                    overwrite = T)
  
  DBI::dbDisconnect(conn = con)
}

###################################
## CREATE MANE TABLE
###################################

create_mane_table <- function() {
  
  con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
  
  ## LOAD THIS FILE
  hg_MANE <- rtracklayer::import(con = "/data/references/MANE/MANE.GRCh38.v1.0.ensembl_genomic.gtf")
  
  hg_MANE_tidy <- hg_MANE %>%
    as_tibble() %>%
    select(-source, -score, -phase, -gene_id, -gene_type, -tag, -protein_id, 
           -db_xref,-transcript_type,-exon_id,-exon_number, -width ) %>%
    mutate(transcript_id = transcript_id %>% str_sub(start = 1, end = 15)) %>%
    drop_na()
  
  
  DBI::dbWriteTable(conn = con,
                    name = "mane",
                    value = hg_MANE_tidy,
                    overwrite = T)
  
  # dbRemoveTable(conn = con, "mane")
  DBI::dbDisconnect(conn = con)
}


###################################
## CREATE GENE TABLE
###################################

create_gene_table <- function() {
  
  # DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=0")
  # dbRemoveTable(conn = con, "gene")
  
  con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
  DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=1")
  
  ## GET ALL DATA FROM MASTER
  query = paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
  
  
  gtf_version <- "105"
  SRA_projects <- (df_metadata$SRA_project %>% unique())
  
  ## Loop through the cluster to obtain all the genes
  gene_ids <- map_df(SRA_projects, function(db) {
    
    # db <- SRA_projects[1]
    
    print(paste0(Sys.time(), " --> Working with '", db, "' DataBase..."))
    base_folder <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", db, "/")
    
    clusters <- df_metadata %>%
      filter(SRA_project == db) %>%
      distinct(cluster) %>%
      pull()
    
    map_df(clusters, function(cluster) { 
      
      # cluster <- clusters[1]
      
      print(paste0(Sys.time(), " --> ", cluster))
      
      if (file.exists(paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
                             cluster, "/v", gtf_version, "/", cluster, "_db_introns.rds"))) {
        
        df_introns <- readRDS(file = paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
                                            cluster, "/v", gtf_version, "/", cluster, "_db_introns.rds"))
        
        genes <- df_introns %>% select(gene_id) %>% unlist(use.names = F, recursive = T)
      
      
        if (file.exists( paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
                                cluster, "/v", gtf_version, "/", cluster, "_db_novel.rds"))){
          df_novel <- readRDS(file = paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
                                            cluster, "/v", gtf_version, "/", cluster, "_db_novel.rds"))
          

          genes <- c(genes, df_novel %>% select(gene_id) %>% unlist(use.names = F, recursive = T))
         
          
        }
        
        return(data.frame(gene_id = genes %>% unique()))
        
      } else {
        
        return(NULL)
        
      }
    })
    
  })
  
  gene_ids <- gene_ids %>% distinct(gene_id)
  gene_ids %>% nrow()
  
  hg38 <- rtracklayer::import(con = "/data/references/ensembl/gtf/v105/Homo_sapiens.GRCh38.105.chr.gtf")
  hg38 %>% as_tibble %>% distinct(gene_id, .keep_all = T) %>% nrow()
  hg38 <- hg38 %>%
    as_tibble() %>%
    select(gene_id, gene_name) %>%
    mutate(gene_id = str_sub(gene_id, start = 1,end = 15)) %>%
    distinct(gene_id, .keep_all = T) %>%
    filter(gene_id %in% gene_ids$gene_id)
  
  
  
  ## CREATE GENE_NAME TABLE ---------------------------------------------------
  
  sql_statement <- paste0("CREATE TABLE IF NOT EXISTS 'gene'", 
                                 "(id INTEGER PRIMARY KEY NOT NULL,
                                 gene_id TEXT NOT NULL,
                                 gene_name TEXT)")
  res <- DBI::dbSendQuery(conn = con, statement = sql_statement)
  DBI::dbClearResult(res)
  
  #sql_statement <- paste0("CREATE UNIQUE INDEX gene_index ON gene(id)");
  #res <- DBI::dbSendQuery(conn = con, statement = sql_statement)
  #DBI::dbClearResult(res)
  
  
  ## POPULATE GENE_NAME TABLE -------------------------------------------------
  
  hg38_genes <- hg38 %>% 
    as_tibble() %>%
    #drop_na() %>%
    #dplyr::rename(name = value) %>%
    tibble::rowid_to_column("id")
  
  hg38_genes %>%
    filter(str_detect(gene_name, pattern = "c\\("))
  any(str_detect(hg38_genes$gene_name, pattern = "c\\("))
  
  DBI::dbAppendTable(conn = con,
                     name = "gene", 
                     value = hg38_genes)
  
  
  
  DBI::dbDisconnect(conn = con)
}


###################################
## CREATE INTRON TABLE
##################################

create_intron_table <- function() {
  
  ## Connect to the DB
  con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
  DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=1")
  # dbRemoveTable(conn = con, "intron")
  
  ## GET ALL GENES
  query = paste0("SELECT * FROM 'gene'")
  all_genes <- dbGetQuery(con, query)
  
  # GET INFO FROM MASTER
  query = paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
  
  
  ###################################
  ## CREATE INTRONS TABLE
  ###################################
  
  # DECLARE VARIABLES
  gtf_version <- "105"
  SRA_projects <- (df_metadata$SRA_project %>% unique())
  
  
  ## LOOP THROUGH PROJECTS
  df_all_introns <- map_df(SRA_projects, function(db) {
    
    # db <- SRA_projects[1]
    
    print(paste0(Sys.time(), " --> Working with '", db, "' DataBase..."))
    base_folder <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", db, "/")
    
    clusters <- df_metadata %>%
      filter(SRA_project == db) %>%
      distinct(cluster) %>%
      pull()

    map_df(clusters, function(cluster) {
      
      # cluster <- clusters[1]
      
      print(paste0(Sys.time(), " --> ", cluster))
      
      if (file.exists(paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
                             cluster, "/v", gtf_version, "/", cluster, "_db_introns.rds"))) {
        
        df_introns <- readRDS(file = paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
                                            cluster, "/v", gtf_version, "/", cluster, "_db_introns.rds"))
        
        df_introns_tidy <- df_introns %>% 
          dplyr::select(-tx_id_junction) %>%
          dplyr::mutate(#gene_id = gene_id %>% as.character(),
                        ref_mean_counts = round((ref_sum_counts / ref_n_individuals), digits = 2)) %>%
          filter(u2_intron == T) %>%
          distinct(ref_junID, .keep_all = T) %>%
          select(ref_junID,
                 ref_ss5score, ref_ss3score,
                 clinvar_type, 
                 MANE,
                 gene_id)
        
        return(df_introns_tidy)
        
      } else {
        
        return(NULL)
      }
    })  
  })  
  
  ## Check that all introns, albeit repeated, have been linked to the same gene
  ## through the different tissues
  # df_all_introns %>% nrow()
  # df_all_introns %>% head()
  # df_all_introns %>% group_by(ref_junID) %>% distinct(gene_id) 
  # df_all_introns %>% distinct(ref_junID)


  ## Flatten each gene list element internally
  ## Only keep introns belonging to one single gene
  df_all_introns_tidy <- df_all_introns %>%
    distinct(ref_junID, .keep_all = T) %>% 
    #rowwise() %>%
    #filter(gene_id %>% unlist() %>% length() == 1) %>%
    mutate_if(is.list, simplify_all) %>%    
    unnest(gene_id)
  
  ## Add the GENE ID for the foreign key
  df_all_introns_tidy <- merge(x = df_all_introns_tidy %>% as.data.table(),
                               y = all_genes %>% as.data.table(),
                               by = "gene_id",
                               all.x = T) %>%
    select(-gene_name, -gene_id) %>%
    dplyr::rename(gene_id = id)
  
  
  ## Extra QC
  ## Same introns should have been assigned same scores across tissues
  # df_all_introns %>%
  #   dplyr::group_by(ref_junID, gene_id) %>%
  #   distinct(ref_ss5score, .keep_all = T) %>% nrow()
  # df_all_introns %>%
  #   dplyr::group_by(ref_junID, gene_id) %>%
  #   distinct(ref_ss3score, .keep_all = T) %>% nrow()
  # df_all_introns %>%
  #   dplyr::group_by(ref_junID) %>%
  #   distinct(ref_ss3score, .keep_all = T) %>% nrow()
  # ## Each intron should only belong to the same gene
  # df_all_introns %>%
  #   dplyr::count(ref_junID, gene_id) 
  # 
  # which(duplicated(df_all_introns_tidy$ref_junID) == T)
  # df_all_introns_tidy %>%
  #   filter(ref_junID == (df_all_introns_tidy[651,1]$ref_junID))
  
  ## CREATE INTRON TABLE ------------------------------------------------------
  
  query <- paste0("CREATE TABLE IF NOT EXISTS 'intron'",
  "(ref_junID NUMERIC PRIMARY KEY NOT NULL,
  ref_coordinates TEXT NOT NULL, 
  ref_ss5score DOUBLE NOT NULL, 
  ref_ss3score DOUBLE NOT NULL, 
  clinvar TEXT NOT NULL, 
  MANE BOOL NOT NULL,
  gene_id INTEGER NOT NULL,
  FOREIGN KEY (gene_id) REFERENCES 'gene'(id))")
  
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  
  # query <- paste0("CREATE UNIQUE INDEX 'index_intron' ON 'intron'(ref_junID)");
  # res <- DBI::dbSendQuery(conn = con, statement = query)
  # DBI::dbClearResult(res)
  
  
  ## POPULATE INTRON TABLE ----------------------------------------------------
  df_all_introns_tidy <- df_all_introns_tidy %>% 
    as_tibble() %>%
    dplyr::rename(ref_coordinates = ref_junID,
                  clinvar = clinvar_type) %>%
    distinct(ref_coordinates, .keep_all = T) %>%
    tibble::rowid_to_column("ref_junID")
  
  any(df_all_introns_tidy$gene_id %>% is.na())
  
  DBI::dbAppendTable(conn = con,
                     name = "intron", 
                     value = df_all_introns_tidy)
  
  DBI::dbDisconnect(conn = con)
}


###################################
## CREATE NOVEL JUNCTION TABLE
##################################

create_novel_table <- function() {
  
  ## Establish connection to the DB
  con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
  DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=1")
  
  # dbRemoveTable(conn = con, "novel")
  
  ## GET GENES
  query = paste0("SELECT id, gene_id FROM 'gene'")
  all_genes <- dbGetQuery(con, query)
  
  
  ## GET REFERENCE INTRONS
  query = paste0("SELECT ref_junID, ref_coordinates FROM 'intron'")
  df_all_introns <- dbGetQuery(con, query) 

  
  # GET MASTER INFO
  query = paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
  
  
  ###################################
  ## CREATE NOVEL TABLE 
  ###################################
  
  SRA_projects <- (df_metadata$SRA_project %>% unique())
  gtf_version <- "105"
  
  
  df_all_novels <- map_df(SRA_projects, function(db) {
    
    # db <- SRA_projects[1]
    print(paste0(Sys.time(), " --> Working with '", db, "' DataBase..."))
    base_folder <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", db, "/")
    
    clusters <- df_metadata %>%
      filter(SRA_project == db) %>%
      distinct(cluster) %>%
      pull()
    
    df_all_novels <- map_df(clusters, function(cluster) {
      # cluster <- clusters[1]
      
      print(paste0(Sys.time(), " --> ", cluster))
      
      if (file.exists( paste0(base_folder, "results/pipeline3/missplicing-ratio/",
                              cluster, "/v", gtf_version, "/", cluster, "_db_novel.rds"))){
        
        df_novel <- readRDS(file = paste0(base_folder, "results/pipeline3/missplicing-ratio/",
                                          cluster, "/v", gtf_version, "/", cluster, "_db_novel.rds"))
        
        df_novel_tidy <- df_novel %>% 
          dplyr::mutate(novel_type = novel_type %>% as.character(),
                        #gene_id = gene_id %>% as.character(),
                        novel_mean_counts = round((novel_sum_counts / novel_n_individuals), digits = 2)) %>%
          select(ref_junID,
                 novel_junID,
                 novel_ss5score, 
                 novel_ss3score, 
                 novel_type,
                 distance, 
                 gene_id)
        
        return(df_novel_tidy)
        
      } else {
        return(NULL)
      }
    })  
  })
  
  ## DO the data merge with other tables -------------------------------------
  
  ## All novel junctions should have been assigned the same gene regardless the
  ## tissue
  # df_all_novels %>% nrow()
  # df_all_novels %>% group_by(novel_junID) %>% distinct(gene_id) %>% nrow()
  # df_all_novels %>% distinct(novel_junID) %>% nrow()
  # df_all_novels %>% group_by(novel_junID) %>% distinct(gene_id) %>% count(gene_id) %>% filter(n>1)
  
  ## Flatten each gene list element internally
  df_all_novels_tidy <- df_all_novels %>%
    distinct(novel_junID, .keep_all = T) %>%
    #rowwise() %>%
    #filter(gene_id %>% unlist() %>% length() == 1) %>%
    mutate_if(is.list, simplify_all) %>%    
    unnest(gene_id)
  
  
  ## Add the GENE ID for the foreign key
  df_all_novels_tidy <- merge(x = df_all_novels_tidy %>% as.data.table(),
                               y = all_genes %>% as.data.table(),
                               by = "gene_id",
                               all.x = T) %>%
    select(-gene_id) %>%
    dplyr::rename(gene_id = id)
  # df_all_novels_tidy %>% nrow()
  # df_all_novels_tidy %>% as_tibble()
  
  
  
  ## ADD INTRON ID info
  df_all_novels_tidy <- merge(x = df_all_novels_tidy %>% as.data.table(),
                              y = df_all_introns %>% as.data.table(),
                              by.x = "ref_junID",
                              by.y = "ref_coordinates",
                              all.x = T) %>%
    select(-ref_junID) %>%
    dplyr::rename(ref_junID = ref_junID.y)
  
  # df_all_novels_tidy %>% nrow()
  # df_all_novels_tidy %>% as_tibble()
  # 
  # 
  # any(df_all_novels_tidy$ref_junID %>% is.na)
  # any(df_all_novels_tidy$gene_id %>% is.na)
  
  ## Extra QC
  ## Same Novel Junctions should have been assigned same scores across genes
  ## and tissues
  # df_all_novels_tidy %>%
  #   dplyr::group_by(novel_junID, gene_id) %>%
  #   distinct(novel_ss5score) %>% nrow()
  # df_all_novels_tidy %>%
  #   dplyr::group_by(novel_junID, gene_id) %>%
  #   distinct(novel_ss3score) %>% nrow()
  # df_all_novels_tidy %>%
  #   dplyr::group_by(novel_junID) %>%
  #   distinct(novel_ss3score) %>% nrow()
  # df_all_novels_tidy %>%
  #   distinct(novel_junID) %>% nrow()
  
  df_all_novels_tidy <- df_all_novels_tidy %>% as_tibble()
  
  
  # Same novels should have assigned the same intron across tissues
  
  df_ambiguous <- df_all_novels_tidy %>%
    dplyr::group_by(novel_junID) %>%
    distinct(ref_junID) %>%
    dplyr::count() %>% 
    filter(n > 1)
    
  ## If there are ambigous junctions, we remove them
  if (df_ambiguous %>% nrow() > 1) {
    df_all_novels_tidy <- df_all_novels_tidy %>% 
    filter(!(novel_junID %in% df_ambiguous$novel_junID))
  }
  
  
  
  ## CREATE NOVEL JUNCTION TABLE -------------------------------------------

  query <- paste0("CREATE TABLE IF NOT EXISTS 'novel'",
                  "(novel_junID NUMERIC NOT NULL,
                  ref_junID NUMERIC NOT NULL,
                  novel_coordinates TEXT NOT NULL, 
                  novel_ss5score DOUBLE NOT NULL, 
                  novel_ss3score DOUBLE NOT NULL, 
                  novel_type TEXT NOT NULL, 
                  distance INTEGER NOT NULL,
                  PRIMARY KEY (ref_junID, novel_junID),
                  FOREIGN KEY (ref_junID) REFERENCES 'intron'(ref_junID))")
  
  res <- DBI::dbSendQuery(conn = con, statement = query)
  DBI::dbClearResult(res)
  
  # query <- paste0("CREATE UNIQUE INDEX 'index_novel' ON 'novel'(novel_junID)");
  # res <- DBI::dbSendQuery(conn = con, statement = query)
  # DBI::dbClearResult(res)
  
  
  
  ## POPULATE NOVEL JUNCTION TABLE  -----------------------------------------
  df_all_novels_tidy <- df_all_novels_tidy %>% 
    as_tibble() %>%
    dplyr::rename(novel_coordinates = novel_junID ) %>%
    distinct(novel_coordinates, .keep_all = T) %>%
    tibble::rowid_to_column("novel_junID")
  
  
  # any(df_all_novels_tidy$gene_id %>% is.na())
  # 
  # 
  # df_all_novels_tidy %>% distinct(novel_coordinates) %>% nrow()
   
  
  DBI::dbAppendTable(conn = con,
                     name = "novel", 
                     value = df_all_novels_tidy %>%
                       distinct(novel_coordinates, .keep_all = T) %>%
                       select(-gene_id))
  DBI::dbDisconnect(conn = con)
}


###################################
## CREATE AND POPULATE TABLES PER EACH PROJECT
###################################

create_cluster_tables <- function() {
    
  con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
  DBI::dbListTables(conn = con)
  
  DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=1")
  
  gtf_version <- "105"
  
  ## GET FROM MASTER TABLE
  query = paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
  
  ## GET FROM INTRON TABLE
  query = paste0("SELECT * FROM 'intron'")
  df_intron <- dbGetQuery(con, query) 
  
  ## GET FROM NOVEL JUNCTION TABLE
  query = paste0("SELECT * FROM 'novel'")
  df_novel <- dbGetQuery(con, query) 
  
  ## GET FROM GENE TABLE
  query = paste0("SELECT * FROM 'gene'")
  df_gene <- dbGetQuery(con, query)
  
  SRA_projects <- (df_metadata$SRA_project %>% unique())
  
  for (db in SRA_projects) {
    
    # db <- SRA_projects[1]
    
    print(paste0(Sys.time(), " --> Working with '", db, "' DataBase..."))
    base_folder <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount3-projects/", db, "/")
    
    clusters <- df_metadata %>%
      filter(SRA_project == db) %>%
      distinct(cluster) %>%
      pull()
  
    for (cluster in clusters) { 
      
      print(paste0(Sys.time(), " --> ", cluster))
      # cluster <- clusters[11]
      # cluster <- clusters[1]
      # cluster <- clusters[2]
      
      ###############################
      ## CREATE INTRON TABLE
      ###############################
      
      # dbRemoveTable(conn = con, paste0(cluster, "_", db))
      query <- paste0("CREATE TABLE IF NOT EXISTS '", paste0(cluster, "_", db, "_misspliced"), "'", 
                                       "(ref_junID INTEGER NOT NULL,
                                       novel_junID INTEGER NOT NULL,
                                       novel_n_individuals INTEGER NOT NULL, 
                                       novel_mean_counts DOUBLE NOT NULL, 
                                       ref_n_individuals INTEGER NOT NULL, 
                                       ref_mean_counts DOUBLE NOT NULL, 
                                       MSR_D DOUBLE NOT NULL, 
                                       MSR_A DOUBLE NOT NULL, 
                                       ref_type TEXT NOT NULL, 
                                       gene_id INTEGER NOT NULL,
                                       FOREIGN KEY (ref_junID, novel_junID) REFERENCES novel (ref_junID, novel_junID),
                      FOREIGN KEY (gene_id) REFERENCES 'gene'(id))")
      
      res <- DBI::dbSendQuery(conn = con, statement = query)
      DBI::dbClearResult(res)
      
      query <- paste0("CREATE TABLE IF NOT EXISTS '", paste0(cluster, "_", db, "_nevermisspliced"), "'", 
                      "(ref_junID INTEGER NOT NULL,
                                       ref_n_individuals INTEGER NOT NULL, 
                                       ref_mean_counts DOUBLE NOT NULL, 
                                       MSR_D DOUBLE NOT NULL, 
                                       MSR_A DOUBLE NOT NULL, 
                                       ref_type TEXT NOT NULL, 
                                       gene_id INTEGER NOT NULL,
                                       FOREIGN KEY (ref_junID) REFERENCES intron (ref_junID),
                      FOREIGN KEY (gene_id) REFERENCES 'gene'(id))")
      
      res <- DBI::dbSendQuery(conn = con, statement = query)
      DBI::dbClearResult(res)
      
      
      ###############################
      ## INSERT DATA
      ###############################
      
      
      
      if (file.exists(paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
                        cluster, "/v", gtf_version, "/", cluster, "_db_introns.rds")) && 
          file.exists( paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
                              cluster, "/v", gtf_version, "/", cluster, "_db_novel.rds"))) {
        
        
        ## INTRONS -------------------------------------------------------------------------------
        
        df_introns_gr <- readRDS(file = paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
                                            cluster, "/v", gtf_version, "/", cluster, "_db_introns.rds"))
        
        df_novel_gr <- readRDS(file = paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
                                             cluster, "/v", gtf_version, "/", cluster, "_db_novel.rds"))
        
        
        ## Tidy data 
        
        df_introns_tidy <- df_introns_gr %>% 
          dplyr::select(ref_junID, 
                        ref_n_individuals, 
                        ref_sum_counts,
                        MSR_D = ref_missplicing_ratio_tissue_ND,
                        MSR_A = ref_missplicing_ratio_tissue_NA,
                        ref_type) %>%
          dplyr::mutate(ref_mean_counts = round((ref_sum_counts / ref_n_individuals), digits = 2))
        
        df_novel_tidy <- df_novel_gr %>% 
          dplyr::mutate(novel_mean_counts = round((novel_sum_counts / novel_n_individuals), digits = 2)) %>%
          select(ref_junID,
                 novel_junID,
                 novel_n_individuals, 
                 novel_mean_counts)
        
        ## QC
        if ((df_introns_tidy %>% filter(ref_type != "never") %>% distinct(ref_junID) %>% nrow() == 
             df_novel_tidy$ref_junID %>% unique %>% length()) && 
            (df_novel_tidy$ref_junID %>% unique %>% length() ==
             intersect(df_introns_tidy$ref_junID, df_novel_tidy$ref_junID) %>% unique %>% length())) {
          print("OK!")
        }
        
        ## JOIN LOCAL INTRON AND NOVELS
        df_all <- merge(x = df_introns_tidy %>% as.data.table(),
                        y = df_novel_tidy %>% as.data.table(),
                        by = "ref_junID",
                        all.y = T) %>%
          relocate(ref_mean_counts, .after = ref_n_individuals)
        
        
        ## JOIN data with parent NOVEL table
        df_all <- merge(x = df_all %>% as.data.table(),
                        y = df_novel %>% select(novel_junID, novel_coordinates) %>% as.data.table(),
                        by.x = "novel_junID",
                        by.y = "novel_coordinates",
                        all.x = T) %>%
          select(-novel_junID) %>%
          dplyr::rename(novel_junID = novel_junID.y)

        
        ## JOIN data with parent INTRON table
        df_all <- merge(x = df_all %>% as.data.table(),
                         y = df_intron %>% select(ref_junID, ref_coordinates, gene_id) %>% as.data.table(),
                         by.x = "ref_junID",
                         by.y = "ref_coordinates",
                         all.x = T) %>%
          select(-ref_junID) %>%
          dplyr::rename(ref_junID = ref_junID.y)
        
        
        ## PREPARE DATA PRIOR POPULATING THE TABLE
        df_all <- df_all %>%
          relocate(ref_junID, novel_junID) %>%
          select(-ref_sum_counts)
        
        df_all %>% as_tibble()
        
        which(df_all$ref_junID %>% is.na())
       
        # query <- paste0("SELECT * FROM novel WHERE ref_junID IN ('", paste(df_all$ref_junID, collapse = "','"),
        #                 "') AND novel_junID IN ('", paste(df_all$novel_junID, collapse = "','"),"')")
        # 
        # df_novel <- dbGetQuery(con, query) 
        # df_novel %>% nrow()
        # df_all %>% nrow()
        
        ###################################################################
        ## CHECK INTEGRITY WITH PARENT TABLE
        
        master_novel <- df_novel %>%
          filter(novel_junID %in% (df_all %>%
                   filter(ref_type != "never") %>%
                   pull(novel_junID))) %>% 
          select(novel_junID, ref_junID) %>% 
          as.data.table()
        
        diff <- setdiff(df_all %>%
                          filter(ref_type != "never") %>% 
                          arrange(desc(ref_junID)) %>% 
                          select(novel_junID, ref_junID),
                        master_novel %>% 
                          arrange(desc(ref_junID)))
        
        if (diff %>% nrow() > 0) {
          df_all <- df_all %>%
            filter(!(novel_junID %in% diff$novel_junID)) 
          print(diff)
        } else {
          print("good")
        }
        
        
        ###################################################################
        
        DBI::dbAppendTable(conn = con,
                           name = paste0(cluster, "_", db, "_misspliced"), 
                           value = df_all)
        
        
        ###################################################################
        ## NEVER MISSPLICED
        
        df_never <- merge(x = df_introns_tidy %>% as.data.table(),
                          y = df_novel_tidy %>% as.data.table(),
                          by = "ref_junID",
                          all.x = T) %>%
          filter(is.na(novel_junID)) %>%
          relocate(ref_mean_counts, .after = ref_n_individuals) %>%
          select(-novel_junID,-novel_n_individuals,-novel_mean_counts)
        
        df_never <- merge(x = df_never %>% as.data.table(),
                          y = df_intron %>% select(ref_junID, ref_coordinates, gene_id) %>% as.data.table(),
                          by.x = "ref_junID",
                          by.y = "ref_coordinates",
                          all.x = T) %>%
          select(-ref_junID, -ref_sum_counts) %>%
          dplyr::rename(ref_junID = ref_junID.y) %>%
          relocate(ref_junID)
        
        DBI::dbAppendTable(conn = con,
                           name = paste0(cluster, "_", db, "_nevermisspliced"), 
                           value = df_never)
        
        
        print(paste0(Sys.time(), ". '", paste0(cluster, "_", db, "_nevermisspliced"), "' table populated!"))
        
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
  
  #query <- paste0("SELECT * FROM 'Brain - Putamen (basal ganglia)_BRAIN' LIMIT 10")
  #dbGetQuery(con, query)

  
  
             
  
  dbListTables(con)
  DBI::dbDisconnect(conn = con)

}

###################################
## QC
###################################

create_intron_table()
create_novel_table()
create_cluster_tables()

# con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
# query <- paste0("SELECT * FROM 'Brain - Hippocampus_BRAIN_db_intron'")
# all_introns <- dbGetQuery(con, query)
# 
# query <- paste0("SELECT * FROM 'Brain - Hippocampus_BRAIN_db_novel'")
# all_novel <- dbGetQuery(con, query)
# 
# all_novel %>% nrow()
# all_novel %>%
#   filter(ref_junID %in% all_introns$ref_junID) %>% nrow()
# all_introns
