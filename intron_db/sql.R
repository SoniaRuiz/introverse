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
DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=0")
tables <- dbListTables(con)
for (table in tables) {
  #if (table != "master") {
    dbRemoveTable(conn = con, table)
    print(paste0("Table: '", table, "' has been removed!"))
  #}
}
dbListTables(con)



###################################
## CREATE MASTER TABLE
###################################
SRA_projects <- c("SRP058181", "SRP051844")

df_all_projects_metadata <- map_df(SRA_projects, function(project) { 
  
  print(paste0(Sys.time(), " - getting info from ", project))
  # project <- SRA_projects[1]
  # project <- SRA_projects[2]
  
  metadata <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount2-projects/", project, "/raw_data/samples_metadata.rds"))
  
  ## EXTRACT METADATA
  
  df_project_metadata <- map_df(metadata$characteristics %>% as.vector, function(characteristic) {
    
    # characteristic <- (metadata$characteristics %>% as.vector)[[1]]
    print(characteristic)
    
    ## TISSUE
    ind <- str_detect(characteristic, pattern = "tissue")
    if (any(ind))
      donor_tissue <- str_sub(string = characteristic[ind], 
                              start = str_locate(characteristic[ind], pattern = ": ")[[2]] + 1,
                              end = characteristic[ind] %>% str_count())
    else
      donor_tissue <- NA
    
    ## GENDER
    ind <- str_detect(characteristic, pattern = "gender")
    if (any(ind))
      donor_gender <- str_sub(string = characteristic[ind], 
                              start = str_locate(characteristic[ind], pattern = ": ")[[2]] + 1,
                              end = characteristic[ind] %>% str_count())
    else
      donor_gender <- NA
    
    ## RIN
    ind <- str_detect(characteristic, pattern = "rin")
    if (any(ind))
      donor_rin <- str_sub(string = characteristic[ind], 
                           start = str_locate(characteristic[ind], pattern = ": ")[[2]] + 1,
                           end = characteristic[ind] %>% str_count())
    else
      donor_rin <- NA
    
    ## AGE
    ind <- str_detect(characteristic, pattern = "death")
    if (any(ind))
      donor_age <- str_sub(string = characteristic[ind], 
                           start = str_locate(characteristic[ind], pattern = ": ")[[2]] + 1,
                           end = characteristic[ind] %>% str_count())
    else
      donor_age <- NA
    
    ind <- str_detect(characteristic, pattern = "diagnosis")
    if (any(ind))
      donor_diagnosis <- str_sub(string = characteristic[ind], 
                                 start = str_locate(characteristic[ind], pattern = ": ")[[2]] + 1,
                                 end = characteristic[ind] %>% str_count())
    else
      donor_diagnosis <- NA
    
    
    
    
    data.frame(age = donor_age,
               rin = donor_rin,
               gender = donor_gender,
               tissue = donor_tissue,
               diagnosis = donor_diagnosis) %>% 
      return()
    
    
  })
  
  ## CLUSTER
  group <- str_sub(string = metadata$title, 
                   start = 1,
                   end = str_locate(metadata$title, pattern = "_")[[1]]) 
  
  
  df_project_metadata <- df_project_metadata %>%
    mutate(cluster = group,
           avg_read_length = metadata$avg_read_length,
           mapped_read_count = metadata$mapped_read_count,
           SRA_project = metadata$project)
  
  
  if (project == "SRP058181") {
    if (any(df_project_metadata %>%
            filter(SRA_project == project) %>%
            pull(diagnosis) %>% is.na())) {
      df_project_metadata[df_project_metadata$cluster == "C_", "diagnosis"] <- "Neurologically normal"
      df_project_metadata[df_project_metadata$cluster == "P_", "diagnosis"] <- "Parkinson's disease"
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
  }
  
  df_project_metadata %>% 
    return()
  
})

saveRDS(object = df_all_projects_metadata,
        file = "./dependencies/df_all_projects_metadata.rds")


DBI::dbWriteTable(conn = con,
                  name = "master",
                  value = df_all_projects_metadata)




###################################
## CREATE AND POPULATE TABLES
###################################

DBI::dbExecute(conn = con, statement = "PRAGMA foreign_keys=1")

gtf_version <- "105"

# Query to the DB
query = paste0("SELECT * FROM 'master'")
con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
df_all_projects_metadata <- dbGetQuery(con, query) 


SRA_projects <- df_all_projects_metadata$SRA_project %>% unique()

for (db in SRA_projects) {
  
  # db <- "GTEx"
  
  print(paste0(Sys.time(), " --> Working with '", db, "' DataBase..."))
  
  clusters <- df_all_projects_metadata %>%
    filter(SRA_project == db) %>%
    distinct(cluster) %>%
    pull()

  if (db == "GTEx") {
    clusters <-  readRDS(file = "./dependencies/all_tissues_used.rda")[7:17]
    base_folder <- "/home/sruiz/PROJECTS/splicing-project/"
    gtf_version <- "104"
    
  } 
  
  base_folder <- paste0("/home/sruiz/PROJECTS/splicing-project/splicing-recount2-projects/", db, "/")
  
  for (cluster in clusters) { 
    
    print(paste0(Sys.time(), " --> ", cluster))
    # cluster <- clusters[11]
    # cluster <- clusters[1]
    # cluster <- clusters[2]
    
    ###############################
    ## CREATE INTRON TABLE
    ###############################
    
    intron_table_statement <- paste0("CREATE TABLE IF NOT EXISTS '", paste0(cluster, "_", db, "_db_introns"), "'", 
                                     "(ref_junID INTEGER PRIMARY KEY NOT NULL, 
                                     seqnames VARCHAR(2) NOT NULL, start INTEGER NOT NULL, end INTEGER NOT NULL, width INTEGER NOT NULL, strand VARCHAR(1) NOT NULL, 
                                     ref_ss5score DOUBLE NOT NULL, ref_ss3score DOUBLE NOT NULL, 
                                     ref_type VARCHAR(10) NOT NULL, ref_n_individuals INTEGER NOT NULL, ref_mean_counts DOUBLE NOT NULL, 
                                     ref_missplicing_ratio_tissue_ND DOUBLE NOT NULL, ref_missplicing_ratio_tissue_NA DOUBLE NOT NULL, 
                                     clinvar_type TEXT NOT NULL, 
                                     MANE BOOL NOT NULL,
                                     gene_name TEXT)")
    
    res <- DBI::dbSendQuery(conn = con, statement = intron_table_statement)
    DBI::dbClearResult(res)
    
    
    ###############################
    ## CREATE NOVEL TABLE
    ###############################
    
    # Create the novel junctions table
    novel_table_statement <- paste0("CREATE TABLE IF NOT EXISTS '", paste0(cluster, "_", db, "_db_novel"), "'",
                                    "(ref_junID INTEGER NOT NULL, 
                                    novel_junID INTEGER NOT NULL, 
                                    seqnames TEXT NOT NULL, start INTEGER NOT NULL, end INTEGER NOT NULL, width INTEGER NOT NULL, strand TEXT NOT NULL, 
                                    novel_ss5score DOUBLE NOT NULL, novel_ss3score DOUBLE NOT NULL, novel_type TEXT NOT NULL, 
                                    novel_n_individuals INTEGER NOT NULL, novel_mean_counts DOUBLE NOT NULL, 
                                    distance INTEGER NOT NULL,
                                    gene_name TEXT, 
                                    PRIMARY KEY (ref_junID, novel_junID),
                                    FOREIGN KEY (ref_junID) REFERENCES '", paste0(cluster, "_", db, "_db_introns"), "' (ref_junID))")
    
    res <- DBI::dbSendQuery(conn = con, statement = novel_table_statement)
    DBI::dbClearResult(res)
    
    
    
    ###############################
    ## INSERT DATA
    ###############################
    
    ## INTRON TABLE
    
    df_introns <- readRDS(file = paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
                                        cluster, "/v", gtf_version, "/", cluster, "_db_introns.rds"))
    
    if (db == "GTEx") {
      df_introns_tidy <- df_introns %>%
        dplyr::mutate(ref_mean_counts = ref_sum_counts / ref_n_individuals)
    }
    
    df_introns_tidy <- df_introns %>% 
      dplyr::select(-tx_id_junction) %>%
      dplyr::mutate(strand = strand %>% as.character(),
                    seqnames = seqnames %>% as.character(),
                    start = start %>% as.integer(),
                    end = end %>% as.integer(),
                    ref_junID = ref_junID %>% as.integer(),
                    width = width %>% as.integer(),
                    #gene_id = gene_id %>% as.character(),
                    gene_name = gene_name %>% as.character()) %>%
      filter(u2_intron == T) %>%
      distinct(ref_junID, .keep_all = T) %>%
      select(ref_junID,
             seqnames, start, end, width, strand,
             ref_ss5score, ref_ss3score,
             ref_type, ref_n_individuals, ref_mean_counts, 
             ref_missplicing_ratio_tissue_ND, ref_missplicing_ratio_tissue_NA,
             clinvar_type, 
             MANE,
             #protein_coding, 
             gene_name)
    
    
    
    DBI::dbAppendTable(conn = con,
                       name = paste0(cluster, "_", db, "_db_introns"), 
                       value = df_introns_tidy)
    
    
    DBI::dbReadTable(conn = con, name = paste0(cluster, "_", db, "_db_introns"))
    
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
    
    df_novel_gr <- readRDS(file = paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
                                         cluster, "/v", gtf_version, "/", cluster, "_db_novel.rds"))
    
    if (db == "GTEx") {
      df_novel_gr <- df_novel_gr %>%
        dplyr::mutate(novel_mean_counts = novel_sum_counts / novel_n_individuals)
    }
    
    df_novel_tidy <- df_novel_gr %>% 
      dplyr::mutate(ref_junID = ref_junID %>% as.integer(),
                    novel_junID = novel_junID %>% as.integer(),
                    novel_type = novel_type %>% as.character(),
                    strand = strand %>% as.character(),
                    width = width %>% as.integer(),
                    gene_name = gene_name %>% as.character(),
                    #gene_id = gene_id %>% as.character(),
                    seqnames = seqnames %>% as.character(),
                    start = start %>% as.integer(),
                    end = end %>% as.integer()) %>%
      select(ref_junID,
             novel_junID,
             seqnames,start, end, width, strand,
             novel_ss5score, novel_ss3score, novel_type,
             novel_n_individuals, novel_mean_counts,
             distance, 
             #protein_coding, 
             #gene_id, 
             gene_name)
    
    DBI::dbAppendTable(conn = con,
                       name = paste0(cluster, "_", db, "_db_novel"), 
                       value = df_novel_tidy)
    
    
    # ## QC
    # df_novel_tidy_test <- df_novel_tidy[1,]
    # df_novel_tidy_test[, "ref_junID"] <- NULL
    # df_novel_tidy_test[, "novel_junID"] <- NA
    # DBI::dbAppendTable(conn = con,
    #                    name = paste0(cluster, "_db_novel"), 
    #                    value = df_novel_tidy_test)
    # 
    # DBI::dbReadTable(conn = con, name = paste0(cluster, "_db_novel")) %>%
    #   filter(is.na(ref_junID))
    
  }
}

dbListTables(con)

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


## CREATE DB TABLES PER GTEx TISSUE






# dbListTables(con)
# dbRemoveTable(con, "Brain-FrontalCortex_BA9_db_introns")
# dbRemoveTable(con, "Brain-FrontalCortex_BA9_db_novel")
# dbRemoveTable(con, "Brain-Hippocampus_db_introns")
# dbRemoveTable(con, "Brain-Hippocampus_db_novel")
# dbRemoveTable(con, "Brain-Substantianigra_db_introns")
# dbRemoveTable(con, "Brain-Substantianigra_db_novel")
# dbListTables(con)

## PD - CONTROL - SRP058181 ------------------------------------------------------

GTEx <- F
clusters <- c("PD", "control")
project_id <- "SRP058181"
base_folder <- "/home/sruiz/PROJECTS/splicing-project/splicing-recount2-projects/SRP058181/"

dbListTables(con)
dbRemoveTable(con, "control_SRP058181_db_introns")
dbRemoveTable(con, "PD_SRP058181_db_introns")
dbRemoveTable(con, "control_SRP058181_db_novel")
dbRemoveTable(con, "PD_SRP058181_db_novel")
dbRemoveTable(con, "control_SRP058181_db_lm")
dbRemoveTable(con, "PD_SRP058181_db_lm")



## PD - CONTROL - SRP049203 -------------------------------------------------------

GTEx <- F
clusters <- c("PD", "control")
project_id <- "SRP049203"
base_folder <- "/home/sruiz/PROJECTS/splicing-project/splicing-recount2-projects/SRP049203/"

dbRemoveTable(con, "control_SRP049203_db_introns")
dbRemoveTable(con, "PD_SRP049203_db_introns")
dbRemoveTable(con, "control_SRP049203_db_novel")
dbRemoveTable(con, "PD_SRP049203_db_novel")






###################################
## ADDING TABLES LOOP 
###################################

for (cluster in clusters) { 
  # cluster <- clusters[1]
  # cluster <- clusters[2]
  # cluster <- clusters[11]
  
  print(paste0(cluster))
  
  #if (str_detect(string = cluster %>% tolower(), pattern = "brain")) {
  
  df_introns <- readRDS(file = paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
                                        cluster, "/v104/", cluster, "_db_introns.rds"))
    
  df_introns %>% head
  # df_introns[1,]
  # df_introns %>% 
  #   filter(ref_missplicing_ratio_tissue_ND == 0, ref_missplicing_ratio_tissue_NA == 0)
  
  df_introns_tidy <- df_introns %>% 
    dplyr::select(-tx_id_junction) %>%
    dplyr::mutate(strand = strand %>% as.character(),
                  seqnames = seqnames %>% as.character(),
                  start = start %>% as.integer(),
                  end = end %>% as.integer(),
                  ref_junID = ref_junID %>% as.integer(),
                  width = width %>% as.integer(),
                  #gene_id = gene_id %>% as.character(),
                  gene_name = gene_name %>% as.character()) %>%
    filter(u2_intron == T) %>%
    distinct(ref_junID, .keep_all = T) %>%
    select(ref_mean_counts, 
           ref_missplicing_ratio_tissue_ND,
           ref_missplicing_ratio_tissue_NA,
           seqnames,
           start,
           end,
           strand,
           ref_type,
           ref_junID,
           width,
           ref_ss5score,
           ref_ss3score,
           ref_n_individuals,
           clinvar_type,
           #gene_id,
           gene_name)
  
  sapply(df_introns_tidy, class)
  
  if (GTEx) {
    # df_introns_tidy <- df_introns_tidy %>% 
    #   dplyr::select(-transcript_id_start,-transcript_id_end) 
    dbWriteTable(con, paste0(cluster, "_db_introns"), df_introns_tidy, overwrite = T)
  } else {
    dbWriteTable(con, paste0(cluster, "_", project_id, "_db_introns"), df_introns_tidy, overwrite = T)
  }
  
  print(paste0("Table '", cluster, "_db_introns' added!"))
  
  
  
  
  # df_introns_details <- readRDS(file = paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
  #                                            cluster, "/", cluster, "_db_introns_details.rds"))
  # sapply(df_introns_details, class)
  # dbWriteTable(con, paste0(cluster, "_db_introns_details"), df_introns_details, overwrite = T)
  # print(paste0("Table '", cluster, "_db_introns_details' added!"))
  
  
  

  df_novel_gr <- readRDS(file = paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
                                       cluster, "/v104/", cluster, "_db_novel.rds"))

  
  df_novel_gr %>% head()
  df_novel_gr %>% names()
  df_novel_tidy <- df_novel_gr %>% 
    #dplyr::select(-transcript_id_start,-transcript_id_end) %>%
    #dplyr::select(-tx_id_junction) %>%
    dplyr::mutate(ref_junID = ref_junID %>% as.integer(),
                  novel_junID = novel_junID %>% as.integer(),
                  novel_type = novel_type %>% as.character(),
                  strand = strand %>% as.character(),
                  width = width %>% as.integer(),
                  #gene_id = gene_id %>% as.character(),
                  gene_name = gene_name %>% as.character(),
                  seqnames = seqnames %>% as.character(),
                  start = start %>% as.integer(),
                  end = end %>% as.integer()) %>%
    select(seqnames,
           start,
           end,
           strand,
           novel_type,
           novel_junID,
           ref_junID,
           width,
           novel_ss5score,
           novel_ss3score,
           distance,
           novel_sum_counts,
           novel_n_individuals,
           #gene_id,
           gene_name)
  sapply(df_novel_tidy, class)
  if (GTEx) {
    # df_novel_tidy <- df_novel_tidy %>% 
    #   dplyr::select(-transcript_id_start,-transcript_id_end)
    dbWriteTable(con, paste0(cluster, "_db_novel"), df_novel_tidy, overwrite = T)
  } else {
    dbWriteTable(con, paste0(cluster, "_", project_id, "_db_novel"), df_novel_tidy, overwrite = T)
  }
  
  print(paste0("Table '", cluster, "_db_novel' added!"))
  
  
  # df_novel_details_gr <- readRDS(file = paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
  #                                              cluster, "/", cluster, "_db_novel_details.rds"))
  # 
  # df_novel_details_tidy <- df_novel_details_gr %>% 
  #   dplyr::mutate(novel_junID = novel_junID %>% as.character())
  #   
  # sapply(df_novel_details_tidy, class)
  # dbWriteTable(con, paste0(cluster, "_db_novel_details"), df_novel_details_tidy, overwrite = T)
  # print(paste0("Table '", cluster, "_db_novel_details' added!"))
  
  # dbListTables(con) %>% print()
  #}
  
}



### GENERATE SQL TABLES FOR DISEASE SAMPLES ======================================

for (cluster in clusters) {
  
  # cluster <- clusters[1]
  # cluster <- clusters[2]
  print(cluster)
  
  df_gr <- readRDS(file = paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
                                 cluster, "/v104/", cluster, "_db_introns.rds"))
  
  if (str_detect(string = cluster, pattern = "PD")) {

    ## Add 'disease_state' column to dataset
    df_gr <- df_gr %>%
      mutate(disease_state = "PD") 
 

    df_gr_control <- readRDS(file = paste0(base_folder, "results/pipeline3/missplicing-ratio/control/v104/control_db_introns.rds"))
    
    df_gr_control <- df_gr_control %>%
      mutate(disease_state = "control") 
    
    
    ## rbind the two datasets
    df_gr <- rbind(df_gr_control,
                   df_gr)
    
    # df_gr <- df_gr %>%
    #   mutate(disease_state = disease_state %>% as.factor()) %>%
    #   mutate(disease_state = relevel(disease_state, ref = "control"))
    
    df_stats <- df_gr %>%
      dplyr::rename(intron_length = width,
                    intron_5ss_score = ref_ss5score,
                    intron_3ss_score = ref_ss3score,
                    gene_length = gene_width,
                    gene_tpm = tpm,
                    gene_num_transcripts = transcript_number
      )%>%
      select(-transcript_id_start,-transcript_id_end)
    
    
    # ## Store "df_stats" within the database
    sapply(df_stats, class)
    con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
    dbWriteTable(con, paste0(cluster, "_", project_id, "_db_lm"), df_stats, overwrite = T)
    #dbDisconnect(con)
    
  } else {


    ## Add 'disease_state' column to dataset
    df_gr <- df_gr %>%
      mutate(disease_state = "control")
    
    
    df_gr_PD <- readRDS(file = paste0(base_folder, "results/pipeline3/missplicing-ratio/PD/v104/PD_db_introns.rds"))
    
    df_gr_PD <- df_gr_PD %>%
      mutate(disease_state = "PD")
    
    
    
    ## rbind the two datasets
    df_gr <- rbind(df_gr,
                   df_gr_PD)
    
    # df_gr <- df_gr %>%
    #   mutate(disease_state = disease_state %>% as.factor()) %>%
    #   mutate(disease_state = relevel(disease_state, ref = "PD"))
    
    ## RENAME COLUMNS
    df_stats <- df_gr %>%
      dplyr::rename(intron_length = width,
                    intron_5ss_score = ref_ss5score,
                    intron_3ss_score = ref_ss3score,
                    #disease_stateControl = disease_state,
                    gene_tpm = tpm,
                    gene_length = gene_width,
                    gene_num_transcripts = transcript_number
      )%>%
      select(-transcript_id_start,-transcript_id_end)
    
    ## Store "df_stats" within the database
    sapply(df_stats, class)
    con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
    dbWriteTable(con, paste0(cluster, "_", project_id, "_db_lm"), df_stats, overwrite = T)
    #dbDisconnect(con)


# # Query to the DB
# query = paste0("SELECT * FROM '", paste0(tissue, "_db_lm"),"'")
# 
# 
# con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
# df_stats <- dbGetQuery(con, query) 
# dbDisconnect(con)
# 
# fit_donor <- lm(ref_missplicing_ratio_tissue_ND ~ intron_length + 
#                   intron_5ss_score * intron_3ss_score +
#                   disease_stateControl +
#                 #gene_tpm + 
#                 gene_length + 
#                 gene_num_transcripts, 
#                 #u2_intron ,
#                 data = df_stats)
# 
# fit_acceptor <- lm(ref_missplicing_ratio_tissue_NA ~ intron_length + 
#                      intron_5ss_score * intron_3ss_score +
#                      disease_stateControl +
#                    #gene_tpm + 
#                    gene_length + 
#                    gene_num_transcripts,
#                    #protein_coding,
#                    data = df_stats)
 }

}
dbListTables(con)
dbDisconnect(con)
