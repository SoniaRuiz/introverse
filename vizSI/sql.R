library(DBI)

# Create an ephemeral in-memory RSQLite database
setwd("/home/sruiz/PROJECTS/splicing-project-app/vizSI/")
con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")

# Clear the result
# dbClearResult(res)
# 
# # Disconnect from the database
# dbDisconnect(con)

dbListTables(con)



# ########################################################
# ########   ALL TISSUES  ################################
# ########################################################
# 
# 
# # 'Both' datatable
# both_misspliced_junc_gr <- readRDS(file = "./dependencies/all_tissues/alltissues_both_junctions.rds") %>% as.data.table()
# both_misspliced_junc_dt <- both_misspliced_junc_gr %>%
#   as.data.frame() %>%
#   mutate(seqnames = seqnames %>% as.character(),
#          novel_junID = novel_junID %>% as.character(),
#          strand = strand %>% as.character(),
#          gene_name_start = gene_name_start %>% as.character(),
#          novel_type = novel_type %>% as.character())
# 
# 
# sapply(both_misspliced_junc_dt,class)
# dbWriteTable(con, "both", both_misspliced_junc_dt)
# dbGetQuery(con, "SELECT start FROM both")
# 
# # dbRemoveTable(con, "both")
# 
# # 'Donor' datatable
# donor_misspliced_junc_gr <- readRDS(file = "./dependencies/all_tissues/alltissues_donor_junctions.rds") %>% as.data.table()
# donor_misspliced_junc_dt <- donor_misspliced_junc_gr %>%
#   as.data.frame() %>%
#   mutate(seqnames = seqnames %>% as.character(),
#          novel_junID = novel_junID %>% as.character(),
#          strand = strand %>% as.character(),
#          gene_name_start = gene_name_start %>% as.character(),
#          novel_type = novel_type %>% as.character())
# 
# 
# sapply(donor_misspliced_junc_dt,class)
# dbWriteTable(con, "donor", donor_misspliced_junc_dt)
# dbGetQuery(con, "SELECT start FROM donor")
# 
# # dbRemoveTable(con, "donor")
# 
# # 'Acceptor' datatable
# acceptor_misspliced_junc_gr <- readRDS(file = "./dependencies/all_tissues/alltissues_acceptor_junctions.rds") %>% as.data.table()
# acceptor_misspliced_junc_dt <- acceptor_misspliced_junc_gr %>%
#   as.data.frame() %>%
#   mutate(seqnames = seqnames %>% as.character(),
#          novel_junID = novel_junID %>% as.character(),
#          strand = strand %>% as.character(),
#          gene_name_start = gene_name_start %>% as.character(),
#          novel_type = novel_type %>% as.character())
# 
# 
# sapply(acceptor_misspliced_junc_dt,class)
# dbWriteTable(con, "acceptor", acceptor_misspliced_junc_dt)
# dbGetQuery(con, "SELECT start FROM acceptor")
# 
# # dbRemoveTable(con, "acceptor")
# 
# 
# # 'Never' datatable
# never_misspliced_junctions_gr <- readRDS(file = "./dependencies/all_tissues/alltissues_never_junctions.rds") %>% as.data.table()
# never_misspliced_junctions_dt <- never_misspliced_junctions_gr %>%
#   as.data.frame() %>%
#   mutate(seqnames = seqnames %>% as.character(),
#          novel_junID = novel_junID %>% as.character(),
#          strand = strand %>% as.character(),
#          gene_name_start = gene_name_start %>% as.character(),
#          novel_type = novel_type %>% as.character())
# 
# 
# sapply(never_misspliced_junctions_dt,class)
# dbWriteTable(con, "never", never_misspliced_junctions_dt)
# dbGetQuery(con, "SELECT start FROM never")
# 
# # dbRemoveTable(con, "never")
# 
# dbListTables(con)
# 
# 
# for (tissue in gtex_tissues) {
#   
#   dbRemoveTable(con, tissue)
#   
# }
# 
# dbListTables(con)
# dbDisconnect(con)

# dbRemoveTable(con,paste0(cluster, "_SRP049203_db_novel"))
# dbRemoveTable(con,"PD_db_introns_details")
# dbRemoveTable(con,"control_db_introns_details")





########################################################

library(DBI)

setwd("/home/sruiz/PROJECTS/splicing-project-app/vizSI/")
con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")

tables <- dbListTables(con)
for (table in tables) {
  dbRemoveTable(conn = con, table)
}
dbListTables(con)


## GTEX  -------------------------------------------------------------

GTEx <- T
gtex_tissues <-  readRDS(file = "./dependencies/all_tissues_used.rda")
clusters <- gtex_tissues
base_folder <- "/home/sruiz/PROJECTS/splicing-project/"

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


## ADDING TABLES LOOP -------------------------------------------------
for (cluster in clusters) { 
  # cluster <- clusters[1]
  # cluster <- clusters[2]
  
  print(paste0(cluster))
  
  #if (str_detect(string = cluster %>% tolower(), pattern = "brain")) {
  
  df_introns <- readRDS(file = paste0(base_folder, "results/pipeline3/missplicing-ratio/", 
                                        cluster, "/v104/", cluster, "_db_introns.rds"))
    
  
  df_introns[1,]
  df_introns %>% 
    filter(ref_missplicing_ratio_tissue_ND == 0, ref_missplicing_ratio_tissue_NA == 0)
  
  df_introns_tidy <- df_introns %>% 
    dplyr::select(-tx_id_junction) %>%
    dplyr::mutate(strand = strand %>% as.character(),
                  seqnames = seqnames %>% as.character(),
                  gene_name = gene_name %>% as.character(),
                  gene_id = gene_id %>% as.character()) %>%
    filter(u12_intron == T | u2_intron == T) %>%
    distinct(ref_junID, .keep_all = T)
  
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
  df_novel_tidy <- df_novel_gr %>% 
    #dplyr::select(-transcript_id_start,-transcript_id_end) %>%
    #dplyr::select(-tx_id_junction) %>%
    dplyr::mutate(novel_junID = novel_junID %>% as.character(),
                  novel_type = novel_type %>% as.character(),
                  strand = strand %>% as.character(),
                  gene_name = gene_name %>% as.character(),
                  seqnames = seqnames %>% as.character(),
                  gene_id = gene_id %>% as.character())
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
