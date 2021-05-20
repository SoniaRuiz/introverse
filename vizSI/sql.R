library(DBI)

# Create an ephemeral in-memory RSQLite database
# setwd("/home/sruiz/PROJECTS/splicing-project-app/vizSI/")
con <- dbConnect(RSQLite::SQLite(), "/home/sruiz/PROJECTS/splicing-project-app/vizSI/dependencies/splicing.sqlite")

# Clear the result
# dbClearResult(res)
# 
# # Disconnect from the database
# dbDisconnect(con)

dbListTables(con)




# 'Both' datatable
both_misspliced_junc_gr <- readRDS(file = "./dependencies/all_tissues/alltissues_both_junctions.rds") %>% as.data.table()
both_misspliced_junc_dt <- both_misspliced_junc_gr %>%
  as.data.frame() %>%
  mutate(seqnames = seqnames %>% as.character(),
         novel_junID = novel_junID %>% as.character(),
         strand = strand %>% as.character(),
         gene_name_start = gene_name_start %>% as.character(),
         novel_type = novel_type %>% as.character())


sapply(both_misspliced_junc_dt,class)
dbWriteTable(con, "both", both_misspliced_junc_dt)
dbGetQuery(con, "SELECT start FROM both")

# dbRemoveTable(con, "both")

# 'Donor' datatable
donor_misspliced_junc_gr <- readRDS(file = "./dependencies/all_tissues/alltissues_donor_junctions.rds") %>% as.data.table()
donor_misspliced_junc_dt <- donor_misspliced_junc_gr %>%
  as.data.frame() %>%
  mutate(seqnames = seqnames %>% as.character(),
         novel_junID = novel_junID %>% as.character(),
         strand = strand %>% as.character(),
         gene_name_start = gene_name_start %>% as.character(),
         novel_type = novel_type %>% as.character())


sapply(donor_misspliced_junc_dt,class)
dbWriteTable(con, "donor", donor_misspliced_junc_dt)
dbGetQuery(con, "SELECT start FROM donor")

# dbRemoveTable(con, "donor")

# 'Acceptor' datatable
acceptor_misspliced_junc_gr <- readRDS(file = "./dependencies/all_tissues/alltissues_acceptor_junctions.rds") %>% as.data.table()
acceptor_misspliced_junc_dt <- acceptor_misspliced_junc_gr %>%
  as.data.frame() %>%
  mutate(seqnames = seqnames %>% as.character(),
         novel_junID = novel_junID %>% as.character(),
         strand = strand %>% as.character(),
         gene_name_start = gene_name_start %>% as.character(),
         novel_type = novel_type %>% as.character())


sapply(acceptor_misspliced_junc_dt,class)
dbWriteTable(con, "acceptor", acceptor_misspliced_junc_dt)
dbGetQuery(con, "SELECT start FROM acceptor")

# dbRemoveTable(con, "acceptor")


# 'Never' datatable
never_misspliced_junctions_gr <- readRDS(file = "./dependencies/all_tissues/alltissues_never_junctions.rds") %>% as.data.table()
never_misspliced_junctions_dt <- never_misspliced_junctions_gr %>%
  as.data.frame() %>%
  mutate(seqnames = seqnames %>% as.character(),
         novel_junID = novel_junID %>% as.character(),
         strand = strand %>% as.character(),
         gene_name_start = gene_name_start %>% as.character(),
         novel_type = novel_type %>% as.character())


sapply(never_misspliced_junctions_dt,class)
dbWriteTable(con, "never", never_misspliced_junctions_dt)
dbGetQuery(con, "SELECT start FROM never")

# dbRemoveTable(con, "never")

dbListTables(con)


for (tissue in gtex_tissues) {
  
  dbRemoveTable(con, tissue)
  
}

dbListTables(con)
dbDisconnect(con)
################################################

gtex_tissues <-  readRDS(file = "/home/sruiz/PROJECTS/splicing-project/results/base_data/all_tissues_used.rda")
con <- dbConnect(RSQLite::SQLite(), "/home/sruiz/PROJECTS/splicing-project-app/vizSI/dependencies/splicing.sqlite")
dbListTables(con)
for (tissue in gtex_tissues) { # tissue <- gtex_tissues[1]
  
  df_introns <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/", tissue, "/", tissue, "_db_introns.rds"))
  df_introns_tidy <- df_introns %>% 
    dplyr::mutate(strand = strand %>% as.character(),
                  seqnames = seqnames %>% as.character(),
                  gene_id = gene_id %>% as.character(),
                  gene_name = gene_name %>% as.character())
  sapply(df_introns_tidy,class)
  dbWriteTable(con, paste0(tissue, "_db_introns"), df_introns_tidy, overwrite = T)
  
  
  df_introns_details <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/", tissue, "/", tissue, "_db_introns_details.rds"))
  df_introns_details_tidy <- df_introns_details %>%
    dplyr::rename(ref_junID = juncID)
  sapply(df_introns_details_tidy, class)
  dbWriteTable(con, paste0(tissue, "_db_introns_details"), df_introns_details_tidy, overwrite = T)
  
  
  df_novel_gr <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/", tissue, "/", tissue, "_db_novel.rds"))
  df_novel_tidy <- df_novel_gr %>% 
    dplyr::mutate(novel_junID = novel_junID %>% as.character(),
                  novel_type = novel_type %>% as.character(),
                  strand = strand %>% as.character(),
                  seqnames = seqnames %>% as.character(),
                  gene_id = gene_id %>% as.character(),
                  gene_name = gene_name %>% as.character())
  sapply(df_novel_tidy, class)
  dbWriteTable(con, paste0(tissue, "_db_novel"), df_novel_tidy, overwrite = T)
  
  
  df_novel_details_gr <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/", tissue, "/", tissue, "_db_novel_details.rds"))
  df_novel_details_tidy <- df_novel_details_gr %>% 
    dplyr::mutate(novel_junID = novel_junID %>% as.character(),
                  novel_type = novel_type %>% as.character())
  sapply(df_novel_details_tidy, class)
  dbWriteTable(con, paste0(tissue, "_db_novel_details"), df_novel_details_tidy, overwrite = T)
  
  
  dbListTables(con)
}

dbListTables(con)
dbDisconnect(con)
