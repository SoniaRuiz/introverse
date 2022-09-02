################################
## CONNECT TO THE DATABASE
################################

con <- dbConnect(RSQLite::SQLite(), "./database/splicing.sqlite")
dbListTables(con)
query <- paste0("SELECT * FROM 'master'")
df_metadata <- dbGetQuery(con, query)


###############################
## GET DATA FOR FRONTAL CORTEX
###############################

project_id <- "BRAIN"
cluster_id <- "Brain - Frontal Cortex (BA9)"

query <- paste0("SELECT ref_junID, MSR_D,  MSR_A, ref_type FROM '", cluster_id, "_", project_id, "_nevermisspliced'")

query <- paste0("SELECT * FROM '", cluster_id, "_", project_id, "_misspliced'")
introns <- dbGetQuery(con, query) %>% as_tibble()


introns %>% filter(ref_junID == "17")

query <- paste0("SELECT * FROM 'novel' WHERE ref_junID = 17 AND novel_junID IN (", paste(introns %>% filter(ref_junID == "17") %>% pull(novel_junID), collapse = ","),")")
novel <- dbGetQuery(con, query) %>% as_tibble()

df_merged <- merge(x = introns,
                   y = novel,
                   by = "novel_junID")
