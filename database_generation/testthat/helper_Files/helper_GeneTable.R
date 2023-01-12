print(paste0(Sys.time(), " - loading reference transcriptome"))
hg38 <- rtracklayer::import(reference_transcriptome_path)

## Load the necessary tables from the database
print(paste0(Sys.time(), " - connecting to the database and getting info from 'gene' table"))
con <- DBI::dbConnect(RSQLite::SQLite(), database_path)
df_gene <- dplyr::tbl(con, "gene") %>% dplyr::collect()
DBI::dbDisconnect(con)
