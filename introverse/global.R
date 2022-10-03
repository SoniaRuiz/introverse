library(DBI)
library(tidyverse)


con <- DBI::dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
DBI::dbListTables(con)
