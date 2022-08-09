
###################################################
################# FUNCTIONS #######################
###################################################



# # intron=6028640
# # coordinates=1:155235845-155236244:-
# # gene=GBA
# # type=acceptor
# # clinvar=-
# # length=400
# # tissue="Adipose-Subcutaneous"
# get_intron_details <- function(intron_id = NULL,
#                                tissue = NULL) {
# 
#   query = paste0("SELECT * FROM '", paste0(tissue, "_db_intron_details"), "' WHERE ref_junID == '", intron_id, "'")
#   
#   
#   con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
#   df_gr <- dbGetQuery(con, query) 
#   dbDisconnect(con)
#   
#   
#   if (df_gr %>% nrow() > 0) {
#     
#     df_gr %>%
#       as.data.frame() %>%
#       #dplyr::mutate(count = count %>% integer()) %>%
#       dplyr::filter(!is.na(count), count > 0) %>%
#       select(Intron_ID = ref_junID,
#              Sample_ID = sample,
#              Counts = count,
#              Age = age,
#              Sex = sex,
#              Read_Count = mapped_read_count) %>% 
#       return()
#     
#   } else {
#     
#     data.frame(Message = paste0("Intron not found in '", 
#                                 names(tissue_GTEx_choices_alphabetical)[tissue_GTEx_choices_alphabetical == tissue], 
#                                 "' tissue data.")) %>% return()
#   }
#   
#   # } 
#   
# }

# # novel_id = "67197814"
# # tissue <- "control"
# get_novel_details <- function(novel_id = NULL,
#                               db = NULL,
#                               cluster = NULL) {
#   
#   
#     
#     # if (intronType == "never") {
#     # 
#     #   data.frame(Message = paste0("No novel junctions for intron '", intronID,"' in '",
#     #                               names(tissue_GTEx_choices_alphabetical)[tissue_GTEx_choices_alphabetical == tissue], "' tissue data.")) %>% return()
#     # 
#     # } else {
#       
#       # ## Load core shared junctions across tissues files  
#       query = paste0("SELECT * FROM '", paste0(tissue, "_db_novel_details"), "' WHERE novel_junID == '", novel_id, "'")
#       
#      
#       con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
#       df_gr <- dbGetQuery(con, query) 
#       dbDisconnect(con)
#       
#       
#       if (df_gr %>% nrow() > 0) {
#         
#         df_gr %>%
#           as.data.frame() %>%
#           #mutate(coordinates = paste0(seqnames,":",start,"-",end,":",strand)) %>%
#           dplyr::filter(!is.na(novel_counts), novel_counts > 0) %>%
#           select(Novel_ID = novel_junID,
#                  Sample_ID = sample,
#                  Counts = novel_counts,
#                  Age = age,
#                  Sex = sex,
#                  Read_Count = mapped_read_count) %>% 
#           return()
#         
#       } else {
#         
#         data.frame(Message = paste0("Intron not found in '", 
#                                     names(tissue_GTEx_choices_alphabetical)[tissue_GTEx_choices_alphabetical == tissue], 
#                                     "' tissue data.")) %>% return()
#       }
#       
#     # } 
#   
# }


# #INTRON SEARCH
# type ="introns"
# chr ="19"
# start = 44905842
# end = 44906601
# strand = "+"
# gene = "ENSG00000171862"
# threshold=-1#70
# search_type= "radio_bygene_tab1"
# data_bases= "BRAIN"
# clusters= "Brain - Hippocampus"
# mane= TRUE
# clinvar= T
# all_data_bases=F
# 
# 
# search_type= "radio_bygenelist_tab1"
# gene = "NULL"
# df <- data.frame(genes = c("SNCA", "MAPT", "PTEN", "APOE"))
# 
# # NOVEL SEARCH
# 
# type ="novel"
# chr ="19"
# start = 44906263
# end = 44906601
# strand = "+"
# gene = "ENSG00000171862"
# threshold= 1
# search_type= "radio_bygene_tab2"
# data_bases= "BLOOD"
# clusters= "Whole Blood"
# mane= FALSE
# clinvar= T



main_IDB_search <- function(type,
                            chr = NULL,
                            start = NULL,
                            end = NULL,
                            strand = NULL,
                            gene = NULL,
                            gene_file = NULL,
                            threshold,
                            search_type,
                            all_data_bases,
                            data_bases,
                            clusters,
                            mane, 
                            clinvar) {

  print(type)
  print(chr)
  print(start)
  print(end)
  print(strand)
  print(gene)
  print(threshold)
  print(search_type)
  print(all_data_bases)
  print(data_bases)
  print(clusters)
  print(mane)
  print(clinvar)
  print("##########################")

  do_next <- F
  # setwd("/home/sruiz/PROJECTS/splicing-project-app/intron_db/")
  con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
  
  # Query to the DB
  query = paste0("SELECT * FROM 'master'")
  df_all_projects_metadata <- dbGetQuery(con, query) 
  #DBI::dbDisconnect(conn = con)

  if (all_data_bases) {
    data_bases <- df_all_projects_metadata$SRA_project %>% unique()
    clusters <- df_all_projects_metadata %>%
      distinct(cluster) %>%
      pull()
  }
  # data_bases <- df_all_projects_metadata$SRA_project %>% unique()
  
  df_gr <- map_df(data_bases, function(db_IDB) {

    # db_IDB <- data_bases[1]
    # print(db_IDB)
    
    details <- df_all_projects_metadata %>%
      filter(SRA_project == db_IDB)
        
    
    map_df(clusters,  function(clust) {

      # clust <- clusters[1]
 
      if (any(details$cluster == clust)) {
        
        ## Get the metadata for the current cluster
        
        cluster_tidy <- details %>%
          filter(cluster == clust) %>% 
          distinct(cluster_tidy)%>%
          pull()
    
        db_tidy <- df_all_projects_metadata %>%
          filter(SRA_project == db_IDB) %>%
          distinct(SRA_project_tidy) %>%
          pull()
        
        all_samples <- details %>% 
          filter(cluster == clust) %>% 
          nrow()
        
        
      } else {
        do_next = T
      }
     
      #print(all_samples)
      
      if (!do_next) {
       
        ## Start building the query
        if (str_detect(search_type, pattern = "bycoordinates")) {
          
          ID <- paste0("chr", chr, ":", start, "-", end, ":", strand)
          
          if (type == "introns") {
            query <- paste0("SELECT distinct(intron.ref_coordinates), gene.gene_name,
            intron.ref_ss5score, intron.ref_ss3score, intron.clinvar, 
            intron.ref_cons5score, intron.ref_cons3score, intron.ref_CDTS5score, intron.ref_CDTS3score,
            tissue.ref_type, tissue.ref_n_individuals, tissue.ref_mean_counts, tissue.MSR_D, tissue.MSR_A,
            tissue.ref_junID
                            FROM '", clust, "_", db_IDB, "_misspliced' AS tissue
                            INNER JOIN 'intron' ON intron.ref_junID=tissue.ref_junID
                            INNER JOIN 'gene' ON gene.id=intron.gene_id
                            WHERE intron.ref_coordinates == '", ID, "'")
            query_never <- paste0("SELECT intron.ref_coordinates, gene.gene_name,
            intron.ref_ss5score, intron.ref_ss3score, intron.clinvar,
            intron.ref_cons5score, intron.ref_cons3score, intron.ref_CDTS5score, intron.ref_CDTS3score,
            tissue.ref_type, tissue.ref_n_individuals, tissue.ref_mean_counts
                            FROM '", clust, "_", db_IDB, "_nevermisspliced' AS tissue
                            INNER JOIN 'intron' ON intron.ref_junID=tissue.ref_junID
                            INNER JOIN 'gene' ON gene.id=intron.gene_id
                            WHERE intron.ref_coordinates == '", ID, "'")
          } else {
            query <- paste0("SELECT * 
                            FROM '", clust, "_", db_IDB, "_misspliced' AS tissue
                            INNER JOIN intron ON intron.ref_junID=tissue.ref_junID
                            INNER JOIN novel ON novel.novel_junID=tissue.novel_junID
                            WHERE novel.novel_coordinates == '", ID, "'")
          }
          
        } else if (str_detect(search_type, pattern = "bygene_")) {
          
          if (str_detect(string = gene, pattern = "ENSG")) {
            gene_query <- paste0("gene.gene_id == '", gene, "'")
          } else {
            gene_query <- paste0("gene.gene_name == '", gene, "'")
          }
          
          
          if (type == "introns") {
            query <- paste0("SELECT distinct(intron.ref_coordinates), gene.gene_name,
            intron.ref_ss5score, intron.ref_ss3score, intron.clinvar, 
            intron.ref_cons5score, intron.ref_cons3score, intron.ref_CDTS5score, intron.ref_CDTS3score,
            tissue.ref_type, tissue.ref_n_individuals, tissue.ref_mean_counts, tissue.MSR_D, tissue.MSR_A,
            tissue.ref_junID
                            FROM '", clust, "_", db_IDB, "_misspliced' AS tissue
                            INNER JOIN 'intron' ON intron.ref_junID=tissue.ref_junID
                            INNER JOIN 'gene' ON gene.id=intron.gene_id
                            WHERE ", gene_query, "")
            query_never <- paste0("SELECT intron.ref_coordinates, gene.gene_name,
            intron.ref_ss5score, intron.ref_ss3score, intron.clinvar,
            intron.ref_cons5score, intron.ref_cons3score, intron.ref_CDTS5score, intron.ref_CDTS3score,
            tissue.ref_type, tissue.ref_n_individuals, tissue.ref_mean_counts
                            FROM '", clust, "_", db_IDB, "_nevermisspliced' AS tissue
                            INNER JOIN 'intron' ON intron.ref_junID=tissue.ref_junID
                            INNER JOIN 'gene' ON gene.id=intron.gene_id
                            WHERE ", gene_query, "")
            
          } else {
            # query <- paste0("SELECT *  FROM  '", clust, "_", db_IDB, "_misspliced' AS tissue 
            #                 INNER JOIN 'novel' ON novel.novel_junID=tissue.novel_junID
            #                 LIMIT 10")
            # dbGetQuery(con, query) 
            
            query <- paste0("SELECT 
            tissue.novel_n_individuals, tissue.novel_mean_counts,
            novel.novel_ss5score, novel.novel_ss3score, novel.novel_coordinates, novel.novel_type, novel.distance,
            intron.clinvar, gene.gene_name
                            FROM '", clust, "_", db_IDB, "_misspliced' AS tissue
                            INNER JOIN 'intron' ON intron.ref_junID=tissue.ref_junID
                            INNER JOIN 'novel' ON novel.novel_junID=tissue.novel_junID
                            INNER JOIN 'gene' ON gene.id=intron.gene_id
                            WHERE ", gene_query, "")
          }
        } else if (str_detect(search_type, pattern = "bygenelist")) {
          
          tryCatch(
            {
              df <- read.csv(gene_file$datapath)
            },
            error = function(e) {
              # return a safeError if a parsing error occurs
              stop(safeError(e))
            }
          )
          
          #print(paste(df[,1]))
        
          #if (str_detect(string = gene, pattern = "ENSG")) {
          
          gene_query <- paste0("LOWER(gene.gene_name) IN ('", paste(df[,1] %>% tolower(), collapse = "','" ), "') OR LOWER(gene.gene_id) IN ('", 
                               paste(df[,1] %>% tolower(), collapse = "','" ), "')")
          #showNotification("This is a notification.")
          
          #} else {
          #  gene_query <- paste0("gene.gene_name == '", gene, "'")
          #}
          print(paste0(gene_query))
          query <- paste0("SELECT distinct(intron.ref_coordinates), gene.gene_name,
            intron.ref_ss5score, intron.ref_ss3score, intron.clinvar, 
            intron.ref_cons5score, intron.ref_cons3score, intron.ref_CDTS5score, intron.ref_CDTS3score,
            tissue.ref_type, tissue.ref_n_individuals, tissue.ref_mean_counts, tissue.MSR_D, tissue.MSR_A,
            tissue.ref_junID
                            FROM '", clust, "_", db_IDB, "_misspliced' AS tissue
                            INNER JOIN 'intron' ON intron.ref_junID=tissue.ref_junID
                            INNER JOIN 'gene' ON gene.id=intron.gene_id
                            WHERE (", gene_query, ")")
          query_never <- paste0("SELECT intron.ref_coordinates, gene.gene_name,
            intron.ref_ss5score, intron.ref_ss3score, intron.clinvar,
            intron.ref_cons5score, intron.ref_cons3score, intron.ref_CDTS5score, intron.ref_CDTS3score,
            tissue.ref_type, tissue.ref_n_individuals, tissue.ref_mean_counts
                            FROM '", clust, "_", db_IDB, "_nevermisspliced' AS tissue
                            INNER JOIN 'intron' ON intron.ref_junID=tissue.ref_junID
                            INNER JOIN 'gene' ON gene.id=intron.gene_id
                            WHERE (", gene_query, ")")
      
        }

        
        if (mane) {
          query <- paste0(query, " AND intron.MANE == ", mane)
          query_never <- paste0(query_never, " AND intron.MANE == ", mane)
        }
        
        if (clinvar) {
          query <- paste0(query, " AND intron.clinvar != '-'")
          query_never <- paste0(query_never, " AND intron.clinvar != \"-\"")
        }
        
        if (type == "introns") {
          df_gr <- dbGetQuery(con, query) %>%
            mutate("p_ref_ind" = round(x = (ref_n_individuals * 100) / all_samples))
        } else {
          df_gr <- dbGetQuery(con, query) %>%
            mutate("p_novel_ind" = ifelse(round(x = (novel_n_individuals * 100) / all_samples) == 0, 1, 
                                          round(x = (novel_n_individuals * 100) / all_samples)))
        }
        
        ## PERCENTAGE OF INDIVIDUALS FILTERING
        if (threshold > -1 ) {
          if (type == "introns") {
            
            ## Query the novel junctions attached to the introns and retrieve
            query <- paste0("SELECT *
                            FROM '", clust, "_", db_IDB, "_misspliced' AS tissue
                            WHERE tissue.ref_junID IN ('", paste(df_gr$ref_junID %>% unique(), collapse = "','"), "')") 
                            
            
            df_novel <- dbGetQuery(con, query) %>%
              mutate("p_novel_ind" = ifelse(round(x = (novel_n_individuals * 100) / all_samples) == 0, 1, 
                                            round(x = (novel_n_individuals * 100) / all_samples))) %>%
              filter(p_novel_ind >= threshold)
            
            df_gr <- df_gr %>%
              filter(ref_junID %in% df_novel$ref_junID)
            
          } else {
            df_gr <- df_gr %>%
              filter(p_novel_ind >= threshold)
          }
        } else {
          if (type == "introns") {
            # Get never mis-spliced introns too
            df_nevermispliced <- dbGetQuery(con, query_never) 
            df_gr <- plyr::rbind.fill(df_gr %>%
                                        mutate(MSR_D = formatC(x = MSR_D, format = "e", digits = 3),
                                               MSR_A = formatC(x = MSR_A, format = "e", digits = 3)), df_nevermispliced)
            #df_gr[is.na(df_gr)] <- ""
          }
        }
          
        
        if (df_gr %>% nrow() >= 1) {

          df_gr %>%
            mutate(DB = db_tidy, 
                   all_samples = all_samples,
                   Cluster = cluster_tidy) %>%
            return()
          
        }
        
        
      } else {
        return(NULL)
      }
    
    })
  })
  
  #DBI::dbDisconnect(con)
  
  
  if (df_gr %>% nrow() == 0) {
    
    data.frame("Message" = "No data found within the IDB using the criteria selected.") %>% 
      return()
    
  } else {
    
    # all_people_tissue <- readRDS(file = "./dependencies/all_people_used_tissue.rda")[[tissue]] %>% length()
    
    if (any(df_gr %>% names() == "ref_width")) {
      df_gr <- df_gr %>%
        dplyr::rename(width = ref_width)
    }
   
    
    # query <- paste0("SELECT * 
    # FROM '", clust, "_", db_IDB, "_misspliced", "' AS tissue
    # INNER JOIN 'intron' ON intron.ref_junID=tissue.ref_junID
    # INNER JOIN 'novel' ON novel.novel_junID=tissue.novel_junID
    # WHERE intron.gene_id == ", gene)
    # df <- dbGetQuery(con, query) 
    # 

    
    if (type == "introns") {
      
      indx <- str_locate_all(string = df_gr$ref_coordinates, pattern = ":")
      start_positions <- NULL
      end_positions <- NULL
      
      for (i in c(1:(indx %>% length()))) {
       # print(i)
        start_j <- df_gr$ref_coordinates[i] %>%
          str_sub(start = str_locate_all(string = df_gr$ref_coordinates[i], pattern = ":")[[1]][1,2]+1,
                  end = str_locate_all(string = df_gr$ref_coordinates[i], pattern = "-")[[1]][1,2]-1) %>% as.integer()
        end_j <- df_gr$ref_coordinates[i] %>%
          str_sub(start = str_locate_all(string = df_gr$ref_coordinates[i], pattern = "-")[[1]][1,2]+1,
                  end = str_locate_all(string = df_gr$ref_coordinates[i], pattern = ":")[[1]][2,2]-1) %>% as.integer()
        
        start_positions <- c(start_positions, start_j)
        end_positions <- c(end_positions, end_j)
      }
      # start_j <- df_gr$ref_coordinates %>%
      #   str_sub(start = str_locate_all(string = df_gr$ref_coordinates, pattern = ":")[[1]][1,2]+1,
      #           end = str_locate_all(string = df_gr$ref_coordinates, pattern = "-")[[1]][1,2]-1) %>% as.integer()
      # end_j <- df_gr$ref_coordinates %>%
      #   str_sub(start = str_locate_all(string = df_gr$ref_coordinates, pattern = "-")[[1]][1,2]+1,
      #           end = str_locate_all(string = df_gr$ref_coordinates, pattern = ":")[[1]][2,2]-1) %>% as.integer()
      
      
      df_gr %>%
        #rename(width = ref_junID) %>%
        mutate(MeanCounts = round(x = ref_mean_counts, digits = 2),
               MSR_D = formatC(x = MSR_D, format = "e", digits = 3),
               MSR_A = formatC(x = MSR_A, format = "e", digits = 3),
               MANE = "T",
               p_ref_ind = round((ref_n_individuals/all_samples)*100),
               p_ref_ind = paste0(p_ref_ind, "% (", ref_n_individuals, "/", all_samples, ")"),
               Width = abs(start_positions-end_positions)+1,
               ref_junID = paste0(ref_junID, "#", Cluster),
               "Cons_5ss" = ref_cons5score %>% round(digits = 2),
               "Cons_3ss" = ref_cons3score %>% round(digits = 2),
               "CDTS_5ss" = ref_CDTS5score %>% round(digits = 2),
               "CDTS_3ss" = ref_CDTS3score %>% round(digits = 2),
               ref_type = ref_type %>% as.factor(),
               gene_name = gene_name %>% as.factor()) %>%
        #,ifelse(MANE == 0, "F", "T"),
               #coordinates = paste0(seqnames,":",start,"-",end,":",strand)) %>%
        #"% Individuals" = round(x = (ref_n_individuals * 100) / all_people_tissue)) %>%
        dplyr::select(ID = ref_coordinates,
                      #Coordinates = coordinates,
                      "Mis-spliced site" = ref_type,
                      "Intron Length (bp)" = Width, 
                      "MES_5ss" = ref_ss5score,
                      "MES_3ss" = ref_ss3score,
                      "Cons_5ss",
                      "Cons_3ss",
                      "CDTS_5ss",
                      "CDTS_3ss",
                      MSR_D,
                      MSR_A,
                      "Avg. Read Count" = MeanCounts,
                      "% Individuals" = p_ref_ind,
                      #"Total Individuals" = ref_n_individuals,
                      ClinVar = clinvar,
                      #MANE, 
                      Gene = gene_name,
                      Samples = Cluster, 
                      "Body region" = DB,
                      ref_junID) %>% 
        return()
    } else {
      
      ##
      
      
      # start_j <- df_gr$novel_coordinates %>%
      #   str_sub(start = str_locate_all(string = df_gr$novel_coordinates, pattern = ":")[[1]][1,2]+1,
      #           end = str_locate_all(string = df_gr$novel_coordinates, pattern = "-")[[1]][1,2]-1) %>% as.integer()
      # end_j <- df_gr$novel_coordinates %>%
      #   str_sub(start = str_locate_all(string = df_gr$novel_coordinates, pattern = "-")[[1]][1,2]+1,
      #           end = str_locate_all(string = df_gr$novel_coordinates, pattern = ":")[[1]][2,2]-1) %>% as.integer()
      
      indx <- str_locate_all(string = df_gr$ref_coordinates, pattern = ":")
      start_positions <- NULL
      end_positions <- NULL
      
      for (i in c(1:(indx %>% length()))) {
        # print(i)
        start_j <- df_gr$ref_coordinates[i] %>%
          str_sub(start = str_locate_all(string = df_gr$ref_coordinates[i], pattern = ":")[[1]][1,2]+1,
                  end = str_locate_all(string = df_gr$ref_coordinates[i], pattern = "-")[[1]][1,2]-1) %>% as.integer()
        end_j <- df_gr$ref_coordinates[i] %>%
          str_sub(start = str_locate_all(string = df_gr$ref_coordinates[i], pattern = "-")[[1]][1,2]+1,
                  end = str_locate_all(string = df_gr$ref_coordinates[i], pattern = ":")[[1]][2,2]-1) %>% as.integer()
        
        start_positions <- c(start_positions, start_j)
        end_positions <- c(end_positions, end_j)
      }
      
      df_gr %>%
        mutate(novel_type = str_replace_all(string = novel_type, pattern = "_", replacement = " ")) %>%
        mutate("% Individuals" = paste0(p_novel_ind, "% (", novel_n_individuals, "/", all_samples, ")"),
               "Width" = abs(start_positions-end_positions)+1) %>%
        select(ID = novel_coordinates,
               NovelType = novel_type,
               #"Ref.Intron" = ref_coordinates,
               NovelWidth = Width,
               "MES_5ss" = novel_ss5score,
               "MES_3ss" = novel_ss3score, 
               Distance = distance,
               "Avg. Read Count" = novel_mean_counts,
               "% Individuals",
               #"Total Individuals" = novel_n_individuals,
               Gene = gene_name,
               Samples = Cluster,
               "Body region" = DB) %>%
        mutate(Modulo3 = abs(Distance) %% 3) %>%
        relocate(Modulo3, .after = Distance) %>%
        return()
    }
   
    
  }

}


# intron_id <- "139704"
# db <- "BRAIN"
# sample_group <- "Brain - Hippocampus"

get_novel_data_from_intron <- function(intron_id = NULL,
                                       db = NULL,
                                       sample_group = NULL) {
  
  intron_id <- str_split(string = intron_id,pattern = "#")[[1]][1]
  print(paste0("'get_novel_data_from_intron()' function called! ",
               intron_id %>% print(), " - ", 
               db %>% print(), " - ", 
               sample_group %>% print()))
  
  # setwd("/home/sruiz/PROJECTS/splicing-project-app/intron_db/")
  con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
  
  # Query to the DB
  query = paste0("SELECT * FROM 'master'")
  df_all_projects_metadata <- dbGetQuery(con, query) %>%
    filter(SRA_project_tidy == db) 
  
  
  db_IDB <- df_all_projects_metadata %>%
    distinct(SRA_project) %>%
    pull()
  
  all_samples <- df_all_projects_metadata %>%
    filter(cluster == sample_group) %>%
    nrow()
  
  cluster_tidy <-  df_all_projects_metadata %>%
    filter(cluster == sample_group) %>%
    distinct(cluster_tidy) %>% 
    pull()
  
  
  query = paste0("SELECT * 
                 FROM '", paste0(sample_group, "_", db_IDB, "_misspliced"), 
                 "' WHERE ref_junID == '", intron_id, "'")
  df_gr <- dbGetQuery(con, query) 
  
  if (df_gr %>% nrow() > 0) {
    
    query <- paste0("SELECT * 
                   FROM 'novel' 
                   WHERE novel_junID IN ('", paste(df_gr$novel_junID, collapse="','"), "')")
    df_novel <- dbGetQuery(con, query)
    
    df_novel <- merge(x = df_novel,
                      y = df_gr,
                      by = "novel_junID")
    start_j <- df_novel$novel_coordinates %>%
      str_sub(start = str_locate_all(string = df_novel$novel_coordinates, pattern = ":")[[1]][1,2]+1,
              end = str_locate_all(string = df_novel$novel_coordinates, pattern = "-")[[1]][1,2]-1) %>% as.integer()
    end_j <- df_novel$novel_coordinates %>%
      str_sub(start = str_locate_all(string = df_novel$novel_coordinates, pattern = "-")[[1]][1,2]+1,
              end = str_locate_all(string = df_novel$novel_coordinates, pattern = ":")[[1]][2,2]-1) %>% as.integer()
    
    
    df_novel %>%
      as.data.frame() %>%
      #group_by(ref_junID, novel_junID) %>%
      ##distinct(novel_junID, .keep_all = T) %>%
      #ungroup() %>%
      dplyr::mutate(MeanCounts = round(x = novel_mean_counts, digits = 2),
                    p_indi = ifelse(round((novel_n_individuals * 100) / all_samples) == 0, 1,
                                    round((novel_n_individuals * 100) / all_samples)),
                    "% Individuals" = paste0(p_indi, "% (", novel_n_individuals, "/", all_samples, ")"),
                    #Mean_MSR = formatC(x = novel_missplicing_ratio_tissue, format = "e", digits = 3),
                    #coordinates = paste0(seqnames,":",start,"-",end,":",strand),
                    Width = abs(start_j - end_j) + 1,
                    Samples = cluster_tidy,
                    Project = db,
                    novel_type = str_replace(string = novel_type, pattern = "_", replacement = " "),
                    Modulo3 = abs(distance) %% 3,
                    "Frameshift?" = ifelse(Modulo3 == 0, "N", "Y")) %>%
      dplyr::select(NovelID = novel_coordinates,#Coordinates = coordinates,
                    id = novel_junID,
                    "Novel Type" = novel_type,
                    #RefID = ref_junID,#Coordinates = coordinates,
                    Length = Width,
                    "MES_5ss" = novel_ss5score,
                    "MES_3ss" = novel_ss3score,
                    "Distance (bp)" = distance,
                    "Frameshift?",
                    "Avg. Read Count" = MeanCounts,
                    "% Individuals",
                    #"Total Individuals" = novel_n_individuals,
                    Samples,
                    "Body region" = Project) %>% return()
    
  } else {
    
    data.frame(Time = Sys.time(),
               Message = paste0("Intron not found.")) %>% return()
  }
  
}

# novel_id <- "chr4:89725315-89729193:-"
get_novel_data_across_idb <- function(novel_id) {
  
  # Query to the DB
  query = paste0("SELECT * FROM 'master'")
  con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
  df_all_projects_metadata <- dbGetQuery(con, query) 
  
  data_bases <- df_all_projects_metadata %>%
    distinct(SRA_project) %>%
    pull()
  
  
  df_gr <- map_df(data_bases, function(db_IDB) {
    
    # db_IDB <- data_bases[1]
    ## SET VARIABLE REGARDING DATABASE CHOSEN
    
    df_tidy_name <- df_all_projects_metadata %>% 
      filter(SRA_project == db_IDB) %>%
      distinct(SRA_project_tidy) %>%
      pull()
    
    clusters <- df_all_projects_metadata %>% 
      filter(SRA_project == db_IDB) %>%
      distinct(cluster) %>%
      pull()
    
    
    
    
    
    
    map_df(clusters,  function(clustr) {
      # clustr <- clusters[1]
      
      
      cluster_tidy_name <- df_all_projects_metadata %>%
        filter(SRA_project == db_IDB,
               cluster == clustr) %>% 
        distinct(cluster_tidy) %>%
        pull()
      
      all_samples <- df_all_projects_metadata %>% 
        filter(SRA_project == db_IDB,
               cluster == clustr) %>%
        nrow()
      
      # query <- paste0("EXPLAIN QUERY PLAN SELECT * FROM '", paste0(clustr, "_", db_IDB, "_db_novel"), "' WHERE novel_junID == '", novel_id, "'")
      
      query <- paste0("SELECT * FROM '", paste0(clustr, "_", db_IDB, "_db_novel"), "' WHERE novel_junID == '", novel_id, "'")
      df_novel_gr <- dbGetQuery(con, query)
      
      
      if (df_novel_gr %>% nrow() > 0) {
        
        
        query = paste0("SELECT gene_name FROM 'gene_name' WHERE gene_id == ", (df_novel_gr$gene_name %>% unique()))
        gene_n <- dbGetQuery(con, query)[[1]]
        
        df_novel_gr %>%
          #mutate("Coordinates" = paste0(seqnames, ":", start, "-", end, ":", strand)) %>%
          mutate(novel_type = str_replace_all(string = novel_type, pattern = "_", replacement = " ")) %>%
          mutate(novel_mean_counts = round(novel_mean_counts, digits = 2)) %>%
          mutate("% Individuals" = round(novel_n_individuals * 100/all_samples),
                 Samples = cluster_tidy_name,
                 Project = df_tidy_name,
                 gene_name = gene_n) %>%
          select(ID = novel_junID,
                 NovelType = novel_type,
                 "Ref.Intron" = ref_junID,
                 #Coordinates,
                 #"Width" = width,
                 "MES_5ss" = novel_ss5score,
                 "MES_3ss" = novel_ss3score, 
                 Distance = distance,
                 "Avg. Read Count" = novel_mean_counts,
                 "% Individuals" = paste0(`% Individuals`, "% (", novel_n_individuals, " out of ", all_samples, ")"),
                 Gene = gene_name,
                 Samples,
                 "Body region" = Project) %>%
          mutate(Modulo3 = abs(Distance) %% 3) %>%
          relocate(Modulo3, .after = Distance) %>%
          return()
      } else {
        
        return(NULL)
      }
      
      
    })
    
    
  })
  
  DBI::dbDisconnect(conn = con)
  df_gr %>% return()
  
}


# setwd("introverse/")
# novel_id = 2959056
# db <- "BRAIN"
# clust <- "Brain - Cerebellum"

visualise_transcript <- function(novel_id = NULL,
                                 intron_id = NULL,
                                 clust,
                                 db) {
  
  print(novel_id)
  print(intron_id)
  print(clust)
  print(db)
  # intron_id <- "chr10:87894110-87925512:+"
  # db <- "GTEXv8 - BRAIN"
  # clust <- "Brain - Hippocampus"
  
  library(ggplot2)
  library(ggtranscript)
  
  ##############################
  ## GET THE TABLE NAME
  ##############################
  
  query = paste0("SELECT * FROM 'master'")
  con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
  df_all_projects_metadata <- dbGetQuery(con, query) 
  db_master_details <- df_all_projects_metadata %>%
    filter(SRA_project_tidy == db,
           cluster_tidy == clust)
  db_name <- db_master_details$SRA_project %>% unique()
  cluster_name <- db_master_details$cluster %>% unique()
  
  
  ###################################
  ## QUERY THE DATABASE
  ###################################
  
  if (!is.null(intron_id)) {
    sql_statement <- paste0("SELECT * 
                            FROM '", cluster_name, "_", db_name, "_misspliced' AS tissue
                            INNER JOIN 'intron' ON intron.ref_junID=tissue.ref_junID
                            INNER JOIN 'novel' ON novel.novel_junID=tissue.novel_junID
                            INNER JOIN 'gene' ON gene.id=intron.gene_id
                            WHERE tissue.ref_junID == '", intron_id, "'")
  } 
  if (!is.null(novel_id)) {
    sql_statement <- paste0("SELECT * 
                            FROM '", cluster_name, "_", db_name, "_misspliced' AS tissue
                            INNER JOIN 'intron' ON intron.ref_junID=tissue.ref_junID
                            INNER JOIN 'novel' ON novel.novel_junID=tissue.novel_junID
                            INNER JOIN 'gene' ON gene.id=intron.gene_id
                            WHERE tissue.novel_junID == '", novel_id, "'")
  } 
  
  
  print(sql_statement)
  df_intron <- dbGetQuery(con, sql_statement)
  
  
  if (any(df_intron$MANE)) {
    
    #######################################
    ## EXTRACT THE COORDINATES FROM THE ID
    #######################################
    
    novel_junctions <- map_df(df_intron$novel_coordinates, function(junction) {
      # junction <- df_intron$novel_coordinates[1]
      chr_junc <- junction %>%
        str_sub(start = 1,
                end = str_locate_all(string = junction, pattern = ":")[[1]][1,2]-1)
      start_junc <- junction %>%
        str_sub(start = str_locate_all(string = junction, pattern = ":")[[1]][1,2]+1,
                end = str_locate_all(string = junction, pattern = "-")[[1]][1,2]-1)
      end_junc <- junction %>%
        str_sub(start = str_locate_all(string = junction, pattern = "-")[[1]][1,2]+1,
                end = str_locate_all(string = junction, pattern = ":")[[1]][2,2]-1)
      strand_junc <- junction %>%
        str_sub(start = str_locate_all(string = junction, pattern = ":")[[1]][2,2]+1,
                end = junction %>% stringr::str_count())
      
      data.frame(ID = junction,
                 seqnames = chr_junc,
                 start = start_junc %>% as.integer(),
                 end = end_junc %>% as.integer(),
                 strand = strand_junc) %>%
        return()
    })
    
    
    novel_junctions <- merge(x = novel_junctions,
                             y = df_intron %>% 
                               select(novel_coordinates, novel_n_individuals, 
                                      novel_mean_counts, novel_type,
                                      ref_n_individuals),
                             by.x = "ID",
                             by.y = "novel_coordinates")  %>%
      mutate(width = abs(start - end))
      
      
    #############################
    ## GET THE MANE DATA
    #############################
    
    sql_statement <- paste0("SELECT * FROM 'mane'
                            WHERE gene_name == '", df_intron$gene_name %>% unique, "'")
    print(sql_statement)
    df_mane <- dbGetQuery(con, sql_statement)
    
    df_mane_cds <- df_mane %>% dplyr::filter(type == "CDS")
    df_mane_utr <- df_mane %>% dplyr::filter(type == "UTR")
    
    
    #############################
    ## FIND THE REF INTRON 
    ## PRIOR THE ZOOM
    #############################
    
    all_introns <- ggtranscript::to_intron(df_mane %>% filter(type == "exon"))
    
    ref_intron <- rbind(all_introns[which.min(abs(all_introns$start - novel_junctions$start)),],
                        all_introns[which.min(abs(all_introns$end - novel_junctions$end)),]) 
      
    ref_exons <- df_mane %>% arrange(end) %>%
      filter(type == "exon") %>%
      filter(start %in% ref_intron$end | end %in% ref_intron$start) 
      
    
    #############################
    ## PLOT
    #############################
    
    df_mane %>%
      filter(type == "exon") %>% 
      ggplot(aes(
        xstart = start,
        xend = end,
        y = df_mane$transcript_id %>% unique()
      )) +
      #ggtranscript::geom_range() +
      geom_range(
        data = df_mane_cds,
        fill = "purple"
      ) +
      geom_range(
        data = df_mane_utr,
        fill = "purple",
        height = 0.25
      ) +
      ggtranscript::geom_intron(
        data = ggtranscript::to_intron(df_mane %>%
                                         filter(type == "exon"), 
                                       "transcript_name"),
        aes(strand = df_mane$strand %>% unique)
      ) + 
      ggtranscript::geom_junction(
        data = novel_junctions,
        aes(colour = novel_type),
        ncp = 100, 
        junction.y.max = 0.5 
      ) +
      scale_colour_manual(breaks = c("novel_acceptor", "novel_donor"),
                          values = c("#C82803FF", "#23C3E4FF")) +
      #viridis::scale_fill_viridis(discrete = T, option = "F")  +
      theme_light() +
      ggforce::facet_zoom(xlim = c(min(ref_exons$start),
                                   max(ref_exons$end)),
                          zoom.data = zoom) + 
      geom_junction_label_repel(
        data = novel_junctions %>% dplyr::mutate( zoom = TRUE ) ,
        aes(label = paste0("seen in ", novel_n_individuals, " of ", ref_n_individuals, " samples")),
        junction.y.max = 0.5
      ) + 
      theme(axis.text.y = element_text(angle = 90, 
                                       hjust = 0.5,
                                       size = "11"),
            axis.text = element_text(size = "11"),
            axis.title = element_text(size = "11"),
            legend.position = "top",
            legend.text = element_text(size = "11"),
            legend.title = element_text(size = "11")) +
      xlab(paste0("Genomic position (", df_mane$seqnames %>% unique ,")")) + 
      ylab("MANE Transcript") +
      guides(color = guide_legend(title = "Novel event type: ")) +
      labs(title = paste0("Excision of the novel event '", 
                          novel_junctions$ID %>% unique, "'"),
           #caption = "UTR sequences are represented in white and CDS are in purple",
           subtitle = paste0("MANE transcript of ", df_mane$gene_name %>% unique(), " gene (", clust, ")")) %>%
      return()
    
  
  } else {
    ggplot() +
      theme_void() +
      geom_text(aes(0,0,label='The selected novel junction doesn´t belong to a MANE transcript.')) + 
      theme(text = element_text(element_text(size = "14"))) %>%
      return()
  }
  
  
  
}

# setwd("introverse/")
# gene_id = "GBA"
# clust <- "Brain - Hippocampus"

visualise_missplicing <- function(gene_id = "SNCA",
                                  clust = "Brain - Hippocampus") {

  print(gene_id)
  print(clust)
  #print(db)
  # intron_id <- "chr10:87894110-87925512:+"
  # db <- "GTEXv8 - BRAIN"
  # clust <- "Brain - Hippocampus"
  
  library(ggplot2)
  library(ggtranscript)
  
  ## GET THE DETAILS OF THE CLUSTER/PROJECT SELECTED
  query = paste0("SELECT * FROM 'master'")
  con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
  df_all_projects_metadata <- dbGetQuery(con, query) 
  db_master_details <- df_all_projects_metadata %>%
    filter(cluster_tidy == clust)
  db_name <- db_master_details$SRA_project %>% unique()
  cluster_name <- db_master_details$cluster %>% unique()
  
  
  ## GET THE GENE_NAME AND MANE INFO FROM THE INTRON TABLE TO GET THE TRANSCRIPT ID
  sql_statement <- paste0("SELECT * 
                            FROM '", cluster_name, "_", db_name, "_misspliced' AS tissue
                            INNER JOIN 'intron' ON intron.ref_junID=tissue.ref_junID
                            INNER JOIN 'novel' ON novel.novel_junID=tissue.novel_junID
                            INNER JOIN 'gene' ON gene.id=intron.gene_id
                            WHERE gene.gene_name == '", gene_id, "' OR gene.gene_id == '", gene_id, "'")
  
  
  print(sql_statement)
  df_gene_splicing <- dbGetQuery(con, sql_statement)
  
  
  if (any(df_gene_splicing$MANE)) {
    
    
    ref_introns <- map_df(df_gene_splicing$ref_coordinates, function(junction) {
      # junction <- df_gene_splicing$novel_coordinates[1]
      chr_junc <- junction %>%
        str_sub(start = 1,
                end = str_locate_all(string = junction, pattern = ":")[[1]][1,2]-1)
      start_junc <- junction %>%
        str_sub(start = str_locate_all(string = junction, pattern = ":")[[1]][1,2]+1,
                end = str_locate_all(string = junction, pattern = "-")[[1]][1,2]-1)
      end_junc <- junction %>%
        str_sub(start = str_locate_all(string = junction, pattern = "-")[[1]][1,2]+1,
                end = str_locate_all(string = junction, pattern = ":")[[1]][2,2]-1)
      strand_junc <- junction %>%
        str_sub(start = str_locate_all(string = junction, pattern = ":")[[1]][2,2]+1,
                end = junction %>% stringr::str_count())
      
      data.frame(ID = junction,
                 seqnames = chr_junc,
                 start = start_junc %>% as.integer(),
                 end = end_junc %>% as.integer(),
                 strand = strand_junc) %>%
        return()
      
    })
    
    ref_introns_MSR <- merge(x = ref_introns %>% distinct(ID, .keep_all = T),
                             y = df_gene_splicing %>% 
                               dplyr::select(ID = ref_coordinates, MANE, MSR_D, MSR_A, gene_name) %>% 
                               dplyr::distinct(ID, .keep_all = T),
                             by = "ID",
                             all.x = T) 
    
    
    ##################
    ## ADD MANE DATA
    ##################
    
    
    sql_statement <- paste0("SELECT * FROM 'mane' WHERE gene_name == '", ref_introns_MSR$gene_name %>% unique, "'")
    print(sql_statement)
    df_mane <- dbGetQuery(con, sql_statement)
    
    
    #########################
    exons <- df_mane %>% filter(type == "exon")
    # exons_rescaled <- shorten_gaps(
    #   exons, 
    #   to_intron(exons, "transcript_name"), 
    #   group_var = "transcript_name"
    # )
    
    
    width_bars <- abs(ref_introns_MSR$start - ref_introns_MSR$end) %>% min() / 2
    ref_introns_MSRD <- ref_introns_MSR %>%
      mutate(end = start + width_bars) %>%
      filter(MANE == 1) %>%
      select(seqnames , strand, start, end, MSR_D )
    ref_introns_MSRA <- ref_introns_MSR %>%
      mutate(start = end - width_bars) %>%
      filter(MANE == 1) %>%
      select(seqnames , strand, start, end, MSR_A )
  
    
    df_mane_cds <- df_mane %>% dplyr::filter(type == "CDS")
    df_mane_utr <- df_mane %>% dplyr::filter(type == "UTR")
    
    ## Generate the plot
    exons %>%
      filter(type == "exon") %>%
      ggplot(aes(
        xstart = start,
        xend = end,
        y = df_mane$transcript_id %>% unique()
      )) +
      ggtranscript::geom_intron(
        data = to_intron(exons, "transcript_name"),
        aes(strand = exons$strand %>% unique)
      ) +    #geom_range()
      geom_half_range(
        range.orientation = "top",
        data = ref_introns_MSRD,
        mapping = aes(height = MSR_D / 4, 
                      fill = "MSR_Donor")
      ) +
      geom_half_range(
        range.orientation = "top",
        data = ref_introns_MSRA,
        mapping = aes(height = MSR_A / 4, 
                      fill = "MSR_Acceptor")
      ) + 
      
      geom_range(
        data = df_mane_cds,
        fill = "purple",
        #height = 0.5
      ) +
      geom_range(
        data = df_mane_utr,
        fill = "purple",
        height = 0.25
      ) +

      scale_fill_manual(breaks = c("MSR_Donor", "MSR_Acceptor"),
                        values = c("MSR_Donor" = "#23C3E4FF", 
                                   "MSR_Acceptor" = "#C82803FF")) +
      #viridis::scale_fill_viridis(discrete = T, option = "F")  +
      #geom_line(data = ref_introns_MSRD,
      #           mapping = aes(x = start,
       #                        y = MSR_D, colour = "MSR_D"))  +
      #geom_line(data = ref_introns_MSR,
      #         mapping = aes(x = end,
      #                       y = MSR_A, colour = "MSR_A"))   
      #scale_colour_manual(breaks = c("novel_acceptor", "novel_donor"),
      #                    values = c("#F8766D", "#00BFC4")) +
      #scale_size_continuous(range = c(0.1, 1)) + 
      # theme_bw() +
      #coord_cartesian(xlim = c(start_intron - 1000, end_intron + 1000)) +

      theme_light() +
      theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
            axis.text = element_text(size = "14"),
            axis.title = element_text(size = "14"),
            plot.title = element_text(size = "16"),
            legend.position = "top",
            legend.text = element_text(size = "13"),
            legend.title = element_text(size = "13")) +
      xlab(paste0("Genomic position (", exons$seqnames %>% unique() ,")")) + 
      ylab("MANE Transcript") +
      guides(fill = guide_legend(element_blank())) +
      labs(title = paste0(clust),
           #caption = "UTR sequences are represented in white and CDS are in purple",
           subtitle = paste0("Mis-splicing activity in the MANE transcript of the ", 
                             df_mane$gene_name %>% unique(), " gene")) %>%
      return()
    
      
      
    
      
  } else {
    ggplot() +
      theme_void() +
      geom_text(aes(0,0,label='The selected gene doesn\'t have a MANE transcript.')) + 
      theme(text = element_text(element_text(size = "14"))) %>%
      return()
  }
  
  
  
}


plot_sample_numbers <- function() {
  # setwd("/home/sruiz/PROJECTS/splicing-project-app/intron_db/")
  con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
  
  # Query to the DB
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
  
  colours <- ifelse(str_detect(string = as.factor(df_metadata %>%
                                                    count(cluster_tidy) %>%
                                                    arrange(desc(n)) %>% 
                                                    pull(cluster_tidy)), pattern = "Brain"), "red", "black")
  
  ggplot(df_metadata %>%
           count(cluster_tidy) %>%
           arrange(n), aes(x=reorder(cluster_tidy,-n,decreasing = F), y = n, fill = cluster_tidy)) + 
    geom_bar(stat="identity") +
    viridis::scale_fill_viridis(discrete = T, option = "F")  +
    #coord_flip() +
    theme_light() +
    xlab("") +
    geom_hline(yintercept= 70,linetype="dotted") +
    scale_y_continuous(name ="Number of Samples", 
                       breaks = c(0,70,200,400,600,800))+
    theme(axis.text.x = element_text(color = colours,
                                     angle = 70, 
                                     hjust = 1,
                                     size = "11"),
          legend.position="none") %>%
    return()
  
}

plot_metadata <- function() {
  
  # setwd("/home/sruiz/PROJECTS/splicing-project-app/intron_db/")
  con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
  
  # Query to the DB
  query <- paste0("SELECT * FROM 'master'")
  df_all_projects_metadata <- dbGetQuery(con, query) 
  
  return(df_all_projects_metadata %>%
           mutate(gender = ifelse(gender == 1, "Male", "Female")) %>%
           dplyr::select(Age = age,
                         RIN = rin,
                         Gender = gender,
                         
                         "Mapped read count" = mapped_read_count,
                         "Avg. read length" = avg_read_length,
                         "Type of nucleic acid isolation batch" = smnabtcht,
                         Tissue = tissue,
                         "Body region" = SRA_project_tidy))
  
  # for (plot in c("samples", "age", "rin")) {
  #   return(paste0("/home/sruiz/PROJECTS/splicing-project-app/intron_db/dependencies/images/", project, "_", plot, ".png"))
  # }
  
}


get_mode <- function(vector) {
  uniqv <- unique(vector)
  uniqv[which.max(tabulate(match(vector, uniqv)))]
}

###################################################
################# VARIABLES #######################
###################################################

# setwd("introverse/")
con <- DBI::dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
DBI::dbListTables(conn = con) %>% print()
chr_choices <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y")
strand_choices <- c("+", "-")


# Query to the DB
query <- paste0("SELECT * FROM 'master'")
df_metadata <- DBI::dbGetQuery(con, query) %>%
  dplyr::arrange(SRA_project_tidy)  

db_choices <- c(df_metadata$SRA_project %>% unique()) %>% as.list()
names(db_choices) <- c(df_metadata %>% 
                         dplyr::distinct(SRA_project, .keep_all =T) %>% 
                         dplyr::pull(SRA_project_tidy)) %>% as.list()


query <- paste0("SELECT * FROM 'gene'")
genes <- DBI::dbGetQuery(con, query)
genes <- genes %>%
  tidyr::drop_na() %>%
  dplyr::arrange(gene_name) 
genes_choices <- c(genes$gene_name, genes$gene_id) %>% as.list()
names(genes_choices) <- c(genes$gene_name, genes$gene_id) %>% as.list()

intronID <- NULL
intronType <- NULL


any(genes_choices %>% names() %>% is.na())
any(genes_choices %>% is.na())
DBI::dbDisconnect(conn = con)


