
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


#INTRON SEARCH
type ="introns"
chr ="19"
start = 44905842
end = 44906601
strand = "+"
gene = "ENSG00000171862"
threshold=70
search_type= "radio_bygene_tab1"
data_bases= "BRAIN"
clusters= "Brain - Hippocampus"
mane= TRUE
clinvar= FALSE

# NOVEL SEARCH

type ="novel"
chr ="19"
start = 44906263
end = 44906601
strand = "+"
gene = "ENSG00000171862"
threshold= 1
search_type= "radio_bygene_tab2"
data_bases= "BLOOD"
clusters= "Whole Blood"
mane= FALSE
clinvar= FALSE



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
            tissue.ref_type, tissue.ref_n_individuals, tissue.ref_mean_counts, tissue.MSR_D, tissue.MSR_A,
            tissue.ref_junID
                            FROM '", clust, "_", db_IDB, "_misspliced' AS tissue
                            INNER JOIN 'intron' ON intron.ref_junID=tissue.ref_junID
                            INNER JOIN 'gene' ON gene.id=intron.gene_id
                            WHERE intron.ref_coordinates == '", ID, "'")
            query_never <- paste0("SELECT intron.ref_coordinates, gene.gene_name,
            intron.ref_ss5score, intron.ref_ss3score, intron.clinvar,
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
            tissue.ref_type, tissue.ref_n_individuals, tissue.ref_mean_counts, tissue.MSR_D, tissue.MSR_A,
            tissue.ref_junID
                            FROM '", clust, "_", db_IDB, "_misspliced' AS tissue
                            INNER JOIN 'intron' ON intron.ref_junID=tissue.ref_junID
                            INNER JOIN 'gene' ON gene.id=intron.gene_id
                            WHERE ", gene_query, "")
            query_never <- paste0("SELECT intron.ref_coordinates, gene.gene_name,
            intron.ref_ss5score, intron.ref_ss3score, intron.clinvar,
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
          
          gene_query <- paste0("LOWER(gene.gene_name) IN ('", paste(df[,1] %>% tolower(), collapse = "','" ), "') OR LOWER(gene.gene_id) IN ('", paste(df[,1] %>% tolower(), collapse = "','" ), "')")
          #showNotification("This is a notification.")
          
          #} else {
          #  gene_query <- paste0("gene.gene_name == '", gene, "'")
          #}
          print(paste0(gene_query))
          query <- paste0("SELECT distinct(intron.ref_coordinates), gene.gene_name,
            intron.ref_ss5score, intron.ref_ss3score, intron.clinvar, 
            tissue.ref_type, tissue.ref_n_individuals, tissue.ref_mean_counts, tissue.MSR_D, tissue.MSR_A,
            tissue.ref_junID
                            FROM '", clust, "_", db_IDB, "_misspliced' AS tissue
                            INNER JOIN 'intron' ON intron.ref_junID=tissue.ref_junID
                            INNER JOIN 'gene' ON gene.id=intron.gene_id
                            WHERE ", gene_query, "")
          query_never <- paste0("SELECT intron.ref_coordinates, gene.gene_name,
            intron.ref_ss5score, intron.ref_ss3score, intron.clinvar,
            tissue.ref_type, tissue.ref_n_individuals, tissue.ref_mean_counts
                            FROM '", clust, "_", db_IDB, "_nevermisspliced' AS tissue
                            INNER JOIN 'intron' ON intron.ref_junID=tissue.ref_junID
                            INNER JOIN 'gene' ON gene.id=intron.gene_id
                            WHERE ", gene_query, "")
      
        }

        
        if (mane) {
          query = paste0(query, " AND intron.MANE == ", mane)
        }
        
        if (clinvar) {
          query = paste0(query, " AND intron.clinvar != '-'")
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
                                               MSR_A = formatC(x = MSR_D, format = "e", digits = 3)), df_nevermispliced)
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
               ref_junID = paste0(ref_junID, "#", Cluster)) %>%
        #,ifelse(MANE == 0, "F", "T"),
               #coordinates = paste0(seqnames,":",start,"-",end,":",strand)) %>%
        #"% Individuals" = round(x = (ref_n_individuals * 100) / all_people_tissue)) %>%
        dplyr::select(ID = ref_coordinates,
                      #Coordinates = coordinates,
                      "Mis-spliced site" = ref_type,
                      "Intron Length (bp)" = Width, 
                      "MES_5ss" = ref_ss5score,
                      "MES_3ss" = ref_ss3score,
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


# setwd("intron_db/")
# intron_id <- "139704"
# db <- "BRAIN"
# clust <- "Brain - Hippocampus"

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
  
  ## GET THE DETAILS OF THE CLUSTER/PROJECT SELECTED
  query = paste0("SELECT * FROM 'master'")
  con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
  df_all_projects_metadata <- dbGetQuery(con, query) 
  db_master_details <- df_all_projects_metadata %>%
    filter(SRA_project_tidy == db,
           cluster_tidy == clust)
  db_name <- db_master_details$SRA_project %>% unique()
  cluster_name <- db_master_details$cluster %>% unique()
  
  
  ## GET THE GENE_NAME AND MANE INFO FROM THE INTRON TABLE TO GET THE TRANSCRIPT ID
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
  ## ADD OTHER DATA
  novel_junctions <- merge(x = novel_junctions,
                           y = df_intron %>% select(novel_coordinates, novel_mean_counts, novel_type),
                           by.x = "ID",
                           by.y = "novel_coordinates") 

  
  sql_statement <- paste0("SELECT * 
                          FROM 'mane'
                          WHERE gene_name == '", df_intron$gene_name %>% unique, "'")
  print(sql_statement)
  df_mane <- dbGetQuery(con, sql_statement)
  
  ## Declare some variables that we will need for the plot
  tx_name <- df_mane$transcript_id %>% unique()
  chr_intron <- (df_intron$ref_coordinates %>% unique() %>%
                     str_sub(start = 1,
                             end = str_locate_all(string = df_intron$ref_coordinates %>% unique(), pattern = ":")[[1]][1,2]-1)) 
  start_intron <- (df_intron$ref_coordinates %>% unique() %>%
                     str_sub(start = str_locate_all(string = df_intron$ref_coordinates %>% unique(), pattern = ":")[[1]][1,2]+1,
                             end = str_locate_all(string = df_intron$ref_coordinates %>% unique(), pattern = "-")[[1]][1,2]-1)) %>% as.integer()
  end_intron <- (df_intron$ref_coordinates %>% unique() %>%
                   str_sub(start = str_locate_all(string = df_intron$ref_coordinates %>% unique(), pattern = "-")[[1]][1,2]+1,
                           end = str_locate_all(string = df_intron$ref_coordinates %>% unique(), pattern = ":")[[1]][2,2]-1)) %>% as.integer()
  strand_intron <- df_intron$ref_coordinates %>% unique() %>%
    str_sub(start = str_locate_all(string = df_intron$ref_coordinates %>% unique(), pattern = ":")[[1]][2,2]+1,
            end = df_intron$ref_coordinates %>% unique() %>% stringr::str_count())
  
  ## Generate the plot
  df_mane %>%
    filter(type == "exon") %>%
    ggplot(aes(
      xstart = start,
      xend = end,
      y = tx_name
    )) +
    ggtranscript::geom_range(
      #fill = "white", 
      height = 0.35
      ) +
    ggtranscript::geom_intron(
      data = ggtranscript::to_intron(df_mane %>%
                         filter(type == "exon"#,
                                #(end == (start_intron - 1) |
                                #   start == (end_intron + 1))
                         ), 
                       "transcript_name"),
      aes(strand = strand_intron)
    ) + 
    ggtranscript::geom_junction(
      data = novel_junctions,
      aes(#size = novel_mean_counts,
          colour = novel_type),
      #junction.orientation = "alternating",
      #angle = 90,
      
      ncp = 100, 
      junction.y.max = 0.5 
    ) +
    scale_colour_manual(breaks = c("novel_acceptor", "novel_donor"),
                        values = c("#F8766D", "#00BFC4")) +
    #scale_size_continuous(range = c(0.1, 1)) + 
   # theme_bw() +
    #coord_cartesian(xlim = c(start_intron - 1000, end_intron + 1000)) +
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
          axis.text = element_text(size = "14"),
          axis.title = element_text(size = "14"),
          legend.position = "top",
          legend.text = element_text(size = "13"),
          legend.title = element_text(size = "13")) +
    xlab(paste0("Genomic position (",chr_intron,")")) + 
    ylab("MANE Transcript") +
    guides(color = guide_legend(title = "Novel type: ")) %>%
    return()
  
  } else {
    ggplot() +
      theme_void() +
      geom_text(aes(0,0,label='The selected novel junction doesn´t belong to a MANE transcript.')) + 
      theme(text = element_text(element_text(size = "14"))) %>%
      return()
  }
  
  
  
}






plot_metadata <- function() {
  
  # setwd("/home/sruiz/PROJECTS/splicing-project-app/intron_db/")
  con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
  
  # Query to the DB
  query <- paste0("SELECT * FROM 'master'")
  df_all_projects_metadata <- dbGetQuery(con, query) 
  

      return(df_all_projects_metadata)
    
  
  
  # for (plot in c("samples", "age", "rin")) {
  #   return(paste0("/home/sruiz/PROJECTS/splicing-project-app/intron_db/dependencies/images/", project, "_", plot, ".png"))
  # }
  
}

# ###################################################
# ############ PLOT FUNCTIONS #######################
# ###################################################
# 
# 
# # tissue <- "PD"
# # limit_bp = 30
# plot_distances <- function(tissue,
#                            limit_bp = 30) {
#   
#   
#   # ## Load core shared junctions across tissues files  
#   negative_bp <- limit_bp * (-1)
#   query = paste0("SELECT * FROM '", paste0(tissue, "_db_novel"), 
#                  "' WHERE (novel_type == 'novel_donor' OR novel_type == 'novel_acceptor') AND distance <= ", 
#                  limit_bp," AND distance >= ", negative_bp)
#   
#   
#   con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
#   df_gr <- dbGetQuery(con, query) 
#   dbDisconnect(con)
#   
# 
#   ## Replace the underscore
#   df_gr <- df_gr %>%
#     mutate(novel_type = str_replace(string = novel_type, pattern = "_", replacement = " "))
#   
#   df_gr$novel_type = factor(df_gr$novel_type, 
#                                levels = c("novel donor", "novel acceptor"))
#   
#   
#   # y_axes <- c(0, (df_gr %>% 
#   #                   filter(novel_type == "novel_acceptor", 
#   #                          distance == (df_gr$distance %>% get_mode())) %>% 
#   #                     nrow()) + 30)
# 
#   
#   ggplot(data = df_gr) + 
#     geom_histogram(aes(x = distance, group = novel_type, fill = novel_type), 
#                    bins = limit_bp * 2) +
#     
#     facet_grid(vars(novel_type)) +
#     xlab("distance (in bp)") +
#     theme_light() +
#     scale_x_continuous(limits = c((limit_bp * -1), limit_bp),
#                        breaks = c((limit_bp * -1), (round(limit_bp / 2) * -1), 0, round(limit_bp / 2), limit_bp)) +
#     scale_fill_manual(values = c("#35B779FF","#440154FF"),
#                       breaks = c("novel donor", "novel acceptor"),
#                       labels = c("novel donor", "novel acceptor")) +
#     guides(fill = guide_legend(title = NULL, #title = "Junction category & Strand",
#                                override.aes = list(size = 3),
#                                ncol = 3 )) +
#     theme(axis.line = element_line(colour = "black"), 
#           axis.text = element_text(colour = "black", size = "14"),
#           axis.title = element_text(colour = "black", size = "14"),
#           strip.text = element_text(colour = "black", size = "16"), 
#           legend.text = element_text(colour = "black", size = "14"),
#           plot.caption = element_text(colour = "black", size = "14"),
#           plot.title = element_text(colour = "black", size = "16"),
#           legend.title = element_text(colour = "black", size = "14"),
#           legend.position = "top")  %>% return()
# }
# 
# 
# plot_modulo <- function(tissue,
#                         limit_bp = 30) {
#   
#   # Query to the DB
#   negative_bp <- limit_bp * (-1)
#   query = paste0("SELECT * FROM '", paste0(tissue, "_db_novel"), 
#                  "' WHERE (novel_type == 'novel_donor' OR novel_type == 'novel_acceptor') AND distance <= ", 
#                  limit_bp," AND distance >= ", 
#                  negative_bp)
#   
#   
#   con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
#   df_gr <- dbGetQuery(con, query) 
#   dbDisconnect(con)
#   
#   
#   df_gr$novel_type = factor(df_gr$novel_type, 
#                             levels = c("novel_donor", "novel_acceptor"))
#   
#   df_gr <- df_gr %>%
#     filter(distance >= limit_bp*(-1), distance <= limit_bp) %>%
#     mutate(type_p = ifelse(distance < 0, paste0(novel_type, "_intron"), paste0(novel_type, "_exon"))) %>% 
#     mutate(module = round(abs(distance) %% 3, digits = 2))
#   
#   
#   ggplot(data = df_gr, aes(x = module, group = type_p, fill = type_p)) + 
#     geom_bar(aes(y = ..prop..), stat = "count") +
#     geom_text(aes( label = scales::percent(..prop.., accuracy = 0.1),
#                    y= ..prop.. ), stat= "count", vjust = -.5) +
#     scale_x_continuous(breaks = c(0, 1, 2)) +
#     scale_y_continuous(labels = scales::percent) +
#     scale_fill_viridis_d() +
#     facet_grid(~type_p) + 
#     ylab("percentage of novel junctions") +
#     xlab("modulo 3") +
#     theme_light() +
#     theme(axis.line = element_line(colour = "black"), 
#           axis.text = element_text(colour = "black", size = "14"),
#           axis.title = element_text(colour = "black", size = "14"),
#           legend.text = element_text(colour = "black",size = "14"),
#           strip.text = element_text(colour = "black", size = "12"), 
#           plot.caption = element_text(colour = "black",size = "14"),
#           legend.title = element_text(colour = "black", size = "14"),
#           legend.position = "none") +
#     guides(fill = guide_legend(title = NULL,
#                                ncol = 2, 
#                                nrow = 2)) %>% return()
# }
# 
# 
# plot_missplicing <- function(tissue) {
#   
#   # Query to the DB
#   query = paste0("SELECT * FROM '", paste0(tissue, "_db_introns"),"'")
#   
#   
#   con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
#   df_gr <- dbGetQuery(con, query) 
#   dbDisconnect(con)
#   
#   df_gr$ref_type %>% unique %>% print()
#   
#   # df_gr$ref_type = factor(df_gr$ref_type, 
#   #                           levels = c("donor", "acceptor"))
#   
#   # y_axes <- c(0, (df_gr %>% 
#   #                   filter(ref_type == "never", 
#   #                          ref_missplicing_ratio_tissue_NA > 0,
#   #                          ref_missplicing_ratio_tissue_NA < 0.003) %>% 
#   #                   nrow()) + 30)
#   
#   ggplot(df_gr) + 
#     geom_density(aes(x = ref_missplicing_ratio_tissue_NA, fill = "#440154FF"), 
#                    alpha = 0.8, bins = 60) +
#     geom_density(aes(x = ref_missplicing_ratio_tissue_ND, fill = "#35B779FF"), 
#                    alpha = 0.8, bins = 60) +
#     xlab("Mis-splicing ratio (MSR)") +
#     theme_light() +
#     scale_fill_manual(values = c("#35B779FF","#440154FF"),
#                       breaks = c("#35B779FF","#440154FF"),
#                       labels = c("MSR_Donor","MSR_Acceptor")) +
#     theme(axis.line = element_line(colour = "black"), 
#           axis.text = element_text(colour = "black", size = "14"),
#           axis.title = element_text(colour = "black", size = "14"),
#           plot.title = element_text(colour = "black", size = "16"),
#           legend.text = element_text(size = "14"),
#           legend.title = element_text(size = "14"),
#           legend.position = "top") +
#     guides(fill = guide_legend(title = NULL,
#                                ncol = 2, 
#                                nrow = 1)) %>% return()
# 
# }
# 
# 
# 
# 
# plot_intron_proportions <- function(clusters) {
#   
#   # clusters <- gtex_tissues
#   # clusters <- c("PD", "control")
#   # folder_root = "/home/sruiz/PROJECTS/splicing-project/splicing-recount2-projects/data/SRP058181/results/pipeline3/missplicing-ratio/"
#   # folder_results = "/home/sruiz/PROJECTS/splicing-project/splicing-recount2-projects/data/SRP058181/results/pipeline3/alltypes/"
# 
#   folder_results <- "./dependencies/"
#   
#   
#   #df_missplicing_all %>% head() %>% print()
#   #df_missplicing_all <- readRDS(file = paste0(folder_results, "/df_", cluster, "_missplicing_all.rds"))
#   
#   
#   if (any(clusters == "PD") || any(clusters == "control")) {
#     df_missplicing_all <- readRDS(file = paste0(folder_results, "/SRP058181/df_missplicing_all.rds"))
#   } else {
#     df_missplicing_all <- readRDS(file = paste0(folder_results, "/df_missplicing_all.rds"))
#   }
#   
#   
#   ## Plot results --------------------------------------------------------------------
#   
#   df <- data.frame(data = df_missplicing_all$prop_both,
#                    tissue = df_missplicing_all$tissue,
#                    type = "both")
#   df <- rbind(df, 
#               data.frame(data = df_missplicing_all$prop_acceptor,
#                          tissue = df_missplicing_all$tissue,
#                          type = "acceptor"))
#   
#   df <- rbind(df, 
#               data.frame(data = df_missplicing_all$prop_donor,
#                          tissue = df_missplicing_all$tissue,
#                          type = "donor"))
#   
#   
#   df <- rbind(df, 
#               data.frame(data = df_missplicing_all$prop_none,
#                          tissue = df_missplicing_all$tissue,
#                          type = "never"))
#   
#   
#   df$tissue <- factor(df$tissue, levels = df$tissue[order(df %>% filter(type == "never") %>% pull(data) %>% dplyr::desc())] %>% unique)
#   
#   
#   
#   breaks <- levels(as.factor(df$tissue[order(df$data %>% dplyr::desc())] %>% unique))
#   colours <- ifelse(str_detect(string = as.factor(df$tissue[order(df %>% filter(type == "never") %>% pull(data) %>% dplyr::desc())] %>% unique), pattern = "Brain"), 
#                     "red", "black")
#   
#   if (any(clusters == "PD") || any(clusters == "control")) {
#     label = round(df$data * 100, digits = 2)
#     xlabel = "sample type"
#   } else {
#     label = ""
#     xlabel = "tissue"
#   }
#   
#   ggplot(data = df, aes(x = tissue, y = data, fill = type)) +
#     geom_bar(stat = "identity")  +
#     geom_text(aes(label = label), position=position_stack(0.5), color = "red") +
#     
#     theme_light() +
#     ylab("proportion of introns") +
#     xlab(xlabel) +
#     scale_fill_viridis_d() +
#     ggtitle("Proportion of intron type.\nOrdered by the 'never mis-spliced' intron category.") +
#     theme(axis.line = element_line(colour = "black"), 
#           axis.text = element_text(colour = "black", size = "12"),
#           axis.text.x = element_text(angle = 70, 
#                                      vjust = 1,
#                                      color = colours,
#                                      hjust = 1),
#           axis.title = element_text(colour = "black", size = "12"),
#           legend.text = element_text(size = "12"),
#           legend.title = element_text(size = "12"),
#           legend.position = "top") +
#     guides(fill = guide_legend(title = "Junction type: ",
#                                ncol = 2, 
#                                nrow = 2)) %>% return()
#   
#   
#   # if (save_result) {
#   #   file_name <- paste0(folder_results, "/proportion_junction_tissue_never_ordered.png")
#   #   ggplot2::ggsave(filename = file_name,
#   #                  width = 183, height = 183, units = "mm", dpi = 300)
#   # }
#   
#   
# }
# 
# # tissue <- "Brain-FrontalCortex_BA9"
# # tissue <- "Brain-Substantianigra"
# # tissue <- "control_SRP058181"
# # tissue <- "PD_SRP058181"
# # tissue <- "control_SRP049203"
# # tissue <- "PD_SRP049203"
# plot_lm <- function(tissue) {
#   
#   if (str_detect(string = tissue,
#                  pattern = "PD") || str_detect(string = tissue,
#                                                pattern = "control")) {
#     
#     # Query to the DB
#     query = paste0("SELECT * FROM '", paste0(tissue, "_db_lm"),"'")
#     
#     
#     con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
#     df_stats <- dbGetQuery(con, query) 
#     dbDisconnect(con)
#     
#     # if (str_detect(string = tissue, pattern = "PD")) {
#     #   df_stats <- df_stats %>%
#     #     mutate(disease_state = disease_state %>% as.factor()) %>%
#     #     mutate(disease_state = relevel(disease_state, ref = "control"))
#     # } else {
#     #   df_stats <- df_stats %>%
#     #     mutate(disease_state = disease_state %>% as.factor()) %>%
#     #     mutate(disease_state = relevel(disease_state, ref = "PD"))
#     # }
#     # 
#     df_stats[1,]
#     
#     df_stats <- df_stats %>%
#       filter(u2_intron == T | u12_intron == T) %>% 
#       distinct(ref_junID, .keep_all = T) %>%
#       dplyr::rename(is_PD = disease_state,
#                     gene_tpm = tpm) %>%
#       mutate(is_control = is_PD) %>%
#       mutate(is_PD = ifelse(is_PD == "PD", T, F)) %>%
#       mutate(is_control = ifelse(is_control == "control", T, F))
# 
#     df_stats %>% nrow()
#     df_stats %>% distinct(ref_junID) %>% nrow()
#     df_stats[1,]
#     
#     fit_donor <- lm(ref_missplicing_ratio_tissue_ND ~ intron_length + 
#                       intron_5ss_score + 
#                       intron_3ss_score  +
#                       disease_state +
#                       protein_coding +
#                       #is_control+
#                       gene_tpm + 
#                       #u12_intron +
#                       u2_intron +
#                       gene_length + 
#                       gene_num_transcripts, 
#                     #u2_intron ,
#                     data = df_stats)
#     
#     fit_acceptor <- lm(ref_missplicing_ratio_tissue_NA ~ intron_length + 
#                          intron_5ss_score * 
#                          intron_3ss_score +
#                          disease_state +
#                          #is_control+
#                          #u12_intron +
#                          u2_intron +
#                          #gene_tpm + 
#                          gene_length + 
#                          #protein_coding +
#                          gene_num_transcripts,
#                        #protein_coding,
#                        data = df_stats)
#     
#     fit_donor %>% summary()
#     fit_acceptor %>% summary()
#     
#   } else {
#     
#     # Query to the DB
#     query = paste0("SELECT * FROM '", paste0(tissue, "_db_introns"),"'")
#     
#     
#     con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
#     df_gr <- dbGetQuery(con, query) 
#     dbDisconnect(con)
#     
#     ## RENAME COLUMNS
#     df_stats <- df_gr %>%
#       dplyr::rename(intron_length = width,
#                     intron_5ss_score = ref_ss5score,
#                     intron_3ss_score = ref_ss3score,
#                     gene_length = gene_width,
#                     gene_tpm = tpm,
#                     gene_num_transcripts = transcript_number) %>%
#       filter(gene_tpm > 0) %>%
#       filter(u2_intron == T | u12_intron == T)
#     
#     
#     
#     
#     fit_donor <- lm(ref_missplicing_ratio_tissue_ND ~ intron_length + 
#                       intron_5ss_score + intron_3ss_score + 
#                       gene_tpm + 
#                       gene_length + 
#                       gene_num_transcripts + 
#                       #u12_intron +
#                       u2_intron +
#                       protein_coding,
#                     data = df_stats)
#     
#     fit_acceptor <- lm(ref_missplicing_ratio_tissue_NA ~ intron_length + 
#                          intron_5ss_score + intron_3ss_score +
#                          gene_tpm + 
#                          gene_length + 
#                          gene_num_transcripts + 
#                          #u12_intron +
#                          u2_intron +
#                          protein_coding,
#                        data = df_stats)
#   }
#   
#   fit_donor %>% summary()
#   fit_acceptor %>% summary()
#   
#   
#   
#   ##########################
#   ## PLOT LINEAR MODELS
#   ##########################
#   
#   jtools::plot_summs(fit_donor, 
#                      fit_acceptor,
#                      scale = TRUE, 
#                      robust = list("HC3", "HC3"),
#                      #inner_ci_level = .75,
#                      n.sd = 2,
#                      legend.title = "Model:",
#                      #plot.distributions = TRUE,
#                      ci_level = 0.95,
#                      colors = c("#35B779FF","#440154FF"),
#                      model.names = c("MSR_Donor", "MSR_Acceptor")) + 
#     
#     theme_minimal() + 
#     theme(axis.line = element_line(colour = "black"), 
#           axis.text.x = element_text(colour = "black", size = "13"),
#           axis.text.y = element_text(colour = "black", size = "13"),
#           axis.title = element_text(colour = "black", size = "13"),
#           legend.text = element_text(colour = "black", size = "13"),
#           legend.title = element_text(colour = "black", size = "13"),
#           legend.position = "top",
#           panel.grid.major.x = element_blank(),
#           panel.grid.major.y = element_blank(),
#           axis.ticks.x = element_line(colour = "black", size = 2)) +  
#     ylab("Predictors") +
#     geom_hline(yintercept = seq(from = -0.5, to = (length((fit_donor$coefficients %>% names)[-1]) + .5)))  %>% return()
#   
#   
# }
get_mode <- function(vector) {
  uniqv <- unique(vector)
  uniqv[which.max(tabulate(match(vector, uniqv)))]
}

###################################################
################# VARIABLES #######################
###################################################

# setwd("intron_db/")
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


