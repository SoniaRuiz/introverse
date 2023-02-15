
###################################################
################# FUNCTIONS #######################
###################################################

#' Title
#' Main function that receives the parameter selection made by the user in the UI and queries the database accordingly
#' @param type whether to search annotated introns or novel junctions. Values: 'introns', 'novel'
#' @param chr should be the query done by intron coordinates, chromosome of the intron.
#' @param start should be the query done by intron coordinates, starting position of the intron within the selected chromosome.
#' @param end should be the query done by intron coordinates, starting position of the intron within the selected chromosome.
#' @param strand should be the query done by intron coordinates, strand position of the intron within the selected chromosome.
#' @param gene should be the query done by gene name, SYMBOL gene name from the ones contained within the dropdown list (e.g. 'PTEN').
#' @param gene_file should be the query done by gene list, uploaded list
#' @param threshold minimum % of samples in which any of the novel junctions to be retrieved should appear. -1 in case this input is not checked.
#' @param search_type type pf search:
#' - by intron coordinates: 'radio_bycoordinates_tab1'
#' - by gene name: 'radio_bygene_tab1'
#' - by gene list: 'radio_bygenelist_tab1'
#' @param all_data_bases whether the query should be done across all tissues. TRUE/FALSE value.
#' @param data_bases GTEx body part/s from which the selected tissue belongs to (i.e. db = "BRAIN")
#' @param clusters GTEx tissue/s in which the database query should be done (i.e. clust = "Brain - Amygdala")
#' @param mane if the intron should be part of a MANE transcript. TRUE/FALSE value.
#' @param clinvar if the intron should contain pathogenic splicing mutations as reported by the ClinVar database. TRUE/FALSE value.
#'
#' @return datatable with the results of the query to the database
#' @export
#'
#' @examples
#' type = "introns"
#' gene = "PTEN"
#' threshold = -1
#' search_type="radio_bygene_tab1"
#' all_data_bases = F
#' data_bases = "BRAIN"
#' clusters = "Brain - Cerebellum"
#' mane = FALSE
#' clinvar = F
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

  #print(type)
  #print(chr)
  #print(start)
  #print(end)
  #print(strand)
  #print(gene)
  #print(threshold)
  #print(search_type)
  #print(all_data_bases)
  #print(data_bases)
  #print(clusters)
  #print(mane)
  #print(clinvar)
  #print("##########################")

  do_next <- F
  # setwd("/home/sruiz/PROJECTS/splicing-project-app/introverse/")
  con <- dbConnect(RSQLite::SQLite(), "./dependencies/introverse.sqlite")
  
  # Query to the DB
  query = paste0("SELECT * FROM 'master'")
  df_all_projects_metadata <- dbGetQuery(con, query) 


  if (all_data_bases) {
    data_bases <- df_all_projects_metadata$SRA_project %>% unique()
    clusters <- df_all_projects_metadata %>%
      distinct(cluster) %>%
      pull()
  }
  
  # data_bases <- df_all_projects_metadata$SRA_project %>% unique()
  
  df_gr <- map_df(data_bases, function(db_IDB) {

    # db_IDB <- data_bases[1]

    details <- df_all_projects_metadata %>%
      dplyr::filter(SRA_project == db_IDB)
        
    map_df(clusters,  function(clust) {

      # clust <- clusters[1]
 
      if (any(details$cluster == clust)) {
        
        ## Get the metadata for the current cluster
        
        cluster_tidy <- details %>%
          dplyr::filter(cluster == clust) %>% 
          distinct(cluster_tidy)%>%
          pull()
    
        db_tidy <- df_all_projects_metadata %>%
          dplyr::filter(SRA_project == db_IDB) %>%
          distinct(SRA_project_tidy) %>%
          pull()
        
        all_samples <- details %>% 
          dplyr::filter(cluster == clust) %>% 
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
            query <- paste0("SELECT distinct(intron.ref_coordinates), 
            gene.gene_name,
            gene.n_transcripts,
            gene.gene_width,
            intron.ref_ss5score, intron.ref_ss3score, intron.clinvar, intron.TSL,
            intron.u12_intron, intron.u2_intron,
            intron.ref_cons5score, intron.ref_cons3score, intron.ref_CDTS5score, intron.ref_CDTS3score,
            tissue.ref_type, tissue.ref_n_individuals, tissue.ref_sum_counts,
            tissue.MSR_D, tissue.MSR_A,
            tissue.gene_tpm,
            tissue.ref_junID
                            FROM '", clust, "_", db_IDB, "_misspliced' AS tissue
                            INNER JOIN 'intron' ON intron.ref_junID=tissue.ref_junID
                            INNER JOIN 'gene' ON gene.id=intron.gene_id
                            WHERE intron.ref_coordinates == '", ID, "'")
            query_never <- paste0("SELECT intron.ref_coordinates, 
            gene.gene_name,
            gene.n_transcripts,
            gene.gene_width,
            intron.ref_ss5score, intron.ref_ss3score, intron.clinvar, intron.TSL,
            intron.u12_intron, intron.u2_intron,
            intron.ref_cons5score, intron.ref_cons3score, intron.ref_CDTS5score, intron.ref_CDTS3score,
            tissue.MSR_D, tissue.MSR_A,
            tissue.gene_tpm,
            tissue.ref_type, tissue.ref_n_individuals,  tissue.ref_sum_counts
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
            query <- paste0("SELECT distinct(intron.ref_coordinates), 
            gene.gene_name,
            gene.n_transcripts,
            gene.gene_width,
            intron.ref_ss5score, intron.ref_ss3score, intron.clinvar, intron.TSL,
            intron.u12_intron, intron.u2_intron,
            intron.ref_cons5score, intron.ref_cons3score, intron.ref_CDTS5score, intron.ref_CDTS3score,
            tissue.ref_type, tissue.ref_n_individuals, tissue.ref_sum_counts, 
            tissue.MSR_D, tissue.MSR_A,
            tissue.gene_tpm,
            tissue.ref_junID
                            FROM '", clust, "_", db_IDB, "_misspliced' AS tissue
                            INNER JOIN 'intron' ON intron.ref_junID=tissue.ref_junID
                            INNER JOIN 'gene' ON gene.id=intron.gene_id
                            WHERE ", gene_query, "")
            query_never <- paste0("SELECT intron.ref_coordinates, 
            gene.gene_name,
            gene.n_transcripts,
            gene.gene_width,
            intron.ref_ss5score, intron.ref_ss3score, intron.clinvar, intron.TSL, 
            intron.u12_intron, intron.u2_intron,
            tissue.MSR_D, tissue.MSR_A,
            tissue.gene_tpm,
            intron.ref_cons5score, intron.ref_cons3score, intron.ref_CDTS5score, intron.ref_CDTS3score,
            tissue.ref_type, tissue.ref_n_individuals, tissue.ref_sum_counts
                            FROM '", clust, "_", db_IDB, "_nevermisspliced' AS tissue
                            INNER JOIN 'intron' ON intron.ref_junID=tissue.ref_junID
                            INNER JOIN 'gene' ON gene.id=intron.gene_id
                            WHERE ", gene_query, "")
            
          } else {
            
            query <- paste0("SELECT 
            tissue.novel_n_individuals, tissue.novel_sum_counts,
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
          
          gene_query <- paste0("LOWER(gene.gene_name) IN ('", paste(df[,1] %>% tolower(), collapse = "','" ), "') OR LOWER(gene.gene_id) IN ('", 
                               paste(df[,1] %>% tolower(), collapse = "','" ), "')")
          

          query <- paste0("SELECT distinct(intron.ref_coordinates), 
            gene.gene_name,
            gene.n_transcripts,
            gene.gene_width,
            intron.ref_ss5score, intron.ref_ss3score, intron.clinvar, intron.TSL,
            intron.u12_intron, intron.u2_intron,
            intron.ref_cons5score, intron.ref_cons3score, intron.ref_CDTS5score, intron.ref_CDTS3score,
            tissue.ref_type, tissue.ref_n_individuals, tissue.ref_sum_counts, 
            tissue.MSR_D, tissue.MSR_A,
            tissue.gene_tpm,
            tissue.ref_junID
                            FROM '", clust, "_", db_IDB, "_misspliced' AS tissue
                            INNER JOIN 'intron' ON intron.ref_junID=tissue.ref_junID
                            INNER JOIN 'gene' ON gene.id=intron.gene_id
                            WHERE (", gene_query, ")")
          query_never <- paste0("SELECT intron.ref_coordinates, 
            gene.gene_name,
            gene.n_transcripts,
            gene.gene_width,
            intron.ref_ss5score, intron.ref_ss3score, intron.clinvar, intron.TSL,
            intron.u12_intron, intron.u2_intron,
            intron.ref_cons5score, intron.ref_cons3score, intron.ref_CDTS5score, intron.ref_CDTS3score,
            tissue.MSR_D, tissue.MSR_A,
            tissue.gene_tpm,
            tissue.ref_type, tissue.ref_n_individuals, tissue.ref_sum_counts
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
              dplyr::filter(p_novel_ind >= threshold)
            
            df_gr <- df_gr %>%
              dplyr::filter(ref_junID %in% df_novel$ref_junID)
            
          } else {
            df_gr <- df_gr %>%
              dplyr::filter(p_novel_ind >= threshold)
          }
        } else {
          if (type == "introns") {
            # Get never mis-spliced introns too
            df_nevermispliced <- dbGetQuery(con, query_never) 
            df_gr <- plyr::rbind.fill(df_gr %>%
                                        mutate(MSR_D = formatC(x = MSR_D, format = "e", digits = 3),
                                               MSR_A = formatC(x = MSR_A, format = "e", digits = 3)), df_nevermispliced)
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
      
      df_gr %>%
        mutate(MeanCounts = round(x = (ref_sum_counts/ref_n_individuals), digits = 2),
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
               "Gene TPM" = gene_tpm %>% round(digits = 2),
               ref_type = ref_type %>% as.factor(),
               gene_name = gene_name %>% as.factor(),
               TSL = ifelse(TSL == 10, "tslNA", TSL),
               u2_intron = ifelse(u2_intron == 0, "N", "Y"),
               u12_intron = ifelse(u12_intron == 0, "N", "Y")) %>%
        
        dplyr::select(ID = ref_coordinates,
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
                      "U2 Intron" = u2_intron,
                      "U12 Intron" = u12_intron,
                      ClinVar = clinvar,
                      TSL,
                      "Gene TPM",
                      "Gene N. Transcripts" = n_transcripts,
                      "Gene Length (bp)" = gene_width,
                      Gene = gene_name,
                      Samples = Cluster, 
                      "Body region" = DB,
                      ref_junID) %>% 
        return()
    } else {

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
        mutate(novel_type = str_replace_all(string = novel_type, pattern = "_", replacement = " "),
               novel_mean_counts = round(x = novel_sum_counts / novel_n_individuals, digits = 2 )) %>%
        mutate("% Individuals" = paste0(p_novel_ind, "% (", novel_n_individuals, "/", all_samples, ")"),
               "Width" = abs(start_positions-end_positions)+1) %>%
        dplyr::select(ID = novel_coordinates,
               NovelType = novel_type,
               NovelWidth = Width,
               "MES_5ss" = novel_ss5score,
               "MES_3ss" = novel_ss3score, 
               Distance = distance,
               "Avg. Read Count" = novel_mean_counts,
               "% Individuals",
               Gene = gene_name,
               Samples = Cluster,
               "Body region" = DB) %>%
        mutate(Modulo3 = abs(Distance) %% 3) %>%
        relocate(Modulo3, .after = Distance) %>%
        return()
    }
   
    
  }

}




#' Title
#' Retrieves from the database all novel junctions attached to the selected intron.
#' @param intron_id Unique integer identifier assigned to the annotated intron within the IntroVerse database (i.e. intron_id = 139704)
#' @param clust GTEx tissue in which the database query should be done (i.e. clust = "Brain - Amygdala")
#' @param db GTEx body part from which the selected tissue belongs to (i.e. db = "BRAIN")
#'
#' @return A dataframe of novel junctions linked to the selected annotated intron across the samples of the tissue selected.
#' @export
#'
#' @examples
get_novel_data_from_intron <- function(intron_id,
                                       clust,
                                       db) {
  
  intron_id <- str_split(string = intron_id,pattern = "#")[[1]][1]

  con <- dbConnect(RSQLite::SQLite(), "./dependencies/introverse.sqlite")
  
  # Query the DB
  query = paste0("SELECT * FROM 'master'")
  df_all_projects_metadata <- dbGetQuery(con, query) %>%
    dplyr::filter(SRA_project_tidy == db) 
  
  
  db_IDB <- df_all_projects_metadata %>%
    distinct(SRA_project) %>%
    pull()
  
  all_samples <- df_all_projects_metadata %>%
    dplyr::filter(cluster == clust) %>%
    nrow()
  
  cluster_tidy <-  df_all_projects_metadata %>%
    dplyr::filter(cluster == clust) %>%
    distinct(cluster_tidy) %>% 
    pull()
  
  
  query = paste0("SELECT * 
                 FROM '", paste0(clust, "_", db_IDB, "_misspliced"), 
                 "' WHERE ref_junID == '", intron_id, "'")
  df_gr <- dbGetQuery(con, query) 
  
  if ( df_gr %>% nrow() > 0 ) {
    
    query <- paste0("SELECT * FROM 'novel' 
                   WHERE novel_junID IN ('", paste(df_gr$novel_junID, collapse="','"), "')")
    df_novel <- dbGetQuery(con, query)
    
    df_novel <- merge(x = df_novel,
                      y = df_gr,
                      by = "novel_junID")
    
    df_novel_coordinates <- get_genomic_coordinates(df_novel$novel_coordinates)
    
    df_novel %>%
      as.data.frame() %>%
      dplyr::mutate(MeanCounts = round(x = (novel_sum_counts / novel_n_individuals), digits = 2),
                    p_indi = ifelse(round((novel_n_individuals * 100) / all_samples) == 0, 1,
                                    round((novel_n_individuals * 100) / all_samples)),
                    "% Individuals" = paste0(p_indi, "% (", novel_n_individuals, "/", all_samples, ")"),
                    Width = abs(df_novel_coordinates$start - df_novel_coordinates$end) + 1,
                    Samples = cluster_tidy,
                    Project = db,
                    novel_type = str_replace(string = novel_type, pattern = "_", replacement = " "),
                    Modulo3 = abs(distance) %% 3,
                    "Frameshift?" = ifelse(Modulo3 == 0, "N", "Y")) %>%
      dplyr::select(NovelID = novel_coordinates,#Coordinates = coordinates,
                    id = novel_junID,
                    "Novel Type" = novel_type,
                    "Junction Length (bp)" = Width,
                    "MES_5ss" = novel_ss5score,
                    "MES_3ss" = novel_ss3score,
                    "Distance (bp)" = distance,
                    "Frameshift?",
                    "Avg. Read Count" = MeanCounts,
                    "% Individuals",
                    Samples,
                    "Body region" = Project) %>% return()
  } else {
    data.frame(Time = Sys.time(), Message = paste0("Intron not found.")) %>% return()
  }
}



#' Title
#' 
#' @param intron_id  Unique integer ID assigned to the annotated intron within the IntroVerse database (i.e. intron_id = 31836)
#' @param novel_id Unique integer ID assigned to the novel donor or novel acceptor junction within the IntroVerse database (i.e. novel_id = 21345)
#' @param clust GTEx tissue in which the database query should be done (i.e. clust = "Brain - Amygdala")
#' @param db GTEx body part from which the selected tissue belongs to (i.e. db = "BRAIN")
#'
#' @return Using ggtranscript, it visualises the mis-splicing activity of the anntotated intron/novel junction selected within the MANE transcript of its gene.
#' @export
#'
#' @examples
visualise_transcript <- function(intron_id = NULL,
                                 novel_id = NULL,
                                 clust,
                                 db) {

  library(ggplot2)
  library(ggtranscript)
  
  ##############################
  ## GET THE TABLE NAME
  ##############################
  
  query = paste0("SELECT * FROM 'master'")
  con <- dbConnect(RSQLite::SQLite(), "./dependencies/introverse.sqlite")
  df_all_projects_metadata <- dbGetQuery(con, query) 
  db_master_details <- df_all_projects_metadata %>%
    dplyr::filter(SRA_project_tidy == db,
           cluster_tidy == clust)
  db_name <- db_master_details$SRA_project %>% unique()
  cluster_name <- db_master_details$cluster %>% unique()
  
  
  ###################################
  ## QUERY THE DATABASE FOR THE CURRENT INTRON/NOVEL JUNCTION
  ###################################
  
  if ( !is.null(intron_id) ) {
    
    sql_statement <- paste0("SELECT * 
                            FROM '", cluster_name, "_", db_name, "_misspliced' AS tissue
                            INNER JOIN 'intron' ON intron.ref_junID=tissue.ref_junID
                            INNER JOIN 'novel' ON novel.novel_junID=tissue.novel_junID
                            INNER JOIN 'gene' ON gene.id=intron.gene_id
                            WHERE tissue.ref_junID == '", intron_id, "'")
  } 
  
  if ( !is.null(novel_id) ) {
    
    sql_statement <- paste0("SELECT * 
                            FROM '", cluster_name, "_", db_name, "_misspliced' AS tissue
                            INNER JOIN 'intron' ON intron.ref_junID=tissue.ref_junID
                            INNER JOIN 'novel' ON novel.novel_junID=tissue.novel_junID
                            INNER JOIN 'gene' ON gene.id=intron.gene_id
                            WHERE tissue.novel_junID == '", novel_id, "'")
  } 
  
  df_intron <- dbGetQuery(con, sql_statement)
  
  
  if (any(df_intron$MANE)) {
    
    #######################################
    ## EXTRACT THE COORDINATES USING THE ID
    #######################################
    
    
    genomic_coordinates <- get_genomic_coordinates(df_intron$novel_coordinates)
    
    novel_junctions <- merge(x = genomic_coordinates,
                             y = df_intron %>% 
                               dplyr::select(novel_coordinates, novel_n_individuals, 
                                      novel_sum_counts, novel_type, distance),
                             by.x = "ID",
                             by.y = "novel_coordinates")  %>%
      mutate(width = abs(start - end))
      
      
    #############################
    ## GET THE MANE DATA
    #############################
    
    sql_statement <- paste0("SELECT * FROM 'mane'
                            WHERE gene_name == '", df_intron$gene_name %>% unique, "'")
    #print(sql_statement)
    df_mane <- dbGetQuery(con, sql_statement)
    
    df_mane_cds <- df_mane %>% dplyr::filter(type == "CDS")
    df_mane_utr <- df_mane %>% dplyr::filter(type == "UTR")
    
    
    #############################
    ## FIND THE REF INTRON 
    ## PRIOR THE ZOOM
    #############################
    
    all_introns <- ggtranscript::to_intron(df_mane %>% dplyr::filter(type == "exon"))
    
    ref_intron <- rbind(all_introns[which.min(abs(all_introns$start - novel_junctions$start)),],
                        all_introns[which.min(abs(all_introns$end - novel_junctions$end)),]) 
      
    ref_exons <- df_mane %>% arrange(end) %>%
      dplyr::filter(type == "exon") %>%
      dplyr::filter(start %in% ref_intron$end | end %in% ref_intron$start) 
      
    
    #############################
    ## PLOT
    #############################
    
    missplicing_plot <- df_mane %>%
      dplyr::filter(type == "exon") %>% 
      ggplot(aes(
        xstart = start,
        xend = end,
        y = df_mane$transcript_id %>% unique()
      )) +
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
        data = ggtranscript::to_intron(df_mane %>% dplyr::filter(type == "exon"), 
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
      theme_light() +
      ggforce::facet_zoom(xlim = c(min(ref_exons$start),
                                   max(ref_exons$end)),
                          zoom.data = zoom)  + 
      geom_junction_label_repel(
        data = novel_junctions %>% dplyr::mutate( zoom = TRUE ) ,
        aes(label = paste0("seen in ", novel_n_individuals, " out of ", db_master_details %>% nrow, " samples")),
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
      labs(title = paste0("Novel event '", novel_junctions$ID %>% unique, "'"),
           #caption = "UTR sequences are represented in white and CDS are in purple",
           subtitle = paste0("MANE transcript of ", df_mane$gene_name %>% unique(), " (", clust, ")")
        )
    

    missplicing_plot %>%
      return()
    
  
  } else {
    ggplot() +
      theme_void(base_size = 14) +
      geom_text(aes(0,0,label='The selected novel junction doesn´t belong to a MANE transcript.')) %>%
      return()
  }
  
  
  
}


#' Title
#' Using ggtranscript, it visualises the mis-splicing activity of all introns from the selected gene across the samples of the selected tissue
#' @param gene_id Gene symbol or ensemblID to plot its mis-splicing activity (i.e. gene_id = "PTEN")
#' @param clust GTEx tissue in which the query should be done (i.e. clust <- "Brain - Hippocampus")
#'
#' @return ggtranscript plot
#' @export
#'
#' @examples
visualise_missplicing <- function(gene_id,
                                  clust) {


  library(ggplot2)
  library(ggtranscript)
  
  
  ## GET THE DETAILS OF THE CLUSTER/PROJECT SELECTED
  query = paste0("SELECT * FROM 'master'")
  con <- dbConnect(RSQLite::SQLite(), "./dependencies/introverse.sqlite")
  df_all_projects_metadata <- dbGetQuery(con, query) 
  db_master_details <- df_all_projects_metadata %>%
    dplyr::filter(cluster_tidy == clust)
  db_name <- db_master_details$SRA_project %>% unique()
  cluster_name <- db_master_details$cluster %>% unique()
  
  ## GET THE GENE_NAME AND MANE INFO FROM THE INTRON TABLE TO GET THE TRANSCRIPT ID
  sql_statement <- paste0("SELECT * 
                            FROM '", cluster_name, "_", db_name, "_misspliced' AS tissue
                            INNER JOIN 'intron' ON intron.ref_junID=tissue.ref_junID
                            INNER JOIN 'novel' ON novel.novel_junID=tissue.novel_junID
                            INNER JOIN 'gene' ON gene.id=intron.gene_id
                            WHERE gene.gene_name == '", gene_id, "' OR gene.gene_id == '", gene_id, "'")
  
  #print(sql_statement)
  df_gene_splicing <- dbGetQuery(con, sql_statement)
  
  if (nrow(df_gene_splicing) == 0) {
    ggplot() +
      theme_void(base_size = 14) +
      geom_text(aes(0,0,
                    label=paste0("The annotated introns from the selected gene have ",
                    "no evidence of mis-splicing in samples from '", clust, "' tissue."))) %>%
      return()
  } else if (any(df_gene_splicing$MANE)) {
    
    
    genomic_coordinates <- get_genomic_coordinates(df_gene_splicing$ref_coordinates)
    
    ref_introns_MSR <- merge(x = genomic_coordinates %>% 
                               distinct(ID, .keep_all = T),
                             y = df_gene_splicing %>% 
                               dplyr::select(ID = ref_coordinates, MANE, MSR_D, MSR_A, gene_name) %>% 
                               dplyr::distinct(ID, .keep_all = T),
                             by = "ID",
                             all.x = T) 
    
    ##################
    ## ADD MANE DATA
    ##################
    
    
    sql_statement <- paste0("SELECT * FROM 'mane' WHERE gene_name == '", ref_introns_MSR$gene_name %>% unique, "'")
    #print(sql_statement)
    df_mane <- dbGetQuery(con, sql_statement)
    
    
    ## Convert exonic coordinates to introns
    exons <- df_mane %>% dplyr::filter(type == "exon")
    introns <- to_intron(exons, "transcript_name")
    
    
    width_bars <- abs(introns$start - introns$end) %>% min() / 2
    
    ref_introns_MSRD <- ref_introns_MSR %>%
      rowwise() %>%
      mutate(end = start + (abs(start - end) / 4)) %>%
      dplyr::filter(MANE == 1) %>%
      dplyr::select(seqnames , strand, start, end, MSR_D ) %>%
      distinct(start, .keep_all = T) %>%
      distinct(end, .keep_all = T)
    ref_introns_MSRA <- ref_introns_MSR %>%
      rowwise() %>%
      mutate(start = end - (abs(start - end) / 4)) %>%
      dplyr::filter(MANE == 1) %>%
      dplyr::select(seqnames , strand, start, end, MSR_A ) %>%
      distinct(start, .keep_all = T) %>%
      distinct(end, .keep_all = T)
  
    
    df_mane_cds <- df_mane %>% dplyr::filter(type == "CDS")
    df_mane_utr <- df_mane %>% dplyr::filter(type == "UTR")
    
    
    ## Generate the plot
    exons %>%
      mutate(type = str_replace(string = type, pattern = "_", replacement = " ")) %>%
      dplyr::filter(type == "exon") %>%
      ggplot(aes(
        xstart = start,
        xend = end,
        y = df_mane$transcript_id %>% unique()
      )) +
      ggtranscript::geom_intron(
        data = to_intron(exons, "transcript_name"),
        aes(strand = exons$strand %>% unique)
      ) +   
      geom_half_range(
        range.orientation = "top",
        data = ref_introns_MSRD,
        mapping = aes(height = MSR_D, 
                      fill = "MSR Donor")
      ) +
      geom_half_range(
        range.orientation = "top",
        data = ref_introns_MSRA,
        mapping = aes(height = MSR_A, 
                      fill = "MSR Acceptor")
      ) + 
      geom_range(
        data = df_mane_cds,
        fill = "purple"
      ) +
      geom_range(
        data = df_mane_utr,
        fill = "purple",
        height = 0.25
      ) +
      geom_junction_label_repel(
        data = rbind(ref_introns_MSRD %>% dplyr::rename("MSRd" = MSR_D) %>% gather(key = type, value = MSR, -seqnames, -strand,-start,-end),
                     ref_introns_MSRA %>% dplyr::rename("MSRa" = MSR_A) %>% gather(key = type, value = MSR, -seqnames, -strand,-start,-end)) %>%
          filter(round(MSR,digits = 2) >= 0.01),
        aes(label = paste0(type, " = ", round(MSR,digits = 2))),
        nudge_y = 0.5,
        nudge_x = 0,
        junction.y.max = 0,
        point.padding = unit(0.5, "lines"),
        segment.color = 'grey50'
      ) +
      scale_fill_manual(breaks = c("MSR Donor", "MSR Acceptor"),
                        values = c("MSR Donor" = "#23C3E4FF", 
                                   "MSR Acceptor" = "#C82803FF")) +
      theme_light() +
      theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
            axis.text = element_text(size = "14"),
            axis.title = element_text(size = "14"),
            plot.title = element_text(size = "16"),
            plot.caption = element_text(size = "11"),
            legend.position = "top",
            legend.text = element_text(size = "13"),
            legend.title = element_text(size = "13")) +
      xlab(paste0("Genomic position (", exons$seqnames %>% unique() ,")")) + 
      ylab("MANE Transcript") +
      guides(fill = guide_legend(element_blank())) +
      labs(title = paste0(clust),
           caption = "*Only MSR values higher than 0.01 are labelled.",
           subtitle = paste0("Mis-splicing activity in the MANE transcript of ", 
                             df_mane$gene_name %>% unique(), "")) %>%
      return()
      
  } else {
    ggplot() +
      theme_void(base_size = 14) +
      geom_text(aes(0,0,
                    label=paste0("The annotated introns from the selected gene have not been found in a MANE transcript in '",
                                 clust, "' tissue."))) %>%
      return()
  }

}


#' Title
#' Using ggplot, it produces a bar plot containing the number of samples used per GTEx tissue.
#' @return ggplot2 plot
#' @export
#'
#' @examples
plot_sample_numbers <- function() {
  
  con <- dbConnect(RSQLite::SQLite(), "./dependencies/introverse.sqlite")
  
  # Query the DB
  query <- paste0("SELECT * FROM 'master'")
  df_metadata <- dbGetQuery(con, query) 
  
  ggplot(df_metadata %>%
           dplyr::count(cluster_tidy) %>%
           arrange(n), 
         aes(x=reorder(cluster_tidy,-n,decreasing = F), y = n, fill = cluster_tidy)) + 
    geom_bar(stat="identity") +
    viridis::scale_fill_viridis(discrete = T, option = "F")  +
    theme_light() +
    xlab("") +
    geom_hline(yintercept= 70,linetype="dotted") +
    scale_y_continuous(name ="Number of Samples", 
                       breaks = c(0,70,200,400,600,800))+
    theme(axis.text.x = element_text(color = "black",
                                     angle = 70, 
                                     hjust = 1,
                                     size = "11"),
          legend.position="none") %>%
    return()
}



#' Title
#' Function to visualise the metadata of the samples used per GTEx project. 
#' @return Table with the metadata details of each sample used
#' @export
#'
#' @examples
plot_metadata <- function() {

  con <- dbConnect(RSQLite::SQLite(), "./dependencies/introverse.sqlite")
  
  # Query the DB
  query <- paste0("SELECT * FROM 'master'")
  df_all_projects_metadata <- dbGetQuery(con, query) 
  
  return(df_all_projects_metadata %>%
           mutate(gender = ifelse(gender == 1, "Male", "Female")) %>%
           dplyr::select(Age = age,
                         RIN = rin,
                         Gender = gender,
                         
                         "Mapped read count" = mapped_read_count,
                         "Avg. read length (bp)" = avg_read_length,
                         "Type of nucleic acid isolation batch" = smnabtcht,
                         Tissue = tissue,
                         "Body region" = SRA_project_tidy))
  
}



#############
## UTILS   ##
#############

#' Title
#' Utils function to extract the genomic coordinates from an annotated intron or novel junction ID
#' @param coordinates Genomic coordinates of an annotated intron or novel junction (i.e. coordinates='chr19:44905842-44906601:+')
#'
#' @return A dataframe containing the genomic coordinates of the intron/annotated junction: e.g.
#'  data.frame(ID = 'chr19:44905842-44906601:+', seqnames = 'chr19', start = 44905842, end = 44906601, strand = '-')
#' @export
#'
#' @examples
get_genomic_coordinates <- function(coordinates) {
  
  map_df(coordinates, function(coordinate) {
    # coordinate <- df_gene_splicing$novel_coordinates[1]
    chr_junc <- coordinate %>%
      str_sub(start = 1,
              end = str_locate_all(string = coordinate, pattern = ":")[[1]][1,2]-1)
    start_junc <- coordinate %>%
      str_sub(start = str_locate_all(string = coordinate, pattern = ":")[[1]][1,2]+1,
              end = str_locate_all(string = coordinate, pattern = "-")[[1]][1,2]-1)
    end_junc <- coordinate %>%
      str_sub(start = str_locate_all(string = coordinate, pattern = "-")[[1]][1,2]+1,
              end = str_locate_all(string = coordinate, pattern = ":")[[1]][2,2]-1)
    strand_junc <- coordinate %>%
      str_sub(start = str_locate_all(string = coordinate, pattern = ":")[[1]][2,2]+1,
              end = coordinate %>% stringr::str_count())
    
    data.frame(ID = coordinate,
               seqnames = chr_junc,
               start = start_junc %>% as.integer(),
               end = end_junc %>% as.integer(),
               strand = strand_junc) %>%
      return()
    
  })
}





################################
## PRE-FILLING INPUTS    #######
################################


# setwd("introverse/")
con <- DBI::dbConnect(RSQLite::SQLite(), "./dependencies/introverse.sqlite")

## Chr and strand lists
chr_choices <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y")
strand_choices <- c("+", "-")

# Query the DB to extract the gene info
query <- paste0("SELECT * FROM 'master'")
df_metadata <- DBI::dbGetQuery(con, query) %>%
  dplyr::arrange(SRA_project_tidy)  

## Gene List
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


