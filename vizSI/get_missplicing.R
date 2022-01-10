
###################################################
################# FUNCTIONS #######################
###################################################


get_missplicing_ratio <- function(intron_coordinates = GRanges(),
                                  tissue = NULL,
                                  all_tissues = F) {
  
  # if (all_tissues) {
  #   
  #   
  #   # ## Load core shared junctions across tissues files  
  #   query = paste0("SELECT * FROM both WHERE seqnames == ", intron_coordinates@seqnames , 
  #                  " AND start == ", intron_coordinates@ranges@start, 
  #                  " AND end == ", intron_coordinates %>% ranges %>% end, 
  #                  " AND strand == '", intron_coordinates@strand, "'")
  #   
  #   con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
  #   both_misspliced_junc_gr <- dbGetQuery(con, query) 
  #   dbDisconnect(con)
  #   
  #   both_x <- findOverlaps(query = both_misspliced_junc_gr %>% GRanges(),
  #                          subject = intron_coordinates,
  #                          ignore.strand = F,
  #                          type = "equal")
  #   
  #   if (both_misspliced_junc_gr[queryHits(both_x),] %>% length() == 0) {
  #     
  #     query = paste0("SELECT * FROM donor WHERE seqnames == ", intron_coordinates@seqnames , 
  #                    " AND start == ", intron_coordinates@ranges@start, 
  #                    " AND end == ", intron_coordinates %>% ranges %>% end, 
  #                    " AND strand == '", intron_coordinates@strand, "'")
  #     
  #     con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
  #     donor_misspliced_junc_gr <- dbGetQuery(con, query) 
  #     dbDisconnect(con)
  #     
  #     donor_x <- findOverlaps(query = donor_misspliced_junc_gr,
  #                             subject = intron_coordinates,
  #                             ignore.strand = F,
  #                             type = "equal")
  #     
  #     if (donor_misspliced_junc_gr[queryHits(donor_x),] %>% length() == 0) {
  #       
  #       query = paste0("SELECT * FROM acceptor WHERE seqnames == ", intron_coordinates@seqnames , 
  #                      " AND start == ", intron_coordinates@ranges@start, 
  #                      " AND end == ", intron_coordinates %>% ranges %>% end, 
  #                      " AND strand == '", intron_coordinates@strand, "'")
  #       
  #       con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
  #       acceptor_misspliced_junc_gr <- dbGetQuery(con, query) 
  #       dbDisconnect(con)
  #       
  #       acceptor_x <- findOverlaps(query = acceptor_misspliced_junc_gr,
  #                                  subject = intron_coordinates,
  #                                  ignore.strand = F,
  #                                  type = "equal")
  #       
  #       if (acceptor_misspliced_junc_gr[queryHits(acceptor_x),] %>% length() == 0) {
  #         
  #         query = paste0("SELECT * FROM never WHERE seqnames == ", intron_coordinates@seqnames , 
  #                        " AND start == ", intron_coordinates@ranges@start, 
  #                        " AND end == ", intron_coordinates %>% ranges %>% end, 
  #                        " AND strand == '", intron_coordinates@strand, "'")
  #         
  #         con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
  #         never_misspliced_junctions_gr <- dbGetQuery(con, query)
  #         dbDisconnect(con)
  #         
  #         never_x <- findOverlaps(query = never_misspliced_junctions_gr,
  #                                 subject = intron_coordinates,
  #                                 ignore.strand = F,
  #                                 type = "equal")
  #         if (never_misspliced_junctions_gr[queryHits(never_x),] %>% length() == 0) {
  #           data.frame(Time = Sys.time(),
  #                      Message = paste0("Intron not found to be mis-spliced across all GTEx tissues. Please, try it again using one single tissue instead.")) %>% return()
  #           
  #         } else {
  #           df_all_gr <- never_misspliced_junc_gr[queryHits(never_x),] 
  #           rm(never_misspliced_junctions_gr)
  #           intronType <<- "never"
  #         }
  #         
  #       } else {
  #         df_all_gr <- acceptor_misspliced_junc_gr[queryHits(acceptor_x),] 
  #         rm(acceptor_misspliced_junc_gr)
  #         intronType <<- "acceptor"
  #       }
  #       
  #     } else {
  #       df_all_gr <- donor_misspliced_junc_gr[queryHits(donor_x),] 
  #       rm(donor_misspliced_junc_gr)
  #       intronType <<- "donor"
  #     }
  #     
  #   } else {
  #     df_all_gr <- both_misspliced_junc_gr[queryHits(both_x),] 
  #     rm(both_misspliced_junc_gr)
  #     intronType <<- "both"
  #   }
  #   
  #   intronID <<- df_all_gr$ref_junID %>% unique()
  #    
  #   df_all_gr %>%
  #     as.data.frame() %>%
  #     group_by(ref_junID) %>%
  #     mutate(MissplicingRatio_Donor = ref_missplicingratio_ND_tissues %>% mean(),
  #            MissplicingRatio_Acceptor = ref_missplicingratio_NA_tissues %>% mean(),
  #            Type = type) %>%
  #     ungroup() %>% 
  #     distinct(ref_junID, .keep_all = T) %>%
  #     mutate(MSR_D = formatC(x = MissplicingRatio_Donor, format = "e", digits = 3),
  #            MSR_A = formatC(x = MissplicingRatio_Acceptor, format = "e", digits = 3),
  #            MSR_D_Var = formatC(x = missplicingratio_ND_tissues_var, format = "e", digits = 3),
  #            MSR_A_Var = formatC(x = missplicingratio_NA_tissues_var, format = "e", digits = 3)) %>%
  #     select(Intron_ID = ref_junID,
  #            Type,
  #            Width = width, 
  #            Ss5score = ref_ss5score,
  #            Ss3score = ref_ss3score,
  #            MSR_D,
  #            MSR_A,
  #            MSR_D_Var,
  #            MSR_A_Var,
  #            Gene = gene_name_start) %>% return()
  #     
  #   
  # } else {
    
    # Create the query  
    query = paste0("SELECT * FROM '", paste0(tissue, "_db_introns"), "' WHERE seqnames == ", intron_coordinates@seqnames , 
                   " AND start == ", intron_coordinates@ranges@start, 
                   " AND end == ", intron_coordinates %>% ranges %>% end, 
                   " AND strand == '", intron_coordinates@strand, "'")
    
    
    # Execute the RSQLite query
    con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
    df_gr <- dbGetQuery(con, query) 
    dbDisconnect(con)
    
    
    if (df_gr %>% nrow() > 0) {
      
      intronID <<- df_gr$ref_junID
      intronType <<- df_gr$ref_type
      
      df_gr %>%
        as.data.frame() %>%
        distinct(ref_junID, .keep_all = T) %>%
        mutate(MeanCounts = round(x = ref_mean_counts, digits = 2),
               MSR_D = formatC(x = ref_missplicing_ratio_tissue_ND, format = "e", digits = 3),
               MSR_A = formatC(x = ref_missplicing_ratio_tissue_NA, format = "e", digits = 3)) %>%
        select(Intron_ID = ref_junID,
               Type = ref_type,
               Chr = seqnames,
               Start = start,
               End = end,
               Strand = strand,
               Width = ref_width, 
               Ss5score = ref_ss5score,
               Ss3score =ref_ss3score,
               MeanCounts,
               Individuals = ref_n_individuals,
               
               MSR_D,
               MSR_A,
               Gene = gene_name) %>% return()
      
      
      
    } else {
      intronID <<- NULL
      intronType <<- NULL
      data.frame(Time = Sys.time(),
                 Message = paste0("Intron not found in '", names(tissue_GTEx_choices_alphabetical)[tissue_GTEx_choices_alphabetical == tissue], "' tissue data.")) %>% return()
    }
    
  # }

}
# intron_ID=6028640
# coordinates=1:155235845-155236244:-
# gene=GBA
# type=acceptor
# clinvar=-
# length=400
# tissue="Adipose-Subcutaneous"
get_novel_data <- function(intron_ID = NULL,
                           tissue = NULL,
                           all_tissues = F) {
  

  all_people_tissue <- readRDS(file = "./dependencies/all_people_used_tissue.rda")[[tissue]] %>% length()
  if (!is.null(intronType) && intronType == "never") {
    
    data.frame(Time = Sys.time(),
               Message = paste0("Intron not found in '", 
                                names(tissue_GTEx_choices_alphabetical)[tissue_GTEx_choices_alphabetical == tissue], 
                                "' tissue data.")) %>% return()
    
  } else {
  
    # intron_ID <- "105665"
    # if (is.null(intron_ID)) {
    # Load core shared junctions across tissues files  
      query = paste0("SELECT * FROM '", paste0(tissue, "_db_novel"), "' WHERE ref_junID == '", intron_ID, "'")
    #} else {
    #  query = paste0("SELECT * FROM '", paste0(tissue, "_db_novel"), "' ")
    #}
    
    
    con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
    df_gr <- dbGetQuery(con, query) 
    dbDisconnect(con)
    
    
    if (df_gr %>% nrow() > 0) {
        
      df_gr %>%
        as.data.frame() %>%
        group_by(ref_junID, novel_junID) %>%
        distinct(novel_junID, .keep_all = T) %>%
        ungroup() %>%
        dplyr::mutate(MeanCounts = round(x = novel_sum_counts/novel_n_individuals, digits = 2),
                      "% Individuals" = ifelse(round((novel_n_individuals * 100) / all_people_tissue) == 0, 1,
                                               round((novel_n_individuals * 100) / all_people_tissue)),
                      #Mean_MSR = formatC(x = novel_missplicing_ratio_tissue, format = "e", digits = 3),
                      coordinates = paste0(seqnames,":",start,"-",end,":",strand)) %>%
        dplyr::select(Novel_ID = coordinates,
                      Type = novel_type,
                      recountID = novel_junID,
                      Width = width,
                      Ss5score = novel_ss5score,
                      Ss3score = novel_ss3score,
                      Distance = distance,
                      MeanCounts,
                      "% Individuals") %>% 
        return()
    } else {
      
      data.frame(Time = Sys.time(),
                 Message = paste0("Intron not found in '", names(tissue_GTEx_choices_alphabetical)[tissue_GTEx_choices_alphabetical == tissue], "' tissue data.")) %>% return()
    }
  }
}

# intron=6028640
# coordinates=1:155235845-155236244:-
# gene=GBA
# type=acceptor
# clinvar=-
# length=400
# tissue="Adipose-Subcutaneous"
get_intron_details <- function(intron_id = NULL,
                               tissue = NULL) {

  query = paste0("SELECT * FROM '", paste0(tissue, "_db_intron_details"), "' WHERE ref_junID == '", intron_id, "'")
  
  
  con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
  df_gr <- dbGetQuery(con, query) 
  dbDisconnect(con)
  
  
  if (df_gr %>% nrow() > 0) {
    
    df_gr %>%
      as.data.frame() %>%
      #dplyr::mutate(count = count %>% integer()) %>%
      dplyr::filter(!is.na(count), count > 0) %>%
      select(Intron_ID = ref_junID,
             Sample_ID = sample,
             Counts = count,
             Age = age,
             Sex = sex,
             Read_Count = mapped_read_count) %>% 
      return()
    
  } else {
    
    data.frame(Message = paste0("Intron not found in '", 
                                names(tissue_GTEx_choices_alphabetical)[tissue_GTEx_choices_alphabetical == tissue], 
                                "' tissue data.")) %>% return()
  }
  
  # } 
  
}

# novel_id = "67197814"
# tissue <- "control"
get_novel_details <- function(novel_id = NULL,
                              tissue = NULL) {
  
  
    
    # if (intronType == "never") {
    # 
    #   data.frame(Message = paste0("No novel junctions for intron '", intronID,"' in '",
    #                               names(tissue_GTEx_choices_alphabetical)[tissue_GTEx_choices_alphabetical == tissue], "' tissue data.")) %>% return()
    # 
    # } else {
      
      # ## Load core shared junctions across tissues files  
      query = paste0("SELECT * FROM '", paste0(tissue, "_db_novel_details"), "' WHERE novel_junID == '", novel_id, "'")
      
     
      con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
      df_gr <- dbGetQuery(con, query) 
      dbDisconnect(con)
      
      
      if (df_gr %>% nrow() > 0) {
        
        df_gr %>%
          as.data.frame() %>%
          #mutate(coordinates = paste0(seqnames,":",start,"-",end,":",strand)) %>%
          dplyr::filter(!is.na(novel_counts), novel_counts > 0) %>%
          select(Novel_ID = novel_junID,
                 Sample_ID = sample,
                 Counts = novel_counts,
                 Age = age,
                 Sex = sex,
                 Read_Count = mapped_read_count) %>% 
          return()
        
      } else {
        
        data.frame(Message = paste0("Intron not found in '", 
                                    names(tissue_GTEx_choices_alphabetical)[tissue_GTEx_choices_alphabetical == tissue], 
                                    "' tissue data.")) %>% return()
      }
      
    # } 
  
}


# gene_id <- "SNCA"
# gene_id <- "ENSG00000145335"
# tissue <- "Brain-FrontalCortex_BA9"
# mane <- TRUE
# enovel <- TRUE
# clinvar <- FALSE

get_gene_intron_data <- function(gene_id, tissue, mane, clinvar, enovel, threshold) {
  
  all_people_tissue <- readRDS(file = "./dependencies/all_people_used_tissue.rda")[[tissue]] %>% length()
  
  if (enovel) {
    ## TODO filter only by novel events present in more than XXX people
    query = paste0("SELECT * FROM '", paste0(tissue, "_db_novel"), "' WHERE gene_name == '", gene_id, "' AND novel_n_individuals >= ", round(x = (threshold * all_people_tissue)/100))
    # Create an ephemeral in-memory RSQLite database
    con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
    df_novel_gr <- dbGetQuery(con, query) 
    dbDisconnect(con)
  }
  
  if (mane) {
    query = paste0("SELECT * FROM '", paste0(tissue, "_db_introns"), "' WHERE gene_name == '", gene_id, "' AND MANE == ", mane)
  } else {
    query = paste0("SELECT * FROM '", paste0(tissue, "_db_introns"), "' WHERE gene_name == '", gene_id, "'")
  }
  
  if (clinvar) {
    query = paste0(query, " AND clinvar_type != '-'")
  }
  if (enovel) {
    query = paste0(query, " AND ref_junID IN ('", paste(df_novel_gr$ref_junID %>% unique, collapse = "','"),"')")
  }
  
  # Create an ephemeral in-memory RSQLite database
  con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
  df_gr <- dbGetQuery(con, query) 
  dbDisconnect(con)
  
  if (any(df_gr %>% names() == "ref_width")) {
    df_gr <- df_gr %>%
      dplyr::rename(width = ref_width)
  }
  
  df_gr %>%
    #rename(width = ref_junID) %>%
    mutate(MeanCounts = round(x = ref_mean_counts, digits = 2),
           MSR_D = formatC(x = ref_missplicing_ratio_tissue_ND, format = "e", digits = 3),
           MSR_A = formatC(x = ref_missplicing_ratio_tissue_NA, format = "e", digits = 3),
           MANE = "T",#ifelse(MANE == 0, "F", "T"),
           coordinates = paste0(seqnames,":",start,"-",end,":",strand),
           "% Individuals" = round(x = (ref_n_individuals * 100) / all_people_tissue)) %>%
    dplyr::select(ID = coordinates,
                  "Mis-spliced site" = ref_type,
                  recountID = ref_junID,
                  Width = width, 
                  Ss5score = ref_ss5score,
                  Ss3score = ref_ss3score,
                  MeanCounts,
                  "% Individuals",
                  MSR_D,
                  MSR_A,
                  ClinVar = clinvar_type,
                  MANE, 
                  Gene = gene_name) %>% return()
  
}


###################################################
############ PLOT FUNCTIONS #######################
###################################################


# tissue <- "PD"
# limit_bp = 30
plot_distances <- function(tissue,
                           limit_bp = 30) {
  
  
  # ## Load core shared junctions across tissues files  
  negative_bp <- limit_bp * (-1)
  query = paste0("SELECT * FROM '", paste0(tissue, "_db_novel"), 
                 "' WHERE (novel_type == 'novel_donor' OR novel_type == 'novel_acceptor') AND distance <= ", 
                 limit_bp," AND distance >= ", negative_bp)
  
  
  con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
  df_gr <- dbGetQuery(con, query) 
  dbDisconnect(con)
  

  ## Replace the underscore
  df_gr <- df_gr %>%
    mutate(novel_type = str_replace(string = novel_type, pattern = "_", replacement = " "))
  
  df_gr$novel_type = factor(df_gr$novel_type, 
                               levels = c("novel donor", "novel acceptor"))
  
  
  # y_axes <- c(0, (df_gr %>% 
  #                   filter(novel_type == "novel_acceptor", 
  #                          distance == (df_gr$distance %>% get_mode())) %>% 
  #                     nrow()) + 30)

  
  ggplot(data = df_gr) + 
    geom_histogram(aes(x = distance, group = novel_type, fill = novel_type), 
                   bins = limit_bp * 2) +
    
    facet_grid(vars(novel_type)) +
    xlab("distance (in bp)") +
    theme_light() +
    scale_x_continuous(limits = c((limit_bp * -1), limit_bp),
                       breaks = c((limit_bp * -1), (round(limit_bp / 2) * -1), 0, round(limit_bp / 2), limit_bp)) +
    scale_fill_manual(values = c("#35B779FF","#440154FF"),
                      breaks = c("novel donor", "novel acceptor"),
                      labels = c("novel donor", "novel acceptor")) +
    guides(fill = guide_legend(title = NULL, #title = "Junction category & Strand",
                               override.aes = list(size = 3),
                               ncol = 3 )) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "14"),
          axis.title = element_text(colour = "black", size = "14"),
          strip.text = element_text(colour = "black", size = "16"), 
          legend.text = element_text(colour = "black", size = "14"),
          plot.caption = element_text(colour = "black", size = "14"),
          plot.title = element_text(colour = "black", size = "16"),
          legend.title = element_text(colour = "black", size = "14"),
          legend.position = "top")  %>% return()
}


plot_modulo <- function(tissue,
                        limit_bp = 30) {
  
  # Query to the DB
  negative_bp <- limit_bp * (-1)
  query = paste0("SELECT * FROM '", paste0(tissue, "_db_novel"), 
                 "' WHERE (novel_type == 'novel_donor' OR novel_type == 'novel_acceptor') AND distance <= ", 
                 limit_bp," AND distance >= ", 
                 negative_bp)
  
  
  con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
  df_gr <- dbGetQuery(con, query) 
  dbDisconnect(con)
  
  
  df_gr$novel_type = factor(df_gr$novel_type, 
                            levels = c("novel_donor", "novel_acceptor"))
  
  df_gr <- df_gr %>%
    filter(distance >= limit_bp*(-1), distance <= limit_bp) %>%
    mutate(type_p = ifelse(distance < 0, paste0(novel_type, "_intron"), paste0(novel_type, "_exon"))) %>% 
    mutate(module = round(abs(distance) %% 3, digits = 2))
  
  
  ggplot(data = df_gr, aes(x = module, group = type_p, fill = type_p)) + 
    geom_bar(aes(y = ..prop..), stat = "count") +
    geom_text(aes( label = scales::percent(..prop.., accuracy = 0.1),
                   y= ..prop.. ), stat= "count", vjust = -.5) +
    scale_x_continuous(breaks = c(0, 1, 2)) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_viridis_d() +
    facet_grid(~type_p) + 
    ylab("percentage of novel junctions") +
    xlab("modulo 3") +
    theme_light() +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "14"),
          axis.title = element_text(colour = "black", size = "14"),
          legend.text = element_text(colour = "black",size = "14"),
          strip.text = element_text(colour = "black", size = "12"), 
          plot.caption = element_text(colour = "black",size = "14"),
          legend.title = element_text(colour = "black", size = "14"),
          legend.position = "none") +
    guides(fill = guide_legend(title = NULL,
                               ncol = 2, 
                               nrow = 2)) %>% return()
}


plot_missplicing <- function(tissue) {
  
  # Query to the DB
  query = paste0("SELECT * FROM '", paste0(tissue, "_db_introns"),"'")
  
  
  con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
  df_gr <- dbGetQuery(con, query) 
  dbDisconnect(con)
  
  df_gr$ref_type %>% unique %>% print()
  
  # df_gr$ref_type = factor(df_gr$ref_type, 
  #                           levels = c("donor", "acceptor"))
  
  # y_axes <- c(0, (df_gr %>% 
  #                   filter(ref_type == "never", 
  #                          ref_missplicing_ratio_tissue_NA > 0,
  #                          ref_missplicing_ratio_tissue_NA < 0.003) %>% 
  #                   nrow()) + 30)
  
  ggplot(df_gr) + 
    geom_density(aes(x = ref_missplicing_ratio_tissue_NA, fill = "#440154FF"), 
                   alpha = 0.8, bins = 60) +
    geom_density(aes(x = ref_missplicing_ratio_tissue_ND, fill = "#35B779FF"), 
                   alpha = 0.8, bins = 60) +
    xlab("Mis-splicing ratio (MSR)") +
    theme_light() +
    scale_fill_manual(values = c("#35B779FF","#440154FF"),
                      breaks = c("#35B779FF","#440154FF"),
                      labels = c("MSR_Donor","MSR_Acceptor")) +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "14"),
          axis.title = element_text(colour = "black", size = "14"),
          plot.title = element_text(colour = "black", size = "16"),
          legend.text = element_text(size = "14"),
          legend.title = element_text(size = "14"),
          legend.position = "top") +
    guides(fill = guide_legend(title = NULL,
                               ncol = 2, 
                               nrow = 1)) %>% return()

}




plot_intron_proportions <- function(clusters) {
  
  # clusters <- gtex_tissues
  # clusters <- c("PD", "control")
  # folder_root = "/home/sruiz/PROJECTS/splicing-project/splicing-recount2-projects/data/SRP058181/results/pipeline3/missplicing-ratio/"
  # folder_results = "/home/sruiz/PROJECTS/splicing-project/splicing-recount2-projects/data/SRP058181/results/pipeline3/alltypes/"

  folder_results <- "./dependencies/"
  
  
  #df_missplicing_all %>% head() %>% print()
  #df_missplicing_all <- readRDS(file = paste0(folder_results, "/df_", cluster, "_missplicing_all.rds"))
  
  
  if (any(clusters == "PD") || any(clusters == "control")) {
    df_missplicing_all <- readRDS(file = paste0(folder_results, "/SRP058181/df_missplicing_all.rds"))
  } else {
    df_missplicing_all <- readRDS(file = paste0(folder_results, "/df_missplicing_all.rds"))
  }
  
  
  ## Plot results --------------------------------------------------------------------
  
  df <- data.frame(data = df_missplicing_all$prop_both,
                   tissue = df_missplicing_all$tissue,
                   type = "both")
  df <- rbind(df, 
              data.frame(data = df_missplicing_all$prop_acceptor,
                         tissue = df_missplicing_all$tissue,
                         type = "acceptor"))
  
  df <- rbind(df, 
              data.frame(data = df_missplicing_all$prop_donor,
                         tissue = df_missplicing_all$tissue,
                         type = "donor"))
  
  
  df <- rbind(df, 
              data.frame(data = df_missplicing_all$prop_none,
                         tissue = df_missplicing_all$tissue,
                         type = "never"))
  
  
  df$tissue <- factor(df$tissue, levels = df$tissue[order(df %>% filter(type == "never") %>% pull(data) %>% dplyr::desc())] %>% unique)
  
  
  
  breaks <- levels(as.factor(df$tissue[order(df$data %>% dplyr::desc())] %>% unique))
  colours <- ifelse(str_detect(string = as.factor(df$tissue[order(df %>% filter(type == "never") %>% pull(data) %>% dplyr::desc())] %>% unique), pattern = "Brain"), 
                    "red", "black")
  
  if (any(clusters == "PD") || any(clusters == "control")) {
    label = round(df$data * 100, digits = 2)
    xlabel = "sample type"
  } else {
    label = ""
    xlabel = "tissue"
  }
  
  ggplot(data = df, aes(x = tissue, y = data, fill = type)) +
    geom_bar(stat = "identity")  +
    geom_text(aes(label = label), position=position_stack(0.5), color = "red") +
    
    theme_light() +
    ylab("proportion of introns") +
    xlab(xlabel) +
    scale_fill_viridis_d() +
    ggtitle("Proportion of intron type.\nOrdered by the 'never mis-spliced' intron category.") +
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black", size = "12"),
          axis.text.x = element_text(angle = 70, 
                                     vjust = 1,
                                     color = colours,
                                     hjust = 1),
          axis.title = element_text(colour = "black", size = "12"),
          legend.text = element_text(size = "12"),
          legend.title = element_text(size = "12"),
          legend.position = "top") +
    guides(fill = guide_legend(title = "Junction type: ",
                               ncol = 2, 
                               nrow = 2)) %>% return()
  
  
  # if (save_result) {
  #   file_name <- paste0(folder_results, "/proportion_junction_tissue_never_ordered.png")
  #   ggplot2::ggsave(filename = file_name,
  #                  width = 183, height = 183, units = "mm", dpi = 300)
  # }
  
  
}

# tissue <- "Brain-FrontalCortex_BA9"
# tissue <- "Brain-Substantianigra"
# tissue <- "control_SRP058181"
# tissue <- "PD_SRP058181"
# tissue <- "control_SRP049203"
# tissue <- "PD_SRP049203"
plot_lm <- function(tissue) {
  
  if (str_detect(string = tissue,
                 pattern = "PD") || str_detect(string = tissue,
                                               pattern = "control")) {
    
    # Query to the DB
    query = paste0("SELECT * FROM '", paste0(tissue, "_db_lm"),"'")
    
    
    con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
    df_stats <- dbGetQuery(con, query) 
    dbDisconnect(con)
    
    # if (str_detect(string = tissue, pattern = "PD")) {
    #   df_stats <- df_stats %>%
    #     mutate(disease_state = disease_state %>% as.factor()) %>%
    #     mutate(disease_state = relevel(disease_state, ref = "control"))
    # } else {
    #   df_stats <- df_stats %>%
    #     mutate(disease_state = disease_state %>% as.factor()) %>%
    #     mutate(disease_state = relevel(disease_state, ref = "PD"))
    # }
    # 
    df_stats[1,]
    
    df_stats <- df_stats %>%
      filter(u2_intron == T | u12_intron == T) %>% 
      distinct(ref_junID, .keep_all = T) %>%
      dplyr::rename(is_PD = disease_state,
                    gene_tpm = tpm) %>%
      mutate(is_control = is_PD) %>%
      mutate(is_PD = ifelse(is_PD == "PD", T, F)) %>%
      mutate(is_control = ifelse(is_control == "control", T, F))

    df_stats %>% nrow()
    df_stats %>% distinct(ref_junID) %>% nrow()
    df_stats[1,]
    
    fit_donor <- lm(ref_missplicing_ratio_tissue_ND ~ intron_length + 
                      intron_5ss_score + 
                      intron_3ss_score  +
                      disease_state +
                      protein_coding +
                      #is_control+
                      gene_tpm + 
                      #u12_intron +
                      u2_intron +
                      gene_length + 
                      gene_num_transcripts, 
                    #u2_intron ,
                    data = df_stats)
    
    fit_acceptor <- lm(ref_missplicing_ratio_tissue_NA ~ intron_length + 
                         intron_5ss_score * 
                         intron_3ss_score +
                         disease_state +
                         #is_control+
                         #u12_intron +
                         u2_intron +
                         #gene_tpm + 
                         gene_length + 
                         #protein_coding +
                         gene_num_transcripts,
                       #protein_coding,
                       data = df_stats)
    
    fit_donor %>% summary()
    fit_acceptor %>% summary()
    
  } else {
    
    # Query to the DB
    query = paste0("SELECT * FROM '", paste0(tissue, "_db_introns"),"'")
    
    
    con <- dbConnect(RSQLite::SQLite(), "./dependencies/splicing.sqlite")
    df_gr <- dbGetQuery(con, query) 
    dbDisconnect(con)
    
    ## RENAME COLUMNS
    df_stats <- df_gr %>%
      dplyr::rename(intron_length = width,
                    intron_5ss_score = ref_ss5score,
                    intron_3ss_score = ref_ss3score,
                    gene_length = gene_width,
                    gene_tpm = tpm,
                    gene_num_transcripts = transcript_number) %>%
      filter(gene_tpm > 0) %>%
      filter(u2_intron == T | u12_intron == T)
    
    
    
    
    fit_donor <- lm(ref_missplicing_ratio_tissue_ND ~ intron_length + 
                      intron_5ss_score + intron_3ss_score + 
                      gene_tpm + 
                      gene_length + 
                      gene_num_transcripts + 
                      #u12_intron +
                      u2_intron +
                      protein_coding,
                    data = df_stats)
    
    fit_acceptor <- lm(ref_missplicing_ratio_tissue_NA ~ intron_length + 
                         intron_5ss_score + intron_3ss_score +
                         gene_tpm + 
                         gene_length + 
                         gene_num_transcripts + 
                         #u12_intron +
                         u2_intron +
                         protein_coding,
                       data = df_stats)
  }
  
  fit_donor %>% summary()
  fit_acceptor %>% summary()
  
  
  
  ##########################
  ## PLOT LINEAR MODELS
  ##########################
  
  jtools::plot_summs(fit_donor, 
                     fit_acceptor,
                     scale = TRUE, 
                     robust = list("HC3", "HC3"),
                     #inner_ci_level = .75,
                     n.sd = 2,
                     legend.title = "Model:",
                     #plot.distributions = TRUE,
                     ci_level = 0.95,
                     colors = c("#35B779FF","#440154FF"),
                     model.names = c("MSR_Donor", "MSR_Acceptor")) + 
    
    theme_minimal() + 
    theme(axis.line = element_line(colour = "black"), 
          axis.text.x = element_text(colour = "black", size = "13"),
          axis.text.y = element_text(colour = "black", size = "13"),
          axis.title = element_text(colour = "black", size = "13"),
          legend.text = element_text(colour = "black", size = "13"),
          legend.title = element_text(colour = "black", size = "13"),
          legend.position = "top",
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.ticks.x = element_line(colour = "black", size = 2)) +  
    ylab("Predictors") +
    geom_hline(yintercept = seq(from = -0.5, to = (length((fit_donor$coefficients %>% names)[-1]) + .5)))  %>% return()
  
  
}

###################################################
################# VARIABLES #######################
###################################################



tissue_GTEx_choices <- c(
  "Adipose - subcutaneous" =	"Adipose-Subcutaneous",
  "Adipose - visceral" =	"Adipose-Visceral_Omentum",
  "Adrenal gland" =	"AdrenalGland",
  "Aorta" =	"Artery-Aorta",
  "Artery - coronary" =	"Artery-Coronary",
  "Artery - tibial" =	"Artery-Tibial",
  "Brain Amygdala" =	"Brain-Amygdala",
  "Brain Anterior cingulate cortex" =	"Brain-Anteriorcingulatecortex_BA24",
  "Brain Caudate" =	"Brain-Caudate_basalganglia",
  "Brain Cerebellar hemisphere" =	"Brain-CerebellarHemisphere",
  "Brain Frontal Cortex" =	"Brain-FrontalCortex_BA9",
  "Brain Substantia nigra" = "Brain-Substantianigra",
  "Brain Hippocampus" =	"Brain-Hippocampus",
#"SRP049203 - PD" =	"PD_SRP049203",
#"SRP049203 - Control" = "control_SRP049203",
#"SRP058181 - PD" =	"PD_SRP058181",
#"SRP058181 - Control" = "control_SRP058181"
"Brain Hypothalamus"	= "Brain-Hypothalamus",
"Brain Nucleus accumbens" =	"Brain-Nucleusaccumbens_basalganglia",
"Brain Putamen" =	"Brain-Putamen_basalganglia",
"Brain Spinal cord" =	"Brain-Spinalcord_cervicalc-1",
"Brain Substantia nigra" = "Brain-Substantianigra",
"Lymphocytes" = "Cells-EBV-transformedlymphocytes",
"Fibroblasts" = "Cells-Transformedfibroblasts",
"Colon Sigmoid" =	"Colon-Sigmoid",
"Colon Transverse" =	"Colon-Transverse",
"Gastroesophageal junction" =	"Esophagus-GastroesophagealJunction",
"Mucosa" =	"Esophagus-Mucosa",
"Muscularis" =	"Esophagus-Muscularis",
"Atrial appendage" =	"Heart-AtrialAppendage",
"Left ventricle" =	"Heart-LeftVentricle",
"Liver" =	"Liver",
"Lung" =	"Lung",
"Minor salivary gland" =	"MinorSalivaryGland",
"Skeletal muscle" =	"Muscle-Skeletal",
"Nerve - tibial" =	"Nerve-Tibial",
"Pancreas" =	"Pancreas",
"Pituitary" =	"Pituitary",
"Skin (suprapubic)" =	"Skin-NotSunExposed_Suprapubic",
"Skin (lower leg)"	= "Skin-SunExposed_Lowerleg",
"Small Intestine" =	"SmallIntestine-TerminalIleum",
"Spleen" =	"Spleen",
"Stomach" =	"Stomach",
"Thyroid" =	"Thyroid",
"Whole blood" =	"WholeBlood")

tissue_GTEx_choices_alphabetical <- tissue_GTEx_choices[names(tissue_GTEx_choices) %>% order()]

chr_choices <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y")
strand_choices <- c("+", "-")


get_mode <- function(vector) {
  uniqv <- unique(vector)
  uniqv[which.max(tabulate(match(vector, uniqv)))]
}

# setwd("./vizSI/")
# genes <- readRDS(file = "./dependencies/all_genes.rds")
# saveRDS(object = gene_names$hgnc_symbol,
#         file = "./dependencies/all_genes_names.rds")
genes <- readRDS(file = "./dependencies/all_genes_names.rds")


intronID <- NULL
intronType <- NULL

# tissue <- "Brain-FrontalCortex_BA9"
# chr <- 10
# start <- 133366994
# end <- 133368922
# strand <- "+"
# intron_coordinates <- GRanges(seqnames = chr,
#                               ranges = IRanges(start = start,
#                                                end = end),
#                               strand = strand)
# 
#  
# 
# DNAJA3_junctions <- read.csv("/home/dzhang/projects/RNA_seq_diag_mito/tmp/DNAJA3_junctions.csv")
# DNAJA3_junctions
# ECHS1_junctions <- read.csv("/home/dzhang/projects/RNA_seq_diag_mito/tmp/ECHS1_junctions.csv")
# ECHS1_junctions
# # ECHS1_junctions@ranges
# # 10:133366994-133368922:-
# 
# # DNAJA3_junctions@ranges
# # DNAJA3_mcols
# # 16:4434518-4437401:+
# 
