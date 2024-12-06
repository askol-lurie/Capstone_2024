load_data_excel <- function(file){
  
  ## values from Andy via Alamute Visual.
  hd <- c(292, 471)
  exons <- c(241, 429, 945)   
  
  data <- readxl::read_excel(file, na=c("","NA"), sheet="Analysis Dataset", col_names = TRUE)
  
  data <- data %>%
    dplyr::rename(PathAnnot = `Pathogenicity Annotation`) %>%
    dplyr::rename(g. = `g. (chr 4, GRCh38)`) %>%
    dplyr::select(-`REVEL Score (UCSC)`) %>%
    dplyr::mutate(PathAnnot = ifelse(PathAnnot == "Disease-causing", "Pathogenic", PathAnnot),
            PathAnnot =
              factor(PathAnnot, levels=c("Benign","Uncertain","Pathogenic")),
            dis_ben = ifelse(PathAnnot == "Uncertain", NA, 1 * (PathAnnot == "Pathogenic")),
            dis_unc = ifelse(PathAnnot == "Benign", NA, 1 * (PathAnnot == "Pathogenic")))
  
  ## pretty up output variables
  names(data) <- gsub(" Score.*", "", names(data))
  names(data) <- tolower(names(data))
  
  ## create genomic and aa positions
  ## and IMPUTE exons for missing values:
  ## exon 1 (1 - 80); exon 2 (81 - 143); exon 3 (>= 144 )
  data <- data %>% 
    dplyr::mutate(pos_g = as.numeric(gsub("[ACGT].+", "", g.)), 
                  pos_c = as.numeric(gsub("[ACGT_].+", "", gsub("c.", "", c.))),
                  homeo = 1 * (pos_c >= 98 & pos_c <= 157))
  
  ## get ref and alt alleles 
  data <- data %>%
    dplyr::mutate(ref = case_when(
        grepl("A>", c.) ~ "T",
        grepl("C>", c.) ~ "G",
        grepl("G>", c.) ~ "C",
        grepl("T>", c.) ~ "A", 
        .default = "."), 
        alt = case_when(
          grepl(">A", c.) ~ "T",
          grepl(">C", c.) ~ "G",
          grepl(">G", c.) ~ "C",
          grepl(">T", c.) ~ "A", 
          .default = "."))
  
  ## Determine exon ##
   data <- data %>% dplyr::mutate(
    exon_inf = NA, 
    exon_inf = ifelse(pos_g >= 41748369 & pos_g <= 41748610, 1, exon_inf),
    exon_inf = ifelse(pos_g >= 41747348 & pos_g <= 41747536, 2, exon_inf), 
    exon_inf = ifelse(pos_g >= 41745806 & pos_g <= 41746322, 3, exon_inf),
    exon = ifelse(is.na(exon), exon_inf, exon) )
   
  ## infer missing exons
  exon_max <- data %>% group_by(exon) %>% select(pos_g, pos_c) %>% 
    summarize(min_posg = min(pos_g), max_posg = max(pos_g),
              min_posc = min(pos_c), max_posc = max(pos_c)) 
  
  data <- data %>% dplyr::mutate(exon_inf = case_when(
      pos_c >= exon_max$min_posc[1] & pos_c <= exon_max$max_posc[1] ~ exon_max$exon[1],
      pos_c >= exon_max$min_posc[2] & pos_c <= exon_max$max_posc[2] ~ exon_max$exon[2],
      pos_c >= exon_max$min_posc[3] & pos_c <= exon_max$max_posc[3] ~ exon_max$exon[3]), 
      exon = ifelse(is.na(exon), exon_inf, exon))
  
  model <- lm(data = data, pos_g ~ pos_c + as.factor(exon))
  data$pred_pos <- predict(model, newdata = data, type = "response")
  data <- data %>% dplyr::mutate(pos_g = ifelse(is.na(pos_g), pred_pos, pos_g))
  
  ## Choose output variables
  data <- data %>% dplyr::select(c., p., g., pos_g, pos_c, ref, alt, cadd, revel, 
                                 alphamissense, bayesdel, pathannot, dis_ben, dis_unc, exon, homeo)
  
  ## note that this was part of the original data
  data$xlsx <- 1

  ## get bayesdel and alphamissense data
  bayesdel <- read_excel(file, na=c("","NA"), sheet="BayesDel", col_names = TRUE)
  alphamissense <- read_excel(file, na=c("","NA"), sheet="AlphaMissense", col_names = TRUE)
  all_snv_pred <- read_excel(file, na=c("","NA"), sheet="All Possible SNV Predictions", col_names = TRUE)
  
  cadd <- read_excel(file, na=c("","NA"), sheet="CADD", col_names = TRUE)
  all_snv_pred <- left_join(all_snv_pred, cadd, by=c("chr", "GRCh38 position" = "pos", 
                                                     "Nt Ref" = "ref", 
                                                     "Nt Alt" = "alt"))

  return(list(data = data, bayesdel=bayesdel, alphamissense=alphamissense, 
              all_snv_pred = all_snv_pred))
}


get_gnomad_data <- function(gnomad_ex, gnomad_ge, gnomad_faf){
  
  
  ## get new gnomad data ##
  gnomad_ex <- fread(gnomad_ex_file, header=T, sep="\t")
  gnomad_ge <- fread(gnomad_ge_file, header=T, sep="\t")
  gnomad_faf <- fread(gnomad_faf_file, header=F, sep="\t")
  
  ## add columns and determine usability: add my_filter for filtering ##
  gnomad_faf <- adjust_gnomad_pass(gnomad_faf)
  
  gnomad <- combine_gnomad(gnomad_ex, gnomad_ge, gnomad_faf, p_path=prop_path )
  
  return(gnomad)
}

adjust_gnomad_pass <- function(gnomad_faf){
  
  names(gnomad_faf) <- c("chr", "pos", "ref", "alt", "filter", "faf95", "faf95max", 
                         "faf95max_group", "p_cmh", "filt_exon", "filt_genom", "ac_exons", "an_exons",
                         "ac_genom", "an_genome") 
  
  ## keep only exonic / non-utr variants ##
  gnomad_faf <- gnomad_faf %>% 
    dplyr::mutate(keep = 0,
                  keep = ifelse(pos >= 41748369 & pos <= 41748610, 1, keep),
                  keep = ifelse(pos >= 41747348 & pos <= 41747536, 1, keep), 
                  keep = ifelse(pos >= 41745806 & pos <= 41746322, 1, keep), 
                  ## only snvs
                  keep = ifelse(nchar(ref) > 1 | nchar(alt) > 1, 0, keep))

  
  ## change char to numeric ##
  gnomad_faf <- gnomad_faf %>%
    dplyr::mutate(faf95 = as.numeric(faf95),
                  faf95max = as.numeric(faf95max),
                  p_cmh = as.numeric(p_cmh), 
                  ac_exons = as.numeric(ac_exons)) 
  
  ## if faf95max is na, use faf95 if ac_genom or ac_exons > 0
  gnomad_faf <- gnomad_faf %>%
    dplyr::mutate(faf95max = ifelse(is.na(faf95max) & (ac_exons+ac_genom > 0), faf95, faf95max))
  
  ## create alternative filter ##
  gnomad_faf <- gnomad_faf %>%
    dplyr::mutate(my_filter = filter, 
                  ac_exons = ifelse(ac_exons == ".", NA, ac_exons), 
                  ac_exons = as.numeric(ac_exons),
                  my_filter = ifelse(filter == "EXOMES_FILTERED" & ac_genom > 0, "PASS", filter),
                  my_filter = ifelse(filter == "GENOMES_LITERED" & ac_exons > 0, "PASS", my_filter)) 

  return(gnomad_faf)
}

combine_gnomad <- function(gnomad_ex, gnomad_ge, gnomad_faf, faf_cutoff=0.000001, p_path){
  
    gex <- gnomad_ex %>% dplyr::mutate(id = paste(chr, pos, ref, alt, sep="_")) %>%
      dplyr::select(id)
      
    gnomad_ge <- gnomad_ge %>% dplyr::mutate(id = paste(chr, pos, ref, alt, sep="_"))
    gnomad_ge_only <- anti_join(gnomad_ge, gex, by="id") %>% dplyr::select(-id) %>%
      dplyr::mutate(source="genome")
    
    gnomad_annot <- bind_rows(gnomad_ex, gnomad_ge_only)
    gnomad <- left_join(gnomad_faf, gnomad_annot, by=c("chr", "pos", "ref", "alt"))

    gnomad <- gnomad %>% dplyr::filter(!is.na(exon)) %>% 
      dplyr::mutate(pathannot_faf = ifelse(faf95max <= faf_cutoff, "Uncertain", "Benign"))
    gnomad$gnomad = 1
    
    ## my_filter = ifelse(pathannot == "Disease-causeing", "PASS", my_filter)
    
    ## create two gnomad data set
    ## 1) All coding exonic SNV variants with pass filters (used to determine proportion of SNVs,
    ##   that are being excluded for not being missense)
    ## 2) All coding exonic SNV variants that are missense (to join with data so we can understand, 
    ##   which variants might need to be removed because of bad filter)
    gnomad_coding_snv_pass <- gnomad %>% 
      dplyr::filter(keep == 1 & my_filter == "PASS" )
    gnomad_coding_snv_missense <- gnomad %>%
      dplyr::filter(keep == 1 & consequence == "missense_variant")
    
    conseq <- gnomad_coding_snv_pass %>%
      dplyr::select(consequence) %>% unlist()
    n_mis = sum(conseq == "missense_variant")
    n_syn = sum(grepl("synonymous", conseq))
    prop_missense <- n_mis/(n_mis + n_syn)
    prob_path_update = p_path * (n_mis + n_syn) / n_mis
                    
    
    return(list(gnomad_coding_snv_pass = gnomad_coding_snv_pass, 
                gnomad_coding_snv_missense = gnomad_coding_snv_missense, 
                prob_path_update = prob_path_update))
}

add_annot <- function(d, bayesdel, alphamissense){
  
  # d <- d %>% dplyr::select(-alphamissense, -bayesdel)
  b <- bayesdel %>% dplyr::select(pos, ref, alt, BayesDel_noAF_score) %>%
    dplyr::rename(bayesdel_new = BayesDel_noAF_score)
  a <- alphamissense %>% dplyr::select(POS, REF, ALT, am_pathogenicity) %>%
    dplyr::rename(alphamissense_new = am_pathogenicity, 
                  pos = POS, 
                  alt = ALT, 
                  ref = REF)
  
  d <- left_join(d, b, by=c("pos", "ref", "alt"))
  d <- left_join(d, a, by=c("pos", "ref", "alt"))
  
  ## this annotation should already exist for the variants in the excel sheet, so it 
  ## only needs to be updated for the gnomad variants
  d <- d %>% dplyr::mutate(revel.y = as.numeric(revel.y))
  d <- d %>% dplyr::mutate(bayesdel = ifelse(gnomad==1 & is.na(xlsx), bayesdel_new, bayesdel),
                           alphamissense = ifelse(gnomad==1 & is.na(xlsx), alphamissense_new, alphamissense),
                           cadd = ifelse(gnomad==1 & is.na(xlsx), cadd_phred, cadd), 
                           revel = ifelse(is.na(revel.y), revel.x, revel.y),
                           pos_g = ifelse(gnomad==1 & is.na(xlsx), pos, pos_g), 
                           exon = ifelse(is.na(exon.x), exon.y, exon.x))
  
  ## add position to pos for xlsx variants
  d <- d %>% dplyr::mutate(pos = ifelse(xlsx == 1 & is.na(gnomad), pos_g, pos))

  return(d)
}

update_pathannot <- function(d){
  
  d <- d %>% 
    dplyr::mutate(pathannot_update = case_when(
                  is.na(pathannot) ~ pathannot_faf ,
                  !is.na(pathannot) & !is.na(pathannot_faf) & pathannot == "Pathogenic" ~ pathannot,
                  !is.na(pathannot) & !is.na(pathannot_faf) & pathannot != "Pathogenic" ~ pathannot_faf,
                  is.na(pathannot_faf) ~ pathannot))
  
  ## update dis_ben
  d <- d %>%
    dplyr::mutate(dis_ben_orig = dis_ben,
                  dis_ben = NA,
                  dis_ben = case_when(
                  pathannot_update == "Pathogenic" ~ 1,
                  pathannot_update == "Benign" ~ 0, 
                  .default = NA))
          
  return(d)
                  
}


remove_dupe_pdots <- function(d){
  
  d <- d %>% group_by(p.) %>%
    dplyr::mutate(pathannot_update = ifelse(any(pathannot_update == "Pathogenic"), 
                                                "Pathogenic", pathannot_update),
                  pathannot_update = ifelse(all(pathannot_update != "Pathogenic") &
                                              any(pathannot_update == "Benign"), 
                                                  "Benign", pathannot_update),
                  alphamissense = mean(alphamissense),
                  revel = mean(revel),
                  cadd = mean(cadd),
                  bayesdel = mean(bayesdel), 
                  n = row_number()) %>%
    ungroup()
  
  cat("Collapsing", sum(d$n==2)*2 + sum(d$n>2), " variants into ", sum(d$n == 2), 
      "unique p. variants\n\n")
  d <- d %>% dplyr::filter(n == 1)
  
  return(d)
}                                 
                                            

make_table1 <- function(data, tbl_file){
  
    mns_raw <- data %>% group_by(pathannot) %>%
        summarize(n=n(), across(CADD:BayesDel, list(mean =~ mean(.x), sd =~ sd(.x))))
    
    ## perform t-test for each
    tbl <- c()
    for (cat in c("Benign", "Uncertain")){
        for (stat in c("CADD","REVEL", "BayesDel", "AlphaMissense")){
            keep = c("pathannot", stat)
            d <- data %>% dplyr::select(all_of(keep)) %>% 
                dplyr::rename(st := !!stat ) %>%
                dplyr::filter(pathannot %in% c(cat, "Pathogenic"))
            t <- t.test(st~pathannot, data=d)
            tbl <- rbind(tbl, c(cat, stat, t$p.value))
        }
    }
    
    tbl <- as.data.frame(tbl, stringsAsFactors=FALSE)
    names(tbl) <- c("pathannot", "annot", "pval")
    tbl$pval <- as.numeric(as.character(tbl$pval))
    tbl <- tbl %>% tidyr::pivot_wider(names_from="annot", values_from="pval",
                                       names_glue = "{annot}_pval")
                    
    mns <- left_join(mns_raw, tbl, by="pathannot")
    mns <- mns %>% dplyr::select(pathannot, starts_with("CADD"), starts_with("REVEL"),
                                 starts_with("BayesDel"), starts_with("AlphaMissense"))

    ## Reshaping the data
    table1 <- mns %>% dplyr::filter(pathannot != "Uncertain") %>%
        dplyr::select(-ends_with("sd")) %>%
        tidyr::pivot_longer(
                   cols = -c(pathannot), 
                   names_to = c("category", ".value"), 
                   names_pattern = "(.+)_(mean|pval)" # Regex pattern to separate the original column names into 'category' and 'mean/sd'
               )

    table1 <- table1 %>%
        tidyr::pivot_wider(names_from="pathannot", values_from=c("mean","pval"),
                           names_glue = "{pathannot}_{.value}") %>%
        dplyr::select(category, starts_with("Benign"), starts_with("Pathogenic"))
    table1 <- table1 %>%
          dplyr::mutate(Benign_mean = signif(Benign_mean, digits=2), 
                    Pathogenic_mean = signif(`Pathogenic_mean`, digits=2)) %>%
          dplyr::select(category, Benign_mean, Pathogenic_mean, Benign_pval)

    
    write.table(file = tbl_file, table1, quote=F, row.names=F, col.names=T, sep="\t")
    cat("Wrote results for table 1 to", tbl_file, "\n")
    
    return(table1 = table1)
}


homeo_summary <- function(data, file, file_supp){
  
  ## data is all cadd data in the phox2b gene
  homeodomain <- c(41747486, 41746281)
  
  
  ## Determine exon ##
  data <- data %>% 
    dplyr::rename(pos = `GRCh38 position`) %>%
    dplyr::mutate(
    keep = 0, 
    keep = ifelse(pos >= 41748369 & pos <= 41748610, 1, keep),
    keep = ifelse(pos >= 41747348 & pos <= 41747536, 1, keep), 
    keep = ifelse(pos >= 41745806 & pos <= 41746322, 1, keep), 
    hd = ifelse(pos >= homeodomain[2] & pos <= homeodomain[1], 1, 0))
  
  ## remove aa change equivalent variants
  d <- data %>% dplyr::filter(keep == 1) %>%
    dplyr::rename(ref = `Nt Ref`,
                  alt= `Nt Alt`,
                  p. = p.nomenclature) %>%
    dplyr::select(pos, ref, alt, p., hd, cadd, REVEL, BayesDel, AlphaMissense) %>%
    group_by(pos, ref, alt, p.) %>% 
      dplyr::mutate(n = row_number()) %>%
      dplyr::filter(n == 1) %>% ungroup() %>%
    dplyr::select(hd, cadd:AlphaMissense)
  
  
  write.table(file = file_supp, d, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
  cat("wrote", file_supp, "\n\n")
  
  d <- d %>% tidyr::pivot_longer(cols = cadd:AlphaMissense, names_to = "annot",
                                 values_to = "value")
  
  tbl_means <-  d %>% group_by(annot, hd) %>% dplyr::summarize(mean = mean(value, na.rm=T), n = n()) %>%
    tidyr::pivot_wider(names_from=c("annot"), values_from=c("mean")) %>%
    dplyr::select(hd, n, cadd, REVEL, BayesDel, AlphaMissense)
  
  write.table(file = file, tbl_means, quote=F, row.names=F, col.names=T, sep="\t")
  cat("Wrote", file, "\n\n")
  
  
}

plot_roc_imputed <- function(d, plotfile = "ROC_curves.pdf"){
  
  pdf(plotfile)
  
  colors <- c('CADD' = 'darkgreen', 'REVEL' = 'darkorange',
              'BayesDel' = 'lightblue', 'AlphaMissense' = 'red')
  annots = names(colors)
  
  p <- ggplot(data=d, aes(d=dis_ben, m=score, group=annot)) +
    plotROC::geom_roc(aes(color = annot), labels=F, size=.8, pointsize=0)
  
  p <- p + 
    style_roc() + 
    ggtitle("Pathogenic vs. Benign (Imputed)") + 
    scale_color_manual(name = 'Annotation',
                       breaks = annots,
                       values = colors[annots]) +
    theme_minimal() +
    theme(legend.title = element_text(size = 20),
          legend.text = element_text(size = 16),
          axis.text = element_text(size=12), 
          axis.title = element_text(size = 12))
   
  
    print(p)
  
  dev.off()
  
  cat("Wrote ROC plot to file", plotfile, "\n")
  
}
 


acmg_table <- function(data, regression_results, all_combinations, adjust, 
                       table_file="acmg_table_weighted.txt"){

    rslts <- c()
    vus_data <- data %>% dplyr::filter(dis_unc == 0)
    cutoffs <- c(2.406, 5.790, 33.53, 1124)
    cutlabel <- c("supporting", "moderate", "strong", "very_strong")
    pD = 0.044
    odds_path = pD / (1 - pD)

    
    
    min_func <- function(x, cutoff, model, odds_path, var_type ){

        if (! var_type %in% c("path","benign")){
            cat("var_type must be path or benign\n")
            stop()
        }
        
        ## calc se of log odds ##
        var_adj <- ifelse(var_type == "path", 1, -1)
        beta <- model$coefficients # Regression coefficients
        V <- vcov(model) # Variance-covariance matrix
        se_log_odds <- sqrt(t(c(1, x)) %*% V %*% c(1, x)) # Standard error
        z <- qnorm(.95)
        logodds_bound <- (beta[1] + beta[2]*x) - var_adj * z * se_log_odds
        odds_bound <- exp( var_adj * logodds_bound)
        odds_path = odds_path^var_adj

        return(cutoff - odds_bound / odds_path)
    }
    
    for (i in 1:4){
        
        glm.fit <- regression_results[[i]]
        cat = all_combinations[[i]]
        int <- glm.fit$coefficients[1] + adjust
        lor <- glm.fit$coefficients[2] 

        rng <- range(data[, cat])
        lower = -50; upper = 50
        if (i == 4){ lower = 0; upper=1} ## alphamissense
        
        for (j in 1:length(cutoffs)){
            cutoff <- cutoffs[j]
            
            #minimize_func_path <- function(x){ return( cutoff - exp(int + x*lor)/odds_path) }
            #minimize_func_ben <- function(x){ return( cutoff -  exp(-(int + x*lor))*odds_path) }

            result <- tryCatch({
              uniroot(min_func, lower = lower, upper = upper, cutoff = cutoff,
                              model = glm.fit, odds_path = odds_path, var_type = "path")
            }, error = function(e){
              cat("No satisfying threshold was able to be found for:\n")
              cat(cat, " using cutoff: ", cutoff, " for ",var_type, " variant\n")
              return(list(root = NA))
            })
            rslts <- rbind(rslts, c(cat, int, lor, "pathogenic", cutlabel[j], result$root, rng))

            result <- tryCatch({
              uniroot(min_func, lower = lower, upper = upper, cutoff = cutoff,
                              model = glm.fit, odds_path = odds_path, var_type = "benign")
            }, error = function(e){
              cat("No satisfying threshold was able to be found for:\n")
              cat(cat, " using cutoff: ", cutoff, " for ",var_type, " variant\n")
              return(list(root = NA))
            })
            rslts <- rbind(rslts, c(cat, int, lor, "benign", cutlabel[j], result$root, rng))
        }       
    }

    rslts <- as.data.frame(rslts, stringsAsFactors=F)
    names(rslts) <- c("annot", "int", "slope", "var_type", "acmg_cat", "annot_value",
                      "min_obs", "max_obs")

    tbl <- rslts %>% tidyr::pivot_wider(names_from=c("var_type","acmg_cat"),
                                        values_from="annot_value",
                                        names_glue = "{var_type}_{acmg_cat}") %>%
        dplyr::select(annot, int, slope, starts_with("pathogenic"), starts_with("benign"),
                      min_obs, max_obs)
    
    
    write.table(file = table_file, tbl, quote=F, row.names=F, col.names=T, sep="\t")
    cat("Wrote threshold for different ACMG support categories to", file, "\n")
    
    return(tbl)
}
        
plot_scores <- function(data, tbls, tick_marks, prefix = "plot_annot_wcutoffs"){
  
  plot_func <- function(d, disease = FALSE, file, pt=2.5){
    
    p <- ggplot(d, aes(x = factor(1), y = score, color = pathannot, shape = factor(exon)))
    
    p <- p + geom_point(position = position_jitter(width = 0.2, height=0.1), 
                        aes(group = exon), size=pt, alpha=.4) +
      facet_wrap(~factor(annot), scales = "free_y") 
    if (disease) {
      p <- p + scale_color_manual(values = c("Benign" = "darkgreen", "Pathogenic"="darkred"))
    }else{
      p <- p + scale_color_manual(values = c("Uncertain" = "darkblue"))
    }
    p <- p + scale_shape_manual(values = c("1" = 15, "2" = 16, "3" = 17, "?" = 18)) + 
      theme_bw() + 
      geom_text(data = tick_marks, aes(x = .5, y = score, label = Label, color=NULL),
                vjust = .5, hjust = 0, show.legend=FALSE, parse=T, ) +
      geom_segment(data = tick_marks, 
                   aes(x = .7, xend = 1.3, y = score, yend = score, color=NULL),
                       color="grey", linewidth = .5) + 
      labs(x = NULL, shape = "Exon Number") # Remove x-axis label, Change title for shape legend
      p <- p + labs(color = "Variant Class")  # Change title for color legend)
    p <- p + guides(color = guide_legend(override.aes = list(shape = 16 ))) + # Remove characters from color legend
      theme(
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(), # Hide x-axis ticks and text
        strip.text = element_text(size = 14, face = "bold"), # Facet label size and bold
        axis.title.y = element_text(size = 14, face = "bold") # Y-axis title size and bold
    )
    
    jpeg(file, width = 6, height = 6, units = "in", res = 300)
    print(p)
    dev.off()
  }
  
  scores <- c("CADD", "REVEL", "BayesDel", "AlphaMissense")
  
  d <- data %>% dplyr::filter(dis_ben %in% c(0,1)) %>%
    dplyr::mutate(exon = ifelse(is.na(exon), "?", exon)) %>%
    dplyr::mutate(exon = factor(exon, levels = c("1", "2", "3", "?")))
  d <- d %>% tidyr::pivot_longer(cols = all_of(scores), names_to = "annot", 
                                 values_to = "score") %>%
    dplyr::mutate(annot = factor(annot, levels = scores))
  
  file <- paste0(prefix, "_observed.jpg")
  plot_func(d, disease=TRUE, file = file)
  cat("Plotted to file", file, "\n")
  
  
  d <- data %>% dplyr::filter(pathannot == "Uncertain") %>%
    dplyr::mutate(exon = ifelse(is.na(exon), "?", exon)) %>%
    dplyr::mutate(exon = factor(exon, levels = c("1", "2", "3", "?")))
  d <- d %>% tidyr::pivot_longer(cols = all_of(scores), names_to = "annot", values_to = "score") %>%
    dplyr::mutate(annot = factor(annot, levels = scores))
  file <- paste0(prefix, "_uncert.jpg")
  plot_func(d, disease=FALSE, file = file, pt=1)
  cat("Plotted to file", file, "\n")
}

make_support_labels <- function(tbls){
  
    tbls <- tbls %>% dplyr::select(annot, matches("path|beni"), min_obs, max_obs) %>% 
    dplyr::rename(p_su = pathogenic_supporting, 
                  p_m = pathogenic_moderate, 
                  p_st = pathogenic_strong, 
                  p_vs = pathogenic_very_strong, 
                  b_su = benign_supporting, 
                  b_m = benign_moderate, 
                  b_st = benign_strong, 
                  b_vs = benign_very_strong) %>%
    pivot_longer(col = p_su:b_vs, names_to = "Label", values_to = "score") %>% 
    dplyr::mutate(Label = paste0(gsub("_", "[", Label), "]"),
                  score = as.numeric(score), 
                  min_obs = as.numeric(min_obs),
                  max_obs = as.numeric(max_obs)) %>%
    dplyr::filter(   (substring(Label, 1,1) == "p" & score < max_obs) | 
                     (substring(Label, 1,1) == "b" & score > min_obs) )
  
  tick_marks <- tbls %>% dplyr::select(annot, Label, score)
  tick_marks <- tick_marks %>% dplyr::mutate(pathannot = "Benign", exon = 1)
  return( tick_marks)
  
}

prob_imputation <- function(data, prop_path){
  
  ## data has the following columns::
  ## 1) dis_ben is a 0/1 or T/F indicator if the variants is classified as disease or benign
  ## 2) CADD, REVEL, BayesDel, AlphaMissense: annotation scores
  ## 3)

  ## prop_path = a priori probably of a variant being pathogenic. An 
  ## upper bound would be the frequency the disease in the popultion (for Dom Disease)
  data$dis_ben <- as.integer(data$dis_ben)
  total_obs <- sum(!is.na(data$dis_ben))
  n_case <- sum(data$dis_ben == 1, na.rm=T)
  n_cont <- sum(data$dis_ben == 0, na.rm=T)

  ## we expect the number of observed pathogenic variants to be 
  ## selection biased in clinvar because we don't actually observe all
  ## pathogentic variants. We need to adjust for that in the logistic regression.
  ## Interpret case = pathogentic / cntl = benign
  ## weight adjustment for a variant being benign/pathogentic
  ## given the observed data and expected data () 
  w_case <- total_obs * prop_path / n_case
  w_cntl <- total_obs * (1 - prop_path) / n_cont

  ## number of imputations used to come up with estimates of logistics regression 
  ## parameters.
  n_imputations = 25
    
  ## all data
  data_all <- data
  
  ## observed data (dis / benign) (remove VUS/uncertain which are listed as NAs)
  d_observed <- data %>% dplyr::filter(!is.na(dis_ben)) %>%
    dplyr::mutate(w = ifelse(dis_ben == 1, w_case, w_cntl))
  
  ## uncertain data ()
  d_unc <- data %>% dplyr::filter(is.na(dis_ben))
  
  ## calculate AUC values and get data to make ROC curves
  tmp <- AUCs <- auc_imputation(d_observed, d_unc)
  AUCs <- tmp$aucs
  roc_data <- tmp$roc_data
  
  model_estimates <- list()
  model_estimates_noweight <- list()
  model_estimates_weight <- list()
  plot_data <- c()
  
  for (annot in c("CADD","REVEL","BayesDel", "AlphaMissense")){
    
    ## ############################
    ## Perform weighted logistic regression to get parameter estimate
    ## and estimate probability that uncertain variants are pathogenic
    ## ############################
    
    d_obs <- d_observed %>% dplyr::select(all_of(c(annot, "dis_ben", "w")))
    formula <- paste("dis_ben ~", annot)
    # Converting the string to an actual formula
    formula <- as.formula(formula)
    model <- glm(formula, data = d_obs, family="binomial", weight=w)
    model_noweight <- glm(formula, data = d_obs, family="binomial")
    ## Estimate probably uncertain variants are pathogenic
    d_unc$prob_disease <- predict(model, newdata = d_unc, type = "response")
  
    ## reweight for the probabilistic imputation
    d_all <- data_all %>% dplyr::mutate(weight = 1)
    d_all <- left_join(d_all, d_unc)
    d_all <- d_all %>% dplyr::mutate(weight = ifelse(is.na(dis_ben), abs(prob_disease - .5)*2, weight))

    models <- vector("list", n_imputations) # Store models
    aucs <- numeric(n_imputations) # Store AUC imputation
    
    ## perform naive estimates
    model_estimates_weight[[annot]] <- list(params = summary(model)$coef[,1],
                                            cov = vcov(model), 
                                            t = summary(model)$coef[,3], 
                                            p = summary(model)$coef[,4])
                                            
    model_estimates_noweight[[annot]] <- list(params = summary(model_noweight)$coef[,1],
                                              cov = vcov(model_noweight), 
                                              t = summary(model_noweight)$coef[,3], 
                                              p = summary(model_noweight)$coef[,4])
  
    
    for (i in 1:n_imputations) {
      # Assuming 'df' has columns: 'predictor', 'classification' (NA for uncertain), 'prob_pathogenic'
      d_imp <- d_all
      
      # Assigning probabilistic classifications
      d_imp <- d_imp %>% rowwise() %>%
        dplyr::mutate(dis_ben = ifelse(is.na(dis_ben), rbinom(1, 1, prob_disease), dis_ben))
      
      # Fit the weighted logistic regression model
      models[[i]] <- glm(formula, family = binomial, data = d_imp, weights = weight)
      
      ## ###############
      ## Collect data for plotting 
      ## ###############
      if (i == 1){
        tmp <- d_imp %>% dplyr::select(all_of(c("dis_ben", annot))) %>%
          dplyr::rename(score = !!sym(annot)) %>% 
          dplyr::mutate(annot = annot)
        if (is.null(plot_data)){
          plot_data = tmp
        }else{
          plot_data = rbind(plot_data, tmp)
        }
      }
    }
    
    ## #####
    ## COMBINE INFORMATION ACCROSS IMPUTATIONS
    ## #####
    
    # Placeholder for aggregated variance-covariance matrices
    vcov_within <- matrix(0, nrow = length(models[[1]]$coefficients), ncol = length(models[[1]]$coefficients))
    vcov_between <- matrix(0, nrow = length(models[[1]]$coefficients), ncol = length(models[[1]]$coefficients))
    
    # Calculate W and B
    for(i in 1:n_imputations) {
      vcov_within <- vcov_within + vcov(models[[i]])
    }
    vcov_within <- vcov_within / n_imputations
    
    beta_bar <- sapply(models, function(x) coef(x)) %>% rowMeans
    for (i in 1:n_imputations) {
      diff <- matrix(coef(models[[i]]) - beta_bar, ncol = length(beta_bar))
      vcov_between <- vcov_between + (t(diff) %*% diff)
    }
    vcov_between <- vcov_between / (n_imputations - 1)
    
    # Calculate T
    vcov_total <- vcov_within + (1 + 1/length(models)) * vcov_between
    
    # T and P-value
    t_stat <- beta_bar / sqrt(diag(vcov_total))
    p_value <- 2 * pt(-abs(t_stat), df = length(models) - 1)
    
    model_estimates[[annot]] <- list(params = beta_bar, cov = vcov_total, t = t_stat, p = p_value)
    
    
  }
   
  return(list(model_estimates = model_estimates, 
              model_estimates_noweight = model_estimates_noweight, 
              model_estimates_weight = model_estimates_weight, 
              aucs =AUCs, 
              roc_data = roc_data, 
              plot_data = plot_data))
}


auc_imputation <- function(d_obs, d_unc, annot, p_train=.8, n_imp=200){
  
  AUCs <- c()
  roc_data = c()
  aucs <- numeric(n_imp) # Store AUC imputation
  n_path = sum(d_obs$dis_ben == 1, na.rm=T)
  n_ben = sum(d_obs$dis_ben == 0, na.rm=T)
  n_unc = nrow(d_unc)
  
  for (annot in c("CADD","REVEL","BayesDel", "AlphaMissense")){
  
    cat("Working on ", annot,"\n")
    
    for (i in 1:n_imp){
    
      i_path <- sample(1:n_path,  ceiling(n_path * p_train))
      i_ben  <- sample(1:n_ben,   ceiling(n_ben * p_train))
      i_unc  <- sample(1:n_unc,   ceiling(n_unc * p_train))
      
      ## ##########
      ## Create sampled tested and training datasets
      ## ##########
    
      ## training sets ##
      d_path_tr = d_obs %>% dplyr::filter(dis_ben == 1) %>% 
        dplyr::mutate(row = row_number()) %>%
        dplyr::filter(row %in% i_path)
      d_ben_tr = d_obs %>% dplyr::filter(dis_ben == 0) %>% 
        dplyr::mutate(row = row_number()) %>%
        dplyr::filter(row %in% i_ben)
      d_unc_tr <- d_unc %>% dplyr::mutate(row = row_number()) %>%
        dplyr::filter(row %in% i_unc)
      
    ## testing sets
      d_path_te = d_obs %>% dplyr::filter(dis_ben == 1) %>% 
        dplyr::mutate(row = row_number()) %>%
        dplyr::filter(row %in% i_path == F)
      d_ben_te = d_obs %>% dplyr::filter(dis_ben == 0) %>% 
        dplyr::mutate(row = row_number()) %>%
        dplyr::filter(row %in% i_ben == F)
      d_unc_te <- d_unc %>% dplyr::mutate(row = row_number()) %>%
        dplyr::filter(row %in% i_unc == F)
      
      ## ###############
      ## Add the prob of disease for uncertain variants
      ## ################
      
      ### Get weights for unc variants ##
      ## Training
     
      formula <- paste("dis_ben ~", annot)
      # Converting the string to an actual formula
      formula <- as.formula(formula)
      
      d_o_tr <- rbind(d_path_tr, d_ben_tr) %>%
        dplyr::select(all_of(c(annot, "dis_ben", "w")))
      model <- glm(formula, data = d_o_tr, family="binomial", weight=w)
      
      d_unc_tr$prob_disease <- predict(model, newdata = d_unc_tr, type = "response") 
      d_unc_tr <- d_unc_tr %>% 
        rowwise() %>%
        dplyr::mutate(dis_ben = rbinom(1, 1, prob_disease), 
                      w = abs(prob_disease - .5)*2 ) %>%
        dplyr::select(all_of(c(annot, "dis_ben", "w"))) %>% ungroup()               
      
      d_o_te <- rbind(d_path_te, d_ben_te) %>%
        dplyr::select(all_of(c(annot, "dis_ben", "w")))
      model <- glm(formula, data = d_o_te, family="binomial", weight=w)
      d_unc_te$prob_disease <- predict(model, newdata = d_unc_te, type = "response")
      d_unc_te <- d_unc_te %>% 
        rowwise() %>%
        dplyr::mutate(dis_ben = rbinom(1, 1, prob_disease), 
                      w = abs(prob_disease - .5)*2 ) %>%
        dplyr::select(all_of(c(annot, "dis_ben", "w")))       
      
      ## ###########
      ## Reweight observed data for the probabilistic imputation 
      ## Combine observed and uncertain
      ## #########
      d_o_tr$w = 1; d_o_te$w = 1
      d_all_tr <- rbind(d_o_tr, d_unc_tr)
      d_all_te <- rbind(d_o_te, d_unc_te)
      
      ## ########
      ## Fit the weighted logistic regression model to test sample
      ## ########
      model <- glm(formula, family = binomial, data = d_all_tr, weights = w)
      preds <- predict(model, d_all_te, type="response")
      
      ## ####
      ## Calculate AUC 
      ## ####
      pred_obj <- ROCR::prediction(preds, d_all_te$dis_ben)  # actual_responses needs to match the model's data
      
      auc <- performance(pred_obj, "auc")
      aucs[i] <- auc@y.values[[1]]  # Extract AUC value
      
      ## ###############
      ## Collect data for plotting 
      ## ###############
      if (i == 1){
        tmp <- d_all_te %>% dplyr::select(all_of(c("dis_ben", annot))) %>%
          dplyr::rename(score = !!sym(annot)) %>% 
          dplyr::mutate(annot = annot)
        if (is.null(roc_data)){
          roc_data = tmp
        }else{
          roc_data = rbind(roc_data, tmp)
        }
      }
    }
    AUCs <- rbind(AUCs, c(annot, mean(aucs), sqrt(var(aucs))))
  }
  
  return(list(aucs = AUCs, roc_data = roc_data))
}

acmg_table_impute <- function(est_perm, data, prop_path, file="acmg_table_impute.txt"){
  
  rslts <- c()
  cutoffs <- c(2.406, 5.790, 33.53, 1124)
  cutlabel <- c("supporting", "moderate", "strong", "very_strong")
  pD = prop_path
  odds_path = pD / (1 - pD)
  
  
  min_func <- function(x, beta, v, cutoff, odds_path, var_type ){
    
    if (! var_type %in% c("path","benign")){
      cat("var_type must be path or benign\n")
      stop()
    }
    
    ## calc se of log odds ##
    var_adj <- ifelse(var_type == "path", 1, -1)
    se_log_odds <- sqrt(t(c(1, x)) %*% v %*% c(1, x)) # Standard error
    z <- qnorm(.95)
    logodds_bound <- (beta[1] + beta[2]*x) - var_adj * z * se_log_odds
    odds_bound <- exp( var_adj * logodds_bound)
    odds_path = odds_path^var_adj
    
    return(cutoff - odds_bound / odds_path)
  }
  
  for (i in 1:4){
    
    betas = est_perm[[i]]$params
    v = est_perm[[i]]$cov
    cat = names(est_perm)[i]
    int <- betas[1]
    logor <- betas[2] 
    p <- est_perm[[i]]$p[2]
    se_int <- sqrt(v[1,1])
    se_logor <- sqrt(v[2,2])
    
    rng <- range(data[, cat])
    lower = -50; upper = 50
    #if (i == 4){ lower = 0; upper=1} ## alphamissense
    
    for (j in 1:length(cutoffs)){
      cutoff <- cutoffs[j]
      
      #minimize_func_path <- function(x){ return( cutoff - exp(int + x*lor)/odds_path) }
      #minimize_func_ben <- function(x){ return( cutoff -  exp(-(int + x*lor))*odds_path) }
      
      result <- tryCatch({
        uniroot(min_func, lower = lower, upper = upper, beta=betas, v = v,  
                cutoff = cutoff, odds_path = odds_path, var_type = "path")
      }, error = function(e){
        cat("No satisfying threshold was able to be found for:\n")
        cat(cat, " using cutoff: ", cutoff, " for pathogenic variant\n")
        return(list(root = NA))
      })
      
      rslts <- rbind(rslts, c(cat, int, logor, se_int, se_logor, p, 
                              "pathogenic", cutlabel[j], result$root, rng))
      
      result <- tryCatch({
        uniroot(min_func, lower = lower, upper = upper, beta=betas, v = v,
                cutoff = cutoff, odds_path = odds_path, var_type = "benign")
      }, error = function(e){
        cat("No satisfying threshold was able to be found for:\n")
        cat(cat, " using cutoff: ", cutoff, " for benign variant\n")
        return(list(root = NA))
      })
      rslts <- rbind(rslts, c(cat, int, logor, se_int, se_logor, p, 
                              "benign", cutlabel[j], result$root, rng))
    }       
  }
  
  rslts <- as.data.frame(rslts, stringsAsFactors=F)
  names(rslts) <- c("annot", "int", "logor", "se_int", "se_logor", "pval", "var_type", 
                    "acmg_cat", "annot_value", "min_obs", "max_obs")
  
  tbl <- rslts %>% tidyr::pivot_wider(names_from=c("var_type","acmg_cat"),
                                      values_from="annot_value",
                                      names_glue = "{var_type}_{acmg_cat}") %>%
    dplyr::select(annot, int, logor, se_int, se_logor, pval, 
                  starts_with("pathogenic"), starts_with("benign"),
                  min_obs, max_obs)
  
  
  write.table(file = file, tbl, quote=F, row.names=F, col.names=T, sep="\t")
  cat("Wrote threshold for different ACMG support categories to", file, "\n")
  
  return(tbl)
}

expected_path_counts <- function(data, tbl, file = "expected_path_counts.txt"){
  
  
  t <- tbl %>% dplyr::select(annot, matches("pathog|benign"))
  t <- t %>% tidyr::pivot_longer(cols = matches("pathog|benign"), 
                                 names_to = "category", 
                                 values_to = "cutoff") %>%
    dplyr::mutate(cutoff = as.numeric(cutoff))
  out <- data.frame()
  
  annots <- c("CADD","REVEL","BayesDel", "AlphaMissense")
  
  for (ann in annots){
  
      for (pb in c("path", "benign")){
        categories <- t %>% dplyr::filter(annot == ann & grepl(pb, category)) 
        if (pb == "path"){
          categories <- categories %>% dplyr::arrange(desc(cutoff))
        }else{
          categories <- categories %>% dplyr::arrange(cutoff)
        }
      
        for (cat in categories$category){
          cutoff <- t %>% filter(category == cat & annot == ann) %>% 
          dplyr::select(cutoff) %>% c() %>% as.numeric()
        
          if (grepl("pathogenic", cat)){
            ind <- which(vals >= cutoff)
          }else{
            ind <- which(vals <= cutoff)
          }
          count = length(ind)
          if (count > 0){
            vals <- vals[-ind]
          }
          out <- rbind(out, c(ann, cat, count))
        }
      }
    out <- rbind(out, c(ann, "unclassified", length(vals)))
  }
  
  names(out) <- c("annot", "support", "n_variants")
  
  out <- out %>% tidyr::pivot_wider(names_from = support, values_from = n_variants)
  out <- out %>% dplyr::select(annot, pathogenic_supporting, pathogenic_moderate, pathogenic_strong,
                               pathogenic_very_strong, benign_supporting, benign_moderate, 
                               benign_strong, benign_very_strong, unclassified)
  
  write.table(file = file, out, quote=F, row.names=F, col.names=T, sep = "\t")
  cat("Wrote table to", paste0(getwd(),"/", file), "\n")
  return(out)
}
  
expected_counts_summary <- function(data, tbl, file = "summary_counts.txt"){
  
  # Function to compute the desired columns
  calculate_differences <- function(values) {
    numb_uniq <- length(table(values)) 
    
    diff_matrix <- abs(outer(values, values, "-"))
    diff_vector <- c(diff_matrix)
    max_diff = max(diff_vector)
    mean_diff <- diff_vector[diff_vector!=0]
    mean_diff <- ifelse(length(mean_diff) != 0, mean(mean_diff), 0)
    return(c(numb_uniq, max_diff, mean_diff))
  }
  
  t <- tbl %>% dplyr::select(annot, matches("pathog|benign"))
  t <- t %>% tidyr::pivot_longer(cols = matches("pathog|benign"), names_to = "category", 
                                 values_to = "cutoff") %>%
              dplyr::mutate(cutoff = as.numeric(cutoff))
  out <- data.frame()
  
  cat_table <- c()
  annots <- c("CADD","REVEL","BayesDel", "AlphaMissense")
  
  for (ann in annots){
    vals <- data %>% dplyr::filter(is.na(dis_ben)) %>%
      dplyr::select(all_of(ann)) %>% unlist() %>% as.numeric()
    cutoffs <- t %>% dplyr::filter(annot == ann) %>% 
      dplyr::arrange(cutoff) 
    co <- c(-1000, cutoffs$cutoff, 1000)
    labels <- c(cutoffs$category[1:4], "unclass", cutoffs$category[5:8])
    categories <- cut(vals, breaks=co, labels = labels)  
    
    cat_table <- cbind(cat_table, as.numeric(categories)) 
  }
  
  ## Uniq_combos contains the number of observed annotation combinations, 
  ## the number of unique classifications observed (n_cats), 
  ## the maximum distance between classifications (max_diff) and 
  ## the mean of the mean differences between classifications (mean_diff)
  
  ## 
  
  
  cat_table <- as.data.frame(cat_table, stringsAsFactors = FALSE)
  names(cat_table) <- annots
  uniq_combos <- cat_table %>% group_by(CADD,REVEL,BayesDel, AlphaMissense) %>% 
    dplyr::mutate(n = n()) %>% ungroup() %>% unique() %>% arrange(desc(n)) 
  result <- t(apply(uniq_combos[,-5], 1, calculate_differences))
  result <- as.data.frame(result)
  names(result) <- c("n_cats", "max_diff", "mean_diff")
  uniq_combos <- cbind(uniq_combos, as.data.frame(result))
  
  ## Figure out the total number of variants where all classifications agree
  match_counts <- uniq_combos %>% group_by(n_cats) %>% 
    dplyr::summarize(variant_count = sum(n), 
                     mean_diff = mean(mean_diff))
  max_diffs <- uniq_combos %>% group_by(n_cats, max_diff) %>% 
    dplyr::summarize(variant_count = sum(n), 
                     mean_diff = mean(mean_diff))
  
  write.table(file = file, uniq_combos, quote=F, row.names=F, col.names=T, sep="\t")
  cat("Wrote ", file, "\n")
  
  new_file = gsub(".txt", "_match_counts.txt", file)
  write.table(file = new_file, match_counts, quote=F, row.names=F, col.names=T, sep="\t")
  cat("Wrote", new_file, "\n")
  
  new_file = gsub(".txt", "_max_diffs.txt", file)
  write.table(file = new_file,max_diffs, quote=F, row.names=F, col.names=T, sep="\t")
  cat("Wrote", new_file, "\n")
  
  return(list(uniq_combos = uniq_combos, match_counts = match_counts, max_diffs = max_diffs))
  
}

homeo_plot <- function(data, prefix = "homeo"){
  
  ## plot annotation scores as a function of AA position
 
  plot_func <- function(d){
    
    ## values from Andy via Alamute Visual.
    hd <- c(292, 471)
    exons <- c(241, 429, 945)   
    p <- ggplot(d, aes(x = pos_c, y = score, color = pathannot))
    
    p <- p + geom_point(size=2, alpha=.4) +
      facet_wrap(~factor(annot), scales = "free_y") 
    p <- p + scale_color_manual(values = c("Benign" = "darkgreen", "Pathogenic"="darkred", 
                                           "Uncertain" = "blue"))
    p <- p + geom_vline(xintercept = hd, linetype = "dashed", color = "red")
    p <- p + geom_vline(xintercept = exons)
    p <- p + theme_bw() + xlab("Position") 
    p <- p + labs(color = "Variant Class") # Change title for color legend)
    p <- p + guides(color = guide_legend(override.aes = list(shape = 16 ))) + # Remove characters from color legend
      theme(
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(), # Hide x-axis ticks and text
        strip.text = element_text(size = 14, face = "bold"), # Facet label size and bold
        axis.title.y = element_text(size = 14, face = "bold", ) # Y-axis title size and bold
      )
    
    file <- paste0(prefix, "_score_x_pos.pdf")
    
    pdf(file, width=11, height = 8.5, )
    print(p)
    dev.off()
    cat("Created plot", file, "\n")
  }
  
  scores <- c("CADD", "REVEL", "AlphaMissense", "BayesDel")
  
  d <- data %>% dplyr::mutate(exon = factor(exon, levels = c("1", "2", "3")))
  d <- d %>% tidyr::pivot_longer(cols = all_of(scores), names_to = "annot", values_to = "score") %>%
    dplyr::mutate(annot = factor(annot, levels = c("CADD","REVEL","BayesDel", "AlphaMissense")))
  
  plot_func(d)
}



## ############################## ##
## END MAIN FUNCTIONS
## KEPT SOME OLDER FUNCTIONS FSAG
## ######################



plot_distributions <- function(data, mns, plotfile = "Distribution_plots.pdf"){
  
  pdf(plotfile)
  
  ## violin plots to show distribution with variant type
  annots <- c("CADD","REVEL","BayesDel", "AlphaMissense")
  for (annot in annots){
    
    d <- data %>% dplyr::select(c("pathannot", annot))
    p <- ggplot(data=d, aes_string(x="pathannot", y=annot, fill="pathannot")) + 
      geom_violin(trim=FALSE, alpha=.3) +
      geom_dotplot(binaxis='y', stackdir='center', dotsize=.35) +
      theme_minimal() + xlab("Variant type") + ylab(annot) + 
      scale_fill_manual(breaks = c("Benign", "Pathogenic", "Uncertain"),
                        values = c("darkgreen","darkred", "darkorange"))
    print(p)
  }
  
  
  for (annot in annots){
    
    txt <- paste0(annot,"_mean")
    annot_means <- mns %>% dplyr::select(!!txt, pathannot)
    ## Create density plots grouped by Path annotation
    p <- ggplot(data, aes_string(x = annot, fill = "pathannot")) + 
      geom_vline(data=annot_means, aes(xintercept=!!txt, color=pathannot), linetype="dashed",
                 size=1) +
      geom_histogram(aes(y=..density..), alpha=.4, bins=16) +
      geom_density(alpha=.2, aes(fill=pathannot)) +
      theme_minimal() + xlab(annot) + ylab("Density") 
    print(p)
  }
  
  dev.off()
  
  cat("Wrote distribution plots to", plotfile, "\n")
}



plot_pairwise <- function(data, plotfile = "Pairwise.pdf"){
  
  annots <- c("CADD","REVEL","BayesDel", "AlphaMissense")
  
  pdf(file = plotfile)
  
  for (i in 1:(length(annots)-1)){
    for (j in (i+1):length(annots)){
      annot1 = annots[i]
      annot2 = annots[j]
      
      
      p <- ggplot(data, aes_string(x=annot1, y=annot2,
                                   color = "pathannot", fill = "pathannot")) +
        geom_point(aes(size=pathannot)) +
        theme_minimal() + xlab(annot1) + ylab(annot2) +
        scale_color_manual(breaks = c("Benign", "Pathogenic", "Uncertain"),
                           values = c("darkgreen","darkred", "darkorange")) +
        scale_size_manual(breaks = c("Benign", "Pathogenic", "Uncertain"),
                          values = c(1.75,1.75,.75)) +
        guides(fill="none")
      print(p)
    }
  }
  dev.off()
  cat("Wrote pairwise plots to", plotfile, "\n")
}

