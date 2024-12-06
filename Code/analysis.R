##
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotROC)
library(ROCR)
library(data.table)

## set working directory ##
home <- "/home/win/askol/Research/phox2b_annotation/"
home_win <- "L:/askol/Research/phox2b_annotation/"
out_dir <- paste0(home_win,"results_gnomad41/")

if (dir.exists(home)){
  setwd(home)
}else{
  setwd(home_win)
}
results_dir <- paste0(getwd(),"/Paper_Results/")

source("Code/analysis_funcs.R")

## important constants
## expected proportion of pathogenic variants among all variants in gnomad
prop_path = 0.044

## ######## ##
## Necessary Data
## ######## ##
excel_file <- "C:/Users/askol/OneDrive - Lurie Childrens/PHOX2B/PHOX2B missense data collection.xlsx"
gnomad_faf_file <- "Resources/phox2b_hg38_faf_selFields.txt"
gnomad_ex_file <- "Resources/phox2b_hg38_gnomad_select_exome.txt"
gnomad_ge_file <- "Resources/phox2b_hg38_gnomad_select_genomes.txt"


## #########
## Load data from Excel sheet from Andy 
## #########
data <- load_data_excel(file = excel_file)
bayesdel <- data$bayesdel
alphamissense <- data$alphamissense
snv_preds <- data$all_snv_pred
data <- data$data

## ###########
## Load and merge gnomad data
## ###########
gnomad <- get_gnomad_data(gnomad_ex_file, gnomad_ge_file, gnomad_faf_file)

gnomad_coding_snv_pass = gnomad$gnomad_coding_snv_pass
gnomad_coding_snv_missense = gnomad$gnomad_coding_snv_missense
prop_path_update <- gnomad$prob_path_update

## ############
## merge gnomad data with xlsx sheet
## ############
data_comb <- full_join(data, gnomad_coding_snv_missense, by=c("p." = "pdot", "ref", "alt"))
data_comb <- update_pathannot(data_comb)
data_comb <- add_annot(data_comb, bayesdel, alphamissense)

## ######
## FINAL DATA FOR ANALYSIS 
## Keep gnomad pass data or variants labelled as "disease-causing"
data_use <- data_comb %>% dplyr::filter((my_filter=="PASS" | is.na(my_filter) | pathannot == "Pathogenic"))
data_use <- remove_dupe_pdots(data_use)

## ##########
## keep reduced set of variables
## Rename annotation to easy plotting
## ###########
data_use <- data_use %>% dplyr::select(p., g., pos, pos_c, c., faf95max, ref, alt, cadd, 
                                       revel, alphamissense, bayesdel, 
                                       pathannot_update, dis_ben, dis_unc, exon) %>%
  dplyr::rename(pathannot=pathannot_update, 
                REVEL = revel,
                CADD = cadd,
                BayesDel = bayesdel,
                AlphaMissense = alphamissense)

## Write supplementary tables out to results directory
file = paste0(results_dir,"Supplementary_Table1.txt")
write.table(file = file, data_use, quote=F, row.names=FALSE, col.names=TRUE, sep="\t")

## ########################### ##
## START ANALYSES FOR PAPER ##
## ########################### ##

data <- data_orig <- data_use ## just in case

## ############
## Perform naive regression (no weights or imputation)
## AUCs will be calculated later and need to be added to table in paper
## ############
tbl_file <- paste0(results_dir, "Table1.txt")
make_table1(data, tbl_file = tbl_file)

## #################
## INCORPORATE UNCERTAIN VARIANTS USING IMPUTATION
## ################
tmp <- prob_imputation(data, prop_path = prop_path_update)

est_perm <- tmp$model_estimates
est_noweight <- tmp$model_estimates_noweight
est_weight <- tmp$model_estimates_weight
aucs <- tmp$aucs
roc_data <- tmp$roc_data
plot_data <- tmp$plot_data

## #########
## Write auc table ##
## #########
write.table(file=paste0(results_dir, "AUC_imputed_80_20.txt"), 
           roc_data, quote=F, row.names=F, col.names=T, sep="\t")

## ##############
## plot ROC curve using an imputed dataset ##
## ##############
plot_file <- paste0(results_dir,"FullSample_ROC_plot.pdf")
plot_roc_imputed(plot_data, plotfile = plot_file)

tbl <- acmg_table_impute(est_perm, data=data, prop_path = prop_path_update,
                         file=paste0(results_dir,"acmg_table_impute.txt"))

## ########
## Expect benign / pathogenic counts
## #######
e_tbl <- expected_path_counts(data, tbl, 
                              file = paste0(results_dir,"expected_path_counts.txt"))

tmp <- expected_counts_summary(snv_preds, tbl)

## ####
## Plot distribution of annotation scores with 
## classification thresholds
## ####
tick_marks <- make_support_labels(tbl)
prefix = paste0(results_dir, "scores_and_cutoffs")
plot_scores(data, tbl, tick_marks, prefix = prefix)

## #######
## Plot scores as a function of gene position with 
## homeodomain noted 
## #######
homeo_plot(data, prefix=paste0(results_dir, "homeo_domain"))

## #####
## Calculate mean of annotations scores in and out of the homeodomain
## #####

file = paste0(results_dir, "Supplementary_Table2.txt")
tbl_file <- paste0(results_dir, "homeo_table.txt")
homeo_summary(snv_preds, file = tbl_file, file_supp = file)













