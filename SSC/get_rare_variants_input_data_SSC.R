rm(list = ls())
# workingDir <- "/Users/shining17wq/Google Drive/MorrowLab/SFARI/SSC/jobs"
workingDir <- "/gpfs/data/emorrow/silvio/Qing_Wu/SFARI/batch_jobs"
setwd(workingDir)
getwd()

library(data.table)
library(stringr)
library(dplyr)
library(plyr)
options(stringsAsFactors = F, scipen = 999)


# generate annovar input --------------------------------------------------

var_df <- read.table(file = "../hail_out/SSC_WES_3/SSC_WES_3.vqsr.vcf_families_merged.pass_QC.rare1pct_variants_het_info.txt", sep = '\t', header = F)
names(var_df) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
head(var_df)

# edit locus
locus_df <- var_df[c(1,2,4,5)]
locus_df$CHROM <- str_replace(locus_df$CHROM, "chr", "")
locus_df$POS <-  as.numeric(locus_df$POS)
head(locus_df)
# locus end
locus_df$END <- locus_df$POS + nchar(locus_df$REF) - 1

avinput_df <- locus_df[c(1,2,5,3,4)]
head(avinput_df)

# write.table(avinput_df, file = "../hail_out/SSC_WES_3/SSC_WES_3.vqsr.vcf_families_merged.pass_QC.rare1pct_variants_het_info.avinput", sep = "\t", quote = F, col.names = F, row.names = F)

# convert to HUGO official gene symbol ------------------------------------

var_gene_anno_df <- fread(file = "../hail_out/SSC_WES_3/SSC_WES_3.vqsr.vcf_families_merged.pass_QC.rare1pct_variants_het_info.hg19_multianno.csv")
head(var_gene_anno_df)

length(unique(var_gene_anno_df$Gene.refGene)) # 22718
length(unique(var_gene_anno_df$Gene.knownGene)) # 23364
length(unique(var_gene_anno_df$Gene.ensGene)) # 26435

grep("^\\.$", var_gene_anno_df$Gene.ensGene, value = T)
grep("^$", var_gene_anno_df$Gene.ensGene, value = T)
# no empty gene_name

gene_name_df <- data.frame(refGene = var_gene_anno_df$Gene.refGene, 
                           knownGene = var_gene_anno_df$Gene.knownGene,
                           ensGene = var_gene_anno_df$Gene.ensGene)

gene_name_df$match_refGene2knownGene <- gene_name_df$refGene == gene_name_df$knownGene
gene_name_df$match_refGene2ensGene <- gene_name_df$refGene == gene_name_df$ensGene
gene_name_df$match_knownGene2ensGene <- gene_name_df$knownGene == gene_name_df$ensGene

gene_name_df$gene_name <- gene_name_df$ensGene
# ref=known=ens, do not need to replace gene name
# ref=known=!ens, replace ensembl gene name
idx_identical_ref_known <- gene_name_df$match_refGene2knownGene & !gene_name_df$match_knownGene2ensGene
gene_name_df$gene_name[idx_identical_ref_known] <- gene_name_df$refGene[idx_identical_ref_known]
# ref=ens!=known, do not need to replace gene name
# ref!=known=ens, do not need to replace gene name
# ref!=known!=ens, keep ensembl gene name

# edit multiple gene name issue
idx_multi_gene_in_refGene <- grepl(";", gene_name_df$refGene)
idx_multi_gene_in_knownGene <- grepl(";", gene_name_df$knownGene)
table(grepl(";", gene_name_df$gene_name))
# multi-gene found in gene_name, but single name in refGene, replace gene_name by refGene
idx_single_refGene <- grepl(";", gene_name_df$gene_name) & !idx_multi_gene_in_refGene
gene_name_df$gene_name[idx_single_refGene] <- gene_name_df$refGene[idx_single_refGene]
# multi-gene found in gene_name, but single name in knownGene, replace gene_name by knownGene
idx_single_knownGene <- grepl(";", gene_name_df$gene_name) & !idx_multi_gene_in_knownGene
gene_name_df$gene_name[idx_single_knownGene] <- gene_name_df$refGene[idx_single_knownGene]

# check how many gene has multiple names
idx_genes_name <- grep(";", gene_name_df$gene_name) # 19567
intra_genes_name <- grep(";", gene_name_df$gene_name, value = T)
# split multiple names by ;
intra_genes_lst <- str_split(intra_genes_name, ";")
intra_genes_df <- plyr::ldply(intra_genes_lst, rbind)
num_unique_gene_name <- apply(intra_genes_df, 1, function(x) length(unique(x)))
idx_identical_genes <- (num_unique_gene_name - 1 == 1)
table(idx_identical_genes) # 297

replaced_gene_name <- intra_genes_df$`1`[idx_identical_genes]
gene_name_df$gene_name[idx_genes_name[idx_identical_genes]] <- replaced_gene_name

length(unique(gene_name_df$gene_name)) # 22746

unique_gene_name <- unique(gene_name_df$gene_name)
# unique_gene_name

variants2gene_df <- data.frame(Chr = var_gene_anno_df$Chr,
                               Start = var_gene_anno_df$Start,
                               End = var_gene_anno_df$End,
                               Ref = var_gene_anno_df$Ref,
                               Alt = var_gene_anno_df$Alt,
                               Gene = gene_name_df$gene_name)
head(variants2gene_df)
write.csv(variants2gene_df, file = "../hail_out/SSC_WES_3/SSC_WES_3.vqsr.vcf_families_merged.pass_QC.rare1pct_variants_het_info_gene_converted.csv")

# annovar output ----------------------------------------------------------

var_anno_df <- fread(file = "../hail_out/SSC_WES_3/SSC_WES_3.vqsr.vcf_families_merged.pass_QC.rare1pct_variants_het_info.hg19_multianno.csv")
head(var_anno_df)

gene_name_df <- read.csv(file = "../hail_out/SSC_WES_3/SSC_WES_3.vqsr.vcf_families_merged.pass_QC.rare1pct_variants_het_info_gene_converted.csv", header = T)
head(gene_name_df)
length(unique(gene_name_df$Gene))
var_anno_df$gene_name<- gene_name_df$Gene

dnv_conf_level <- read.table(file = "../hail_out/SSC_WES_3/SSC_WES_3.vqsr.vcf_families_merged.pass_QC.rare1pct_variants_het.w_dnv.tsv", sep = '\t', header = T)
head(dnv_conf_level)
table(dnv_conf_level$dnv_confidence) # high: 6008; medium: 471; low: 1966

var_anno_df$dnv_conf_level<- dnv_conf_level$dnv_confidence
table(var_anno_df$dnv_conf_level)


# variants to include in following analysis -------------------------------

as.data.frame(table(var_anno_df$Func.refGene))
as.data.frame(table(var_anno_df$ExonicFunc.refGene))

# replace "." in ExonicFunc with "splicing"
idx_splicing <- var_anno_df$Func.refGene == "splicing"
var_anno_df$ExonicFunc.refGene[idx_splicing] <- "splicing"

idx_func_1 <- var_anno_df$Func.refGene == "exonic"
idx_func_2 <- var_anno_df$Func.refGene == "exonic;splicing"
idx_func_3 <- var_anno_df$Func.refGene == "splicing"
idx_exonic_func_1 <- var_anno_df$ExonicFunc.refGene == "synonymous SNV"

# function to write index of which indicate variant row to exclude from original variant by sample matrix
get_idx_var2exclude <- function(idx_var2incl, variant_category_name){
  idx_var2excl <- (!idx_var2incl)*c(1:nrow(var_anno_df))
  idx_var2excl <- idx_var2excl[idx_var2excl != 0]
  print(length(idx_var2excl))
  
  if (variant_category_name == "") {
    out_file_name <- "../hail_out/SSC_WES_3/SSC_WES_3.vqsr.vcf_families_merged.pass_QC.rare1pct_variants_het.idx_var2exclude.txt"
  } else {
    out_file_name <- paste("../hail_out/SSC_WES_3/SSC_WES_3.vqsr.vcf_families_merged.pass_QC.rare1pct_variants_het.idx_var2exclude", variant_category_name, "txt", sep = ".")
  }
  print(out_file_name)
  
  # write.table(idx_var2excl, file = out_file_name, sep = "\t", quote = F, col.names = F, row.names = F)
}

# include exonic and splicing variant, but not synonymous variants
idx_var2incl <- (idx_func_1 | idx_func_2 | idx_func_3) & (!idx_exonic_func_1)
table(idx_var2incl) 
# 683614 567910 (F/T)
get_idx_var2exclude(idx_var2incl, "")

# include PTVs only
idx_var_PTVs <- (var_anno_df$ExonicFunc.refGene == "stopgain" | var_anno_df$ExonicFunc.refGene == "frameshift substitution" | var_anno_df$ExonicFunc.refGene == "splicing")
get_idx_var2exclude(idx_var_PTVs, "PTVs")

# include PTVs and missense MPC > 2
var_anno_df$MPC_score <- as.numeric(var_anno_df$MPC_score)
idx_MPC_above2 <- var_anno_df$MPC_score >= 2 & !is.na(var_anno_df$MPC_score)
idx_var_PTVs_MisA <- idx_var_PTVs | idx_MPC_above2
get_idx_var2exclude(idx_var_PTVs_MisA, "PTVs_MisA")

# include PTVs and missense MPC > 1
idx_MPC_above1 <- var_anno_df$MPC_score >= 1 & !is.na(var_anno_df$MPC_score)
idx_var_PTVs_MisAB <- idx_var_PTVs | idx_MPC_above1
get_idx_var2exclude(idx_var_PTVs_MisAB, "PTVs_MisAB")

# generate VxA matrix -----------------------------------------------------

get_var_anno_by_variant_category <- function(var_anno_df, idx_var2incl, variant_category_name){
  
  var_anno_var2incl_df <- var_anno_df[idx_var2incl,]
  print(dim(var_anno_var2incl_df))
  # 567910     134
  
  # get variants to gene indices list ---------------------------------------
  
  gene_df <- as.data.frame(table(var_anno_var2incl_df$gene_name))
  names(gene_df)[1] <- "Gene"
  
  unique_gene <- unique(var_anno_var2incl_df$gene_name)
  
  tmp_idx <- which(var_anno_var2incl_df$gene_name %in% unique_gene)
  print(is.unsorted(tmp_idx)) # FALSE, indicating all gene name in unique_gene are sorted based on their location in var_anno_var2incl_df
  
  var2gene_lst <- lapply(unique_gene, function(x) which(var_anno_var2incl_df$gene_name %in% x))
  
  # check if all indices of each gene are sorted, by row
  # check if all genes are sorted, by the largest indices of each row should be smaller than the smallest indices of next row
  
  bool_if_sorted_wi_gene = c()
  for (i in 1:length(var2gene_lst)){
    bool_if_sorted_wi_gene = c(bool_if_sorted_wi_gene, is.unsorted(var2gene_lst[[i]]))
  }
  table(bool_if_sorted_wi_gene) # 18565 all false for gene_name
  
  # split indices that is not a consecutive number (multiple range)
  new_var2gene_lst = list()
  new_gene_names = c()
  for (i in 1:(length(var2gene_lst)-1)){
    bool_if_sorted_x_consecutive_genes_i = tail(var2gene_lst[[i]], n=1) > head(var2gene_lst[[i+1]], n=1)
    if (bool_if_sorted_x_consecutive_genes_i){
      # generate ranges of consecutive values of gene i
      ranges_of_unconsecutive_gene_i = unname(tapply(var2gene_lst[[i]], cumsum(c(1, diff(var2gene_lst[[i]])) != 1), range))
      # split the indices list of gene i into multiple indices list based on consecutive range
      for (j in 1:length(ranges_of_unconsecutive_gene_i)){
        new_list = list(seq(ranges_of_unconsecutive_gene_i[[j]][1], ranges_of_unconsecutive_gene_i[[j]][2]))
        new_var2gene_lst = append(new_var2gene_lst, new_list)
        new_gene_names = c(new_gene_names, unique_gene[i])
      }
    } else {
      new_var2gene_lst = append(new_var2gene_lst, list(var2gene_lst[[i]]))
      new_gene_names = c(new_gene_names, unique_gene[i])
    }
  }
  # add the last one to list
  new_var2gene_lst = append(new_var2gene_lst, list(var2gene_lst[[length(var2gene_lst)]]))
  new_gene_names = c(new_gene_names, unique_gene[length(unique_gene)])
  print(paste("Length of new var2gene list", length(new_var2gene_lst), sep = ":")) # 18822
  print(paste("Length of new gene name list", length(new_gene_names), sep = ":")) # 18822
  print(paste("Length of original variant 2 gene list", length(var2gene_lst), sep = ":")) # 18556
  
  # sort based on smallest value of each list range
  smallest_value_lst = c()
  for (i in 1:length(new_var2gene_lst)){
    smallest_value_lst = c(smallest_value_lst, range(new_var2gene_lst[[i]])[1])
  }
  
  new_var2gene_lst_sorted = new_var2gene_lst[order(smallest_value_lst)]
  new_gene_names_sorted = new_gene_names[order(smallest_value_lst)]
  
  # check if all elements are sorted
  smallest_value_lst = c()
  for (i in 1:length(new_var2gene_lst_sorted)){
    smallest_value_lst = c(smallest_value_lst, range(new_var2gene_lst_sorted[[i]])[1])
  }
  print(is.unsorted(order(smallest_value_lst))) # FALSE
  
  var2gene_df <- plyr::ldply(new_var2gene_lst_sorted, rbind)
  var2gene_df[is.na(var2gene_df)] <- ""
  
  var2gene_df <- data.frame(cbind(new_gene_names_sorted, var2gene_df))
  print(dim(var2gene_df)) # 18822  2023
  
  if (variant_category_name == "") {
    out_file_name <- "../hail_out/SSC_WES_3/SSC_WES_3.vqsr.vcf_families_merged.pass_QC.rare1pct_variants_het.idx_var2gene_lst_df.unconsecutive_indices_splited.txt"
  } else {
    out_file_name <- paste("../hail_out/SSC_WES_3/SSC_WES_3.vqsr.vcf_families_merged.pass_QC.rare1pct_variants_het.idx_var2gene_lst_df.unconsecutive_indices_splited", variant_category_name, "txt", sep = ".")
  }
  
  # write.table(var2gene_df, file = out_file_name, sep = "\t", quote = F, col.names = F, row.names = F)
  
  
  # get the full list of variants to gene indices list ----------------------
  
  var2gene_df <- plyr::ldply(var2gene_lst, rbind)
  var2gene_df[is.na(var2gene_df)] <- ""
  print(dim(var2gene_df)) # 18556  2022
  
  var2gene_df <- data.frame(cbind(unique_gene, var2gene_df))
  
  if (variant_category_name == "") {
    out_file_name <- "../hail_out/SSC_WES_3/SSC_WES_3.vqsr.vcf_families_merged.pass_QC.rare1pct_variants_het.idx_var2gene_lst_df.txt"
  } else {
    out_file_name <- paste("../hail_out/SSC_WES_3/SSC_WES_3.vqsr.vcf_families_merged.pass_QC.rare1pct_variants_het.idx_var2gene_lst_df", variant_category_name, "txt", sep = ".")
  }
  
  # write.table(var2gene_df, file = out_file_name, sep = "\t", quote = F, col.names = F, row.names = F)
  
  return(var_anno_var2incl_df)
}

var_anno_var2incl_df <- get_var_anno_by_variant_category(var_anno_df, idx_var2incl, "")
var_anno_var2incl_df <- get_var_anno_by_variant_category(var_anno_df, idx_var_PTVs, "PTVs")
var_anno_var2incl_df <- get_var_anno_by_variant_category(var_anno_df, idx_var_PTVs_MisA, "PTVs_MisA")
var_anno_var2incl_df <- get_var_anno_by_variant_category(var_anno_df, idx_var_PTVs_MisAB, "PTVs_MisAB")


# variants annotations ----------------------------------------------------

var_anno_var2incl_df <- var_anno_df[idx_var2incl,]

names(var_anno_var2incl_df)
dim(var_anno_var2incl_df) # 547949    134

var_anno_2train_df <- var_anno_var2incl_df[, c(1:7, 9, 11, 26, 29, 32, 35, 38, 41, 44, 47, 50, 55, 58, 61, 64, 69, 71, 95, 96, 100, 103, 131:134)]
head(var_anno_2train_df)

# replace "." in ExonicFunc with "splicing"
idx_splicing <- var_anno_2train_df$Func.refGene == "splicing"
var_anno_2train_df$ExonicFunc.refGene[idx_splicing] <- "splicing"

# replace NA in dnv_confidence with "."
idx_dnv_nv <- is.na(var_anno_2train_df$dnv_conf_level)
var_anno_2train_df$dnv_conf_level[idx_dnv_nv] <- "."


# exclude SP id -----------------------------------------------------------

# SFARI gene 
SFARI_gene_lst_df <- read.csv(file = "../reference_data/SFARI_gene/SFARI-Gene_genes_01-11-2022release_03-22-2022export.csv")
head(SFARI_gene_lst_df)

unique_gene <- unique(var_anno_df$gene_name)
idx_SFARI_gene <- unique_gene %in% SFARI_gene_lst_df$gene.symbol
table(var_anno_2train_df$gene_name[idx_SFARI_gene]) # 1024
non_SFARI_gene <- ifelse(idx_SFARI_gene, FALSE, TRUE)
table(non_SFARI_gene) # 1024/19059, T/F


idx <- c(1:length(unique_gene))[non_SFARI_gene]
# write.table(idx, file = "../SSC_WES_3/SPARK/SPARK_vcf_families_merged_rare_variant_idx_nonSFARIgene2exclude_intraV_incl.txt", sep = "\t", quote = F, col.names = F, row.names = F)

# continue with edit_rare_variant_info_table.R

getwd()
head(var_anno_2train_df)

table(var_anno_2train_df$ExonicFunc.refGene)

# generate VxA input ------------------------------------------------------

summary(var_anno_2train_df)

# stopgain, stoploss, splicing, frameshift substitution, nonframeshift substitution, nonsynonymous SNV
exonic_func_PTVs <- ifelse(var_anno_2train_df$ExonicFunc.refGene == "stopgain" | var_anno_2train_df$ExonicFunc.refGene == "frameshift substitution" | var_anno_2train_df$ExonicFunc.refGene == "splicing", 1, 0)
exonic_func_missense <- ifelse(var_anno_2train_df$ExonicFunc.refGene == "nonsynonymous SNV", 1, 0)

table(var_anno_2train_df$dnv_conf_level)
denovo <- ifelse((var_anno_2train_df$dnv_conf_level == "HIGH" | var_anno_2train_df$dnv_conf_level == "MEDIUM") & !is.na(var_anno_2train_df$dnv_conf_level), 1, 0)
table(denovo) # 

# SIFT
table(var_anno_2train_df$SIFT_pred)
sift_D <- ifelse(var_anno_2train_df$SIFT_pred == "D" & !is.na(var_anno_2train_df$SIFT_pred), 1, 0)
table(sift_D)

# pph2
table(var_anno_2train_df$Polyphen2_HDIV_pred)
pph2_D <- ifelse(var_anno_2train_df$Polyphen2_HDIV_pred == "D" & !is.na(var_anno_2train_df$Polyphen2_HDIV_pred), 1, 0)
pph2_P <- ifelse(var_anno_2train_df$Polyphen2_HDIV_pred == "P" & !is.na(var_anno_2train_df$Polyphen2_HDIV_pred), 1, 0)

# mutationTaster & mutationAssessor
table(var_anno_2train_df$MutationTaster_pred)
table(var_anno_2train_df$MutationAssessor_pred)

# FATHMM
table(var_anno_2train_df$FATHMM_pred)
FATHMM_D <- ifelse(var_anno_2train_df$FATHMM_pred == "D" & !is.na(var_anno_2train_df$FATHMM_pred), 1, 0)

# PROVEAN
table(var_anno_2train_df$PROVEAN_pred)
PROVEAN_D <- ifelse(var_anno_2train_df$PROVEAN_pred == "D" & !is.na(var_anno_2train_df$PROVEAN_pred), 1, 0)

# metaSVM
table(var_anno_2train_df$MetaSVM_pred)
metaSVM_D <- ifelse(var_anno_2train_df$MetaSVM_pred == "D" & !is.na(var_anno_2train_df$MetaSVM_pred), 1, 0)

# metaLR
table(var_anno_2train_df$MetaLR_pred)
metaLR_D <- ifelse(var_anno_2train_df$MetaLR_pred == "D" & !is.na(var_anno_2train_df$MetaLR_pred), 1, 0)

# MetaRNN
table(var_anno_2train_df$MetaRNN_pred)
metaRNN_D <- ifelse(var_anno_2train_df$MetaRNN_pred == "D" & !is.na(var_anno_2train_df$MetaRNN_pred), 1, 0)

# metaSVM and metaLR are significantly correlated with cor=0.86

# MCAP
table(var_anno_2train_df$`M-CAP_pred`)
MCAP_D <- ifelse(var_anno_2train_df$`M-CAP_pred` == "D" & !is.na(var_anno_2train_df$`M-CAP_pred`), 1, 0)

# MVP 
var_anno_2train_df$MVP_score <- as.numeric(var_anno_2train_df$MVP_score)
summary(var_anno_2train_df$MVP_score)
MVP_above075 <- ifelse(var_anno_2train_df$MVP_score >= 0.75 & !is.na(var_anno_2train_df$MVP_score), 1, 0)

# based on cor.test, FATHMM, MCAP and MVP has a correlation of 0.45

# MPC
var_anno_2train_df$MPC_score <- as.numeric(var_anno_2train_df$MPC_score)
summary(var_anno_2train_df$MPC_score)
MPC_above2 <- ifelse(var_anno_2train_df$MPC_score >= 2 & !is.na(var_anno_2train_df$MPC_score), 1, 0)
MPC_1to2 <- ifelse(var_anno_2train_df$MPC_score >= 1 & var_anno_2train_df$MPC_score < 2 & !is.na(var_anno_2train_df$MPC_score), 1, 0)

# dbscSNV_ADA
var_anno_2train_df$dbscSNV_ADA_SCORE <- as.numeric(var_anno_2train_df$dbscSNV_ADA_SCORE)
summary(var_anno_2train_df$dbscSNV_ADA_SCORE)
dbscSNV_ADA_above06 <- ifelse(var_anno_2train_df$dbscSNV_ADA_SCORE >= 0.60 & !is.na(var_anno_2train_df$dbscSNV_ADA_SCORE), 1, 0)

# # dbscSNV_RF
# var_anno_2train_df$dbscSNV_RF_SCORE <- as.numeric(var_anno_2train_df$dbscSNV_RF_SCORE)
# summary(var_anno_2train_df$dbscSNV_RF_SCORE)
# dbscSNV_RF_above06 <- ifelse(var_anno_2train_df$dbscSNV_RF_SCORE >= 0.60 & !is.na(var_anno_2train_df$dbscSNV_RF_SCORE), 1, 0)
# # splicing variant prediction ADA and RF are highly correlation with cor=0.9568
# cor.test(dbscSNV_ADA_above06, dbscSNV_RF_above06)

# CADD
var_anno_2train_df$CADD_phred <- as.numeric(var_anno_2train_df$CADD_phred)
summary(var_anno_2train_df$CADD_phred)
CADD_above30 <- ifelse(var_anno_2train_df$CADD_phred >= 30 & !is.na(var_anno_2train_df$CADD_phred), 1, 0)
CADD_25to30 <- ifelse(var_anno_2train_df$CADD_phred >= 25 & var_anno_2train_df$CADD_phred < 30 & !is.na(var_anno_2train_df$CADD_phred), 1, 0)

# pLI
var_anno_2train_df$pLi.refGene <- as.numeric(var_anno_2train_df$pLi.refGene)
summary(var_anno_2train_df$pLi.refGene)
pLi_above0.99 <- ifelse(var_anno_2train_df$pLi.refGene >= 0.99 & !is.na(var_anno_2train_df$pLi.refGene), 1, 0)
pLi_0.5to0.99 <- ifelse(var_anno_2train_df$pLi.refGene < 0.99 & var_anno_2train_df$pLi.refGene >= 0.5 & !is.na(var_anno_2train_df$pLi.refGene), 1, 0)

# SFARI gene 
SFARI_gene_lst_df <- read.csv(file = "../reference_data/SFARI_gene/SFARI-Gene_genes_01-11-2022release_03-22-2022export.csv")
head(SFARI_gene_lst_df)

idx_SFARI_gene_s1 <- SFARI_gene_lst_df$gene.score == 1
SFARI_gene_s1 <- ifelse(var_anno_2train_df$gene_name %in% SFARI_gene_lst_df$gene.symbol[idx_SFARI_gene_s1], 1, 0)
idx_SFARI_gene_s2 <- SFARI_gene_lst_df$gene.score == 2
SFARI_gene_s2 <- ifelse(var_anno_2train_df$gene_name %in% SFARI_gene_lst_df$gene.symbol[idx_SFARI_gene_s2], 1, 0)
idx_SFARI_gene_s3 <- SFARI_gene_lst_df$gene.score == 3
SFARI_gene_s3 <- ifelse(var_anno_2train_df$gene_name %in% SFARI_gene_lst_df$gene.symbol[idx_SFARI_gene_s3], 1, 0)
idx_SFARI_gene_S <- SFARI_gene_lst_df$gene.score == "S"
SFARI_gene_S <- ifelse(var_anno_2train_df$gene_name %in% SFARI_gene_lst_df$gene.symbol[idx_SFARI_gene_S], 1, 0)

head(var_anno_2train_df)
var_anno_dummy_df <- as.data.frame(cbind(var_anno_2train_df[, c(1:5, 31)], 
                                         exonic_func_PTVs, exonic_func_missense, denovo, 
                                         sift_D, pph2_P, pph2_D, FATHMM_D, PROVEAN_D, metaSVM_D,
                                         MCAP_D, MVP_above075, MPC_1to2, MPC_above2, 
                                         CADD_above30, CADD_25to30, pLi_above0.99, pLi_0.5to0.99,
                                         SFARI_gene_s1, SFARI_gene_s2, SFARI_gene_s3, SFARI_gene_S))
head(var_anno_dummy_df)
dim(var_anno_dummy_df)

