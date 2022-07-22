rm(list = ls())
# workingDir <- "/Users/shining17wq/Google Drive/MorrowLab/SFARI/SPARK/job"
workingDir <- "/gpfs/data/emorrow/silvio/Qing_Wu/SFARI/batch_jobs"
setwd(workingDir)
getwd()

library(data.table)
library(stringr)
library(dplyr)
library(plyr)
options(stringsAsFactors = F, scipen = 999)


# de novo variants --------------------------------------------------------

dnv_raw <- read.table(file = "../Regeneron/wes1_wes2_combined.deepvariant.dnvs.tsv", header = T)
head(dnv_raw)

table(dnv_raw$confidence)

# generate annovar input --------------------------------------------------

var_df <- read.table(file = "../hail_out/SPARK/wes1_wes2_combined.gatk.rare1pct_variants_het_info.txt", sep = '\t', header = F)
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

write.table(avinput_df, file = "../hail_out/SPARK/wes1_wes2_combined.gatk.rare1pct_variants_het_info.avinput", sep = "\t", quote = F, col.names = F, row.names = F)


# HUGO genes --------------------------------------------------------------

var_gene_anno_df <- fread(file = "../hail_out/SPARK/wes1_wes2_combined.gatk.rare1pct_variants_het_info.hg38_multianno.csv")
head(var_gene_anno_df)

length(unique(var_gene_anno_df$Gene.refGene)) 
length(unique(var_gene_anno_df$Gene.knownGene)) 
length(unique(var_gene_anno_df$Gene.ensGene)) 

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
idx_genes_name <- grep(";", gene_name_df$gene_name) 
intra_genes_name <- grep(";", gene_name_df$gene_name, value = T)
# split multiple names by ;
intra_genes_lst <- str_split(intra_genes_name, ";")
intra_genes_df <- plyr::ldply(intra_genes_lst, rbind)
num_unique_gene_name <- apply(intra_genes_df, 1, function(x) length(unique(x)))
idx_identical_genes <- (num_unique_gene_name - 1 == 1)
table(idx_identical_genes) # 538

replaced_gene_name <- intra_genes_df$`1`[idx_identical_genes]
gene_name_df$gene_name[idx_genes_name[idx_identical_genes]] <- replaced_gene_name

length(unique(gene_name_df$gene_name)) # 25060

unique_gene_name <- unique(gene_name_df$gene_name)
# unique_gene_name

variants2gene_df <- data.frame(Chr = var_gene_anno_df$Chr,
                               Start = var_gene_anno_df$Start,
                               End = var_gene_anno_df$End,
                               Ref = var_gene_anno_df$Ref,
                               Alt = var_gene_anno_df$Alt,
                               Gene = gene_name_df$gene_name)
head(variants2gene_df)
write.csv(variants2gene_df, file = "../hail_out/SPARK/wes1_wes2_combined.gatk.rare1pct_variants_het_info.gene_converted.csv")



# annovar output ----------------------------------------------------------

var_anno_df <- fread(file = "../hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het_info.hg38_multianno.csv")
head(var_anno_df)

gene_name_df <- read.csv(file = "../hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het_info.gene_converted.csv", header = T)
head(gene_name_df)
length(unique(gene_name_df$Gene))
var_anno_df$gene_name<- gene_name_df$Gene

dnv_conf_level <- read.table(file = "../hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.w_dnv.tsv", sep = '\t', header = T)
head(dnv_conf_level)
table(dnv_conf_level$dnv_confidence) # high: 25509; medium: 1674; low: 1947

var_anno_df$dnv_conf_level<- dnv_conf_level$dnv_confidence
table(var_anno_df$dnv_conf_level)
# 25509   1947   1674

# comparison between two caller -------------------------------------------

var_anno_df_gatk <- fread(file = "../hail_out/SPARK/wes1_wes2_combined.gatk.rare1pct_variants_het_info.hg38_multianno.csv")
head(var_anno_df_gatk)

gene_name_df_gatk <- read.csv(file = "../hail_out/SPARK/wes1_wes2_combined.gatk.rare1pct_variants_het_info.gene_converted.csv", header = T)
head(gene_name_df_gatk)
length(unique(gene_name_df_gatk$Gene))
var_anno_df_gatk$gene_name<- gene_name_df_gatk$Gene

var_deepvariant <- paste(var_anno_df$Chr, var_anno_df$Start, var_anno_df$End, var_anno_df$Ref, var_anno_df$Alt, sep = ":")
var_gatk <- paste(var_anno_df_gatk$Chr, var_anno_df_gatk$Start, var_anno_df_gatk$End, var_anno_df_gatk$Ref, var_anno_df_gatk$Alt, sep = ":")

length(var_gatk) - length(var_deepvariant)

idx_deepvariant_matched_gatk <- var_deepvariant %in% var_gatk
table(idx_deepvariant_matched_gatk)

idx_exonic_deepvariant_matched_gatk <- idx_deepvariant_matched_gatk & idx_var2incl
table(idx_exonic_deepvariant_matched_gatk)

idx_impactful_deepvariant_matched_gatk <- idx_deepvariant_matched_gatk & idx_var_PTVs_MisAB
table(idx_impactful_deepvariant_matched_gatk)

idx_non_impactful_deepvariant_matched_gatk <- idx_deepvariant_matched_gatk & idx_var_non_impactful
table(idx_non_impactful_deepvariant_matched_gatk)

61766-4188
73739-14788

selected_gene_name <- read.table("/users/qwu24/data/silvio/Qing_Wu/SFARI/batch_jobs/python_script/machine_learning_script/tpot_out/wes1_wes2_combined.deepvariant.feature_name_exonic_selPct4.txt", header = T)
idx_selected_genes_deepvariant <- var_anno_df$gene_name %in% selected_gene_name$selected_features[-1]
table(idx_selected_genes_deepvariant)
# DeepVariant: 3752277  105248
idx_selected_genes_gatk <- var_anno_df_gatk$gene_name %in% selected_gene_name$selected_features[-1]
# GATK: 3959943  112279

idx_selected_gene_deepvariant_matched_gatk <- idx_selected_genes_deepvariant & idx_deepvariant_matched_gatk
table(idx_selected_gene_deepvariant_matched_gatk)
# DeepVariant: 3757871   99654 (105248-99654=5594)
# GATK: (112279-99654=12625)

idx_selected_genes_impactful_deepvariant_matched_gatk <- idx_selected_genes_deepvariant & idx_deepvariant_matched_gatk & idx_var_PTVs_MisAB
table(idx_selected_genes_impactful_deepvariant_matched_gatk)
# DeepVariant: 3849001    8524
table(idx_selected_genes_deepvariant &  idx_var_PTVs_MisAB)
# Deepvariant: 3848710    (8815-8524=291)

idx_var_PTVs_gatk <- (var_anno_df_gatk$ExonicFunc.refGene == "stopgain" | var_anno_df_gatk$ExonicFunc.refGene == "frameshift substitution" | var_anno_df_gatk$ExonicFunc.refGene == "splicing")
idx_MPC_above1_gatk <- var_anno_df_gatk$MPC_score >= 1 & !is.na(var_anno_df_gatk$MPC_score)
idx_var_PTVs_MisAB_gatk <- idx_var_PTVs_gatk | idx_MPC_above1_gatk
table(idx_selected_genes_gatk & idx_var_PTVs_MisAB_gatk)
# GATK: 4061265   (10957-8524=2433)

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

# CADD score --------------------------------------------------------------
library(ggplot2)

var_anno_df$CADD_phred <- as.numeric(var_anno_df$CADD_phred)
summary(var_anno_df$CADD_phred)

p <- ggplot(var_anno_df[idx_func_1,], aes(x=CADD_phred)) + 
  geom_histogram(color="black", fill="white") + 
  theme_classic() +
  labs(x = "CADD_phred (exonic)")

pdf("CADD_exonic_ggplot.pdf")
print(p)
dev.off()

p <- ggplot(var_anno_df, aes(x=CADD_phred, color = ExonicFunc.refGene)) +
  geom_histogram(fill="white", position="dodge") +
  theme_classic()

pdf("CADD_exonic_func_ggplot.pdf")
print(p)
dev.off()


# continue variants index generation --------------------------------------

# function to write index of which indicate variant row to exclude from original variant by sample matrix
get_idx_var2exclude <- function(idx_var2incl, variant_category_name){
  idx_var2excl <- (!idx_var2incl)*c(1:nrow(var_anno_df))
  idx_var2excl <- idx_var2excl[idx_var2excl != 0]
  print(length(idx_var2excl))
  
  if (variant_category_name == "") {
    out_file_name <- "../hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2exclude.txt"
  } else {
    out_file_name <- paste("../hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2exclude", variant_category_name, "txt", sep = ".")
  }
  print(out_file_name)
  
  write.table(idx_var2excl, file = out_file_name, sep = "\t", quote = F, col.names = F, row.names = F)
}

# include exonic and splicing variant, but not synonymous variants
idx_var2incl <- (idx_func_1 | idx_func_2 | idx_func_3) & (!idx_exonic_func_1)
table(idx_var2incl) 
# DeepVariant: 2388489+1469036 (F/T)
# GATK: 2486676+1585546 (F/T)
get_idx_var2exclude(idx_var2incl, "")

1585546-1469036

# include PTVs only
idx_var_PTVs <- (var_anno_df$ExonicFunc.refGene == "stopgain" | var_anno_df$ExonicFunc.refGene == "frameshift substitution" | var_anno_df$ExonicFunc.refGene == "splicing")
table(idx_var_PTVs)
# DeepVariant: 3761745+95780
# GATK: 3932170+140052
get_idx_var2exclude(idx_var_PTVs, "PTVs")

# include PTVs and missense MPC > 2
var_anno_df$MPC_score <- as.numeric(var_anno_df$MPC_score)
idx_MPC_above2 <- var_anno_df$MPC_score >= 2 & !is.na(var_anno_df$MPC_score)
idx_var_PTVs_MisA <- idx_var_PTVs | idx_MPC_above2
table(idx_var_PTVs_MisA)
# DeepVaint: 3738592+118933
# GATK: 3905302+166920
get_idx_var2exclude(idx_var_PTVs_MisA, "PTVs_MisA")

# include PTVs and missense MPC > 1
idx_MPC_above1 <- var_anno_df$MPC_score >= 1 & !is.na(var_anno_df$MPC_score)
idx_var_PTVs_MisAB <- idx_var_PTVs | idx_MPC_above1
table(idx_var_PTVs_MisAB)
# DeepVaint: 3582552+274973
# GATK: 3739671+332551
get_idx_var2exclude(idx_var_PTVs_MisAB, "PTVs_MisAB")

# include missense MPC > 2
idx_MPC_above2 <- var_anno_df$MPC_score >= 2 & !is.na(var_anno_df$MPC_score)
table(idx_MPC_above2) 
# DeepVariantL: 3834370+23155 (F/T)
# GATK: 4045351+26871 
get_idx_var2exclude(idx_var2incl = idx_MPC_above2, "MisA")

# include missense MPC between 1 and 2
idx_MPC_bw12 <- (var_anno_df$MPC_score <= 2 & var_anno_df$MPC_score >= 1) & !is.na(var_anno_df$MPC_score)
table(idx_MPC_bw12) 
# DeepVariant: 3701392  156133
# GATK: 3906492  165730
get_idx_var2exclude(idx_var2incl = idx_MPC_bw12, "MisB")

# include missense MPC above 1
idx_MPC_above1 <- var_anno_df$MPC_score >= 1 & !is.na(var_anno_df$MPC_score)
table(idx_MPC_above1) 
# DeepVaraint: 3678310  179215
# GATK: 3879700  192522
get_idx_var2exclude(idx_var2incl = idx_MPC_above1, "MisAB")

# include other variants, not PTVs, not MisAB
idx_var2incl <- (idx_func_1 | idx_func_2 | idx_func_3) & (!idx_exonic_func_1)
idx_var_PTVs <- (var_anno_df$ExonicFunc.refGene == "stopgain" | var_anno_df$ExonicFunc.refGene == "frameshift substitution" | var_anno_df$ExonicFunc.refGene == "splicing")
idx_MPC_above1 <- var_anno_df$MPC_score >= 1 & !is.na(var_anno_df$MPC_score)
idx_var_non_impactful <- idx_var2incl & !(idx_MPC_above1 | idx_var_PTVs)
table(idx_var2incl) # 2388489 1469036
table(idx_MPC_above1 | idx_var_PTVs) # 3582552  274973
table(idx_var_non_impactful) 
# DeepVariant: 2663218 1194307
# GATK: 2818964 1253258
get_idx_var2exclude(idx_var2incl = idx_var_non_impactful, "non_impactful")


# generate VxA matrix -----------------------------------------------------

get_var_anno_by_variant_category <- function(var_anno_df, idx_var2incl, variant_category_name){
  
  var_anno_var2incl_df <- var_anno_df[idx_var2incl,]
  print(dim(var_anno_var2incl_df))

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
  print(paste("Length of new var2gene list", length(new_var2gene_lst), sep = ":"))
  print(paste("Length of new gene name list", length(new_gene_names), sep = ":"))
  print(paste("Length of original variant 2 gene list", length(var2gene_lst), sep = ":"))

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
  print(is.unsorted(order(smallest_value_lst)))

  var2gene_df <- plyr::ldply(new_var2gene_lst_sorted, rbind)
  var2gene_df[is.na(var2gene_df)] <- ""

  var2gene_df <- data.frame(cbind(new_gene_names_sorted, var2gene_df))
  print(dim(var2gene_df))

  if (variant_category_name == "") {
    out_file_name <- "../hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.unconsecutive_indices_splited.txt"
  } else {
    out_file_name <- paste("../hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.unconsecutive_indices_splited", variant_category_name, "txt", sep = ".")
  }
  
  write.table(var2gene_df, file = out_file_name, sep = "\t", quote = F, col.names = F, row.names = F)

  # get the full list of variants to gene indices list ----------------------

  var2gene_df <- plyr::ldply(var2gene_lst, rbind)
  var2gene_df[is.na(var2gene_df)] <- ""
  print(dim(var2gene_df))

  var2gene_df <- data.frame(cbind(unique_gene, var2gene_df))

  if (variant_category_name == "") {
    out_file_name <- "../hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.txt"
  } else {
    out_file_name <- paste("../hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df", variant_category_name, "txt", sep = ".")
  }
  
  write.table(var2gene_df, file = out_file_name, sep = "\t", quote = F, col.names = F, row.names = F)

  return(var_anno_var2incl_df)
}


var_anno_var2incl_df <- get_var_anno_by_variant_category(var_anno_df, idx_var2incl, "")
var_anno_var2incl_df <- get_var_anno_by_variant_category(var_anno_df, idx_var_PTVs, "PTVs")
var_anno_var2incl_df <- get_var_anno_by_variant_category(var_anno_df, idx_var_PTVs_MisA, "PTVs_MisA")
var_anno_var2incl_df <- get_var_anno_by_variant_category(var_anno_df, idx_var_PTVs_MisAB, "PTVs_MisAB")

var2_anno_var2incl_df <- get_var_anno_by_variant_category(var_anno_df, idx_var2incl = idx_MPC_above2, variant_category_name="MisA")
var2_anno_var2incl_df <- get_var_anno_by_variant_category(var_anno_df, idx_var2incl = idx_MPC_bw12, variant_category_name="MisB")
var2_anno_var2incl_df <- get_var_anno_by_variant_category(var_anno_df, idx_var2incl = idx_MPC_above1, variant_category_name="MisAB")

var2_anno_var2incl_df <- get_var_anno_by_variant_category(var_anno_df, idx_var2incl = idx_var_non_impactful, variant_category_name="non_impactful")


# variants annotations ----------------------------------------------------

# idx_var2incl <- idx_var_PTVs
# idx_var2incl <- idx_var_PTVs_MisA
# idx_var2incl <- idx_var_PTVs_MisAB

var_anno_var2incl_df <- var_anno_df[idx_var2incl,]
dim(var_anno_var2incl_df)

names(var_anno_var2incl_df)
# head(var_anno_var2incl_df)

var_anno_2train_df <- var_anno_var2incl_df[, c(1:7, 9, 11, 26, 29, 32, 35, 38, 41, 44, 47, 50, 55, 58, 61, 64, 69, 71, 95, 96, 100, 103, 131:134)]
head(var_anno_2train_df)


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

# save
write.table(var_anno_dummy_df, file = "../hail_out/SSC_WES_3/SSC_WES_3.vqsr.vcf_families_merged.pass_QC.rare1pct_variants_het.var_anno_dummy_df.txt", sep = "\t", quote = F, col.names = F, row.names = F)

# pLi and MPC annotation only ---------------------------------------------

head(var_anno_2train_df)
var_anno_dummy_df_def <- as.data.frame(cbind(var_anno_2train_df[, c(1:5, 31)], 
                                         exonic_func_PTVs, exonic_func_missense, denovo, 
                                         MPC_1to2, MPC_above2, pLi_above0.99, pLi_0.5to0.99,
                                         SFARI_gene_s1, SFARI_gene_s2, SFARI_gene_s3, SFARI_gene_S))
head(var_anno_dummy_df_def)
dim(var_anno_dummy_df_def)

# save
write.table(var_anno_dummy_df_def, file = "../hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.var_anno_dummy_df_pLi_MPC.txt", sep = "\t", quote = F, col.names = F, row.names = F)



# PTVs and MPC2 SFARI gene only
PTVs_in_SFARI_gene <- exonic_func_PTVs * SFARI_gene
MPC2_in_SFARI_gene <- MPC_above2 * SFARI_gene
# var_anno_2train_df[1:4]
var_anno_dummy_df <- as.data.frame(cbind(var_anno_2train_df[, c(1:5, 23)], 
                                         PTVs_in_SFARI_gene, MPC2_in_SFARI_gene))
head(var_anno_dummy_df)

# save
write.csv(var_anno_dummy_df, file = "../hail_out/SPARK/SPARK_vcf_families_merged_rare_var_anno_dummy_df_0216_PTVs_MPC2_in_SFARIgene.csv", col.names = T, row.names = F)


# SFARI gene --------------------------------------------------------------

SFARI_gene_lst_df <- read.csv(file = "../reference_data/SFARI_gene/SFARI-Gene_genes_01-11-2022release_03-22-2022export.csv")
head(SFARI_gene_lst_df)
SFARI_genes <- SFARI_gene_lst_df$gene.symbol # 1031

# DDG2P gene
DDG2P_gene_lst_df <- read.csv(file = "../reference_data/DDG2P/DDG2P_7_4_2022.csv")
head(DDG2P_gene_lst_df)
DDG2P_genes <- DDG2P_gene_lst_df$gene.symbol # 2580

SFARI_genes_01 <- ifelse(var_anno_2train_df$gene_name %in% SFARI_genes, 1, 0)
DDG2P_genes_01 <- ifelse(var_anno_2train_df$gene_name %in% DDG2P_genes, 1, 0)

var_anno_dummy_df$SFARI_genes_01 = SFARI_genes_01
var_anno_dummy_df$DDG2P_genes_01 = DDG2P_genes_01

head(var_anno_dummy_df)

denovo_df <- var_anno_dummy_df[denovo == 1,]
head(denovo_df)

# SFARI_genes and DDG2P genes
table(denovo_df$SFARI_genes_01) # 1/(7425/1189) (0/1) 
table(denovo_df$DDG2P_genes_01) # 1/(6765/1849) (0/1)

table(var_anno_dummy_df$SFARI_genes_01) # 1/(1337945/131091)
table(var_anno_dummy_df$DDG2P_genes_01) # 1/(1232225/236811)

length(unique(denovo_df$gene_name)) # 5811
length(unique(denovo_df$gene_name[denovo_df$SFARI_genes_01 == 1])) # 544
length(unique(denovo_df$gene_name[denovo_df$DDG2P_genes_01 == 1])) # 981

length(unique(var_anno_dummy_df$gene_name)) # 19117
length(unique(var_anno_dummy_df$gene_name[(var_anno_dummy_df$SFARI_genes_01 == 1) & (var_anno_dummy_df$denovo == 0)])) # 1019
length(unique(var_anno_dummy_df$gene_name[(var_anno_dummy_df$DDG2P_genes_01 == 1) & (var_anno_dummy_df$denovo == 0)])) # 2218

# PTV
table(denovo_df$exonic_func_PTVs) # 1/(9460/1422)
table(var_anno_dummy_df$exonic_func_PTVs[(var_anno_dummy_df$denovo == 0)]) # 1/(2386830/273088)

# MPC > 2
table(denovo_df$MPC_above2) # 1/(7577/1037)
table(var_anno_dummy_df$MPC_above2[(var_anno_dummy_df$denovo == 0)]) # 1/(1365679/94743)

# MPC 1 ~ 2
table(denovo_df$MPC_1to2) # 1/(7306/1308)
table(var_anno_dummy_df$MPC_1to2[(var_anno_dummy_df$denovo == 0)]) # 1/(1305889/154533)

# pLi_above0.99
table(denovo_df$pLi_above0.99) # 1/(6749/1865)
table(var_anno_dummy_df$pLi_above0.99[(var_anno_dummy_df$denovo == 0)]) # 1/(1252640/207782)

# pLi_0.5to0.99
table(denovo_df$pLi_0.5to0.99) # 1/(7285/1329)
table(var_anno_dummy_df$pLi_0.5to0.99[(var_anno_dummy_df$denovo == 0)]) # 1/(1271679/188743)







