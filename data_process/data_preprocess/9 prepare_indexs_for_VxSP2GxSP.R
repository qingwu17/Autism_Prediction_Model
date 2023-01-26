rm(list = ls())
workingDir <- "/path_to_working_directory/"
setwd(workingDir)
getwd()

library(data.table)
library(stringr)
library(dplyr)
library(plyr)
options(stringsAsFactors = F, scipen = 999)

# annovar output ----------------------------------------------------------

dir_name <- "/path_to/"
sample_set_name <- "wes1_wes2_combined"

var_anno_filename <- paste(dir_name, sample_set_name, ".deepvariant.rare1pct_variants_het_info.hg38_multianno.csv", sep = "")
var_anno_df <- fread(file = var_anno_filename)
head(var_anno_df)

gene_name_filename <- paste(dir_name, sample_set_name, ".deepvariant.rare1pct_variants_het_info.gene_converted.csv", sep = "")
gene_name_df <- read.csv(file = gene_name_filename, header = T)
head(gene_name_df)
length(unique(gene_name_df$Gene))
var_anno_df$gene_name<- gene_name_df$Gene

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


# index to exclude variants from VxSP matrix ------------------------------

# function to export row index, indicating which variant row to exclude from original variant by sample matrix
get_idx_var2exclude <- function(idx_var2incl, variant_category_name){
  idx_var2excl <- (!idx_var2incl)*c(1:nrow(var_anno_df))
  idx_var2excl <- idx_var2excl[idx_var2excl != 0]
  print(length(idx_var2excl))
  
  if (variant_category_name == "") {
    out_file_name <- paste(dir_name, sample_set_name, ".deepvariant.rare1pct_variants_het.idx_var2exclude.txt", sep = "")
  } else {
    out_file_name <- paste(dir_name, sample_set_name, ".deepvariant.rare1pct_variants_het.idx_var2exclude.", variant_category_name, ".txt", sep = "")
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


# index of variants to gene conversion for VxSP to GxSP -------------------

# function to generate variants index of a gene, some variants are consecutive, some are not
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
    out_file_name <- paste(dir_name, sample_set_name, ".deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.unconsecutive_indices_splited.txt", sep = "")
  } else {
    out_file_name <- paste(dir_name, sample_set_name, ".deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.unconsecutive_indices_splited.", variant_category_name, ".txt", sep = "")
  }
  print(out_file_name)
  
  write.table(var2gene_df, file = out_file_name, sep = "\t", quote = F, col.names = F, row.names = F)
  
  # get the full list of variants to gene indices list ----------------------
  
  var2gene_df <- plyr::ldply(var2gene_lst, rbind)
  var2gene_df[is.na(var2gene_df)] <- ""
  print(dim(var2gene_df))
  
  var2gene_df <- data.frame(cbind(unique_gene, var2gene_df))
  
  if (variant_category_name == "") {
    out_file_name <- paste(dir_name, sample_set_name, ".deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.txt", sep = "")
  } else {
    out_file_name <- paste(dir_name, sample_set_name, ".deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.", variant_category_name, ".txt", sep = "")
  }
  print(out_file_name)
  
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