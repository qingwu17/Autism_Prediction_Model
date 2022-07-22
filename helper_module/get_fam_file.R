rm(list = ls())
workingDir <- "/Users/shining17wq/Google Drive/MorrowLab/SFARI/SPARK/job"
setwd(workingDir)
getwd()

Packages <- c("stringr", "plyr", "dplyr")
lapply(Packages, require, character.only = TRUE)

library(ggplot2)
library(ggpubr)

options(stringsAsFactors = F)

# input file --------------------------------------------------------------

collection_version3_df_0 <- read.csv(file = "../phenotype/SPARK_Collection_Version3/individuals.csv", header = T)
collection_version4_df_0 <- read.csv(file = "../phenotype/SPARK_Collection_Version4/individuals.csv", header = T)
collection_version5_df_0 <- read.csv(file = "../phenotype/SPARK_Collection_Version5/individuals.csv", header = T)

# sample_id_from_vcf_family_df <- read.table(file = "../Regeneron/sample_id.tsv", sep = '\t', header = TRUE)

wes1_27281_exome_sample_id_df <- read.table(file = "../Regeneron/wes1_27281_exome_sample_id.tsv", sep = '\t', header = TRUE)
head(wes1_27281_exome_sample_id_df)
wes2_15995_exome_sample_id_df <- read.table(file = "../Regeneron/wes2_15995_exome_sample_id.tsv", sep = '\t', header = TRUE)
head(wes2_15995_exome_sample_id_df)


# edit collection data frame into a uniform format ------------------------

names(collection_version3_df_0)
names(collection_version4_df_0)
names(collection_version5_df_0)

# cv3
# add biofather id and biomother id to cv3
collection_version3_df_0$biofather_id <- ""
collection_version3_df_0$biomother_id <- ""

fam_id <- unique(collection_version3_df_0$family_id)
idx_bool_fa <- collection_version3_df_0$role == "Father"
idx_bool_mo <- collection_version3_df_0$role == "Mother"
idx_bool_child <- (collection_version3_df_0$role == "Proband") | (collection_version3_df_0$role == "Sibling")

for (i in 1:length(fam_id)){

  idx_fam_i <- collection_version3_df_0$family_id == fam_id[i]

  # assign biomother id
  if (any(idx_fam_i&idx_bool_mo)) {
    collection_version3_df_0$biomother_id[idx_fam_i&idx_bool_child] <- collection_version3_df_0$subject_sp_id[idx_fam_i&idx_bool_mo]
  } 
  # assign biofather id
  if (any(idx_fam_i&idx_bool_fa)) {
    collection_version3_df_0$biofather_id[idx_fam_i&idx_bool_child] <- collection_version3_df_0$subject_sp_id[idx_fam_i&idx_bool_fa]
  }
}

collection_version3_df_0 <- collection_version3_df_0[c(1:2, 31:32, 3:30)]

# cv5
collection_version5_df_0$family_id <- collection_version5_df_0$family_sf_id
collection_version5_df_0$biofather_id <- collection_version5_df_0$biofather_sp_id
collection_version5_df_0$biomother_id <- collection_version5_df_0$biomother_sp_id

# merge pheno with geno ---------------------------------------------------

# merge phenotype v4 v5, because these two data frame has biofather_id and biomother_id
comm_names_v45 <- intersect(names(collection_version4_df_0), names(collection_version5_df_0))
collection_version45_df_0 <- full_join(collection_version4_df_0,
                                       collection_version5_df_0,
                                       by = comm_names_v45)
names(collection_version45_df_0)

# find common names of three basic_medical_screening data frame
comm_col_names <- Reduce(intersect, list(names(collection_version4_df_0),
                                         names(collection_version5_df_0)))
# select columns by common colname
cv4_comm_names_df <- subset(collection_version4_df_0, select = comm_col_names)
cv5_comm_names_df <- subset(collection_version5_df_0, select = comm_col_names)

# rbind two dataframe
combined_comm_names_df <- rbind(cv4_comm_names_df, cv5_comm_names_df)
# find and remove duolicated sp by subject sp id
combined_comm_names_uniq_sp_df <- combined_comm_names_df[!duplicated(combined_comm_names_df$subject_sp_id),]
dim(combined_comm_names_uniq_sp_df)
# 252147     35

names(combined_comm_names_uniq_sp_df)

# sample_id_from_vcf_family <- sample_id_from_vcf_family_df$s

sample_id_wes1_wes2_combined <- unique(c(wes1_27281_exome_sample_id_df$s, wes2_15995_exome_sample_id_df$s))

idx_0 <- match(sample_id_wes1_wes2_combined, combined_comm_names_uniq_sp_df$subject_sp_id)
idx <- idx_0[!is.na(idx_0)]

pheno_SPARK_df_0 <- combined_comm_names_uniq_sp_df[idx,c(1:6)]
head(pheno_SPARK_df_0)
names(pheno_SPARK_df_0)
dim(pheno_SPARK_df_0)

missing_SP_ID <- setdiff(sample_id_wes1_wes2_combined, combined_comm_names_uniq_sp_df$subject_sp_id)
pheno_SPARK_df_rest <- collection_version3_df_0[match(missing_SP_ID, collection_version3_df_0$subject_sp_id),c(1:4, 8, 12)]
names(pheno_SPARK_df_rest)

pheno_SPARK_df <- data.frame(rbind(pheno_SPARK_df_0, pheno_SPARK_df_rest))

# analysis ----------------------------------------------------------------

fam_df <- data.frame(FID = pheno_SPARK_df$family_id, IID = pheno_SPARK_df$subject_sp_id,
                     fa_id = pheno_SPARK_df$biofather_id, mo_id = pheno_SPARK_df$biomother_id, 
                     sex = tolower(pheno_SPARK_df$sex), pheno = pheno_SPARK_df$asd*1+1)

fam_df$sex <- ifelse(fam_df$sex =="male", 1, ifelse(fam_df$sex == "female", 2, 0))
fam_df[is.na(fam_df)] <- 0
fam_df[fam_df == ""] <- 0

# check
table(pheno_SPARK_df$asd) # 21357 21870 (F/T)
table(is.na(fam_df))
head(pheno_SPARK_df$asd)

head(fam_df)

write.table(fam_df, file = "../Regeneron/wes1_wes2_combined.fam", sep = "\t", quote = F, col.names = F, row.names = F)


