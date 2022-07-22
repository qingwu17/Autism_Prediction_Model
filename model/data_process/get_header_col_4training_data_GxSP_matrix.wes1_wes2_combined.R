rm(list = ls())
workingDir <- "/Users/shining17wq/Google Drive/MorrowLab/SFARI/SPARK/job"
setwd(workingDir)
getwd()


library(data.table)
library(stringr)
library(dplyr)
library(plyr)
options(stringsAsFactors = F)

# phenotype sp annotation -------------------------------------------------
sample_id_pass_QC <- read.table(file = "../pub/WES12/wes1_wes2_combined.gatk.sample_id.tsv", header = T)
sp_pass_QC_is_case <- read.table(file = "../pub/WES12/wes1_wes2_combined.gatk.is_female.tsv", header = T)
sp_pass_QC_is_female <- read.table(file = "../pub/WES12/wes1_wes2_combined.gatk.is_case.tsv", header = T)
head(sample_id_pass_QC)
head(sp_pass_QC_is_case)
head(sp_pass_QC_is_female)

sp_pass_QC_is_case$is_case <- as.integer(as.logical(sp_pass_QC_is_case$X__uid_18))
sp_pass_QC_is_female$is_female <- as.integer(as.logical(sp_pass_QC_is_female$X__uid_19))

head(sp_pass_QC_is_case)
head(sp_pass_QC_is_female)

table(sp_pass_QC_is_case$is_case)

identical(sample_id_pass_QC$s, sp_pass_QC_is_case$s)
identical(sample_id_pass_QC$s, sp_pass_QC_is_female$s)
identical(sp_pass_QC_is_case$s, sp_pass_QC_is_female$s)

variant_sample_matrix_sp_id_SSC <- data.frame(sample_id = sample_id_pass_QC$s, is_case = sp_pass_QC_is_case$is_case, is_female = sp_pass_QC_is_female$is_female)
head(variant_sample_matrix_sp_id_SSC)
dim(variant_sample_matrix_sp_id_SSC)
write.table(variant_sample_matrix_sp_id_SSC, file = "../pub/WES12/wes1_wes2_combined.gatk.rare1pct_variants_het_by_sample_matrix_cleaned.sp_sex_case_42651.txt", sep = "\t", quote = F, col.names = F, row.names = F)

table(is.na(sp_pass_QC_is_case))
