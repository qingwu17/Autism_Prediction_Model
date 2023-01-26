rm(list = ls())
workingDir <- "/path_to_working_directory/"
setwd(workingDir)
getwd()

library(data.table)
library(stringr)
library(dplyr)
library(plyr)
options(stringsAsFactors = F, scipen = 999)

# phenotype sp annotation -------------------------------------------------

dir_name <- "/path_to/"
sample_set_name <- "wes1_wes2_combined" # other option

sample_id_pass_QC <- read.table(file = paste(dir_name, sample_set_name, ".deepvariant.sample_id.tsv", sep = ""), header = T)
sp_pass_QC_is_female <- read.table(file = paste(dir_name, sample_set_name, ".deepvariant.is_female.tsv", sep = ""), header = T)
sp_pass_QC_is_case <- read.table(file = paste(dir_name, sample_set_name, ".deepvariant.is_case.tsv", sep = ""), header = T)
head(sample_id_pass_QC)
head(sp_pass_QC_is_case)
head(sp_pass_QC_is_female)

sp_pass_QC_is_case$is_case <- as.integer(as.logical(sp_pass_QC_is_case$X__uid_18))
sp_pass_QC_is_female$is_female <- as.integer(as.logical(sp_pass_QC_is_female$X__uid_19))


identical(sample_id_pass_QC$s, sp_pass_QC_is_case$s)
identical(sample_id_pass_QC$s, sp_pass_QC_is_female$s)
identical(sp_pass_QC_is_case$s, sp_pass_QC_is_female$s)

variant_sample_matrix_sp_id_SSC <- data.frame(sample_id = sample_id_pass_QC$s, is_female = sp_pass_QC_is_female$is_female, is_case = sp_pass_QC_is_case$is_case)
out_filename <- paste(dir_name, sample_set_name, ".deepvariant.rare1pct_variants_by_sample_matrix_cleaned.sp_sex_case.txt", sep = "")

write.table(variant_sample_matrix_sp_id_SSC, file = out_filename, sep = "\t", quote = F, col.names = F, row.names = F)

