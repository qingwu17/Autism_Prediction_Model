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
sample_id_pass_QC <- read.table(file = "../pub/iWES_v1/wes_70487_exome.deepvariant.sample_id.tsv", header = T)
sp_pass_QC_is_case <- read.table(file = "../pub/iWES_v1/wes_70487_exome.deepvariant.is_female.tsv", header = T)
sp_pass_QC_is_female <- read.table(file = "../pub/iWES_v1/wes_70487_exome.deepvariant.is_case.tsv", header = T)
head(sample_id_pass_QC)
head(sp_pass_QC_is_case)
head(sp_pass_QC_is_female)

sp_pass_QC_is_case$is_case <- as.integer(as.logical(sp_pass_QC_is_case$X__uid_59))
sp_pass_QC_is_female$is_female <- as.integer(as.logical(sp_pass_QC_is_female$X__uid_60))

head(sp_pass_QC_is_case)
head(sp_pass_QC_is_female)

identical(sample_id_pass_QC$s, sp_pass_QC_is_case$s)
identical(sample_id_pass_QC$s, sp_pass_QC_is_female$s)
identical(sp_pass_QC_is_case$s, sp_pass_QC_is_female$s)

variant_sample_matrix_sp_id_SSC <- data.frame(sample_id = sample_id_pass_QC$s, is_case = sp_pass_QC_is_case$is_case, is_female = sp_pass_QC_is_female$is_female)
head(variant_sample_matrix_sp_id_SSC)
dim(variant_sample_matrix_sp_id_SSC)

# write.table(variant_sample_matrix_sp_id_SSC, file = "../pub/iWES_v1/wes_70487_exome.deepvariant.rare1pct_variants_het_by_sample_matrix_cleaned.sp_sex_case_70464.txt", sep = "\t", quote = F, col.names = F, row.names = F)

table(is.na(sp_pass_QC_is_case))


# sample not in WES12 -----------------------------------------------------

WES12_sample_id_pass_QC <- read.table(file = "../pub/WES12/wes1_wes2_combined.deepvariant.sample_id.tsv", header = T)

table(sample_id_pass_QC$s %in% WES12_sample_id_pass_QC$s)

# pick p1 and s1
idx_WES12 <- sample_id_pass_QC$s %in% WES12_sample_id_pass_QC$s
idx_sp2exclude <- c(1:nrow(sample_id_pass_QC))[idx_WES12]
write.table(idx_sp2exclude, file = "../pub/iWES_v1/wes_70487_exome.deepvariant.rare1pct_variants_het_by_sample_matrix_idx_sp2exclude.txt", sep = "\t", quote = F, col.names = F, row.names = F)



