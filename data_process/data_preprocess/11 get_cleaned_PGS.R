rm(list = ls())
workingDir <- "/path_to_working_directory/"
setwd(workingDir)
getwd()

library(data.table)
library(stringr)
library(dplyr)
library(plyr)
options(stringsAsFactors = F, scipen = 999)



dir_name <- "/path_to/"
sample_set_name <- "wes1_wes2_combined" # other option

input_filename <- paste(dir_name, sample_set_name, ".deepvariant.rare1pct_variants_by_sample_matrix_cleaned.sp_sex_case.txt", sep = "") # wes12: 43203, iwes1: 70464

header_df <- read.table(file = input_filename)
# head(header_df)
# dim(header_df)

PGS_out <- read.table(file = "../pub/iWES_v1/GWAS_PGS_out/SPARK.iWES_v1.array.2022_02.cleaned.pgs.profile", header = T)
head(PGS_out)

idx_PGS2spid_header <- match(header_df$V1, PGS_out$IID)



spid2excl_noPGS <- c(1:nrow(header_df))[is.na(idx_PGS2spid_header)]
out_filename_spid2excl_noPGS <- paste(dir_name, sample_set_name, ".deepvariant.rare1pct_variants_by_sample_matrix_cleaned.spid2excl_noPGS.txt", sep = "")

write.table(spid2excl_noPGS, file = out_filename_spid2excl_noPGS, sep = "\t", quote = F, col.names = F, row.names = F)


write.table(spid2PGS_df, file = "../pub/iWES_v1/iWES_v1.deepvariant.rare1pct_variants_het_by_sample_matrix_cleaned.spid2PGS_68301.txt",  sep = "\t", quote = F, col.names = F, row.names = F) 


# if subset, iWES_v1 with PGS but not in WES12
WES12_sample_id_pass_QC <- read.table(file = "../pub/WES12/wes1_wes2_combined.deepvariant.sample_id.tsv", header = T)
spid_WES12 <- WES12_sample_id_pass_QC$s

idx_sp2exclude_WES12 <- c(1:nrow(spid2PGS_df))[spid2PGS_df$IID %in% spid_WES12]
length(idx_sp2exclude_WES12) # 40871

write.table(idx_sp2exclude_WES12, file = "../pub/iWES_v1/iWES_v1.deepvariant.rare1pct_variants_het_by_sample_matrix_cleaned.spid2subset_wPGS_27430.txt", sep = "\t", quote = F, col.names = F, row.names = F)


# sample not in WES12 -----------------------------------------------------

sample_id_pass_QC <- read.table(file = "/path_to/wes_70487_exome.deepvariant.sample_id.tsv", header = T)

WES12_sample_id_pass_QC <- read.table(file = "/path_to/wes1_wes2_combined.deepvariant.sample_id.tsv", header = T)
head(WES12_sample_id_pass_QC)

table(sample_id_pass_QC$s %in% WES12_sample_id_pass_QC$s)

# pick p1 and s1
idx_WES12 <- sample_id_pass_QC$s %in% WES12_sample_id_pass_QC$s
idx_sp2exclude <- c(1:nrow(sample_id_pass_QC))[idx_WES12]
write.table(idx_sp2exclude, file = "/path_to/wes_70487_exome.deepvariant.rare1pct_variants_het_by_sample_matrix_idx_sp2exclude.txt", sep = "\t", quote = F, col.names = F, row.names = F)



