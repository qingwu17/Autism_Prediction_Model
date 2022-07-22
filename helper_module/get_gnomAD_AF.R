rm(list = ls())
workingDir <- "/Users/shining17wq/Google Drive/MorrowLab/SFARI/SPARK/job"
setwd(workingDir)

library(stringr)
library(dplyr)
library(plyr)
options(stringsAsFactors = F)

# allele frequency --------------------------------------------------------

gnomAD_exome211_AF_df_0 <- read.table(file = "/Users/shining17wq/Google Drive/MorrowLab/gnomAD/gnomAD_exome211_AF_non_neuro.tsv", header = T)
head(gnomAD_exome211_AF_df_0)

# add a Variant column
alleles_0 <- str_replace_all(gnomAD_exome211_AF_df_0$alleles, "(\\[|\\]|\")", "") 
# replace the first , by :
gnomAD_exome211_AF_df_0$alleles <- str_replace(alleles_0, ",", ":")
gnomAD_exome211_AF_df_0$Variant <- paste(gnomAD_exome211_AF_df_0$locus, gnomAD_exome211_AF_df_0$alleles, sep = ":")
head(gnomAD_exome211_AF_df_0)

# select Variants and gnomAD_exome_flag column
gnomAD_exome211_AF_df <- gnomAD_exome211_AF_df[,c("Variant", "gnomad_freq")]
head(gnomAD_exome211_AF_df)

write.table(gnomAD_exome211_AF_df, file = "/Users/shining17wq/Google Drive/MorrowLab/gnomAD/gnomAD_exome211_AF_non_neuro_na.tsv", sep = "\t", quote = F, col.names = T, row.names = F)

idx_nona <- !is.na(gnomAD_exome211_AF_df_edit$gnomad_freq)
table(!idx_nona)
gnomAD_exome211_AF_df_nona <- gnomAD_exome211_AF_df_edit[idx_nona,]
head(gnomAD_exome211_AF_df_nona)

write.table(gnomAD_exome211_AF_df_nona, file = "/Users/shining17wq/Google Drive/MorrowLab/gnomAD/gnomAD_exome211_AF_non_neuro_nona.tsv", sep = "\t", quote = F, col.names = T, row.names = F)

gnomAD_exome211_AF_df_na0 <- gnomAD_exome211_AF_df
gnomAD_exome211_AF_df_na0[is.na(gnomAD_exome211_AF_df_na0)] <- 0
head(gnomAD_exome211_AF_df_na0)

write.table(gnomAD_exome211_AF_df_na0, file = "/Users/shining17wq/Google Drive/MorrowLab/gnomAD/gnomAD_exome211_AF_non_neuro_na0.tsv", sep = "\t", quote = F, col.names = T, row.names = F)

