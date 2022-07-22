rm(list = ls())
workingDir <- "/Users/shining17wq/Google Drive/MorrowLab/SFARI/SPARK/job"
setwd(workingDir)

library(stringr)
library(dplyr)
library(plyr)
options(stringsAsFactors = F)


# input -------------------------------------------------------------------
SPARK_CSQ <- read.table(file = "../Regeneron/wes2_15995_exome.gatk.individual_1to5.CSQ.tsv", sep = '\t', header = TRUE)
head(SPARK_CSQ)

CSQ_0 <- SPARK_CSQ$X__uid_4
CSQ_list <- str_split(CSQ_0, "\\|", 5)
CSQ_list

SPARK_segdup <- unlist(lapply(CSQ_list, tail, n = 1L))
table(SPARK_segdup)

idx_na <- is.na(SPARK_segdup)
SPARK_segdup[idx_na] <- ""
table(SPARK_segdup)

# add a Variant column
alleles_0 <- str_replace_all(SPARK_CSQ$alleles, "(\\[|\\])", "") 
# replace the first , by :
alleles_sub <- sub(":", "|", alleles_0)
SPARK_CSQ$alleles <- str_replace(alleles_0, ",", ":")
SPARK_CSQ$Variant <- paste(SPARK_CSQ$locus, SPARK_CSQ$alleles, sep = ":")

# add a gnomAD exome flag column
SPARK_CSQ$gnomAD_exome_flag <- SPARK_segdup
head(SPARK_CSQ)

# select Variants and gnomAD_exome_flag column
SPARK_CSQ <- SPARK_CSQ[,c("Variant", "gnomAD_exome_flag")]
head(SPARK_CSQ)

write.table(SPARK_CSQ, file = "../Regeneron/wes2_15995_exome.gatk.segdup.tsv", row.names = F, col.names = T, quote = F, sep = "\t")

# sample id ---------------------------------------------------------------

sample_id <- read.table(file = "../Regeneron/SPARK_vcf_families_merged_mall_sample_id.tsv", sep = '\t', header = TRUE)

table(nchar(sample_id$s))

idx <- (nchar(sample_id$s) > 9)
sample_id$s[idx]

sample_id$sample_id_wdup <- str_replace(sample_id$s, "^\\d+:", "")
sample_id$sample_id_wdup_bool <- str_to_title(idx)
head(sample_id)

write.table(sample_id, file = "../Regeneron/SPARK_vcf_families_merged_mall_sample_id_wdup.tsv", row.names = F, col.names = T, quote = F, sep = "\t")


