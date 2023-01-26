rm(list = ls())
workingDir <- "/path_to_working_directory/"
setwd(workingDir)

library(stringr)
library(dplyr)
library(plyr)
options(stringsAsFactors = F)


# input -------------------------------------------------------------------
SPARK_CSQ <- read.table(file = "/path_to/wes.deepvariant.individual_1to5.CSQ.tsv", sep = '\t', header = TRUE)
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

write.table(SPARK_CSQ, file = "/path_to/wes.deepvariant.segdup.tsv", row.names = F, col.names = T, quote = F, sep = "\t")


