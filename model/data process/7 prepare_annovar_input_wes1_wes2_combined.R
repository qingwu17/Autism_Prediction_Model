rm(list = ls())
workingDir <- "/path_to_working_directory/"
setwd(workingDir)
getwd()

library(data.table)
library(stringr)
library(dplyr)
library(plyr)
options(stringsAsFactors = F, scipen = 999)


# generate annovar input --------------------------------------------------

var_df <- read.table(file = "/path_to/wes1_wes2_combined.deepvariant.rare1pct_variants_info.txt", sep = '\t', header = F)
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

write.table(avinput_df, file = "/path_to/wes1_wes2_combined.deepvariant.rare1pct_variants_info.avinput", sep = "\t", quote = F, col.names = F, row.names = F)


