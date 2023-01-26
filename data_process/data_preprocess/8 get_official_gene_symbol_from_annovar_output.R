rm(list = ls())
workingDir <- "/path_to_working_directory/"
setwd(workingDir)
getwd()

library(data.table)
library(stringr)
library(dplyr)
library(plyr)
options(stringsAsFactors = F, scipen = 999)



# HUGO genes --------------------------------------------------------------

var_gene_anno_df <- fread(file = "/path_to/wes1_wes2_combined.deepvariant.rare1pct_variants_info.hg38_multianno.csv")
head(var_gene_anno_df)

length(unique(var_gene_anno_df$Gene.refGene)) 
length(unique(var_gene_anno_df$Gene.knownGene)) 
length(unique(var_gene_anno_df$Gene.ensGene)) 

grep("^\\.$", var_gene_anno_df$Gene.ensGene, value = T)
grep("^$", var_gene_anno_df$Gene.ensGene, value = T)
# no empty gene_name

gene_name_df <- data.frame(refGene = var_gene_anno_df$Gene.refGene, 
                           knownGene = var_gene_anno_df$Gene.knownGene,
                           ensGene = var_gene_anno_df$Gene.ensGene)

gene_name_df$match_refGene2knownGene <- gene_name_df$refGene == gene_name_df$knownGene
gene_name_df$match_refGene2ensGene <- gene_name_df$refGene == gene_name_df$ensGene
gene_name_df$match_knownGene2ensGene <- gene_name_df$knownGene == gene_name_df$ensGene

gene_name_df$gene_name <- gene_name_df$ensGene
# ref=known=ens, do not need to replace gene name
# ref=known=!ens, replace ensembl gene name
idx_identical_ref_known <- gene_name_df$match_refGene2knownGene & !gene_name_df$match_knownGene2ensGene
gene_name_df$gene_name[idx_identical_ref_known] <- gene_name_df$refGene[idx_identical_ref_known]
# ref=ens!=known, do not need to replace gene name
# ref!=known=ens, do not need to replace gene name
# ref!=known!=ens, keep ensembl gene name

# edit multiple gene name issue
idx_multi_gene_in_refGene <- grepl(";", gene_name_df$refGene)
idx_multi_gene_in_knownGene <- grepl(";", gene_name_df$knownGene)
table(grepl(";", gene_name_df$gene_name))
# multi-gene found in gene_name, but single name in refGene, replace gene_name by refGene
idx_single_refGene <- grepl(";", gene_name_df$gene_name) & !idx_multi_gene_in_refGene
gene_name_df$gene_name[idx_single_refGene] <- gene_name_df$refGene[idx_single_refGene]
# multi-gene found in gene_name, but single name in knownGene, replace gene_name by knownGene
idx_single_knownGene <- grepl(";", gene_name_df$gene_name) & !idx_multi_gene_in_knownGene
gene_name_df$gene_name[idx_single_knownGene] <- gene_name_df$refGene[idx_single_knownGene]

# check how many gene has multiple names
idx_genes_name <- grep(";", gene_name_df$gene_name) 
intra_genes_name <- grep(";", gene_name_df$gene_name, value = T)
# split multiple names by ;
intra_genes_lst <- str_split(intra_genes_name, ";")
intra_genes_df <- plyr::ldply(intra_genes_lst, rbind)
num_unique_gene_name <- apply(intra_genes_df, 1, function(x) length(unique(x)))
idx_identical_genes <- (num_unique_gene_name - 1 == 1)
table(idx_identical_genes) # 538

replaced_gene_name <- intra_genes_df$`1`[idx_identical_genes]
gene_name_df$gene_name[idx_genes_name[idx_identical_genes]] <- replaced_gene_name

length(unique(gene_name_df$gene_name)) # 25060

grep("FOXP", gene_name_df$gene_name)
table(gene_name_df$refGene[grep("FOXP", gene_name_df$gene_name)])


unique_gene_name <- unique(gene_name_df$gene_name)
# unique_gene_name

variants2gene_df <- data.frame(Chr = var_gene_anno_df$Chr,
                               Start = var_gene_anno_df$Start,
                               End = var_gene_anno_df$End,
                               Ref = var_gene_anno_df$Ref,
                               Alt = var_gene_anno_df$Alt,
                               Gene = gene_name_df$gene_name)
head(variants2gene_df)
write.csv(variants2gene_df, file = "/path_to/wes1_wes2_combined.deepvariant.rare1pct_variants_info.gene_converted.csv")






