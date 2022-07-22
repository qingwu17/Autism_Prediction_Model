rm(list = ls())
workingDir <- "/Users/shining17wq/Google Drive/MorrowLab/SFARI/SPARK/job"
# workingDir <- "/gpfs/data/emorrow/silvio/Qing_Wu/SFARI/batch_jobs"
setwd(workingDir)
getwd()

library(data.table)
library(stringr)
library(dplyr)
library(plyr)
library(ggplot2)
library(biomaRt)
library(gprofiler2)
library(clusterProfiler)
options(stringsAsFactors = F, scipen = 999)



# input feature name ------------------------------------------------------

selected_features_tpot <- read.table(file = "../pub/WES12/wes1_wes2_combined.deepvariant.feature_name_exonic_selPct4.txt", header = T)
selected_gene_name <- selected_features_tpot$selected_features[-1]
grep("^GAB", selected_gene_name, value = T)

# SFARI and DDG2P genes ---------------------------------------------------

SFARI_gene_lst_df <- read.csv(file = "../reference_data/SFARI_gene/SFARI-Gene_genes_01-11-2022release_03-22-2022export.csv")
head(SFARI_gene_lst_df)
SFARI_genes <- SFARI_gene_lst_df$gene.symbol # 1031

# DDG2P gene
DDG2P_gene_lst_df <- read.csv(file = "../reference_data/DDG2P_genes/DDG2P_7_4_2022.csv")
head(DDG2P_gene_lst_df)
DDG2P_genes <- DDG2P_gene_lst_df$gene.symbol # 2580

intersect(SFARI_genes, DDG2P_genes)
grep("^CLDN", unique(c(SFARI_genes, DDG2P_genes)), value = T)

# venn diagram ------------------------------------------------------------


sort(intersect(selected_gene_name, SFARI_genes)) # 79
sort(intersect(selected_gene_name, DDG2P_genes)) # 157
sort(intersect(selected_gene_name, intersect(SFARI_genes, DDG2P_genes))) # 58

library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

# Chart
vennDiagram <- venn.diagram(
  x = list(selected_gene_name, SFARI_genes, DDG2P_genes),
  category.names = c("Selected Features", "SFARI Genes", "DDG2P Genes"),
  filename = NULL,
  
  # Output features
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = 1,
  fontface = "bold",
  # fontfamily = "sans",
  
  # Set names
  cat.cex = 1,
  cat.fontface = "bold",
  # cat.default.pos = "outer",
  cat.pos = c(-20, 20, 180),
  # cat.dist = c(0.055, 0.055, 0.055),
  # cat.fontfamily = "sans",
  rotation = 1
)

grid.draw(vennDiagram)
pdf(file = "venn.pdf")
grid.draw(vennDiagram)
dev.off()

# feature importance ------------------------------------------------------

feature_importance_df <- read.csv(file = "../_out/feature_importance.csv", header = T)
head(feature_importance_df)
str(feature_importance_df)

feature_importance_df$feature_name <- factor(feature_importance_df$feature_name, levels = feature_importance_df$feature_name)

feature_importance_df$pos <- (feature_importance_df$coef_ > 0)
feature_importance_df$pos <- ifelse(feature_importance_df$pos == TRUE, "Positive", "Negative")

feature_importance_df$color <- "New"
feature_importance_df$color[feature_importance_df$feature_name %in% SFARI_genes] <- "SFARI"
feature_importance_df$color[feature_importance_df$feature_name %in% DDG2P_genes] <- "DDG2P"
feature_importance_df$color[feature_importance_df$feature_name %in% intersect(SFARI_genes, DDG2P_genes)] <- "DDG2P & SFARI"

head(feature_importance_df)

input_df4ggplot <- feature_importance_df

ggplot(input_df4ggplot, aes(x=feature_name, y=coef_, fill=color)) +
  geom_bar(stat = "identity") +
  # theme_minimal() +
  theme_classic() +
  theme(axis.text.x = element_text(size = 6, angle = 45))



# gene pathway enrichment analysis ----------------------------------------

gp_link <- gost(list("combined_SFARI_genes" = unique(c(selected_gene_name, SFARI_genes)),
                     "SFARI_genes" = SFARI_genes), 
                multi_query = TRUE, significant = TRUE, 
                as_short_link = TRUE)
gp_link

# gene symbol to gene id
library(DOSE)

# example: https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
data(geneList)
de <- names(geneList)[abs(geneList) > 2]
de

library('org.Hs.eg.db')

symbols <- selected_gene_name
geneID <- mapIds(org.Hs.eg.db, keys = symbols, column = 'ENTREZID', keytype = 'SYMBOL')
edo <- enrichDO(geneID)
edo_df <- as.data.frame(edo@result)
edo_df <- as.data.frame(edo)
dotplot(enrichDO(geneID))

X_linked_monogenic_disease_geneID <- c("412", "8573", "64840", "7454", "4952", "2332", "215", "2157")
X_linked_monogenic_disease_gene_name <- mapIds(org.Hs.eg.db, keys = X_linked_monogenic_disease_geneID, column = 'SYMBOL', keytype = 'ENTREZID')

edo <- enrichDGN(geneID)
edox <- setReadable(edo, 'org.Hs.eg.db', "ENTREZID")
cnetplot(edox)
cnetplot(edox, categorySize = "pvalue", node_label="gene")

# enriched pathway in genes
combined_SFARI_gene_lst_gp <- gost(unique(c(selected_gene_name, SFARI_genes)),
                                   multi_query = FALSE, significant = TRUE, correction_method = "fdr", user_threshold = 0.01)
SFARI_gene_lst_gp <- gost(SFARI_genes, significant = TRUE, correction_method = "fdr", user_threshold = 0.01)

combined_DDG2P_gene_lst_gp <- gost(unique(c(selected_gene_name, DDG2P_genes)),
                                   multi_query = FALSE, significant = TRUE, correction_method = "fdr", user_threshold = 0.01)
DDG2P_gene_lst_gp <- gost(DDG2P_genes, significant = TRUE, correction_method = "fdr", user_threshold = 0.01)

selected_gene_lst_gp <- gost(selected_gene_name, significant = TRUE, correction_method = "fdr", user_threshold = 0.01)

# add combined SFARI genes
combined_SFARI_gene_lst_gp_df <- combined_SFARI_gene_lst_gp$result
SFARI_gene_lst_gp_df <- SFARI_gene_lst_gp$result
selected_gene_lst_gp_df <- selected_gene_lst_gp$result

idx_uniq_comp2SFARI <- !(combined_SFARI_gene_lst_gp_df$term_name %in% SFARI_gene_lst_gp_df$term_name)
table(idx_uniq_comp2SFARI)
dim(combined_SFARI_gene_lst_gp_df[idx_uniq_comp2SFARI,])

uniq_combined_SFARI_gene_df <- combined_SFARI_gene_lst_gp_df[idx_uniq_comp2SFARI,]
uniq_combined_SFARI_gene_df <- uniq_combined_SFARI_gene_df[order(uniq_combined_SFARI_gene_df$p_value, decreasing = T),]

enriched_GO <- rbind(tail(uniq_combined_SFARI_gene_df[uniq_combined_SFARI_gene_df$source == "GO:MF",], 10),
                     tail(uniq_combined_SFARI_gene_df[uniq_combined_SFARI_gene_df$source == "GO:CC",], 10),
                     tail(uniq_combined_SFARI_gene_df[uniq_combined_SFARI_gene_df$source == "GO:BP",], 10))


enriched_GO$log_p <- -log10(enriched_GO$p_value)
enriched_GO$term_name <- factor(enriched_GO$term_name, levels = unique(enriched_GO$term_name))
ggplot(enriched_GO, aes(x = term_name, y = log_p)) +
  geom_col() +
  coord_flip() +
  theme_classic()


selected_gene_lst_gp_df <- selected_gene_lst_gp_df[order(selected_gene_lst_gp_df$p_value, decreasing = T), ]
selected_gene_lst_gp_df$log_p <- -log10(selected_gene_lst_gp_df$p_value)
selected_gene_lst_gp_df$term_name <- factor(selected_gene_lst_gp_df$term_name, levels = unique(selected_gene_lst_gp_df$term_name))

ggplot(selected_gene_lst_gp_df[selected_gene_lst_gp_df$source == "GO:BP",], aes(x = term_name, y = log_p)) +
  geom_col() +
  coord_flip() +
  theme_classic()


ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
get_cytoband <- function(gene_lst){
  
  # get cytoband for input gene list
  input_genes_cytoband_df <- getBM(attributes=c('hgnc_symbol', 'chromosome_name', "band"), 
                                   filters=c('hgnc_symbol'), 
                                   values=list(gene_lst), 
                                   mart=ensembl)
  idx_rm <- grepl("^CHR_", input_genes_cytoband_df$chromosome_name)
  input_genes_cytoband_df <- input_genes_cytoband_df[!idx_rm,]
  input_genes_cytoband_df$cyto_band <- paste(input_genes_cytoband_df$chromosome_name, input_genes_cytoband_df$band, sep = "")
  
  out_df <- data.frame(gene_name = gene_lst, 
                       cytoband = input_genes_cytoband_df$cyto_band[match(gene_lst, input_genes_cytoband_df$hgnc_symbol)])
  return(out_df)
}




# feature importance to y -------------------------------------------------

negFeat_to_y_df <- read.csv(file = "/Users/shining17wq/Google Drive/MorrowLab/SFARI/SPARK/_out/wes1_wes2_combined.deepvariant.exonic.negFeat_to_y.csv", header = T)
negFeat_to_y_df
table(negFeat_to_y_df$is_female)
table(negFeat_to_y_df$GAGE1)
table(negFeat_to_y_df$CT45A1)
table(negFeat_to_y_df$CSAG3)
table(negFeat_to_y_df$TLK1)

table(negFeat_to_y_df$y == 1 & negFeat_to_y_df$is_female == 1)
table(negFeat_to_y_df$y == 1 & negFeat_to_y_df$GAGE1 == 1)
table(negFeat_to_y_df$y == 1 & negFeat_to_y_df$CT45A1 == 1)
table(negFeat_to_y_df$y == 1 & negFeat_to_y_df$CSAG3 == 1)
table(negFeat_to_y_df$y == 1 & negFeat_to_y_df$TLK1 == 1) 

posFeat_to_y_df <- read.csv(file = "/Users/shining17wq/Google Drive/MorrowLab/SFARI/SPARK/_out/wes1_wes2_combined.deepvariant.exonic.posFeat_to_y.csv", header = T)
posFeat_to_y_df
table(posFeat_to_y_df$WDR44)
table(posFeat_to_y_df$CLDN34)
table(posFeat_to_y_df$PNPLA4)
table(posFeat_to_y_df$RAB41)
table(posFeat_to_y_df$NOX1)

table(posFeat_to_y_df$y == 1 & posFeat_to_y_df$WDR44 == 1) # 82/97
table(posFeat_to_y_df$y == 1 & posFeat_to_y_df$CLDN34 == 1) # 41/47
table(posFeat_to_y_df$y == 1 & posFeat_to_y_df$PNPLA4 == 1) # 29/32
table(posFeat_to_y_df$y == 1 & posFeat_to_y_df$RAB41 == 1) # 70/85
table(posFeat_to_y_df$y == 1 & posFeat_to_y_df$NOX1 == 1) # 125/151


# probability to sex and y ------------------------------------------------

prob_sex_on_test_selected_genes_df <- read.csv("/Users/shining17wq/Downloads/wes1_wes2_combined.deepvariant.exonic.linearSVC.prob_sex_on_test_selected_genes.csv", header = T)
names(prob_sex_on_test_selected_genes_df) <- c("score", "sex", "ASD")

prob_sex_on_test_selected_genes_df$ASD <- factor(prob_sex_on_test_selected_genes_df$ASD)
prob_sex_on_test_selected_genes_df$sex <- factor(prob_sex_on_test_selected_genes_df$sex, levels = c(1,0))
head(prob_sex_on_test_selected_genes_df)

prob_sex_on_test_selected_genes_df$ASD_sex <- ""
idx_ASD_male <- prob_sex_on_test_selected_genes_df$sex == 0 & prob_sex_on_test_selected_genes_df$ASD == 1
idx_ASD_female <- prob_sex_on_test_selected_genes_df$sex == 1 & prob_sex_on_test_selected_genes_df$ASD == 1
idx_contr_male <- prob_sex_on_test_selected_genes_df$sex == 0 & prob_sex_on_test_selected_genes_df$ASD == 0
idx_contr_female <- prob_sex_on_test_selected_genes_df$sex == 1 & prob_sex_on_test_selected_genes_df$ASD == 0
prob_sex_on_test_selected_genes_df$ASD_sex[idx_ASD_male] <- "ASD_male"
prob_sex_on_test_selected_genes_df$ASD_sex[idx_ASD_female] <- "ASD_female"
prob_sex_on_test_selected_genes_df$ASD_sex[idx_contr_male] <- "control_male"
prob_sex_on_test_selected_genes_df$ASD_sex[idx_contr_female] <- "control_female"

mean(prob_sex_on_test_selected_genes_df$score[prob_sex_on_test_selected_genes_df$ASD_sex == "ASD_male"])
mean(prob_sex_on_test_selected_genes_df$score[prob_sex_on_test_selected_genes_df$ASD_sex == "ASD_female"])
mean(prob_sex_on_test_selected_genes_df$score[prob_sex_on_test_selected_genes_df$ASD_sex == "control_male"])
mean(prob_sex_on_test_selected_genes_df$score[prob_sex_on_test_selected_genes_df$ASD_sex == "control_female"])

head(prob_sex_on_test_selected_genes_df)
stat_df <- ddply(prob_sex_on_test_selected_genes_df, "sex", summarise, grp.mean=mean(score))
stat_df

ggplot(prob_sex_on_test_selected_genes_df, aes(x = score, fill = sex, color = sex)) +
  geom_histogram(alpha=0.5, position = "dodge", binwidth = 0.01) + 
  geom_vline(data = stat_df, aes(xintercept=grp.mean, color = sex), linetype="dashed") +
  labs(x = "score") +
  scale_color_brewer(palette="Set1")+
  theme_classic() +
  theme(legend.position = c(0.9, 0.8)) 

t.test(prob_sex_on_test_selected_genes_df$score[idx_contr_male],
       prob_sex_on_test_selected_genes_df$score[idx_ASD_male], paired = F, alternative = "less", var.equal = F)

table((prob_sex_on_test_selected_genes_df$score > 0) & idx_ASD_female)
table(idx_ASD_female)
529/1264

mean(prob_sex_on_test_selected_genes_df$score) - 2 * sd(prob_sex_on_test_selected_genes_df$score)

42487-16721
23409-11101
