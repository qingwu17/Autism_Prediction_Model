#!/bin/bash
#SBATCH -J table_annovar_wes1_wes2_combined.gatk
#SBATCH -n 16
#SBATCH -t 36:00:00
#SBATCH --mem=64G
#SBATCH --mail-type=END
#SBATCH --mail-user=qing_wu@brown.edu
 
module load annovar/2018Apr16

# annotate variants
table_annovar.pl /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.gatk.rare1pct_variants_het_info.avinput \
/users/qwu24/data/silvio/Qing_Wu/SFARI/annovar/humandb/ -buildver hg38 \
-out /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.gatk.rare1pct_variants_het_info \
-remove -protocol refGene,knownGene,ensGene,dbnsfp42a,dbscsnv11 \
-operation gx,g,g,f,f \
-nastring . -csvout -polish -xref /users/qwu24/data/silvio/Qing_Wu/SFARI/annovar/humandb/gene_xref_pLi.txt
