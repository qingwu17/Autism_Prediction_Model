#!/bin/bash
#SBATCH -J SPARK_wes_vep-dbnsfp-segdup
#SBATCH -n 32
#SBATCH -t 36:00:00
#SBATCH --mem=128G

module load vep/104 perl/5.24.1 mysql/8.0.13

# 
vep --fork 8 --species homo_sapiens --assembly GRCh38 --symbol --gencode_basic --distance 1 \
--cache --offline --format vcf --vcf --force_overwrite --fields "Allele,Consequence,Feature_type,Feature,gnomAD_exomes_flag" \
--dir_cache /gpfs/runtime/opt/vep/104/cache \
--plugin dbNSFP,/path_to/dbNSFP42a/dbNSFP4.2a_hg38.gz,gnomAD_exomes_flag \
--input_file /path_to/wes1_27281_exome.gatk.individual_1to5.vcf.bgz \
--output_file /path_to/wes1_27281_exome.gatk.individual_1to5.vep_anno_segdup.vcf