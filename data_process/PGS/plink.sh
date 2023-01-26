#!/bin/bash
#SBATCH -J plink
#SBATCH -n 16
#SBATCH -t 36:00:00
#SBATCH --mem=64G

module load plink/1.90

plink -bfile /path_to/array/genotype/SPARK.iWES_v1.array.2022_02 \
--maf 0.01 \
--mind 0.1 \
--geno 0.1 \
--hwe 1e-6 \
--indep-pairwise 200 50 0.2 \
--make-bed --out /path_to/GWAS/SPARK.iWES_v1.array.2022_02.cleaned

# --mind excludes individuals with too much missing genotype data
# --maf excludes SNPs with minor allele frequency below 1%
# --geno excludes SNPs with > 10% genotyping data missed
# --hwe excludes SNPs which have Hardy-Weinberg equilibrium exact test p.value below 1e-6.

plink -bfile /path_to/GWAS/SPARK.iWES_v1.array.2022_02.cleaned \
--score /path_to/iPSYCH_PGC_comm_var/daner_iPSYCH-DBS_Aut_PGC1_Eur_ph3.out.snpRes 2 5 8 header sum center \
--out /path_to/GWAS/SPARK.iWES_v1.array.2022_02.cleaned.pgs

# Before QC, 69592 people (39578 males, 30014 females) loaded from .fam
# setting --mind 0.01, 15516 people pass filters and QC
# setting --mind 0.05, 25801 people pass filters and QC
# setting --mind 0.1, 68316 people pass filters and QC

# output: /path_to/GWAS/SPARK.iWES_v1.array.2022_02.cleaned.pgs.profile
