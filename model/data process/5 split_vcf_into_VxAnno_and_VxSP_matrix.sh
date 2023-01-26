#!/bin/bash
#SBATCH -J split_vcf_into_VxAnno_and_VxSP_matrix
#SBATCH -n 8
#SBATCH -t 96:00:00
#SBATCH --mem=128G


module load bcftools/1.13 samtools/1.13 tabix/0.2.6 gsl/2.5 gcc/8.3 perl/5.24.1


# gunzip the bgz file
gunzip -c /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.vcf.bgz \
> /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.vcf

# remove the vcf header lines
bcftools view -H /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.vcf \
-Ov -o /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants_by_sample_table.txt


# substitute ./X with 0/X
echo "substitute ./X with 0/X"
sed -i -e 's/\.\//0\//g' /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants_by_sample_table.txt
# substitute X/. with X/0
echo "substitute X/. with X/0"
sed -i -e 's/\/\./\/0/g' /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants_by_sample_table.txt
# substitute | with /
echo "substitute | with /"
sed -i -e 's/|/\//g' /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants_by_sample_table.txt


# substitute 0/0 with 0
echo "substitute 0/0 with 0"
sed -i -e 's/0\/0/0/g' /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants_by_sample_table.txt
# other unchanged variants format: 0/1 1/1

# substitute 0/1 with 1
echo "substitute 0/1 with 1"
sed -i -e 's/0\/1/1/g' /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants_by_sample_table.txt
# other unchanged variants format: 1/1

# substitute 1/1 with 2
echo "substitute 1/1 with 2"
sed -i -e 's/1\/1/2/g' /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants_by_sample_table.txt
# no variants should have format other than: 0, 1, 2


# split the raw variant by sample table into 1) variant info table and 2) variant by sample data frame 
cut -f 1-9 /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants_by_sample_table.txt \
> /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants_info.txt
cut -f 10- /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants_by_sample_table.txt \
> /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants_by_sample_matrix.txt

