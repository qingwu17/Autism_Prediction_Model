#!/bin/bash
#SBATCH -J edit_VxSp_file
#SBATCH -n 8
#SBATCH -t 72:00:00
#SBATCH --mem=256G
#SBATCH --mail-type=END
#SBATCH --mail-user=qing_wu@brown.edu

module load bcftools/1.13 samtools/1.13 tabix/0.2.6 gsl/2.5 gcc/8.3 perl/5.24.1

# gunzip the bgz file
gunzip -c /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_0216.vcf.bgz \
> /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_0216.vcf

# remove the vcf header lines
bcftools view -H /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_0216.vcf \
-Ov -o /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_sample_table_0216.txt

# substitute ./X with 0/X
sed -i -e 's/\.\//0\//g' /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_sample_table_0216.txt
# substitute X/. with X/0
sed -i -e 's/\/\./\/0/g' /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_sample_table_0216.txt

# substitute | with /
sed -i -e 's/|/\//g' /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_sample_table_0216.txt

# substitute 0/0 with 0
echo "substitute 0/0 with 0"
sed -i -e 's/0\/0/0/g' /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_sample_table_0216.txt
# other unchanged variants format: 0/1 1/1

# substitute 0/1 with 1
echo "substitute 0/1 with 1"
sed -i -e 's/0\/1/1/g' /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_sample_table_0216.txt
# other unchanged variants format: 1/1

# substitute 1/1 with 2
echo "substitute 1/1 with 2"
sed -i -e 's/1\/1/2/g' /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_sample_table_0216.txt
# no variants should have format other than: 0, 1, 2

# split the raw variant by sample table into 1) variant info and 2) variant by sample data only
cut -f 1-9 /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_sample_table_0216.txt \
> /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_info_0216.txt
cut -f 10- /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_sample_table_0216.txt \
> /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_sample_matrix_0216.txt

# remove variants by index (including exonic and splicing variants, exclude synonymous and the rest)
sed "$(sed 's/$/d/' /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_idx_var2exclude_0216.txt)" \
/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_sample_matrix_0216.txt \
> /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_sample_matrix_cleaned_0216.txt

# add "d" to the end of each character
sed 's/$/d/' /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_idx_var2exclude_0216.txt \
> /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_idx_var2exclude_0216_2.sed
# substitute /n by ;
sed -z -i 's/\n/;/g' /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_idx_var2exclude_0216_2.sed
sed -i "s/;$//g" /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_idx_var2exclude_0216_2.sed
# filter 
sed -f /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_idx_var2exclude_0216_2.sed \
/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_sample_matrix_0216.txt \
> /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_sample_matrix_cleaned_0216_2.txt


# remove variants by index (stopgain)
# add "d" to the end of each character
sed 's/$/d/' /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_idx_var2exclude_stopgain.txt \
> /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_idx_var2exclude_stopgain.sed
# substitute /n by ;
sed -z -i 's/\n/;/g' /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_idx_var2exclude_stopgain.sed
sed -i "s/;$//g" /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_idx_var2exclude_stopgain.sed
# filter 
sed -f /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_idx_var2exclude_stopgain.sed \
/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_sample_matrix_cleaned.txt \
> /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_sample_matrix_stopgain.txt

# remove variants by index (stoploss)
sed "$(sed 's/$/d/' /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_idx_var2exclude_stoploss.txt)" \
/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_sample_matrix_cleaned.txt \
> /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_sample_matrix_stoploss.txt

# remove variants by index (fs_sub)
# add "d" to the end of each character
sed 's/$/d/' /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_idx_var2exclude_fs_sub.txt \
> /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_idx_var2exclude_fs_sub.sed
# substitute /n by ;
sed -zi 's/\n/;/g' /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_idx_var2exclude_fs_sub.sed
sed -i "s/;$//g" /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_idx_var2exclude_fs_sub.sed
# filter 
sed -f /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_idx_var2exclude_fs_sub.sed \
/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_sample_matrix_cleaned.txt \
> /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_sample_matrix_fs_sub.txt

# remove variants by index (nfs_sub)
sed "$(sed 's/$/d/' /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_idx_var2exclude_nfs_sub.txt)" \
/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_sample_matrix_cleaned.txt \
> /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_sample_matrix_nfs_sub.txt

# remove variants by index (splicing)
# add "d" to the end of each character
sed 's/$/d/' /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_idx_var2exclude_splicing.txt \
> /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_idx_var2exclude_splicing.sed
# substitute /n by ;
sed -z -i 's/\n/;/g' /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_idx_var2exclude_splicing.sed
sed -i "s/;$//g" /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_idx_var2exclude_splicing.sed
# filter 
sed -f /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_idx_var2exclude_splicing.sed \
/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_sample_matrix_cleaned.txt \
> /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_sample_matrix_splicing.txt

# remove variants by index (missense MPC>1)
sed 's/$/d/' /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_idx_var2exclude_MPC1.txt \
> /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_idx_var2exclude_MPC1.sed
# substitute /n by ;
sed -z -i 's/\n/;/g' /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_idx_var2exclude_MPC1.sed
sed -i "s/;$//g" /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_idx_var2exclude_MPC1.sed
# filter 
sed -f /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_idx_var2exclude_MPC1.sed \
/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_sample_matrix_cleaned.txt \
> /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_sample_matrix_MPC1.txt


# remove variants by index (stopgain_fs_sub_splicing)
# add "d" to the end of each character
sed 's/$/d/' /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_idx_var2exclude_stopgain_fs_sub_splicing.txt \
> /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_idx_var2exclude_stopgain_fs_sub_splicing.sed
# substitute /n by ;
sed -z -i 's/\n/;/g' /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_idx_var2exclude_stopgain_fs_sub_splicing.sed
sed -i "s/;$//g" /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_idx_var2exclude_stopgain_fs_sub_splicing.sed
# filter 
sed -f /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_idx_var2exclude_stopgain_fs_sub_splicing.sed \
/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_sample_matrix_cleaned.txt \
> /users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_sample_matrix_stopgain_fs_sub_splicing.txt





