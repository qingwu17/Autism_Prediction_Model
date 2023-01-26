#!/usr/bin/env python
# -*- coding: utf-8 -*-

# PYSPARK_SUBMIT_ARGS="--driver-memory 200g --executor-memory 125g --executor-cores 2 --num-executors 8 pyspark-shell" ipython
# PYSPARK_SUBMIT_ARGS="--driver-memory 200g --executor-memory 350g pyspark-shell" ipython

import hail as hl
import os

hl.init(spark_conf=dict({'spark.driver.memory': '200g',
                         'spark.executor.memory': '350g'}))

# region import wes1 and wes2 matrix table and merge them
wes1_mt = hl.read_matrix_table('/path_to/wes1.deepvariant.pass_QC.mt')
print("Shape of wes1: ", wes1_mt.count())
# (4276376, 27266)

wes2_mt = hl.read_matrix_table('/path_to/wes2.deepvariant.pass_QC.mt')
print("Shape of wes2: ", wes2_mt.count())
# (3875705, 15986)

combined_mt = wes1_mt.union_cols(wes2_mt, row_join_type='outer')
# print("Shape of combined mt: ", combined_mt.count())
# # (6281625, 43252)

# remove duplicates
combined_mt = combined_mt.distinct_by_col()
# print("After remove duplicated samples across wes1 and wes2, sample size: ", combined_mt.count_cols())
# # After remove duplicated samples across wes1 and wes2, sample size: 43203

# endregion



# region annotate variants by gnomAD non-neuro population allele frequency

# download
if not os.path.isfile("/path_to/gnomad.exomes.r2.1.1.sites.liftover_grch38.ht"):
    gsutil -m cp -r gs://gcp-public-data--gnomad/release/2.1.1/liftover_grch38/ht/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.ht .
else:
    print("gnomAD allele frequency table has been downloaded. Continue ...")

# import gnomAD hail table
gnomAD_exome211_ht = hl.read_table('/path_to/gnomad.exomes.r2.1.1.sites.liftover_grch38.ht')
# gnomAD_exome211_ht.count() # 17201296

# check which population represents non-neuronal population in gnomAD 2.1.1 ==> index:102
gnomAD_exome211_freq_index_dict_txt_filename = "/path_to/gnomAD_exome211_hg38_ht_freq_index_dict.txt"

if 'non_neuro' in hl.eval(gnomAD_exome211_ht.freq_index_dict):
    print("Has non_neuro population.")

    if os.path.isfile(gnomAD_exome211_freq_index_dict_txt_filename):
        print(f"Check the index of non_neuro population in {gnomAD_exome211_freq_index_dict_txt_filename}")

    else:
        print(f"Export {gnomAD_exome211_freq_index_dict_txt_filename}.")
        gnomAD_exome211_ht.freq_index_dict.export('/path_to/gnomAD_exome211_hg38_ht_freq_index_dict.txt')

# annotate variants by gnomAD non-neuro population AF
combined_mt = combined_mt.annotate_rows(gnomAD_non_neuro_pop_AF=gnomAD_exome211_ht[combined_mt.row_key].freq[102].AF)
# handle missing data from gnomAD
combined_mt = combined_mt.annotate_rows(non_neuro_pop_AF=hl.coalesce(combined_mt.gnomAD_non_neuro_pop_AF, hl.float(0)))

# endregion



# region generate rare variants and export

rare_variants_maf = 0.01

# region heterozygous or hemizygous
combined_mt_het = combined_mt.filter_entries((combined_mt.GT.is_het()))

fraction_filtered_ab = combined_mt.aggregate_entries(hl.agg.fraction(~(combined_mt.GT.is_het())))
# print(f'Filtering {fraction_filtered_ab * 100:.2f}% entries out of downstream analysis.')
# # Filtering 64.88% entries out of downstream analysis.（verified 01213023）

# confidence variants call
combined_mt_het = combined_mt_het.filter_rows(combined_mt_het.variant_qc.call_rate >= 0.90)
# print('Removed variants with low confidence, %d variants remain.' % combined_mt_het.count_rows())
# # Removed variants with low confidence, 3963595 variants remain.（verified 01213023）

# rare variants: remove variants with SPARK dataset MAF of >= 1% 
rare_variant_wi_SPARK = combined_mt_het.filter_rows(combined_mt_het.variant_qc.AF[1] > rare_variants_maf, keep=False) 
# print(f"Removed variants with MAF >= {rare_variants_maf}, {rare_variant_wi_SPARK.count_rows()} variants remain.")
# # Removed variants with MAF >= 0.01, 3874914 variants remain.（verified 01213023）

# rare variants: remove variants with non-neuro gnomAD population AF >= 1% 
rare_variants_by_gnomAD = rare_variant_wi_SPARK.filter_rows(rare_variant_wi_SPARK.non_neuro_pop_AF > rare_variants_maf, keep=False)
# print(f"Removed variants with AF >= {rare_variants_maf} in non-neuro population, {rare_variants_by_gnomAD.count_rows()} variants remain.")
# # Removed variants with AF >= 0.01 in non-neuro population, 3857525 variants remain.（verified 01213023）

# select genotype to export
mt_out = rare_variants_by_gnomAD.select_cols().select_rows().select_entries('GT')
# print("Shape of output matrix table: ", mt_out.count())
# # Shape of output matrix table: (3857525, 43203)

# export to vcf.bgz
hl.export_vcf(mt_out, '/path_to/wes1_wes2_combined.deepvariant.rare1pct_variants_het.vcf.bgz')

# endregion



# region heterozygous/homozygous ref/homozygous alt 

# confident variants call
combined_mt_all = combined_mt.filter_rows(combined_mt.variant_qc.call_rate >= 0.90)
# print('Removed variants with low confidence, %d variants remain.' % combined_mt.count_rows())
# # Removed variants with low confidence, 6281625 variants remain.（verified 01213023）

# rare variants: remove variants with SPARK dataset MAF of >= 0.1% 
rare_variant_wi_SPARK = combined_mt_all.filter_rows(combined_mt_all.variant_qc.AF[1] > rare_variants_maf, keep=False) 
# print(f"Removed variants with MAF >= {rare_variants_maf}, {rare_variant_wi_SPARK.count_rows()} variants remain.")
# # Removed variants with MAF >= 0.01, 3874914 variants remain. (verified 01202023)

# rare variants: remove variants with non-neuro gnomAD population AF >= 0.1% 
rare_variants_by_gnomAD = rare_variant_wi_SPARK.filter_rows(rare_variant_wi_SPARK.non_neuro_pop_AF > rare_variants_maf, keep=False)
# print(f"Removed variants with AF >= {rare_variants_maf} in non-neuro population, {rare_variants_by_gnomAD.count_rows()} variants remain.")
# # Removed variants with AF >= 0.01 in non-neuro population, 3857525 variants remain.（verified 01213023）

# select genotype to export
mt_out = rare_variants_by_gnomAD.select_cols().select_rows().select_entries('GT')
# print("Shape of output matrix table: ", mt_out.count())
# # Shape of output matrix table: (3857525, 43203)

# export to vcf.bgz
hl.export_vcf(mt_out, '/path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.vcf.bgz')

# endregion


# endregion



# region export sample id in trio families

pedigree = hl.Pedigree.read('/path_to/wes1_wes2_combined.fam')
trio_mt = hl.trio_matrix(combined_mt, pedigree, complete_trios=True)
print("Shape of trio families matrix table:", trio_mt.count())
# Shape of trio families matrix table: (6281625, 12121)

print("Count the case trio families and control trio families:", trio_mt.aggregate_cols(hl.agg.counter(trio_mt.proband.pheno.is_case)))
# Count the case trio families and control trio families: frozendict({False: 3250, True: 8871})

# export
combined_mt.s.export("/path_to/wes1_wes2_combined.deepvariant.sample_id.tsv")
trio_mt.proband.s.export("/path_to/wes1_wes2_combined.deepvariant.trio_proband_id.tsv")
trio_mt.father.s.export("/path_to/wes1_wes2_combined.deepvariant.trio_father_id.tsv")
trio_mt.mother.s.export("/path_to/wes1_wes2_combined.deepvariant.trio_mother_id.tsv")

# endregion

# region export sample id, is_female, is_case

combined_mt.s.export("/path_to/wes1_wes2_combined.deepvariant.sample_id.tsv")
combined_mt.pheno.is_female.export("/path_to/wes1_wes2_combined.deepvariant.is_female.tsv")
combined_mt.pheno.is_case.export("/path_to/wes1_wes2_combined.deepvariant.is_case.tsv")

# endregion

