#!/usr/bin/env python
# -*- coding: utf-8 -*-

# PYSPARK_SUBMIT_ARGS="--driver-memory 200g --executor-memory 125g --executor-cores 2 --num-executors 8 pyspark-shell" ipython
# PYSPARK_SUBMIT_ARGS="--driver-memory 200g --executor-memory 350g pyspark-shell" ipython

import hail as hl
from hail.plot import output_notebook, show
import bokeh.plotting

hl.init(spark_conf=dict({'spark.driver.memory': '200g',
                         'spark.executor.memory': '350g'}))

# region: read wes1 and wes2, merge, and annotate variants by AF from gnomAD
wes1_mt = hl.read_matrix_table('/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/pub/wes1_27281_exome.deepvariant.pass_QC.mt')
print("Shape of wes1: ", wes1_mt.count())
# GATK (4650932, 26851) -> (4276376, 27266)

wes2_mt = hl.read_matrix_table('/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/pub/wes2_15995_exome.deepvariant.pass_QC.mt')
print("Shape of wes2: ", wes2_mt.count())
# GATK (4617039, 15846) -> (3875705, 15986)

combined_mt = wes1_mt.union_cols(wes2_mt, row_join_type='outer')
print("Shape of combined mt: ", combined_mt.count())
# GATK (6906370, 43161) -> DeepVariant (6281625, 43252)

# remove duplicates
combined_mt = combined_mt.distinct_by_col()
print("After remove duplicated samples across wes1 and wes2, sample size: ", combined_mt.count_cols())
# 43112 -> 43203

# region: gnomAD non-neuro population AF annotation

# # download
# gsutil -m cp -r gs://gcp-public-data--gnomad/release/2.1.1/liftover_grch38/ht/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.ht .

# import gnomAD hail table
gnomAD_exome211_ht = hl.read_table('/users/qwu24/data/silvio/Qing_Wu/index_files/gnomAD/gnomad.exomes.r2.1.1.sites.liftover_grch38.ht')
# gnomAD_exome211_ht.count() # 17201296


# gnomAD_exome211_ht.freq_index_dict.export('/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/gnomAD_exome211_ht_freq_index_dict.txt')
# gnomAD_exome211_ht.select(gnomad_freq=gnomAD_exome211_ht.freq[102].AF).export("/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/gnomAD_exome211_AF_non_neuro.tsv")

# annotate variants by gnomAD non-neuro population AF
combined_mt = combined_mt.annotate_rows(gnomAD_non_neuro_pop_AF=gnomAD_exome211_ht[combined_mt.row_key].freq[102].AF)
# handle missing data from gnomAD
combined_mt = combined_mt.annotate_rows(non_neuro_pop_AF=hl.coalesce(combined_mt.gnomAD_non_neuro_pop_AF, hl.float(0)))

# endregion

# region: de novo 

pedigree = hl.Pedigree.read('/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.fam')

de_novo_ht = hl.de_novo(combined_mt, pedigree, pop_frequency_prior=combined_mt.non_neuro_pop_AF, \
	min_gq=20, min_p=0.05, max_parent_ab=0.03, min_child_ab=0.25, min_dp_ratio=0.3, ignore_in_sample_allele_frequency=False)

# # check the number of de novo variants with H/M/L confidence level
# print("de novo: ", de_novo_ht.aggregate(hl.agg.counter(de_novo_ht.confidence)))
# # de novo:  frozendict({'HIGH': 24017, 'LOW': 7292, 'MEDIUM': 3214})6 + 1) / 4197]

# # export
# de_novo_ht.select(
#     is_case = de_novo_ht.proband.pheno.is_case,
#     is_female = de_novo_ht.proband.pheno.is_female,
# 	  confidence=de_novo_ht.confidence
#     ).export('/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.dnvs.tsv')

# region: de novo using default parameters setting

# de_novo_def_ht = hl.de_novo(combined_mt, pedigree, pop_frequency_prior=combined_mt.non_neuro_pop_AF)

# # check the number of de novo variants with H/M/L confidence level
# print("de novo with default parameters : ", de_novo_def_ht.aggregate(hl.agg.counter(de_novo_def_ht.confidence)))
# # de novo with default parameters :  frozendict({'HIGH': 24461, 'LOW': 7729, 'MEDIUM': 3541})

# # export
# de_novo_def_ht.select(
#     is_case=de_novo_def_ht.proband.pheno.is_case,
#     is_female=de_novo_def_ht.proband.pheno.is_female,
# 	  confidence=de_novo_def_ht.confidence
#     ).export('/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.dnvs_def.tsv')

# endregion

# endregion

# region: rare variants

# region: heterozygous or hemizygous
combined_mt_het = combined_mt.filter_entries((combined_mt.GT.is_het()))

fraction_filtered_ab = combined_mt.aggregate_entries(hl.agg.fraction(~(combined_mt.GT.is_het())))
print(f'Filtering {fraction_filtered_ab * 100:.2f}% entries out of downstream analysis.')
# Filtering 64.88% entries out of downstream analysis.

# confidence variants call
combined_mt_het = combined_mt_het.filter_rows(combined_mt_het.variant_qc.call_rate >= 0.90)
print('Removed variants with low confidence, %d variants remain.' % combined_mt_het.count_rows())
# Removed variants with low confidence, 3963595 variants remain.

# rare variants: remove variants with SPARK dataset MAF of >= 1% 
rare_variant_wi_SPARK = combined_mt_het.filter_rows(combined_mt_het.variant_qc.AF[1] > 0.01, keep=False) 
# print('Removed variants with MAF >= 0.01, %d variants remain.' % rare_variant_wi_SPARK.count_rows())
# # Removed variants with MAF >= 0.01, 3874914 variants remain.

# rare variants: remove variants with non-neuro gnomAD population AF >= 1% 
rare_variants_by_gnomAD = rare_variant_wi_SPARK.filter_rows(rare_variant_wi_SPARK.non_neuro_pop_AF > 0.01, keep=False)
# print('Removed variants with AF >= 0.01 in non-neuro population, %d variants remain.' % rare_variants_by_gnomAD.count_rows())
# # Removed variants with AF >= 0.001 in non-neuro population, 3857525 variants remain.

# select genotype to export
mt_out = rare_variants_by_gnomAD.select_cols().select_rows().select_entries('GT')
# print("Shape of output matrix table: ", mt_out.count())
# # Shape of output matrix table: (3857525, 43203)

# # export to vcf.bgz
# hl.export_vcf(mt_out, '/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.vcf.bgz')

# annotate rare variants by dnv confidence level
de_novo_ht_tmp = de_novo_ht.key_by('locus', 'alleles')
mt_out_w_dnv = mt_out.annotate_rows(dnv_confidence = de_novo_ht_tmp[mt_out.locus, mt_out.alleles].confidence)
# print("Count of dnvs in rare mt output: ", mt_out_w_dnv.filter_rows(hl.is_defined(mt_out_w_dnv.dnv_confidence), keep = True).count_rows() )
# # Count of dnvs in rare mt output:  19916

# # export
# mt_out_w_dnv.dnv_confidence.export('/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.w_dnv.tsv')

# endregion

# region: heterozygous/homozygous ref/homozygous alt 

# confident variants call
combined_mt_all = combined_mt.filter_rows(combined_mt.variant_qc.call_rate >= 0.90)
# print('Removed variants with low confidence, %d variants remain.' % combined_mt.count_rows())
# # Removed variants with low confidence, 6281625 variants remain.

# rare variants: remove variants with SPARK dataset MAF of >= 0.1% 
rare_variant_wi_SPARK = combined_mt_all.filter_rows(combined_mt_all.variant_qc.AF[1] > 0.01, keep=False) 
# print('Removed variants with MAF >= 0.01, %d variants remain.' % rare_variant_wi_SPARK.count_rows())
# # Removed variants with MAF >= 0.01, 3874914 variants remain.

# rare variants: remove variants with non-neuro gnomAD population AF >= 0.1% 
rare_variants_by_gnomAD = rare_variant_wi_SPARK.filter_rows(rare_variant_wi_SPARK.non_neuro_pop_AF > 0.01, keep=False)
# print('Removed variants with AF >= 0.01 in non-neuro population, %d variants remain.' % rare_variants_by_gnomAD.count_rows())
# # Removed variants with AF >= 0.001 in non-neuro population, 3857525 variants remain.

# select genotype to export
mt_out = rare_variants_by_gnomAD.select_cols().select_rows().select_entries('GT')
# print("Shape of output matrix table: ", mt_out.count())
# # Shape of output matrix table: (3857525, 43203)

# # export to vcf.bgz
# hl.export_vcf(mt_out, '/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants.vcf.bgz')

# annotate rare variants by dnv confidence level
de_novo_ht_tmp = de_novo_ht.key_by('locus', 'alleles')
mt_out_w_dnv = mt_out.annotate_rows(dnv_confidence = de_novo_ht_tmp[mt_out.locus, mt_out.alleles].confidence)
# print("Count of dnvs in rare mt output: ", mt_out_w_dnv.filter_rows(hl.is_defined(mt_out_w_dnv.dnv_confidence), keep = True).count_rows() )
# # Count of dnvs in rare mt output:  19916

# export
mt_out_w_dnv.dnv_confidence.export('/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants.w_dnv.tsv')


# endregion

# endregion

# region: seperate trio matrix and case-control matrix

pedigree = hl.Pedigree.read('/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.fam')
trio_mt = hl.trio_matrix(combined_mt, pedigree, complete_trios=True)
trio_mt.count()
# (6281625, 12121)

trio_mt.aggregate_cols(hl.agg.counter(trio_mt.proband.pheno.is_case))
# frozendict({False: 3250, True: 8871})

# combined_mt.s.export("/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.sample_id.tsv")
# trio_mt.proband.s.export("/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.trio_proband_id.tsv")
# trio_mt.father.s.export("/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.trio_father_id.tsv")
# trio_mt.mother.s.export("/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.trio_mother_id.tsv")

trio_sample_id = hl.import_table('/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.trio_id_flag.tsv').key_by('s')
trio_sample_id.show()

combined_mt = combined_mt.annotate_cols(in_trio = trio_sample_id[combined_mt.s].in_trio)
cc_mt = combined_mt.filter_cols(combined_mt.in_trio == "FALSE", keep=True)
# cc_mt.count()
# # 6909001, 16753

# endregion

# region: export sample id, is_female, is_case

combined_mt.s.export("/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.sample_id.tsv")
combined_mt.pheno.is_female.export("/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.is_female.tsv")
combined_mt.pheno.is_case.export("/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.is_case.tsv")

# endregion

