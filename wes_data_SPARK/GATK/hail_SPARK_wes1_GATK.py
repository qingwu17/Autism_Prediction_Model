#!/usr/bin/env python
# -*- coding: utf-8 -*-

import hail as hl
from hail.plot import output_notebook, show
import bokeh.plotting

# PYSPARK_SUBMIT_ARGS="--driver-memory 200g --executor-memory 350g pyspark-shell" ipython

hl.init(spark_conf=dict({'spark.driver.memory': '200g',
                         'spark.executor.memory': '350g'}))

# region: import merged vcf file

# hl.import_vcf('/users/qwu24/data/silvio/SPARK_Data_Aug_2021/pub/WES1/Variants/GATK/wes1_27281_exome.gatk.vcf.gz', reference_genome='GRCh38', \
# 	contig_recoding={'1': 'chr1', '2': 'chr2', '3': 'chr3', '4': 'chr4', '5': 'chr5', '6': 'chr6', '7': 'chr7', '8': 'chr8', '9': 'chr9', '10': 'chr10', \
# 					'11': 'chr11', '12': 'chr12', '13': 'chr13', '14': 'chr14', '15': 'chr15', '16': 'chr16', '17': 'chr17', '18': 'chr18', '19': 'chr19', '20': 'chr20', \
# 					'21': 'chr21', '22': 'chr22', 'X': 'chrX', 'Y': 'chrY'}, \
# 	force_bgz=True, array_elements_required=False).write('/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_27281_exome.gatk.mt', overwrite=True)
mt = hl.read_matrix_table('/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_27281_exome.gatk.mt')
# shape of the matrix table
mt.count()
# (5203288, 27281)

# endregion

# region: data overview

# mt.distinct_by_col().count_cols() # 27281
# # export sample id
# mt.s.export("/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_27281_exome_sample_id.tsv")

# import fam file 
SPARK_fam = hl.import_fam('/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.fam')
# annotate sample with sex and case/control status as pheno to column field
mt = mt.annotate_cols(pheno = SPARK_fam[mt.s])

# region: vep segment duplication or low complexity region

# # select certain sample, export their variants data
# mt.s.show(5)
# samples_to_keep = {'SP0000002', 'SP0000003', 'SP0000006', 'SP0000007', 'SP0000009'}
# set_to_keep = hl.literal(samples_to_keep)
# mt_1to5 = mt.filter_cols(set_to_keep.contains(mt['s']))
# hl.export_vcf(mt_1to5, '/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_27281_exome.gatk.individual_1to5.vcf.bgz')

# # vep annotated 1to5
# # import
# hl.import_vcf('/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_27281_exome.gatk.individual_1to5.vep_anno_segdup.vcf', reference_genome='GRCh38', \
# 	contig_recoding={'1': 'chr1', '2': 'chr2', '3': 'chr3', '4': 'chr4', '5': 'chr5', '6': 'chr6', '7': 'chr7', '8': 'chr8', '9': 'chr9', '10': 'chr10', \
# 					'11': 'chr11', '12': 'chr12', '13': 'chr13', '14': 'chr14', '15': 'chr15', '16': 'chr16', '17': 'chr17', '18': 'chr18', '19': 'chr19', '20': 'chr20', \
# 					'21': 'chr21', '22': 'chr22', 'X': 'chrX', 'Y': 'chrY'}, \
# 	force_bgz=True, find_replace=('nul', '.'), array_elements_required=False).write('/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_27281_exome.gatk.individual_1to5.vep_anno_segdup.mt', overwrite=True)
# mt_1to5_vep = hl.read_matrix_table('/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_27281_exome.gatk.individual_1to5.vep_anno_segdup.mt')
# # shape of the matrix table
# print("shaoe of mt_1to5: ", mt_1to5_vep.count())
# # (4681309, 5)

# # export CSQ and filters to pick select info to R
# mt_1to5_vep.describe()
# # export csq from vcf files
# mt_1to5_vep.info.CSQ[0].show(5)
# mt_1to5_vep.info.CSQ[0].export('/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_27281_exome.gatk.individual_1to5.CSQ.tsv')

# edit CSQ column using R (gnomAD_exome_flag.R)

# endregion

# endregion

# region: variant QC1

# import the variants vep annotation on segdup 
vep_anno = hl.import_table('/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_27281_exome.gatk.segdup.tsv', impute=True)
vep_anno = vep_anno.transmute(**hl.parse_variant(vep_anno.Variant, reference_genome='GRCh38')).key_by('locus', 'alleles')
mt_vep_segdup = mt.annotate_rows(info = mt.info.annotate(segdup_flag = vep_anno[mt.row_key].gnomAD_exome_flag))

mt_pass_vep_segdup = mt_vep_segdup.filter_rows(hl.len(mt_vep_segdup.info.segdup_flag) == 0)
# print('Removed variants on segment duplicated region, %d variants remain.' % mt_pass_vep_segdup.count_rows())
# # Removed variants on segment duplicated region, 5062186 variants remain.

# # filter not monoallelic
# print(mt_pass_vep_segdup.aggregate_rows(hl.agg.counter(mt_pass_vep_segdup.filters)))
# # frozendict({frozenset({'MONOALLELIC'}): 51761, None: 5010425})

mt_pass_VQSR = mt_pass_vep_segdup.filter_rows(hl.is_defined(mt_pass_vep_segdup.filters), keep=False)

# endregion

# region: genotype QC

mt_geno_QC = mt_pass_VQSR

# exclude chromosome Y for female participants
mt_geno_QC = mt_geno_QC.annotate_rows(is_y = mt_geno_QC.locus.in_y_nonpar())
mt_geno_QC = mt_geno_QC.filter_entries(mt_geno_QC.is_y & mt_geno_QC.pheno.is_female, keep=False)

# fraction_filtered_Y_female = mt_geno_QC.aggregate_entries(hl.agg.fraction(mt_geno_QC.is_y & mt_geno_QC.pheno.is_female))
# print(f'Filtering {fraction_filtered_Y_female * 100:.2f}% entries out of downstream analysis.')
# # Filtering 0.00% entries out of downstream analysis.

# DP: Approximate read depth; some reads may have been filtered
# High coverage: DP 10-1,000
filter_condition_DP = ((mt_geno_QC.DP >= 10) & (mt_geno_QC.DP <= 1000))
mt_geno_QC_DP = mt_geno_QC.filter_entries(filter_condition_DP)

# fraction_filtered_DP = mt_geno_QC.aggregate_entries(hl.agg.fraction(~filter_condition_DP))
# print(f'Filtering {fraction_filtered_DP * 100:.2f}% entries out of downstream analysis.')
# # Filtering 3.01% entries out of downstream analysis.

# AD: Allelic depths for the ref and alt alleles in the order listed" 
# allele balance
ab = mt_geno_QC_DP.AD[1] / hl.sum(mt_geno_QC_DP.AD)
# if we find the following, it is likely to be an error.
# 1) a genotype called homozygous reference with >10% alternate reads, 
# 2) a genotype called homozygous alternate with >10% reference reads, 
# 3) a genotype called heterozygote alt / ref balance < 0.25 
filter_condition_ab = ((mt_geno_QC_DP.GT.is_hom_ref() & (ab <= 0.1) & (mt_geno_QC_DP.GQ >= 25)) | (mt_geno_QC_DP.GT.is_hom_var() & (ab >= 0.9) & (mt_geno_QC_DP.PL[0] > 25)) | (mt_geno_QC_DP.GT.is_het() & (ab >= 0.25) & (mt_geno_QC_DP.PL[0] > 25)))
mt_pass_GTQC = mt_geno_QC_DP.filter_entries(filter_condition_ab)

# fraction_filtered_ab = mt_geno_QC_DP.aggregate_entries(hl.agg.fraction(~filter_condition_ab))
# print(f'Filtering {fraction_filtered_ab * 100:.2f}% entries out of downstream analysis.')
# # Filtering 1.73% entries out of downstream analysis.

# endregion

# region: sample QC

mt_sample_QC = hl.sample_qc(mt_pass_GTQC)
# print("Sample QC stats:", mt_sample_QC.aggregate_cols(hl.agg.stats(mt_sample_QC.sample_qc.call_rate)))
# # Sample QC stats: Struct(mean=0.9516898777503763, stdev=0.021166068652723913, min=0.8398509108508759, max=0.9902048229441615, n=27281, sum=25963.051554908016)

# remove duplicate samples
mt_sample_QC = mt_sample_QC.distinct_by_col()
# print('Removed duplicate samples, %d samples remain.' % mt_sample_QC.count_cols())
# # Removed duplicate samples, 27281 samples remain.

# # impute sex, check impute sex with reported sex
# imputed_sex = hl.impute_sex(mt_sample_QC.GT)
# mt_sample_QC.filter_cols(imputed_sex[mt_sample_QC.sample_id].is_female != mt_sample_QC.pheno.is_female, keep=False).count_cols()
# # 8580 samples remaining

# filter_cols removes entire columns from the matrix table
mt_pass_sampleQC = mt_sample_QC.filter_cols((mt_sample_QC.sample_qc.dp_stats.mean >= 10) & (mt_sample_QC.sample_qc.gq_stats.mean >= 20) & (mt_sample_QC.sample_qc.call_rate >= 0.90))
# print('After filter, %d samples remain.' % mt_pass_sampleQC.count_cols())
# # After filter, 26851 samples remain.

# endregion

# region: variant QC2

mt_var_QC = hl.variant_qc(mt_pass_sampleQC)

mt_pass_varQC = mt_var_QC.filter_rows((mt_var_QC.variant_qc.call_rate >= 0.10) & (mt_var_QC.variant_qc.p_value_hwe >= 1e-12))
# print('After filter, %d variants remain.' % mt_pass_varQC.count_rows())
# # After filter, 4650932 variants remain.

# export to new matrix table
mt_pass_varQC.write('/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_27281_exome.gatk.pass_QC.mt', overwrite=True)

# endregion


#region: annotation and de novo identification

mt_pass_QC = hl.read_matrix_table('/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_27281_exome.gatk.pass_QC.mt')
# shape of the matrix table
mt_pass_QC.count()


#region gnomAD non-neuro population AF annotation

# # download
# gsutil -m cp -r gs://gcp-public-data--gnomad/release/2.1.1/liftover_grch38/ht/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.ht .

# import gnomAD hail table
gnomAD_exome211_ht = hl.read_table('/users/qwu24/data/silvio/Qing_Wu/index_files/gnomAD/gnomad.exomes.r2.1.1.sites.liftover_grch38.ht')
# gnomAD_exome211_ht.count()
# # 17201296

# gnomAD_exome211_ht.freq_index_dict.export('/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/gnomAD_exome211_ht_freq_index_dict.txt')
# gnomAD_exome211_ht.select(gnomad_freq=gnomAD_exome211_ht.freq[102].AF).export("/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/gnomAD_exome211_AF_non_neuro.tsv")

# annotate variants by gnomAD non-neuro population AF
mt_pass_QC = mt_pass_QC.annotate_rows(gnomAD_non_neuro_pop_AF=gnomAD_exome211_ht[mt_pass_QC.row_key].freq[102].AF)
# handle missing data from gnomAD
mt_pass_QC = mt_pass_QC.annotate_rows(non_neuro_pop_AF=hl.coalesce(mt_pass_QC.gnomAD_non_neuro_pop_AF, hl.float(0)))

# endregion

# region: de novo

# de novo detection
pedigree = hl.Pedigree.read('/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.fam')
# priors = hl.import_table('/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/gnomAD_exome211_AF_hg38_non_neuro_na.tsv', impute=True)
# priors = priors.transmute(**hl.parse_variant(priors.Variant, reference_genome='GRCh38')).key_by('locus', 'alleles')

de_novo_ht = hl.de_novo(mt_pass_QC, pedigree, pop_frequency_prior=mt_pass_QC.non_neuro_pop_AF)

# check the number of de novo variants with H/M/L confidence level
print("de novo: ", de_novo_ht.aggregate(hl.agg.counter(de_novo_ht.confidence)))
# de novo:  frozendict({'HIGH': 17, 'LOW': 49, 'MEDIUM': 1})

# # export
# de_novo_ht.select(
#     is_case=de_novo_ht.proband.pheno.is_case, 
#     is_female=de_novo_ht.proband.pheno.is_female, 
# 	confidence=de_novo_ht.confidence
# 	).export('/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_27281_exome.gatk.dnvs.tsv')

# endregion

# endregion
