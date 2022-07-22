#!/usr/bin/env python
# -*- coding: utf-8 -*-

import hail as hl
from hail.plot import output_notebook, show
import bokeh.plotting

hl.init(spark_conf={'spark.driver.memory': '350g'})

# region: import merged vcf file

# hl.import_vcf('/users/qwu24/data/silvio/SPARK_Data_Aug_2021/pub/WES2/Variants/DeepVariant/wes2_15995_exome.deepvariant.vcf.gz', reference_genome='GRCh38', \
# 	contig_recoding={'1': 'chr1', '2': 'chr2', '3': 'chr3', '4': 'chr4', '5': 'chr5', '6': 'chr6', '7': 'chr7', '8': 'chr8', '9': 'chr9', '10': 'chr10', \
# 					'11': 'chr11', '12': 'chr12', '13': 'chr13', '14': 'chr14', '15': 'chr15', '16': 'chr16', '17': 'chr17', '18': 'chr18', '19': 'chr19', '20': 'chr20', \
# 					'21': 'chr21', '22': 'chr22', 'X': 'chrX', 'Y': 'chrY'}, \
# 	force_bgz=True, array_elements_required=False).write('/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes2_15995_exome.deepvariant.mt', overwrite=True)
mt = hl.read_matrix_table('/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes2_15995_exome.deepvariant.mt')
# shape of the matrix table
mt.count()
# (4212994, 15995)

# endregion

# region: data overview

# mt.distinct_by_col().count_cols() # 15995
# # export sample id
# mt.s.export("/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes2_15995_exome.deepvariant.sample_id.tsv")

# import fam file 
SPARK_fam = hl.import_fam('/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.fam')
# annotate sample with sex and case/control status as pheno to column field
mt = mt.annotate_cols(pheno = SPARK_fam[mt.s])

# # count of is_case 
# # mt.aggregate_rows(hl.agg.counter(mt.pheno.is_case)) # not working
# mt.filter_cols(mt.pheno.is_case == True).count_cols() # 12109
# mt.filter_cols(mt.pheno.is_case == False).count_cols() # 3886
# mt.filter_cols(hl.is_defined(mt.pheno.is_case)).count_cols() # 15995

# # count of is_female 
# mt.filter_cols(mt.pheno.is_female == True).count_cols() # 4902
# mt.filter_cols(mt.pheno.is_female == False).count_cols() # 11093
# mt.filter_cols(hl.is_defined(mt.pheno.is_female)).count_cols() # 15995

# region: vep segment duplication or low complexity region

# # select certain sample, export their variants data
# mt.s.show(5)
# samples_to_keep = {'SP0000017', 'SP0000104', 'SP0000143', 'SP0000216', 'SP0000229'}
# set_to_keep = hl.literal(samples_to_keep)
# mt_1to5 = mt.filter_cols(set_to_keep.contains(mt['s']))
# hl.export_vcf(mt_1to5, '/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes2_15995_exome.deepvariant.individual_1to5.vcf.bgz')

# # vep annotated 1to5
# # import
# hl.import_vcf('/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes2_15995_exome.deepvariant.individual_1to5.vep_anno_segdup.vcf', reference_genome='GRCh38', \
# 	contig_recoding={'1': 'chr1', '2': 'chr2', '3': 'chr3', '4': 'chr4', '5': 'chr5', '6': 'chr6', '7': 'chr7', '8': 'chr8', '9': 'chr9', '10': 'chr10', \
# 					'11': 'chr11', '12': 'chr12', '13': 'chr13', '14': 'chr14', '15': 'chr15', '16': 'chr16', '17': 'chr17', '18': 'chr18', '19': 'chr19', '20': 'chr20', \
# 					'21': 'chr21', '22': 'chr22', 'X': 'chrX', 'Y': 'chrY'}, \
# 	force_bgz=True, find_replace=('nul', '.'), array_elements_required=False).write('/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes2_15995_exome.deepvariant.individual_1to5.vep_anno_segdup.mt', overwrite=True)
# mt_1to5_vep = hl.read_matrix_table('/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes2_15995_exome.deepvariant.individual_1to5.vep_anno_segdup.mt')
# # shape of the matrix table
# print("shaoe of mt_1to5: ", mt_1to5_vep.count())
# # (4212994, 5)

# # export CSQ and filters to pick select info to R
# mt_1to5_vep.describe()
# # export csq from vcf files
# mt_1to5_vep.info.CSQ[0].show(5)
# mt_1to5_vep.info.CSQ[0].export('/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes2_15995_exome.deepvariant.individual_1to5.CSQ.tsv')

# edit CSQ column using R (gnomAD_exome_flag.R)

# endregion

# endregion

#region: variant QC1

# import the segdup table
vep_anno = hl.import_table('/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes2_15995_exome.deepvariant.segdup.tsv', impute=True)
vep_anno = vep_anno.transmute(**hl.parse_variant(vep_anno.Variant, reference_genome='GRCh38')).key_by('locus', 'alleles')
mt_vep_segdup = mt.annotate_rows(info = mt.info.annotate(segdup_flag = vep_anno[mt.row_key].gnomAD_exome_flag))

mt_pass_vep_segdup = mt_vep_segdup.filter_rows(hl.len(mt_vep_segdup.info.segdup_flag) == 0)
# print('Removed variants on segment duplicated region, %d variants remain.' % mt_pass_vep_segdup.count_rows())
# # Removed variants on segment duplicated region, 4094344 variants remain.

# # filter not monoallelic
# print(mt_pass_vep_segdup.aggregate_rows(hl.agg.counter(mt_pass_vep_segdup.filters)))
# # frozendict({frozenset({'MONOALLELIC'}): 8835, None: 4085509})

# filter not monoallelic
mt_pass_VQSR = mt_pass_vep_segdup.filter_rows(hl.is_defined(mt_pass_vep_segdup.filters), keep=False)
# print('Removed variants did not pass filter, %d variants remain.' % mt_pass_VQSR.count_rows())
# # Removed variants did not pass filter, 4085509 variants remain.

#endregion

#region: genotype QC

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
# # Filtering 1.82% entries out of downstream analysis.

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
# # Filtering 0.44% entries out of downstream analysis.

#endregion

#region: sample QC

mt_sample_QC = hl.sample_qc(mt_pass_GTQC)
# print("sample QC stats:", mt_sample_QC.aggregate_cols(hl.agg.stats(mt_sample_QC.sample_qc.call_rate)))
# # sample QC stats: Struct(mean=0.9773431694125141, stdev=0.015463212797472384, min=0.8801319492871023, max=0.9964183165426879, n=15995, sum=15632.603994753163)

# remove duplicate samples
mt_sample_QC = mt_sample_QC.distinct_by_col()
# print('Removed duplicate samples, %d samples remain.' % mt_sample_QC.count_cols())
# # Removed duplicate samples, 15995 samples remain.

# # impute sex, check impute sex with reported sex
# imputed_sex = hl.impute_sex(mt_sample_QC.GT)
# print("Samples with same sex as imputed sex: ", mt_sample_QC.filter_cols(imputed_sex[mt_sample_QC.s].is_female != mt_sample_QC.pheno.is_female, keep=False).count_cols())
# #  samples remaining

# filter_cols removes entire columns from the matrix table
mt_pass_sampleQC = mt_sample_QC.filter_cols((mt_sample_QC.sample_qc.dp_stats.mean >= 10) & (mt_sample_QC.sample_qc.gq_stats.mean >= 20) & (mt_sample_QC.sample_qc.call_rate >= 0.90))
# print('After filter, %d samples remain.' % mt_pass_sampleQC.count_cols())
# After filter, 15986 samples remain.

#endregion

#region: variant QC2

mt_var_QC = hl.variant_qc(mt_pass_sampleQC)

mt_pass_varQC = mt_var_QC.filter_rows((mt_var_QC.variant_qc.call_rate >= 0.10) & (mt_var_QC.variant_qc.p_value_hwe >= 1e-12))
# print('After filter, %d variants remain.' % mt_pass_varQC.count_rows())
# After filter, 3875705 variants remain.

# export to new matrix table
mt_pass_varQC.write('/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes2_15995_exome.deepvariant.pass_QC.mt', overwrite=True)

#endregion

