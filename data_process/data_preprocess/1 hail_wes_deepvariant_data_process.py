#!/usr/bin/env python
# -*- coding: utf-8 -*-

import hail as hl
import os

hl.init(spark_conf={'spark.driver.memory': '256g'})

# region import vcf file and converted to hail matrix table

wes_deepvariant_mt_filename = '/path_to/wes.deepvariant.mt'

if os.path.isfile(wes_deepvariant_mt_filename):

    print("Converting wes_deepvariant vcf file to hail matrix table ...")

    hl.import_vcf('/path_to/wes.deepvariant.vcf.gz', reference_genome='GRCh38', \
        contig_recoding={'1': 'chr1', '2': 'chr2', '3': 'chr3', '4': 'chr4', '5': 'chr5', '6': 'chr6', '7': 'chr7', '8': 'chr8', '9': 'chr9', '10': 'chr10', \
                        '11': 'chr11', '12': 'chr12', '13': 'chr13', '14': 'chr14', '15': 'chr15', '16': 'chr16', '17': 'chr17', '18': 'chr18', '19': 'chr19', '20': 'chr20', \
                        '21': 'chr21', '22': 'chr22', 'X': 'chrX', 'Y': 'chrY'}, \
        force_bgz=True, array_elements_required=False).write('/path_to/wes.deepvariant.mt', overwrite=True)

else:

    print(f"wes_deepvariant vcf file has been converted to hail matrix table {wes_deepvariant_mt_filename}.")


mt = hl.read_matrix_table('/path_to/wes.deepvariant.mt')
print("Shape of wes_deepvariant hail matrix table:", mt.count())
# wes1: (4681309, 27281)
# wes2: (4212994, 15995)
# iwes1: (8336937, 70487)

# endregion

# region import phenotypic data 

# import fam file and annotate sample with sex and case/control status as pheno to column field
SPARK_fam = hl.import_fam('/path_to/wes_deepvariant.fam')
mt = mt.annotate_cols(pheno = SPARK_fam[mt.s])

# count number of cases and controls in genotypic data
num_tot_participants = mt.distinct_by_col().count_cols() 
num_tot_participants_by_disease = mt.filter_cols(hl.is_defined(mt.pheno.is_case)).count_cols() 
num_tot_participants_by_gender = mt.filter_cols(hl.is_defined(mt.pheno.is_female)).count_cols() 

if num_tot_participants == num_tot_participants_by_disease:    
    print("All participants has disease status.")

else:   
    print("{} participants do not have disease status.", format(num_tot_participants - num_tot_participants_by_disease))


if num_tot_participants == num_tot_participants_by_gender:    
    print("All participants has gender record.")

else:
    print("{} participants do not have gender record.", format(num_tot_participants - num_tot_participants_by_gender))


num_male_cases = mt.filter_cols((mt.pheno.is_case == True) & (mt.pheno.is_female == False)).count_cols()
num_male_controls = mt.filter_cols((mt.pheno.is_case == False) & (mt.pheno.is_female == False)).count_cols()
num_female_cases = mt.filter_cols((mt.pheno.is_case == True) & (mt.pheno.is_female == True)).count_cols() 
num_female_controls = mt.filter_cols((mt.pheno.is_case == False) & (mt.pheno.is_female == True)).count_cols()

print("Total number of {} participants in wes1 dataset. Among them, {} are male cases, {} are male controls, {} are female cases, and {} are female controls.", 
        format(num_tot_participants, num_male_cases, num_male_controls, num_female_cases, num_female_controls)) 


# export sample id
sample_id_tsv_filename = "/path_to/wes.deepvariant.sample_id.tsv"
mt.s.export(sample_id_tsv_filename)

# endregion

# region vep annotation to variants on segment duplication or low complexity region

participants_1to5_vcf_filename = "/path_to/wes.deepvariant.individual_1to5.vcf.bgz"

if os.path.isfile(participants_1to5_vcf_filename):

    print(f"The file {participants_1to5_vcf_filename} exist. Continue ...")

else:
    # select certain sample, export their variants data
    print(f"Generating {participants_1to5_vcf_filename} for first 5 participants in the record ...")
    print(mt.s.show(5))
    samples_to_keep = {'SP0000002', 'SP0000003', 'SP0000006', 'SP0000007', 'SP0000009'}
    set_to_keep = hl.literal(samples_to_keep)
    mt_1to5 = mt.filter_cols(set_to_keep.contains(mt['s']))
    hl.export_vcf(mt_1to5, '/path_to/wes.deepvariant.individual_1to5.vcf.bgz')
    print(f"The file {participants_1to5_vcf_filename} for first 5 participants in the record has been generated. VEP annotation outside in bash ...")


# vep annotated by vep_annotation_on_segment_duplication.sh

segment_duplication_vep_annotation_CSQ_tsv_filename = "/path_to/wes.deepvariant.individual_1to5.CSQ.tsv"
segment_duplication_vep_annotation_tsv_filename = "/path_to/wes.deepvariant.segdup.tsv"

if os.path.isfile(segment_duplication_vep_annotation_tsv_filename):

    print(f"The file {segment_duplication_vep_annotation_tsv_filename} exist. Continue ...")

elif os.path.isfile(segment_duplication_vep_annotation_CSQ_tsv_filename):

    print(f"The file {segment_duplication_vep_annotation_CSQ_tsv_filename} exist. Need to edit CSQ column using R (get_gnomAD_exome_flag.R).")

else:

    print("Import VEP annotation and embed to matrix table ... ")
    # import vep annotated 1to5 vcf file and converted to matrix table
    hl.import_vcf('/path_to/wes.deepvariant.individual_1to5.vep_anno_segdup.vcf', reference_genome='GRCh38', \
        contig_recoding={'1': 'chr1', '2': 'chr2', '3': 'chr3', '4': 'chr4', '5': 'chr5', '6': 'chr6', '7': 'chr7', '8': 'chr8', '9': 'chr9', '10': 'chr10', \
                        '11': 'chr11', '12': 'chr12', '13': 'chr13', '14': 'chr14', '15': 'chr15', '16': 'chr16', '17': 'chr17', '18': 'chr18', '19': 'chr19', '20': 'chr20', \
                        '21': 'chr21', '22': 'chr22', 'X': 'chrX', 'Y': 'chrY'}, \
        force_bgz=True, find_replace=('nul', '.'), array_elements_required=False).write('/path_to/wes.deepvariant.individual_1to5.vep_anno_segdup.mt', overwrite=True)
    
    mt_1to5_vep = hl.read_matrix_table('/path_to/wes.deepvariant.individual_1to5.vep_anno_segdup.mt')

    # shape of the matrix table
    print("Number of variants been vep annotated: ", mt_1to5_vep.count_rows()) 

    mt_1to5_vep.info.CSQ[0].export('/path_to/wes.deepvariant.individual_1to5.CSQ.tsv')

    print("Need to edit CSQ column using R (gnomAD_exome_flag.R).")
   
# edit CSQ column using R (get_gnomAD_exome_flag.R)

# endregion


# region variant quality control, round 1

# import the variants vep annotation on segdup 
vep_anno = hl.import_table('/path_to/wes.deepvariant.segdup.tsv', impute=True)
vep_anno = vep_anno.transmute(**hl.parse_variant(vep_anno.Variant, reference_genome='GRCh38')).key_by('locus', 'alleles')
mt_vep_segdup = mt.annotate_rows(info = mt.info.annotate(segdup_flag = vep_anno[mt.row_key].gnomAD_exome_flag))

mt_pass_vep_segdup = mt_vep_segdup.filter_rows(hl.len(mt_vep_segdup.info.segdup_flag) == 0)
print('Removed variants on segment duplicated region, %d variants remain.' % mt_pass_vep_segdup.count_rows())
# wes1: Removed variants on segment duplicated region, 4549169 variants remain.
# wes2: Removed variants on segment duplicated region, 4094344 variants remain.
# iwes1: Removed variants on segment duplicated region, 8144172 variants remain.


# filter not monoallelic
mt_pass_VQSR = mt_pass_vep_segdup.filter_rows(hl.is_defined(mt_pass_vep_segdup.filters), keep=False)
print('Removed variants did not pass filter, %d variants remain.' % mt_pass_VQSR.count_rows())
# wes1: Removed variants did not pass filter, 4537198 variants remain.
# wes2: Removed variants did not pass filter, 4085509 variants remain.
# iwes1: Removed variants did not pass filter,  8110073 variants remain.


# endregion

# region genotype quality control

mt_geno_QC = mt_pass_VQSR

# exclude chromosome Y for female participants
mt_geno_QC = mt_geno_QC.annotate_rows(is_y = mt_geno_QC.locus.in_y_nonpar())
mt_geno_QC = mt_geno_QC.filter_entries(mt_geno_QC.is_y & mt_geno_QC.pheno.is_female, keep=False)

fraction_filtered_Y_female = mt_geno_QC.aggregate_entries(hl.agg.fraction(mt_geno_QC.is_y & mt_geno_QC.pheno.is_female))
print(f'Filtering {fraction_filtered_Y_female * 100:.2f}% entries out of downstream analysis.')
# wes1: Filtering 0.00% entries out of downstream analysis.
# wes2: Filtering 0.00% entries out of downstream analysis.
# iwes1: Filtering 0.00% entries out of downstream analysis.

# DP: Approximate read depth; some reads may have been filtered
# High coverage: DP 10-1,000
filter_condition_DP = ((mt_geno_QC.DP >= 10) & (mt_geno_QC.DP <= 1000))
mt_geno_QC_DP = mt_geno_QC.filter_entries(filter_condition_DP)

fraction_filtered_DP = mt_geno_QC.aggregate_entries(hl.agg.fraction(~filter_condition_DP))
print(f'Filtering {fraction_filtered_DP * 100:.2f}% entries out of downstream analysis.')
# wes1: Filtering 2.54% entries out of downstream analysis.
# wes2: Filtering 1.82% entries out of downstream analysis.
# iwes1: Filtering 1.91% entries out of downstream analysis.

# AD: Allelic depths for the ref and alt alleles in the order listed" 
# allele balance
ab = mt_geno_QC_DP.AD[1] / hl.sum(mt_geno_QC_DP.AD)
# if we find the following, it is likely to be an error.
# 1) a genotype called homozygous reference with >10% alternate reads, 
# 2) a genotype called homozygous alternate with >10% reference reads, 
# 3) a genotype called heterozygote alt / ref balance < 0.25 
filter_condition_ab = ((mt_geno_QC_DP.GT.is_hom_ref() & (ab <= 0.1) & (mt_geno_QC_DP.GQ >= 25)) | (mt_geno_QC_DP.GT.is_hom_var() & (ab >= 0.9) & (mt_geno_QC_DP.PL[0] > 25)) | (mt_geno_QC_DP.GT.is_het() & (ab >= 0.25) & (mt_geno_QC_DP.PL[0] > 25)))
mt_pass_GTQC = mt_geno_QC_DP.filter_entries(filter_condition_ab)

fraction_filtered_ab = mt_geno_QC_DP.aggregate_entries(hl.agg.fraction(~filter_condition_ab))
print(f'Filtering {fraction_filtered_ab * 100:.2f}% entries out of downstream analysis.')
# wes1: Filtering 0.49% entries out of downstream analysis.
# wes2: Filtering 0.44% entries out of downstream analysis.
# iwes1: Filtering 0.37% entries out of downstream analysis.

# endregion

# region sample quality control

mt_sample_QC = hl.sample_qc(mt_pass_GTQC)
print("Sample QC stats:", mt_sample_QC.aggregate_cols(hl.agg.stats(mt_sample_QC.sample_qc.call_rate)))
# wes1: Sample QC stats: Struct(mean=0.9696631762716488, stdev=0.016718991048648806, min=0.8759487683808377, max=0.9963686839322419, n=27281, sum=26453.38111186685)
# wes2: Sample QC stats: Struct(mean=0.9773431694125141, stdev=0.015463212797472384, min=0.8801319492871023, max=0.9964183165426879, n=15995, sum=15632.603994753163)
# iwes1: Sample QC stats: Struct(mean=0.9771233940397143, stdev=0.01623572327110824, min=0.8599428143248525, max=0.9975589615531204, n=70487, sum=68874.49667567734)

# remove duplicate samples
mt_sample_QC = mt_sample_QC.distinct_by_col()
print('Removed duplicate samples, %d samples remain.' % mt_sample_QC.count_cols())
# wes1: Removed duplicate samples, 27281 samples remain.
# wes2: Removed duplicate samples, 15995 samples remain.
# iwes1: Removed duplicate samples, 70487 samples remain.

# filter_cols removes entire columns from the matrix table
mt_pass_sampleQC = mt_sample_QC.filter_cols((mt_sample_QC.sample_qc.dp_stats.mean >= 10) & (mt_sample_QC.sample_qc.gq_stats.mean >= 20) & (mt_sample_QC.sample_qc.call_rate >= 0.90))
print('After filter, %d samples remain.' % mt_pass_sampleQC.count_cols())
# wes1: After filter, 27266 samples remain.
# wes2: After filter, 15986 samples remain.
# iwes1: After filter, 70464 samples remain.

# endregion

# region variant quality control round 2

mt_var_QC = hl.variant_qc(mt_pass_sampleQC)

mt_pass_varQC = mt_var_QC.filter_rows((mt_var_QC.variant_qc.call_rate >= 0.10) & (mt_var_QC.variant_qc.p_value_hwe >= 1e-12))
print('After filter, %d variants remain.' % mt_pass_varQC.count_rows())
# wes1: After filter, 4276376 variants remain.
# wes2: After filter, 3875705 variants remain.
# iwes1: After filter, 7390523 variants remain.

# export to new matrix table
mt_pass_varQC.write('/path_to/wes.deepvariant.pass_QC.mt', overwrite=True)

# endregion

