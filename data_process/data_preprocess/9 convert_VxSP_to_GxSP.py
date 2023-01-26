#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
   This file read and write data.
"""

import numpy as np

# 
sample_dir = "/path_to/"
sample_name = "wes1_wes2_combined"
variant_type = "" # other option including "PTVs.", "PTVs_MisA."", "PTVs_MisAB.", "MisAB.", "MisA.", "MisB.", "non_impactful."

# read input file of variant to gene indices 
V2G_idx_file = sample_dir + sample_name + ".deepvariant.rare1pct_variants.idx_var2gene_lst_df.unconsecutive_indices_splited." + variant_type + "txt"
# read input VxSP matrix, and export to GxSP matrix
VxSP_matrix_file = sample_dir + sample_name + ".deepvariant.rare1pct_variants_by_sample_matrix_cleaned." + variant_type + "txt"
GxSP_matrix_file = sample_dir + sample_name + ".deepvariant.rare1pct_variants_by_sample_matrix_cleaned_2GxSP." + variant_type + "txt"


# read variant to gene index file, each line represent the index of variant to one gene
V2G_idx_lst = []
with open(V2G_idx_file) as tmp_file:
    for line in tmp_file:
        var2gene_idx_w_gene = line.rstrip().split('\t')
        print(var2gene_idx_w_gene[0])
        var2gene_idx = [int(x) for x in var2gene_idx_w_gene[1:]]
        # print(var2gene_idx)
        V2G_idx_lst.append(var2gene_idx)

print(sample_name, variant_type, len(V2G_idx_lst))
# 19626 wes1_wes2.deepvariant rare het exonic only (19117 unique)
# 16832 wes1_wes2.deepvariant rare het (PTVs) exonic only (16832 unique)
# 17680 wes1_wes2.deepvariant rare het (PTVs_MisA) exonic only ( unique)
# 18528 wes1_wes2.deepvariant rare het (PTVs_MisAB) exonic only ( unique)
# 3813 wes1_wes2.deepvariant rare het (MisA) exonic only ( unique)
# 9624 wes1_wes2.deepvariant rare het (MisAB) exonic only ( unique)
# 9495 wes1_wes2.deepvariant rare het (MisB) exonic only ( unique)
# 19181 wes1_wes2.deepvariant rare het (MisC) exonic only ( unique)

# 895 wes1_wes2.deepvariant rare de novo PTVs exonic only
# 1438 wes1_wes2.deepvariant rare de novo MisAB exonic only
# 4241 wes1_wes2.deepvariant rare de novo MisC exonic only
# 16754 wes1_wes2.deepvariant rare inherited PTVs exonic only
# 9615 wes1_wes2.deepvariant rare inherited MisAB exonic only

# 19877 wes_70487_exome.deepvariant rare het exonic only
# 18385 wes_70487_exome.deepvariant rare het (PTVs) exonic only
# 18830 wes_70487_exome.deepvariant rare het (PTVs_MisA) exonic only
# 19121 wes_70487_exome.deepvariant rare het (PTVs_MisAB) exonic only
# XXXX wes_70487_exome.deepvariant rare het (MisA) exonic only
# 9997 wes_70487_exome.deepvariant rare het (MisAB) exonic only
# XXXX wes_70487_exome.deepvariant rare het (MisB) exonic only
# 19435 wes_70487_exome.deepvariant rare het (MisC) exonic only


# count the number of variant to each gene by line
num_var2gene = [len(x) for x in V2G_idx_lst]

output_file = open(GxSP_matrix_file, 'a')

VxSP_lst = []
idx_gene = 0
num_var = num_var2gene[idx_gene]
with open(VxSP_matrix_file, "r", buffering = 200*(1024**2)) as f:
    for line in f:
        
        if num_var > 0:
            VxSP_lst.append(line.rstrip().split('\t'))
        num_var -= 1
        
        if num_var == 0:
            # sum up the VxSP 2D array by row to GxSP 1D array
            VxSP_2D_array = np.array(VxSP_lst).astype(int)
            GxSP_1D_array = np.sum(VxSP_2D_array, axis = 0)
            GxSP_1D_array = GxSP_1D_array.astype(str)
            # write append to the output file
            output_file.write("\t".join(GxSP_1D_array) + "\n")
            # restart the idx of variant to gene lst and 
            VxSP_lst = []

            # idx_gene
            idx_gene += 1
            print(idx_gene)

            if idx_gene == len(V2G_idx_lst):
                break

            num_var = num_var2gene[idx_gene]

output_file.close()





