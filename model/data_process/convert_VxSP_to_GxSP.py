#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
   This file read and write data.
"""

import numpy as np

# read variant to gene index file, each line represent the index of variant to one gene
V2G_idx_file = '/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.unconsecutive_indices_splited.non_impactful.txt'
V2G_idx_lst = []
with open(V2G_idx_file) as tmp_file:
    for line in tmp_file:
        var2gene_idx_w_gene = line.rstrip().split('\t')
        # print(var2gene_idx_w_gene[0])
        var2gene_idx = [int(x) for x in var2gene_idx_w_gene[1:]]
        # print(var2gene_idx)
        V2G_idx_lst.append(var2gene_idx)

len(V2G_idx_lst) 
# 19626 wes1_wes2.deepvariant rare het exonic only (19117 unique)
# 16832 wes1_wes2.deepvariant rare het (PTVs) exonic only (16832 unique)
# 17680 wes1_wes2.deepvariant rare het (PTVs_MisA) exonic only ( unique)
# 18528 wes1_wes2.deepvariant rare het (PTVs_MisAB) exonic only ( unique)
# 3813 wes1_wes2.deepvariant rare het (MisA) exonic only ( unique)
# 9624 wes1_wes2.deepvariant rare het (MisAB) exonic only ( unique)
# 9495 wes1_wes2.deepvariant rare het (MisB) exonic only ( unique)
# 19181 wes1_wes2.deepvariant rare het (non_impactful) exonic only ( unique)

# 19877 wes_70487_exome.deepvariant rare het exonic only
# 18385 wes_70487_exome.deepvariant rare het (PTVs) exonic only
# 18830 wes_70487_exome.deepvariant rare het (PTVs_MisA) exonic only
# 19121 wes_70487_exome.deepvariant rare het (PTVs_MisB) exonic only

# count the number of variant to each gene by line
num_var2gene = [len(x) for x in V2G_idx_lst]

VxSP_matrix_file = "/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het_by_sample_matrix_cleaned.non_impactful.txt"
GxSP_matrix_file = "/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het_by_sample_matrix_cleaned_2GxSP.non_impactful.txt"

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





