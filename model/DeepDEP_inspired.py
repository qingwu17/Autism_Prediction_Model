#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow import keras
from matplotlib import pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.svm import LinearSVC
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, roc_auc_score, roc_curve, precision_recall_curve, auc, confusion_matrix

import sys
sys.path.insert(0, '/users/qwu24/data/silvio/Qing_Wu/SFARI/batch_jobs/python_script/machine_learning_script/')
from get_SPxG_input_data_exonic_wes1_wes2_combined import generate_input_df_SPARK_wes12
from get_SPxG_input_data_exonic_wes_70487_exome import generate_input_df_SPARK_iwes1, data_frame_to_scipy_sparse_matrix
from get_SPxG_input_data_exonic_SSC import generate_input_df_SSC


# region: input gene name files, including SFARI gene, DDG2P genes and selected gene feature from top 4% select Percentile()
SFARI_gene_df = pd.read_csv("/users/qwu24/data/silvio/Qing_Wu/SFARI/reference_data/SFARI_gene/SFARI-Gene_genes_01-11-2022release_03-22-2022export.csv")
SFARI_gene_lst = SFARI_gene_df['gene-symbol'].tolist() # 1031

DDG2P_gene_df = pd.read_csv("/users/qwu24/data/silvio/Qing_Wu/SFARI/reference_data/DDG2P/DDG2P_7_4_2022.csv")
DDG2P_gene_lst = DDG2P_gene_df['gene symbol'].tolist() # 2580

selFeat_names_wes12 = pd.read_table("/users/qwu24/data/silvio/Qing_Wu/SFARI/batch_jobs/python_script/machine_learning_script/tpot_out/wes1_wes2_combined.deepvariant.feature_name_exonic_selPct4.txt")
selFeat_names_wes12_lst = selFeat_names_wes12['selected_features'].tolist() # 765 

selFeat_names_PTVs_wes12 = pd.read_table("/users/qwu24/data/silvio/Qing_Wu/SFARI/batch_jobs/python_script/machine_learning_script/tpot_out/wes1_wes2_combined.deepvariant.feature_name_exonic_PTVs_selPct007.txt")
selFeat_names_PTVs_wes12_lst = selFeat_names_PTVs_wes12['selected_features'].tolist() # 129 

selFeat_names_PTVs_MisA_wes12 = pd.read_table("/users/qwu24/data/silvio/Qing_Wu/SFARI/batch_jobs/python_script/machine_learning_script/tpot_out/wes1_wes2_combined.deepvariant.feature_name_exonic_PTVs_MisA_selPct1.txt")
selFeat_names_PTVs_MisA_wes12_lst = selFeat_names_PTVs_MisA_wes12['selected_features'].tolist() # 175 

selFeat_names_PTVs_MisAB_wes12 = pd.read_table("/users/qwu24/data/silvio/Qing_Wu/SFARI/batch_jobs/python_script/machine_learning_script/tpot_out/wes1_wes2_combined.deepvariant.feature_name_exonic_PTVs_MisAB_selPct2.txt")
selFeat_names_PTVs_MisAB_wes12_lst = selFeat_names_PTVs_MisAB_wes12['selected_features'].tolist() # 347 

selFeat_names_iwes1 = pd.read_table("/users/qwu24/data/silvio/Qing_Wu/SFARI/batch_jobs/python_script/machine_learning_script/tpot_out/wes_70487_exome.deepvariant.feature_name_exonic_selPct6.txt")
# selFeat_names_iwes1 = pd.read_table("/users/qwu24/data/silvio/Qing_Wu/SFARI/batch_jobs/python_script/machine_learning_script/tpot_out/wes_70487_exome.deepvariant.feature_name_exonic_selPct0462.txt")
selFeat_names_iwes1_lst = selFeat_names_iwes1['selected_features'].tolist() # 1153

selFeat_names_PTVs_iwes1 = pd.read_table("/users/qwu24/data/silvio/Qing_Wu/SFARI/batch_jobs/python_script/machine_learning_script/tpot_out/wes_70487_exome.deepvariant.feature_name_exonic_PTVs_selPct5.txt")
selFeat_names_PTVs_iwes1_lst = selFeat_names_PTVs_iwes1['selected_features'].tolist() # 900

selFeat_names_PTVs_MisA_iwes1 = pd.read_table("/users/qwu24/data/silvio/Qing_Wu/SFARI/batch_jobs/python_script/machine_learning_script/tpot_out/wes_70487_exome.deepvariant.feature_name_exonic_PTVs_MisA_selPct3.txt")
selFeat_names_PTVs_MisA_iwes1_lst = selFeat_names_PTVs_MisA_iwes1['selected_features'].tolist() # 547

# get full variant annotation matrix
VxA_anno_file_exonic = '/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.var_anno_dummy_df.txt'

# endregion

# region: get SPxG input data file from wes12
GxSP_matrix_file_exonic_wes12 = "/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het_by_sample_matrix_cleaned_2GxSP.txt"
V2G_unconsec_lst_file_exonic_wes12 ='/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.unconsecutive_indices_splited.txt'
V2G_lst_file_exonic_wes12 ='/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.txt'

# PTVs
GxSP_matrix_file_exonic_PTVs_wes12 = "/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het_by_sample_matrix_cleaned_2GxSP.PTVs.txt"
V2G_unconsec_lst_file_exonic_PTVs_wes12 ='/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.unconsecutive_indices_splited.PTVs.txt'
V2G_lst_file_exonic_PTVs_wes12 ='/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.PTVs.txt'

# PTVs_MisA
GxSP_matrix_file_exonic_PTVs_MisA_wes12 = "/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het_by_sample_matrix_cleaned_2GxSP.PTVs_MisA.txt"
V2G_unconsec_lst_file_exonic_PTVs_MisA_wes12 ='/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.unconsecutive_indices_splited.PTVs_MisA.txt'
V2G_lst_file_exonic_PTVs_MisA_wes12 ='/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.PTVs_MisA.txt'

# PTVs_MisAB
GxSP_matrix_file_exonic_PTVs_MisAB_wes12 = "/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het_by_sample_matrix_cleaned_2GxSP.PTVs_MisAB.txt"
V2G_unconsec_lst_file_exonic_PTVs_MisAB_wes12 ='/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.unconsecutive_indices_splited.PTVs_MisAB.txt'
V2G_lst_file_exonic_PTVs_MisAB_wes12 ='/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.PTVs_MisAB.txt'

# MisA
GxSP_matrix_file_exonic_MisA_wes12 = "/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het_by_sample_matrix_cleaned_2GxSP.MisA.txt"
V2G_unconsec_lst_file_exonic_MisA_wes12 ='/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.unconsecutive_indices_splited.MisA.txt'
V2G_lst_file_exonic_MisA_wes12 ='/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.MisA.txt'

# MisB
GxSP_matrix_file_exonic_MisB_wes12 = "/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het_by_sample_matrix_cleaned_2GxSP.MisB.txt"
V2G_unconsec_lst_file_exonic_MisB_wes12 ='/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.unconsecutive_indices_splited.MisB.txt'
V2G_lst_file_exonic_MisB_wes12 ='/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.MisB.txt'

# MisAB
GxSP_matrix_file_exonic_MisAB_wes12 = "/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het_by_sample_matrix_cleaned_2GxSP.MisAB.txt"
V2G_unconsec_lst_file_exonic_MisAB_wes12 ='/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.unconsecutive_indices_splited.MisAB.txt'
V2G_lst_file_exonic_MisAB_wes12 ='/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.MisAB.txt'

# nonimpactful
GxSP_matrix_file_exonic_nif_wes12 = "/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het_by_sample_matrix_cleaned_2GxSP.non_impactful.txt"
V2G_unconsec_lst_file_exonic_nif_wes12 ='/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.unconsecutive_indices_splited.non_impactful.txt'
V2G_lst_file_exonic_nif_wes12 ='/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.non_impactful.txt'

# endregion

# region: read input file X,y from wes12

X_exonic_wes12, y_exonic_wes12, _, _ = generate_input_df_SPARK_wes12(GxSP_matrix_file=GxSP_matrix_file_exonic_wes12,
                                                                     var2gene_lst_file=V2G_lst_file_exonic_wes12, var2gene_unconsec_lst_file=V2G_unconsec_lst_file_exonic_wes12,
                                                                     VxA_anno_file=VxA_anno_file_exonic)
print(X_exonic_wes12.shape, len(y_exonic_wes12)) # exonic: (43203, 19118) 43203 

# PTVs
X_exonic_PTVs_wes12, y_exonic_PTVs_wes12, _, _ = generate_input_df_SPARK_wes12(GxSP_matrix_file=GxSP_matrix_file_exonic_PTVs_wes12, 
                                                                               var2gene_lst_file=V2G_lst_file_exonic_PTVs_wes12, var2gene_unconsec_lst_file=V2G_unconsec_lst_file_exonic_PTVs_wes12,
                                                                               VxA_anno_file=VxA_anno_file_exonic)
print(X_exonic_PTVs_wes12.shape, len(y_exonic_PTVs_wes12)) # PTVs: (43203, 16655) 43203 

# # PTVs_MisA
# X_exonic_PTVs_MisA_wes12, y_exonic_PTVs_MisA_wes12, _, _ = generate_input_df_SPARK_wes12(GxSP_matrix_file=GxSP_matrix_file_exonic_PTVs_MisA_wes12, 
#                                                                                          var2gene_lst_file=V2G_lst_file_exonic_PTVs_MisA_wes12, var2gene_unconsec_lst_file=V2G_unconsec_lst_file_exonic_PTVs_MisA_wes12,
#                                                                                          VxA_anno_file=VxA_anno_file_exonic)
# print(X_exonic_PTVs_MisA_wes12.shape, len(y_exonic_PTVs_MisA_wes12)) # PTVs_MisA: (43203, 17484) 43203 

# # PTVs_MisAB
# X_exonic_PTVs_MisAB_wes12, y_exonic_PTVs_MisAB_wes12, _, _ = generate_input_df_SPARK_wes12(GxSP_matrix_file=GxSP_matrix_file_exonic_PTVs_MisAB_wes12, 
#                                                                                            var2gene_lst_file=V2G_lst_file_exonic_PTVs_MisAB_wes12, var2gene_unconsec_lst_file=V2G_unconsec_lst_file_exonic_PTVs_MisAB_wes12,
#                                                                                            VxA_anno_file=VxA_anno_file_exonic)
# print(X_exonic_PTVs_MisAB_wes12.shape, len(y_exonic_PTVs_MisAB_wes12)) # PTVs_MisAB: (43203, 18285) 43203 

# # MisA
# X_exonic_MisA_wes12, y_exonic_MisA_wes12, _, _ = generate_input_df_SPARK_wes12(GxSP_matrix_file=GxSP_matrix_file_exonic_MisA_wes12, 
#                                                                                var2gene_lst_file=V2G_lst_file_exonic_MisA_wes12, var2gene_unconsec_lst_file=V2G_unconsec_lst_file_exonic_MisA_wes12,
#                                                                                VxA_anno_file=VxA_anno_file_exonic)
# print(X_exonic_MisA_wes12.shape, len(y_exonic_MisA_wes12)) # MisA: (43203, 3805) 43203

# # MisB
# X_exonic_MisB_wes12, y_exonic_MisB_wes12, _, _ = generate_input_df_SPARK_wes12(GxSP_matrix_file=GxSP_matrix_file_exonic_MisB_wes12, 
#                                                                                var2gene_lst_file=V2G_lst_file_exonic_MisB_wes12, var2gene_unconsec_lst_file=V2G_unconsec_lst_file_exonic_MisB_wes12,
#                                                                                VxA_anno_file=VxA_anno_file_exonic)
# print(X_exonic_MisB_wes12.shape, len(y_exonic_MisB_wes12)) # MisB: (43203, 9461) 43203

# MisAB
X_exonic_MisAB_wes12, y_exonic_MisAB_wes12, _, _ = generate_input_df_SPARK_wes12(GxSP_matrix_file=GxSP_matrix_file_exonic_MisAB_wes12, 
                                                                                 var2gene_lst_file=V2G_lst_file_exonic_MisAB_wes12, var2gene_unconsec_lst_file=V2G_unconsec_lst_file_exonic_MisAB_wes12,
                                                                                 VxA_anno_file=VxA_anno_file_exonic)
print(X_exonic_MisAB_wes12.shape, len(y_exonic_MisAB_wes12)) # MisAB: (43203, 9582) 43203

# non_impactful
X_exonic_nif_wes12, y_exonic_nif_wes12, _, _ = generate_input_df_SPARK_wes12(GxSP_matrix_file=GxSP_matrix_file_exonic_nif_wes12, 
                                                                                 var2gene_lst_file=V2G_lst_file_exonic_nif_wes12, var2gene_unconsec_lst_file=V2G_unconsec_lst_file_exonic_nif_wes12,
                                                                                 VxA_anno_file=VxA_anno_file_exonic)
print(X_exonic_nif_wes12.shape, len(y_exonic_nif_wes12)) # non_impactful: (43203, 18761) 43203

# endregion

# prepare sample indices for training and testing
_, _, _, _, id_train, id_test = train_test_split(X_exonic_wes12, y_exonic_wes12, range(len(y_exonic_wes12)), test_size=0.1, random_state=17)
print("Training/validation on %d samples and testing on %d samples." % (len(id_train), len(id_test)))
_, _, _, _, id_train, id_val = train_test_split(X_exonic_wes12.iloc[id_train], y_exonic_wes12[id_train], id_train, test_size=0.1, random_state=17)
print("Training on %d samples and validation on %d samples." % (len(id_train), len(id_val)))

# helper function
def get_X_by_selected_feature(set_name, feature_lst, X, y):
    
    if set_name == "selFeat":
        # select genes by select percentile 4
        selected_columns = [g for g in feature_lst if g in X.columns.tolist()] # 765
        X_filtered = X[selected_columns]
    elif set_name == "SFARI_genes":
        # select SFARI genes column
        selected_columns = [g for g in SFARI_gene_lst if g in X.columns.tolist()] # 1019
        selected_columns.insert(0, 'is_female')
        X_filtered = X[selected_columns]
    elif set_name == "DDG2P_genes":
        # select DDG2P genes column
        selected_columns = [g for g in DDG2P_gene_lst if g in X.columns.tolist()] # 2547
        selected_columns.insert(0, 'is_female')
        X_filtered = X[selected_columns]
    elif set_name == "selFeat_SFARI_genes":
        selFeat_SFARI_gene_lst = list(set(feature_lst + SFARI_gene_lst)) # 1717
        # select genes by select percentile 4 and SFARI genes
        selected_columns = [g for g in selFeat_SFARI_gene_lst  if g in X.columns.tolist()] 
        X_filtered = X[selected_columns]
    elif set_name == "selFeat_DDG2P_genes":
        selFeat_DDG2P_gene_lst = list(set(feature_lst + DDG2P_gene_lst)) # 2857
        # select genes by select percentile 4 and DDG2P genes
        selected_columns = [g for g in selFeat_DDG2P_gene_lst  if g in X.columns.tolist()] 
        X_filtered = X[selected_columns]
    elif set_name == "SFARI_DDG2P_genes":
        SFARI_DDG2P_gene_lst = list(set(SFARI_gene_lst + DDG2P_gene_lst)) # 2837
        # select genes by select percentile 4 and DDG2P genes
        selected_columns = [g for g in SFARI_DDG2P_gene_lst  if g in X.columns.tolist()] 
        selected_columns.insert(0, 'is_female')
        X_filtered = X[selected_columns]
    elif set_name == "selFeat_SFARI_DDG2P_genes":
        selFeat_SFARI_DDG2P_gene_lst = list(set(feature_lst + SFARI_gene_lst + DDG2P_gene_lst)) # 3424
        # select genes by select percentile 4, SFARI genes and DDG2P genes
        selected_columns = [g for g in selFeat_SFARI_DDG2P_gene_lst  if g in X.columns.tolist()]  
        X_filtered = X[selected_columns]
    elif set_name == "selFeat_SFARI_DDG2P_genes_int":
        selFeat_SFARI_DDG2P_int_gene_lst = list(set(feature_lst).intersection(set(SFARI_gene_lst), set(DDG2P_gene_lst)))
        # 
        selected_columns = [g for g in selFeat_SFARI_DDG2P_int_gene_lst if g in X.columns.tolist()]
        X_filtered = X[selected_columns]
    else:
        raise ValueError("Invalid set name.")

    # # check if all participants have mutations in listed genes
    # sum_row = X_filtered.iloc[:,1::].sum(axis=1) 
    # # print((sum_row > 0).value_counts()) # all participants have selected type of mutations in listed genes
    # # if not all participants have mutations in selected features, remove the participants
    # row_filter = sum_row > 0
    # if any(row_filter) > 0:
    #     X_filtered = X_filtered.loc[row_filter,:]
    #     y_filtered = y[row_filter]
    # # print(X_filtered.shape, len(y_filtered))
    # # print(y_filtered.value_counts())

    # check if all genes have mutations in certain people among the whole dataset
    sum_col = X_filtered.sum(axis=0)
    # print((sum_col > 0).value_counts()) # all listed genes get mutated in certain participants
    # if not all geen have mutations in all participants, remove the gene columns
    col_filter = sum_col > 0
    if any(col_filter) > 0:
        X_filtered = X_filtered.loc[:,col_filter]
    # print(X_filtered.shape)
    # print(len(X_filtered.columns.tolist()))

    y_filtered = y.astype(int)

    print(X_filtered.shape, len(y_filtered))

    return X_filtered, y_filtered

X_cleaned_wes12, y_filtered_wes12 = get_X_by_selected_feature(set_name="selFeat_SFARI_DDG2P_genes", feature_lst=selFeat_names_wes12_lst, X=X_exonic_wes12, y=y_exonic_wes12)
X_cleaned_PTVs_wes12, y_filtered_PTVs_wes12 = get_X_by_selected_feature(set_name="selFeat_SFARI_DDG2P_genes", feature_lst=selFeat_names_wes12_lst, X=X_exonic_PTVs_wes12, y=y_exonic_PTVs_wes12)
# X_cleaned_PTVs_MisA_wes12, y_filtered_PTVs_MisA_wes12 = get_X_by_selected_feature(set_name="selFeat_SFARI_DDG2P_genes", feature_lst=selFeat_names_wes12_lst, X=X_exonic_PTVs_MisA_wes12, y=y_exonic_PTVs_MisA_wes12)
# X_cleaned_PTVs_MisAB_wes12, y_filtered_PTVs_MisAB_wes12 = get_X_by_selected_feature(set_name="selFeat_SFARI_DDG2P_genes", feature_lst=selFeat_names_wes12_lst, X=X_exonic_PTVs_MisAB_wes12, y=y_exonic_PTVs_MisAB_wes12)
# X_cleaned_MisA_wes12, y_filtered_MisA_wes12 = get_X_by_selected_feature(set_name="selFeat_SFARI_DDG2P_genes", feature_lst=selFeat_names_wes12_lst, X=X_exonic_MisA_wes12, y=y_exonic_MisA_wes12)
# X_cleaned_MisB_wes12, y_filtered_MisB_wes12 = get_X_by_selected_feature(set_name="selFeat_SFARI_DDG2P_genes", feature_lst=selFeat_names_wes12_lst, X=X_exonic_MisB_wes12, y=y_exonic_MisB_wes12)
X_cleaned_MisAB_wes12, y_filtered_MisAB_wes12 = get_X_by_selected_feature(set_name="selFeat_SFARI_DDG2P_genes", feature_lst=selFeat_names_wes12_lst, X=X_exonic_MisAB_wes12, y=y_exonic_MisAB_wes12)
X_cleaned_nif_wes12, y_filtered_nif_wes12 = get_X_by_selected_feature(set_name="selFeat_SFARI_DDG2P_genes", feature_lst=selFeat_names_wes12_lst, X=X_exonic_nif_wes12, y=y_exonic_nif_wes12)


print(X_cleaned_PTVs_wes12.shape, X_cleaned_MisAB_wes12.shape, X_cleaned_nif_wes12.shape, len(y_filtered_PTVs_wes12), len(y_filtered_MisAB_wes12), len(y_filtered_wes12), len(y_filtered_nif_wes12))


# region: make model using multiple source of inputs

METRICS = [
    tf.keras.metrics.TruePositives(name='tp'),
    tf.keras.metrics.TrueNegatives(name='tn'),
    tf.keras.metrics.FalsePositives(name='fp'),
    tf.keras.metrics.FalseNegatives(name='fn'), 
    tf.keras.metrics.BinaryAccuracy(name='accuracy'),
    tf.keras.metrics.Precision(name='precision'),
    tf.keras.metrics.Recall(name='recall'),
    tf.keras.metrics.AUC(name='auc'),
    tf.keras.metrics.AUC(name='prc', curve='PR'), # precision-recall curve
]

early_stopping = tf.keras.callbacks.EarlyStopping(
    monitor='val_loss', 
    verbose=1,
    patience=20,
    mode='min',
    restore_best_weights=True)

def scheduled_lr_decay(epoch, lr):
    if epoch < 100:
        return lr
    else:
        decay = 1e-4/30
        return float(lr * 1/(1+decay*epoch))

lr_decay = tf.keras.callbacks.LearningRateScheduler(scheduled_lr_decay, verbose=0)

for i in range(10):
    print(i)
    # subnetwork of PTVs
    input_PTVs = tf.keras.Input(shape = (X_cleaned_PTVs_wes12.shape[-1],))
    model_PTVs = tf.keras.layers.Dense(units = 1, activation="sigmoid", use_bias=True, bias_initializer=tf.keras.initializers.GlorotUniform())(input_PTVs)
    # model_PTVs = tf.keras.layers.Dropout(0.5)(model_PTVs)
    # model_PTVs = tf.keras.layers.Dense(units = 128, activation="relu", use_bias=True, bias_initializer=tf.keras.initializers.GlorotUniform())(model_PTVs)
    model_PTVs = tf.keras.Model(inputs=input_PTVs, outputs=model_PTVs)
    # model_PTVs.summary()

    # subnetwork of MisAB
    input_MisAB = tf.keras.Input(shape = (X_cleaned_MisAB_wes12.shape[-1],))
    model_MisAB = tf.keras.layers.Dense(units = 1, activation="sigmoid", use_bias=True, bias_initializer=tf.keras.initializers.GlorotUniform())(input_MisAB)
    # model_MisAB = tf.keras.layers.Dropout(0.5)(model_MisAB)
    # model_MisAB = tf.keras.layers.Dense(units = 128, activation="relu", use_bias=True, bias_initializer=tf.keras.initializers.GlorotUniform())(model_MisAB)
    model_MisAB = tf.keras.Model(inputs=input_MisAB, outputs=model_MisAB)
    # model_MisAB.summary()

    # subnetwork of non_impactful variants
    input_nif = tf.keras.Input(shape = (X_cleaned_nif_wes12.shape[-1],))
    model_nif = tf.keras.layers.Dense(units = 1, activation="sigmoid", use_bias=True, bias_initializer=tf.keras.initializers.GlorotUniform())(input_nif)
    # model_nif = tf.keras.layers.Dropout(0.5)(model_nif)
    # model_nif = tf.keras.layers.Dense(units = 128, activation="relu", use_bias=True, bias_initializer=tf.keras.initializers.GlorotUniform())(model_nif)
    model_nif = tf.keras.Model(inputs=input_nif, outputs=model_nif)
    # model_nif.summary()

    # merged PTVs, MisAB, non_impactful
    # combined = tf.keras.layers.concatenate([model_PTVs.output, model_MisAB.output, model_nif.output])
    combined = tf.keras.layers.concatenate([model_nif.output])
    # out = tf.keras.layers.Dense(512, activation="relu", use_bias=True, bias_initializer=tf.keras.initializers.GlorotUniform())(combined)
    # out = tf.keras.layers.Dropout(0.5)(out)
    # out = tf.keras.layers.Dense(128, activation="relu", use_bias=True, bias_initializer=tf.keras.initializers.GlorotUniform())(combined)
    # out = tf.keras.layers.Dropout(0.2)(out)
    out = tf.keras.layers.Dense(1, activation="sigmoid", use_bias=True, bias_initializer=tf.keras.initializers.GlorotUniform())(combined)

    model = tf.keras.Model(inputs=[model_PTVs.input, model_MisAB.input, model_nif.input], outputs=out)

    model.compile(
        optimizer=tf.keras.optimizers.Adam(learning_rate=1e-4),
        loss=tf.keras.losses.binary_crossentropy,
        metrics=METRICS)

    EPOCHS = 2000
    BATCH_SIZE = 128

    history = model.fit([X_cleaned_PTVs_wes12.iloc[id_train], X_cleaned_MisAB_wes12.iloc[id_train], X_cleaned_nif_wes12.iloc[id_train]], y_filtered_wes12[id_train],
                        batch_size=BATCH_SIZE,
                        epochs=EPOCHS,
                        callbacks=[early_stopping, lr_decay],
                        validation_data=([X_cleaned_PTVs_wes12.iloc[id_val], X_cleaned_MisAB_wes12.iloc[id_val], X_cleaned_nif_wes12.iloc[id_val]], y_filtered_wes12[id_val]),
                        verbose=0)

    results = model.evaluate([X_cleaned_PTVs_wes12.iloc[id_test], X_cleaned_MisAB_wes12.iloc[id_test], X_cleaned_nif_wes12.iloc[id_test]], y_filtered_wes12[id_test], verbose=0)
    for name, value in zip(model.metrics_names, results):
        print(name, ': ', value)
    print()

# endregion

# region: plot roc_auc
plt.figure()
lw=2
plt.plot(
    fpr_lst[0], tpr_lst[0], label="Selected_features (%.4f)" % roc_auc_lst[0], color="b", lw=lw
)
plt.plot(
    fpr_lst[1], tpr_lst[1], label="SFARI genes (%.4f)" % roc_auc_lst[1], color="g", linestyle=":", lw=lw
)
plt.plot(
    fpr_lst[2], tpr_lst[2], label="DDG2P genes (%.4f)" % roc_auc_lst[2], color="r", linestyle=":", lw=lw
)
plt.plot([0, 1], [0, 1], "k--", lw=2)
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Receiver operating characteristic")
plt.legend(loc="lower right")
plt.savefig("/users/qwu24/data/silvio/Qing_Wu/SFARI/batch_jobs/python_script/machine_learning_script/tpot_out/wes1_wes2_combined.deepvariant.exonic.DeepDEP.roc_auc.png")
plt.close

# endregion

# region: plot precision-recall
plt.figure()
lw=2
plt.plot(
    recall_lst[0], precision_lst[0], label="Exonic", color="b", lw=lw
)
plt.plot(
    recall_lst[29], precision_lst[29], label="PTVs", color="g", lw=lw
)
plt.plot(
    recall_lst[30], precision_lst[30], label="PTVs_MisA", color="r", linestyle=":", lw=lw
)
plt.plot(
    recall_lst[31], precision_lst[31], label="PTVs_MisAB", color="y", linestyle=":", lw=lw
)
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision Recall Tradeoff')
plt.legend(loc="lower right")
plt.savefig("/users/qwu24/data/silvio/Qing_Wu/SFARI/batch_jobs/python_script/machine_learning_script/tpot_out/wes1_wes2_combined.deepvariant.exonic.DeepDEP.precisio_recall.png")
plt.close

# endregion

