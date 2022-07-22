#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
   This file read data and run the TPOT default classifier,
   and print results to stdout.
"""

import numpy as np
import pandas as pd

from scipy.sparse import lil_matrix
from sklearn.model_selection import train_test_split


def get_SPxis_case_is_female_header_2D_arr(SP_col_header_file):
    # read SP column header file, columns: sp_id, is_case, is_female
    SP_col_header_lst = []
    with open(SP_col_header_file) as tmp_file:
        for line in tmp_file:
            SP_col_header_lst.append(line.rstrip().split("\t"))
    SP_col_header_2D_array = np.array(SP_col_header_lst)
    SP_col_header_df = pd.DataFrame(SP_col_header_2D_array)
    # add column names to SPxV matrix
    col_names = ['SP_ID', 'is_case', 'is_female']
    SP_col_header_df.columns = col_names
    print("Shape of SP_id_header numpy array:", SP_col_header_2D_array.shape) # (70464, 3)
    return SP_col_header_2D_array

def get_GxSP_2D_arr(GxSP_matrix_file):
    # read gene by sp matrix file
    GxSP_matrix_file = GxSP_matrix_file
    GxSP_matrix_lst = []
    with open(GxSP_matrix_file) as tmp_file:
        for line in tmp_file:
            GxSP_matrix_lst.append(line.rstrip().split('\t'))
    GxSP_matrix = np.array(GxSP_matrix_lst)
    print("Shape of GxSP matrix numpy array:", GxSP_matrix.shape) # (19877, 70464)
    return GxSP_matrix

def get_gene_lst_of_GxSP(var2gene_lst_file):
    # read gene list, read variant to gene index file, each line represent the index of variant to one gene
    # the first vector of each index line is the gene name
    V2G_idx_file = var2gene_lst_file
    gene_lst = []
    V2G_idx_lst = []
    with open(V2G_idx_file) as tmp_file:
        for line in tmp_file:
            var2gene_idx_w_gene = line.rstrip().split('\t')
            gene_name = var2gene_idx_w_gene[0]
            gene_lst.append(gene_name)
            # get the idx of variant to the gene of that row
            var2gene_idx = [int(x)-1 for x in var2gene_idx_w_gene[1:]]
            # var2gene_idx = var2gene_idx[1:]
            V2G_idx_lst.append(var2gene_idx)
    print("Length of gene converted by rare exonic variants:", len(gene_lst)) # 18825
    return gene_lst, V2G_idx_lst

def get_exclude_sp_id_idx(idx_sp2exclude_file):
    # read sp idx file 
    idx_sp2exclude_lst = []
    with open(idx_sp2exclude_file) as tmp_file:
        for line in tmp_file:
            idx_sp2exclude_lst.append(int(line.rstrip())-1)
    print("Length of participants with missing phenotype information to exclude:", len(idx_sp2exclude_lst)) # 640
    return idx_sp2exclude_lst

def get_exonic_variant_annotation_dummy_df(VxA_anno_file):
    # read variant annotation file, each line represent dummy variables from 9 annotation tools of a variant
    VxA_file = VxA_anno_file
    VxA_lst = []
    with open(VxA_file) as tmp_file:
        for line in tmp_file:
            VxA_lst.append(line.rstrip().split('\t'))
    VxA_2D_array = np.array(VxA_lst)
    VxA_df = pd.DataFrame(VxA_2D_array)
    VxA_df = VxA_df.iloc[:, 6:]
    VxA_df = VxA_df.astype(np.float32)
    print("Shape of exonic variants annotation matrix:", VxA_df.shape) # (274973, 21)
    return VxA_df

def generate_input_df_SPARK_iwes1(GxSP_matrix_file, var2gene_lst_file, var2gene_unconsec_lst_file, VxA_anno_file, subset=False):
    
    print("in generate_input_df()")

    SPxis_case_is_female_header_matrix = get_SPxis_case_is_female_header_2D_arr(SP_col_header_file='/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes_70487_exome.deepvariant.rare1pct_variants_het_by_sample_matrix_cleaned.sp_sex_case_70464.txt') # (70464, 3)
    GxSP_matrix = get_GxSP_2D_arr(GxSP_matrix_file) # (19877, 70464)
    
    if subset == True:
        idx_sp2exclude_lst = get_exclude_sp_id_idx(idx_sp2exclude_file='/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes_70487_exome.deepvariant.rare1pct_variants_het_by_sample_matrix_idx_sp2exclude.txt') # 640
        # exclude selected sp from VxSP matrix
        GxSP_matrix_cleaned = np.delete(GxSP_matrix, idx_sp2exclude_lst, axis=1) 
        SPxis_case_is_female_header_matrix_cleaned = np.delete(SPxis_case_is_female_header_matrix, idx_sp2exclude_lst, axis=0)
    else:
        GxSP_matrix_cleaned = GxSP_matrix
        SPxis_case_is_female_header_matrix_cleaned = SPxis_case_is_female_header_matrix

    # transpose VxSP matrix into SPxV matrix
    SPxG_matrix = GxSP_matrix_cleaned.T
    print("Shape of SPxG matrix numpy array:", SPxG_matrix.shape) # (70464, 19877)

    # paste sp_id, is_case, is_female to SPxV matrix
    SPxG_matrix_full = np.hstack((SPxis_case_is_female_header_matrix_cleaned, SPxG_matrix))
    # change numpy 2D array to pd.DataFrame
    SPxG_matrix_full = pd.DataFrame(SPxG_matrix_full)
    # add column names to SPxV matrix
    gene_lst_w_unconsecutive_indices_splited, _ = get_gene_lst_of_GxSP(var2gene_lst_file=var2gene_unconsec_lst_file) # 19877
    col_names = [['SP_ID', 'is_case', 'is_female'], gene_lst_w_unconsecutive_indices_splited]
    SPxG_matrix_full.columns = sum(col_names, [])
    print("Shape of SPxG with header column:", SPxG_matrix_full.shape) # (70464, 19880)

    # combined same gene name columns into one
    new_SPxG_matrix_full = SPxG_matrix_full.groupby(level=0, axis=1).sum()
    # print(new_SPxG_matrix_full) 

    # order columns by original gene name list
    gene_lst, VxG_idx_lst = get_gene_lst_of_GxSP(var2gene_lst_file=var2gene_lst_file) # 19233
    gene_lst.insert(0, 'is_female')
    gene_lst.insert(0, 'is_case')
    gene_lst.insert(0, 'SP_ID')
    SPxG_matrix_full = new_SPxG_matrix_full[gene_lst]
    # print(SPxG_matrix_full) 

    '''
    # exclude selected sp from risk score data frame
    SPxrisk_score_header_df = get_SPxrisk_score_header_df() # (8906, 2)
    # exclude selected sp from SPxrisk_score df
    SPxrisk_score_header_df_cleaned = SPxrisk_score_header_df.drop(idx_sp2exclude_lst)
    # print(SPxrisk_score_header_df_cleaned.shape) # (8266, 2)

    # paste sp_id_cluster_membership to SPxV matrix
    SPxG_matrix_full.insert(loc=2, column="risk_score_comm_vars", value=SPxrisk_score_header_df_cleaned.loc[:,1])
    # print(SPxG_matrix_full) 
    '''

    # prepare input file
    # drop SP_ID column and change all values into float
    input_df = SPxG_matrix_full.drop(['SP_ID'], axis=1)
    input_df = input_df.astype(np.float32)
    # print(input_df.dtypes) # float32
    print("Shape of input data wo SP_id:", input_df.shape) # (70464, 19235)

    # split into input_df_nolabel (X) and label (Y)
    label_column = 'is_case'
    label_Y = input_df[label_column]
    input_df_nolabel = input_df.drop([label_column], axis=1)

    # read variants function annotation file
    VxA_dummy_df = get_exonic_variant_annotation_dummy_df(VxA_anno_file) # (274973, 21)

    return input_df_nolabel, label_Y, VxA_dummy_df, VxG_idx_lst

# helper function
def get_input_df_by_selected_feature(selected_gene_lst, X, y, V2G_idx_lst):

    # select genes by select percentile 4
    selected_col_names = [g for g in selected_gene_lst if g in X.columns.tolist()] # 765
    selected_col_names.insert(0, 'is_female')
    print(len(selected_col_names))
    X_filtered = X[selected_col_names]

    # check if all participants have mutations in listed genes
    sum_row = X_filtered.iloc[:,1::].sum(axis=1) 
    # print((sum_row > 0).value_counts()) # all participants have selected type of mutations in listed genes
    # if not all participants have mutations in selected features, remove the participants
    row_filter = sum_row > 0
    if any(row_filter) > 0:
        X_filtered = X_filtered.loc[row_filter,:]
        y_filtered = y[row_filter]
    # print(X_filtered.shape, len(y_filtered))
    # print(y_filtered.value_counts())

    # check if all genes have mutations in certain people among the whole dataset
    sum_col = X_filtered.sum(axis=0)
    # print((sum_col > 0).value_counts()) # all listed genes get mutated in certain participants
    # if not all geen have mutations in all participants, remove the gene columns
    col_filter = sum_col > 0
    if any(col_filter) > 0:
        X_filtered = X_filtered.loc[:,col_filter]
    # print(X_filtered.shape)
    # print(len(X_filtered.columns.tolist()))

    from itertools import compress
    V2G_idx_lst_filtered = list(compress(V2G_idx_lst, np.isin(X.columns, selected_col_names)[1:]))

    return X_filtered, y_filtered, V2G_idx_lst_filtered
    

# region: convert sparse matrix to scipy sparse matrix
def data_frame_to_scipy_sparse_matrix(df):
    """
    Converts a sparse pandas data frame to sparse scipy csr_matrix.
    :param df: pandas data frame
    :return: csr_matrix
    """
    arr = lil_matrix(df.shape, dtype=np.float32)
    for i, col in enumerate(df.columns):
        ix = df[col] != 0
        arr[np.where(ix), i] = 1

    return arr.tocsr()

def print_memory_usage_of_data_frame(df, BYTES_TO_MB_DIV = 0.000001):
    mem = round(df.memory_usage().sum() * BYTES_TO_MB_DIV, 3) 
    print("Memory usage is " + str(mem) + " MB")

def print_memory_usage_of_csr(matrix, BYTES_TO_MB_DIV = 0.000001):
    mem = (matrix.data.nbytes + matrix.indptr.nbytes + matrix.indices.nbytes) * BYTES_TO_MB_DIV
    print("Memory usage is " + str(mem) + " MB")

# endregion

def main():

    GxSP_matrix_file_exonic = "/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes_70487_exome.deepvariant.rare1pct_variants_het_by_sample_matrix_cleaned_2GxSP.txt"
    var2gene_unconsec_lst_file_exonic ='/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes_70487_exome.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.unconsecutive_indices_splited.txt'
    var2gene_lst_file_exonic ='/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes_70487_exome.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.txt'
    VxA_anno_file_exonic = '/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.var_anno_dummy_df.txt'

    X, y, VxA_dummy_df, V2G_idx_lst = generate_input_df_SPARK_iwes1(GxSP_matrix_file=GxSP_matrix_file_exonic, 
                                                        var2gene_lst_file=var2gene_lst_file_exonic, var2gene_unconsec_lst_file=var2gene_unconsec_lst_file_exonic,
                                                        VxA_anno_file=VxA_anno_file_exonic)
    print(X.shape, len(y), VxA_dummy_df.shape, len(V2G_idx_lst))
    # (70464, 19234) 70464 (1469036, 21) 19233

    # region: select columns by feature name
    # SFARI gene, DDG2P genes and selected gene feature from top 4% select Percentile()
    SFARI_gene_df = pd.read_csv("/users/qwu24/data/silvio/Qing_Wu/SFARI/reference_data/SFARI_gene/SFARI-Gene_genes_01-11-2022release_03-22-2022export.csv")
    SFARI_gene_lst = SFARI_gene_df['gene-symbol'].tolist() # 1031

    DDG2P_gene_df = pd.read_csv("/users/qwu24/data/silvio/Qing_Wu/SFARI/reference_data/DDG2P/DDG2P_7_4_2022.csv")
    DDG2P_gene_lst = DDG2P_gene_df['gene symbol'].tolist() # 2580

    selected_feature_names = pd.read_table("/users/qwu24/data/silvio/Qing_Wu/SFARI/batch_jobs/python_script/machine_learning_script/tpot_out/wes1_wes2_combined.deepvariant.feature_name_exonic_selPct4.txt")
    selected_feature_names_lst = selected_feature_names['selected_features'].tolist() # 765

    X_filtered, y_filtered, V2G_idx_lst_filtered = get_input_df_by_selected_feature(selected_feature_names_lst[1:], X, y, V2G_idx_lst)
    print(X_filtered.shape, len(y), len(V2G_idx_lst_filtered))


    # seperate the dataset into polygenic risk score with sex and sex with mutation matrix from rare variants
    SPxG_X_nolabel = X.iloc[:,1::]
    print(SPxG_X_nolabel) 
    SPxis_female_risk_score_X_nolabel = X.iloc[:,[0,1]]
    print(SPxis_female_risk_score_X_nolabel)

    SPxG_X_nolabel_csr = data_frame_to_scipy_sparse_matrix(X)
    print_memory_usage_of_data_frame(X)
    print_memory_usage_of_csr(SPxG_X_nolabel_csr)


    # X_train, X_test, Y_train, Y_test = train_test_split(input_df_nolabel, label_Y, train_size=0.80, test_size=0.20, random_state=17)
    

    # print(X_train.shape) # (34580, 2)
    # print(Y_train.shape) # (34580,)
    # print(X_test.shape) # (8645, 2)
    # print(Y_test.shape) # (8645,)



if __name__=="__main__":
    main()

