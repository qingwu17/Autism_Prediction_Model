#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
   This file read data and run the TPOT classifier,
   and print results to stdout.
"""

import sys
import time
import numpy as np
import pandas as pd

from scipy.sparse import lil_matrix
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from tpot import TPOTClassifier
from tpot.config import classifier_config_dict

# adjust TPOTClassifier parameters
NUM_GEN = int(sys.argv[1])
POP_SIZE = int(sys.argv[2])
CV = int(sys.argv[3])
METRIC = sys.argv[4]

print(NUM_GEN, POP_SIZE, CV, METRIC, sep=" / ") # 100, 100, 10, "roc_auc"
print(type(NUM_GEN), type(POP_SIZE), type(CV), type(METRIC), sep=" / ")

# read SP column header file, columns: sp_id, is_case, is_female
SP_col_header_file = '/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_sample_matrix_spid_case_sex.txt'
SP_col_header_lst = []
with open(SP_col_header_file) as tmp_file:
    for line in tmp_file:
        SP_col_header_lst.append(line.rstrip().split("\t"))
SP_col_header_2D_array = np.array(SP_col_header_lst)
print(SP_col_header_2D_array.shape) # (43225, 3)

# read gene list, read variant to gene index file, each line represent the index of variant to one gene
# the first vector of each index line is the gene name
V2G_idx_file = '/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_idx_var2gene_lst_df_0216.txt'
gene_lst = []
with open(V2G_idx_file) as tmp_file:
    for line in tmp_file:
        var2gene_idx = line.rstrip().split('\t')
        gene_name = var2gene_idx[0]
        gene_lst.append(gene_name)
print(len(gene_lst)) # 19430

# read gene by sp matrix file
# G_SP_matrix_file = '/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_sample_matrix_cleaned_2GxSP.txt'
G_SP_matrix_file = "/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_sample_matrix_cleaned_2GxSP_0216.txt"
G_SP_matrix_lst = []
with open(G_SP_matrix_file) as tmp_file:
    for line in tmp_file:
        G_SP_matrix_lst.append(line.rstrip().split('\t'))
G_SP_matrix = np.array(G_SP_matrix_lst)
print(G_SP_matrix.shape) # (19430, 43227)

# read sp idx file 
idx_sp2exclude_file = '/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_sample_matrix_idx_sp2exclude.txt'
idx_sp2exclude_lst = []
with open(idx_sp2exclude_file) as tmp_file:
    for line in tmp_file:
        idx_sp2exclude_lst.append(int(line.rstrip())-1)
print(idx_sp2exclude_lst) # [20435, 21290]
# exclude selected sp from VxSP matrix
G_SP_matrix_cleaned = np.delete(G_SP_matrix, idx_sp2exclude_lst, axis = 1) 

# transpose VxSP matrix into SPxV matrix
SP_G_matrix = G_SP_matrix_cleaned.T
print(SP_G_matrix.shape) # (43225, 19430)

# # export SPxV matrix
# SP_G_matrix.to_csv('/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/SPARK_vcf_families_merged_rare_variant_sample_matrix_cleaned_2GxSP_T.csv', sep='\t', header = False, index=False)

# paste sp_id, is_case, is_female to SPxV matrix
SP_G_matrix_full = np.hstack((SP_col_header_2D_array, SP_G_matrix))

# change numpy 2D array to pd.DataFrame
SP_G_matrix_full = pd.DataFrame(SP_G_matrix_full)
# add column names to SPxV matrix
col_names = [['SP_ID', 'is_case', 'is_female'], gene_lst]
SP_G_matrix_full.columns = sum(col_names, [])
print(SP_G_matrix_full.shape) # (43225, 19433)

# drop SP_ID column and change all values into int
input_df = SP_G_matrix_full.drop(['SP_ID'], axis=1)
input_df = input_df.astype(float)
print(input_df.dtypes) # int64
print(input_df.shape) # (43225, 19432)

# split into input_df_nolabel (X) and label (Y)
label_column = 'is_case'
label_Y = input_df[label_column]
input_df_nolabel = input_df.drop([label_column], axis=1)

# convert sparse matrix to scipy sparse matrix
BYTES_TO_MB_DIV = 0.000001
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

def print_memory_usage_of_data_frame(df):
    mem = round(df.memory_usage().sum() * BYTES_TO_MB_DIV, 3) 
    print("Memory usage is " + str(mem) + " MB")
def get_csr_memory_usage(matrix):
    mem = (matrix.data.nbytes + matrix.indptr.nbytes + matrix.indices.nbytes) * BYTES_TO_MB_DIV
    print("Memory usage is " + str(mem) + " MB")

input_df_nolabel_csr = data_frame_to_scipy_sparse_matrix(input_df_nolabel)
print_memory_usage_of_data_frame(input_df_nolabel)
get_csr_memory_usage(input_df_nolabel_csr)

X_train, X_test, Y_train, Y_test = train_test_split(input_df_nolabel_csr, label_Y, train_size=0.80, test_size=0.20, random_state=17)

print(X_train.shape) # (34580, 19431)
print(Y_train.shape) # (34580,)
print(X_test.shape) # (8645, 19431)
print(Y_test.shape) # (8645,)

# config
classifier_config_sparse = {
    'tpot.builtins.OneHotEncoder': {
        'minimum_fraction': [0.05, 0.1, 0.15, 0.2, 0.25]
    },

    'sklearn.feature_selection.SelectFwe': {
        'alpha': np.arange(0, 0.05, 0.001),
        'score_func': {
            'sklearn.feature_selection.f_classif': None
        }
    },

    'sklearn.feature_selection.SelectPercentile': {
        'percentile': range(1, 100),
        'score_func': {
            'sklearn.feature_selection.f_classif': None
        }
    },

    'sklearn.feature_selection.VarianceThreshold': {
        'threshold': np.arange(0.05, 1.01, 0.05)
    },

    'sklearn.feature_selection.RFE': {
        'step': np.arange(0.05, 1.01, 0.05),
        'estimator': {
            'sklearn.ensemble.ExtraTreesClassifier': {
                'n_estimators': [100],
                'criterion': ['gini', 'entropy'],
                'max_features': np.arange(0.05, 1.01, 0.05)
            }
        }
    },

    'sklearn.feature_selection.SelectFromModel': {
        'threshold': np.arange(0, 1.01, 0.05),
        'estimator': {
            'sklearn.ensemble.ExtraTreesClassifier': {
                'n_estimators': [100],
                'criterion': ['gini', 'entropy'],
                'max_features': np.arange(0.05, 1.01, 0.05)
            }
        }
    },

    'sklearn.linear_model.LogisticRegression': {
        'penalty': ["l1", "l2"],
        'C': [1e-4, 1e-3, 1e-2, 1e-1, 0.5, 1., 5., 10., 15., 20., 25.],
        'dual': [True, False]
    },

    'sklearn.svm.LinearSVC': {
        'penalty': ["l1", "l2"],
        'loss': ["hinge", "squared_hinge"],
        'dual': [True, False],
        'tol': [1e-5, 1e-4, 1e-3, 1e-2, 1e-1],
        'C': [1e-4, 1e-3, 1e-2, 1e-1, 0.5, 1., 5., 10., 15., 20., 25.]
    },

    'sklearn.ensemble.RandomForestClassifier': {
        'n_estimators': [100],
        'criterion': ["gini", "entropy"],
        'max_features': np.arange(0.05, 1.01, 0.05),
        'min_samples_split': range(2, 21),
        'min_samples_leaf':  range(1, 21),
        'bootstrap': [True, False]
    },

    'xgboost.XGBClassifier': {
        'n_estimators': [100],
        'max_depth': range(1, 11),
        'learning_rate': [1e-3, 1e-2, 1e-1, 0.5, 1.],
        'subsample': np.arange(0.05, 1.01, 0.05),
        'min_child_weight': range(1, 21),
        'n_jobs': [-1],
        'verbosity': [0]
    }
}


tpot = TPOTClassifier(generations=NUM_GEN, population_size = POP_SIZE, cv = CV, scoring = METRIC, 
                      n_jobs = -1, config_dict=classifier_config_sparse,
                      verbosity=2, random_state = 17)

tpot.fit(X_train, Y_train)
print(tpot.score(X_test, Y_test))

preds = tpot.predict(X_test)
print(accuracy_score(Y_test, preds))

tpot.export('/users/qwu24/data/silvio/Qing_Wu/SFARI/batch_jobs/python_script/machine_learning_script/tpot_out/tpot_GxSP_woA_sparse_for_100_100_10_roc_auc_pipeline.py')
print(tpot.export())




