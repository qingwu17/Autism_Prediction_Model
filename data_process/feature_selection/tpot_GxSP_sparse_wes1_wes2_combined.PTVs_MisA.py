#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
   This file read data and run the TPOT classifier,
   and print results to stdout.
"""

import sys
import numpy as np
import pandas as pd

from scipy.sparse import lil_matrix
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from tpot import TPOTClassifier
from tpot.config import classifier_config_dict

import sys
sys.path.insert(0, '/PATH/SFARI/batch_jobs/python_script/machine_learning_script/')
print("sys.path:", sys.path)
from get_SPxG_input_data_exonic_wes1_wes2_combined import generate_input_df_SPARK_wes12, data_frame_to_scipy_sparse_matrix

# adjust TPOTClassifier parameters
NUM_GEN = int(sys.argv[1])
POP_SIZE = int(sys.argv[2])
CV = int(sys.argv[3])
METRIC = sys.argv[4]

output_file_prefix = "/PATH/SFARI/batch_jobs/python_script/machine_learning_script/tpot_out/wes1_wes2_combined.deepvariant.tpot_GxSP_woA_sparse_exonic"
output_file_surfix = "pipeline.PTVs_MisA.py"

output_pipeline_file_name = "_".join(str(i) for i in [output_file_prefix, NUM_GEN, POP_SIZE, CV, METRIC, output_file_surfix])
print(NUM_GEN, POP_SIZE, CV, METRIC, sep=" / ") # 100, 100, 10, "roc_auc"
print(type(NUM_GEN), type(POP_SIZE), type(CV), type(METRIC), sep=" / ")
print(output_pipeline_file_name)

# get SPxG input data, variant annotation matrix and variant to gene index list
VxA_anno_file_exonic = '/PATH/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.var_anno_dummy_df.txt'

GxSP_matrix_file_exonic_PTVs_MisA_wes12 = "/PATH/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het_by_sample_matrix_cleaned_2GxSP.PTVs_MisA.txt"
V2G_unconsec_lst_file_exonic_PTVs_MisA_wes12 ='/PATH/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.unconsecutive_indices_splited.PTVs_MisA.txt'
V2G_lst_file_exonic_PTVs_MisA_wes12 ='/PATH/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.PTVs_MisA.txt'

X_exonic_PTVs_MisA_wes12, y_exonic_PTVs_MisA_wes12, _, _ = generate_input_df_SPARK_wes12(GxSP_matrix_file=GxSP_matrix_file_exonic_PTVs_MisA_wes12, 
                                                                                        var2gene_lst_file=V2G_lst_file_exonic_PTVs_MisA_wes12, var2gene_unconsec_lst_file=V2G_unconsec_lst_file_exonic_PTVs_MisA_wes12,
                                                                                        VxA_anno_file=VxA_anno_file_exonic)
print(X_exonic_PTVs_MisA_wes12.shape, len(y_exonic_PTVs_MisA_wes12)) # PTVs_MisA: (43203, 17484) 43203 


X = X_exonic_PTVs_MisA_wes12
y = y_exonic_PTVs_MisA_wes12

# check if all participants have mutations in listed genes
sum_row = X.iloc[:,1::].sum(axis=1) 
# print((sum_row > 0).value_counts()) # all participants have selected type of mutations in listed genes
# if not all participants have mutations in selected features, remove the participants
row_filter = sum_row > 0
if any(row_filter) > 0:
    X_filtered = X.loc[row_filter,:]
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

# convert to scipy sparse matrix
X_filtered_csr = data_frame_to_scipy_sparse_matrix(X_filtered)
# split training and test data 
X_train, X_test, y_train, y_test = train_test_split(X_filtered_csr, y_filtered, train_size=0.80, test_size=0.20, random_state=17)

print(X_train.shape, y_train.shape) # (34562, 19118)
print(X_test.shape, y_test.shape) # (8641, 19118)

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

    'sklearn.ensemble.RandomForestClassifier': {
        'n_estimators': [100],
        'criterion': ["gini", "entropy"],
        'max_features': np.arange(0.05, 1.01, 0.05),
        'min_samples_split': range(2, 21),
        'min_samples_leaf':  range(1, 21),
        'bootstrap': [True, False]
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
        'tol': [1e-6, 1e-5, 1e-4, 1e-3],
        'C': [1e-4, 1e-3, 1e-2, 1e-1, 0.5, 1., 5., 10., 15., 20., 25.]
    },

    'xgboost.XGBClassifier': {
        'n_estimators': [100],
        'max_depth': range(1, 21),
        'learning_rate': [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.5],
        'subsample': np.arange(0.05, 1.01, 0.05),
        'min_child_weight': range(1, 21),
        'n_jobs': [-1],
        'verbosity': [0]
    }
}

# # adjust TPOTClassifier parameters
# NUM_GEN = int(50)
# POP_SIZE = int(100)
# CV = int(10)
# METRIC = "roc_auc"

# print(NUM_GEN, POP_SIZE, CV, METRIC, sep=" / ") # 50, 100, 10, "roc_auc"
# print(type(NUM_GEN), type(POP_SIZE), type(CV), type(METRIC), sep=" / ")

tpot = TPOTClassifier(generations=NUM_GEN, population_size = POP_SIZE, cv = CV, scoring = METRIC, 
                      n_jobs = -1, config_dict=classifier_config_sparse,
                      verbosity=2, random_state = 17)

tpot.fit(X_train, y_train)
print(tpot.score(X_test, y_test))

preds = tpot.predict(X_test)
print(accuracy_score(y_test, preds))

tpot.export(output_pipeline_file_name)
print(tpot.export())


