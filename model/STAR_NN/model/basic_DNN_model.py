#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow import keras
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc

from model.model_utils import get_X_by_selected_feature, evaluate

import sys
sys.path.insert(0, '/users/qwu24/data/silvio/Qing_Wu/SFARI/batch_jobs/python_script/machine_learning_script/')
from get_SPxG_input_data_exonic_wes1_wes2_combined import generate_input_df_SPARK_wes12

# SFARI gene, DDG2P genes and selected gene feature from top 4% select Percentile()
SFARI_gene_df = pd.read_csv("/users/qwu24/data/silvio/Qing_Wu/SFARI/reference_data/SFARI_gene/SFARI-Gene_genes_01-11-2022release_03-22-2022export.csv")
SFARI_gene_lst = SFARI_gene_df['gene-symbol'].tolist() # 1031

DDG2P_gene_df = pd.read_csv("/users/qwu24/data/silvio/Qing_Wu/SFARI/reference_data/DDG2P/DDG2P_7_4_2022.csv")
DDG2P_gene_lst = DDG2P_gene_df['gene symbol'].tolist() # 2580

selFeat_names_wes12 = pd.read_csv("/users/qwu24/data/silvio/Qing_Wu/SFARI/batch_jobs/python_script/machine_learning_script/tpot_out/wes12.deepvariant.feature_name_union.csv")
selFeat_names_wes12_lst = selFeat_names_wes12['selected_features'].tolist() # 1489 

# region: get SPxG input data file from wes12
GxSP_matrix_file_exonic_wes12 = "/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/pub/WES12/DeepVariant/wes1_wes2_combined.deepvariant.rare1pct_variants_het_by_sample_matrix_cleaned_2GxSP.txt"
V2G_unconsec_lst_file_exonic_wes12 ='/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/pub/WES12/DeepVariant/wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.unconsecutive_indices_splited.txt'
V2G_lst_file_exonic_wes12 ='/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/pub/WES12/DeepVariant/wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.txt'
# endregion

class Basic_DNN():
    def __init__(self, data_params):

        self.data_params = data_params

    def run(self):

        print("In Basic_DNN.run():")

        # import data
        X_exonic_wes12, y_exonic_wes12, _, _ = generate_input_df_SPARK_wes12(
            GxSP_matrix_file=GxSP_matrix_file_exonic_wes12,
            var2gene_lst_file=V2G_lst_file_exonic_wes12, 
            var2gene_unconsec_lst_file=V2G_unconsec_lst_file_exonic_wes12)
        # filter dataset by selected gene features
        X_cleaned_wes12, y_filtered_wes12 = get_X_by_selected_feature(set_name=self.data_params['id'], feature_lst=selFeat_names_wes12_lst, X=X_exonic_wes12, y=y_exonic_wes12, filter_y=False)

        # split data into train:val:test=8:1:1
        _, _, _, _, id_train, id_test = train_test_split(X_cleaned_wes12, y_filtered_wes12, range(len(y_filtered_wes12)), test_size=0.1, random_state=1)
        print("Training/validation on %d samples and testing on %d samples." % (len(id_train), len(id_test)))
        _, _, _, _, id_train, id_val = train_test_split(X_cleaned_wes12.iloc[id_train], y_filtered_wes12.iloc[id_train], id_train, test_size=0.1, random_state=1)
        print("Training on %d samples and validation on %d samples." % (len(id_train), len(id_val)))
        # Training/validation on 36783 samples and testing on 4088 samples.
        # Training on 33104 samples and validation on 3679 samples.


        # model settings
        METRICS = [
            tf.keras.metrics.BinaryAccuracy(name='accuracy'),
            tf.keras.metrics.Precision(name='precision'),
            tf.keras.metrics.Recall(name='recall'),
            tf.keras.metrics.AUC(name='auc'),
            # tf.keras.metrics.AUC(name='prc', curve='PR'), # precision-recall curve
        ]

        early_stopping = tf.keras.callbacks.EarlyStopping(
            monitor='val_loss', 
            verbose=1,
            patience=30,
            mode='min',
            restore_best_weights=True)

        def scheduled_lr_decay(epoch, lr):
            if epoch < 100:
                return lr
            else:
                decay = 1e-5/30
                return float(lr * 1/(1+decay*epoch))

        lr_decay = tf.keras.callbacks.LearningRateScheduler(scheduled_lr_decay, verbose=0)

        # start training and testing
        test_pred_score_lst = list()
        for i in range(10):
            
            # random seed setting
            print(i)
            
            
            tf.keras.backend.clear_session()

            tf.random.set_seed(i)

            # neg, pos = np.bincount(y_filtered_wes12)
            # initial_bias = np.log([pos/neg])

            input = tf.keras.Input(shape = (X_cleaned_wes12.shape[-1],))
            model = tf.keras.layers.Dense(units = 1, activation="sigmoid", use_bias=True, bias_initializer=tf.keras.initializers.GlorotUniform())(input)
            model = tf.keras.Model(inputs=input, outputs=model)

            model.compile(
                optimizer=tf.keras.optimizers.Adam(learning_rate=3*1e-5),
                loss=tf.keras.losses.binary_crossentropy,
                metrics=METRICS)

            EPOCHS = 2000
            BATCH_SIZE = 128

            history = model.fit(X_cleaned_wes12.iloc[id_train], y_filtered_wes12.iloc[id_train],
                                batch_size=BATCH_SIZE,
                                epochs=EPOCHS,
                                callbacks=[early_stopping, lr_decay],
                                validation_data=(X_cleaned_wes12.iloc[id_val], y_filtered_wes12.iloc[id_val]),
                                verbose=0)

            results = model.evaluate(X_cleaned_wes12.iloc[id_test], y_filtered_wes12.iloc[id_test], verbose=0)
            for name, value in zip(model.metrics_names, results):
                print(name, ': ', value)
            print()

            # print metrics
            test_pred_scores = model.predict(X_cleaned_wes12.iloc[id_test]).ravel()
            test_pred_bool = (test_pred_scores > 0.5).astype(int)
            metrics_out = evaluate(y_filtered_wes12.iloc[id_test], test_pred_bool, test_pred_scores)

            test_pred_score_lst.append(test_pred_scores)

        test_pred_score_avg = np.mean(test_pred_score_lst, axis=0)

        fpr, tpr, _ = roc_curve(y_filtered_wes12.iloc[id_test], test_pred_score_avg)
        auc_score = auc(fpr, tpr)
        print(self.data_params['id'], "\nauc:", auc_score)

        return fpr, tpr, auc_score
