from os import makedirs
import random, yaml
import scipy.sparse

import numpy as np
import pandas as pd

from os.path import join, exists, dirname, realpath
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc

from model.builders.SPxG_matrix_input_data_reader import Combined_SPxG_from_variant_type
from model import STAR_NN_model
from model.model_utils import evaluate

# from model.builders.SPxG_matrix_input_data_reader_PTVs_MisAB import Combined_SPxG_from_variant_type
# from model import STAR_NN_model_PTVs_MisAB
# from model.model_utils import evaluate


class STAR_NN_Pipeline():
    def __init__(self, data_params, model_params):

        self.data_params = data_params
        self.model_params = model_params
        self.prepare_saving_dir()

    def run(self):

        print("In STAR_NN_Pipeline.run():")
        # print(self.data_params)
        # print(self.data_params['params'])

        # import data

        data = Combined_SPxG_from_variant_type(
            variant_type_lst=self.data_params['params']['variant_type_lst'], 
            selected_genes=self.data_params['params']['selected_genes'],
            training=True)

        Combined_X, labels, samples_SPID, col_name_w_variant_type = data.get_data()
        print("col_name_w_variant_type: ", col_name_w_variant_type.shape) # 57375, (57375/3 = 19125 genes, each with three category: 'PTVs', 'MisAB', 'MisC')

        # whether use validation dataset
        if self.data_params['params']['validation']:
            iwes1_data = Combined_SPxG_from_variant_type(
                variant_type_lst=self.data_params['params']['variant_type_lst'], 
                selected_genes=self.data_params['params']['selected_genes'],
                training=False)

            Combined_X_iwes1, labels_iwes1, samples_SPID_iwes1, col_name_w_variant_type_iwes1 = iwes1_data.get_data()
            print(col_name_w_variant_type_iwes1.shape) 

            # remove feature columns          
            gene_name_wes12 = [g[0] for g in col_name_w_variant_type]
            uniq_gene_name_wes12 = list(dict.fromkeys(gene_name_wes12))

            gene_name_iwes1 = [g[0] for g in col_name_w_variant_type_iwes1]
            uniq_gene_name_iwes1 = list(dict.fromkeys(gene_name_iwes1))

            comm_feature_names = list(set(uniq_gene_name_wes12).intersection(set(uniq_gene_name_iwes1)))

            X_cleaned_wes12 = Combined_X.reindex(columns=comm_feature_names, level=0) # [40871 rows x 4461 columns]
            X_cleaned_iwes1 = Combined_X_iwes1.reindex(columns=comm_feature_names, level=0) # [68301 rows x 4461 columns]

            Combined_X = X_cleaned_wes12.values
            Combined_X_iwes1 = X_cleaned_iwes1.values

        # loop to train/test/validation
        test_pred_score_lst = list()
        for i in range(10):
            
            # i=1
            print(i)
            # random seed setting

            # split data into train:val:test=8:1:1
            X_train, X_test, y_train, y_test, samples_SPID_train, samples_SPID_test = train_test_split(Combined_X, labels, samples_SPID, test_size=0.1, stratify=labels.astype(str), random_state=1)
            X_train, X_val, y_train, y_val, samples_SPID_train, samples_SPID_val = train_test_split(X_train, y_train, samples_SPID_train, test_size=len(y_test), stratify=y_train.astype(str), random_state=1)
            print(X_train.shape, y_train.shape, samples_SPID_train.shape) # (32695, ) (32695,) (32695,)
            print(X_val.shape, y_val.shape, samples_SPID_val.shape) # (4088, ) (4088,) (4088,)
            print(X_test.shape, y_test.shape, samples_SPID_test.shape) # (4088, ) (4088,) (4088,)

            # get model, fit/train, predict, save prediction
            model = STAR_NN_model.STAR_NN(self.model_params['params'])
            # print("\nCall model/sparse_mlp_model.py ... fit ...")
            # history = model.fit(X_train, y_train, X_val, y_val, X_test, y_test, col_name_w_variant_type=X_cleaned_wes12.columns, random_seed=i)
            history = model.fit(X_train, y_train, X_val, y_val, X_test, y_test, col_name_w_variant_type=col_name_w_variant_type, random_seed=i)

            # predict on test data
            test_pred_scores = model.predict(X_test)[:,0]
            test_pred_bool = (test_pred_scores > 0.5).astype(int)
            metrics_out = evaluate(y_test, test_pred_bool, test_pred_scores)

            if self.data_params['params']['validation']:
                # predict on validation data
                iwes1_pred_scores = model.predict(Combined_X_iwes1)[:,0]
                iwes1_pred_bool = (iwes1_pred_scores > 0.5).astype(int)
                iwes1_metrics_out = evaluate(labels_iwes1, iwes1_pred_bool, iwes1_pred_scores)
                print("iwes_v1:\n", iwes1_metrics_out)

            test_pred_score_lst.append(test_pred_scores)

        test_pred_score_avg = np.mean(test_pred_score_lst, axis=0)

        test_pred_bool = (test_pred_score_avg > 0.5).astype(int)
        print(evaluate(y_test, test_pred_bool, test_pred_score_avg))

        fpr, tpr, _ = roc_curve(y_test, test_pred_score_avg)
        auc_score = auc(fpr, tpr)

        # self.save4roc_plot(y_test, test_pred_score_avg, file_name="/users/qwu24/data/silvio/Qing_Wu/SFARI/batch_jobs/python_script/net/wes12.deepvariant.selFeat1489.auc_" + self.data_params['id'] + ".csv")
        # self.save_prediction(samples_SPID_test, X_test, test_pred_score_avg, y_test, file_name="/users/qwu24/data/silvio/Qing_Wu/SFARI/batch_jobs/python_script/net/wes12.deepvariant.selFeat1489.y_test_pred_" + self.data_params['id'] + ".csv")
        print(self.data_params['id'], "\nauc:", auc_score)

        return fpr, tpr, auc_score

 
    def save_prediction(self, samples_SPID_test, X_test, y_pred_score, y_test, file_name):

        testing_out_df = pd.DataFrame(list(zip(samples_SPID_test, y_test, X_test[:,3], X_test[:,0], y_pred_score)), columns=['SP_ID_index', 'is_case', 'is_female', 'PGS', 'probabilities'])
        testing_out_df.to_csv(file_name, index=False)

    def save4roc_plot(self, y_test, y_pred_score, file_name):
        # input for roc_auc
        fpr, tpr, _ = roc_curve(y_test, y_pred_score)
        auc_score = auc(fpr, tpr)

        auc_out_df = pd.DataFrame(list(zip(fpr, tpr)), columns=['fpr', 'tpr'])
        auc_out_df.to_csv(file_name, index=False)

    def prepare_saving_dir(self):
        self.saving_dir = "../_saving_dir"
        if not exists(self.saving_dir):
            makedirs(self.saving_dir)