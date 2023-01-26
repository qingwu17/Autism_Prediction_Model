import numpy as np
import pandas as pd
from random import random
from matplotlib import pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC, LinearSVC
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from xgboost import XGBClassifier
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, roc_auc_score, roc_curve, auc, precision_recall_curve, average_precision_score, confusion_matrix

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

class ML_models():
    def __init__(self, data_params):

        self.data_params = data_params

    def run(self):

        print("In traditional_ML_models.run():")

        # import data
        X_exonic_wes12, y_exonic_wes12, _, _ = generate_input_df_SPARK_wes12(
            GxSP_matrix_file=GxSP_matrix_file_exonic_wes12,
            var2gene_lst_file=V2G_lst_file_exonic_wes12, 
            var2gene_unconsec_lst_file=V2G_unconsec_lst_file_exonic_wes12)
        # filter dataset by selected gene features
        X_cleaned_wes12, y_filtered_wes12 = get_X_by_selected_feature(set_name=self.data_params['id'], feature_lst=selFeat_names_wes12_lst, X=X_exonic_wes12, y=y_exonic_wes12, filter_y=False)
        
        # split data for training and testing
        training_X, testing_X, training_y, testing_y = train_test_split(X_cleaned_wes12, y_filtered_wes12, random_state=1)

        ML_model_name_set = ["Decision Trees", "Random Forest", "XGBoost", "L1 Logistic Regression", "L2 Logistic Regression", "Linear Support Vector Classifier"]
        
        fpr_lst, tpr_lst, auc_lst, names = list(), list(), list(), list()
        for mn in ML_model_name_set:

            test_pred_score_lst = list()
            for i in range(10):

                if mn == "Decision Trees":
                    model = make_pipeline(DecisionTreeClassifier(random_state=i))
                elif mn == "Random Forest":
                    model = make_pipeline(RandomForestClassifier(random_state=i))
                elif mn == "XGBoost":
                    model = XGBClassifier(use_label_encoder =False, random_state=i, verbosity=0)
                elif mn == "L1 Logistic Regression":
                    model = make_pipeline(LogisticRegression(C=0.001, penalty="l1", solver="saga", n_jobs=-1, max_iter=2000, random_state=i))
                elif mn == "L2 Logistic Regression":
                    model = make_pipeline(LogisticRegression(C=0.001, penalty="l2", solver='saga', n_jobs=-1, max_iter=2000, random_state=i))
                elif mn == "Multilayer Perceptron":
                    model = make_pipeline(MLPClassifier(random_state=i))
                elif mn == "RBF Support Vector Classifier":
                    model = make_pipeline(SVC(C=0.001, kernel='rbf', probability=True, random_state=i))
                elif mn == "Linear Support Vector Classifier":
                    model = make_pipeline(LinearSVC(C=0.001, penalty="l2", dual=False, random_state=i))

                model.fit(training_X, training_y)
                prediction = model.predict(testing_X)

                if hasattr(model, "predict_proba"):
                    probabilities = model.predict_proba(testing_X)[:,1]
                else:
                    probabilities = model.decision_function(testing_X)

                test_pred_score_lst.append(probabilities)

            # metrics
            test_pred_score_avg = np.mean(test_pred_score_lst, axis=0)
            test_pred_bool = (test_pred_score_avg > 0.5).astype(int)

            fpr, tpr, _ = roc_curve(testing_y, test_pred_score_avg)
            auc_score = auc(fpr, tpr)

            print("Model name: ", mn)
            print(self.data_params['id'], " auc:", auc_score)
            print(evaluate(testing_y, test_pred_bool, test_pred_score_avg))

            fpr_lst.append(fpr)
            tpr_lst.append(tpr)
            auc_lst.append(auc_score)
            names.append(mn)



        return fpr_lst, tpr_lst, auc_lst, names
