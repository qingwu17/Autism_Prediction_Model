
from random import random
import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from sklearn import svm, linear_model
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.svm import LinearSVC
from sklearn.linear_model import LogisticRegression, Lasso, Ridge, ElasticNet, SGDClassifier, RidgeClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.neural_network import MLPClassifier
from xgboost import XGBClassifier

from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, roc_auc_score, roc_curve, precision_recall_curve, average_precision_score, confusion_matrix

import sys
sys.path.insert(0, '/users/qwu24/data/silvio/Qing_Wu/SFARI/batch_jobs/python_script/machine_learning_script/')
from get_SPxG_input_data_exonic_wes1_wes2_combined import generate_input_df_SPARK_wes12
from get_SPxG_input_data_exonic_wes_70487_exome import generate_input_df_SPARK_iwes1
from get_SPxG_input_data_exonic_SSC import generate_input_df_SSC

# SFARI gene, DDG2P genes and selected gene feature from top 4% select Percentile()
SFARI_gene_df = pd.read_csv("/users/qwu24/data/silvio/Qing_Wu/SFARI/reference_data/SFARI_gene/SFARI-Gene_genes_01-11-2022release_03-22-2022export.csv")
SFARI_gene_lst = SFARI_gene_df['gene-symbol'].tolist() # 1031

DDG2P_gene_df = pd.read_csv("/users/qwu24/data/silvio/Qing_Wu/SFARI/reference_data/DDG2P/DDG2P_7_4_2022.csv")
DDG2P_gene_lst = DDG2P_gene_df['gene symbol'].tolist() # 2580

selFeat_names_wes12 = pd.read_table("/users/qwu24/data/silvio/Qing_Wu/SFARI/batch_jobs/python_script/machine_learning_script/tpot_out/wes1_wes2_combined.deepvariant.feature_name_exonic_selPct4.txt")
selFeat_names_wes12_lst = selFeat_names_wes12['selected_features'].tolist() # 765 

selFeat_names_iwes1 = pd.read_table("/users/qwu24/data/silvio/Qing_Wu/SFARI/batch_jobs/python_script/machine_learning_script/tpot_out/wes_70487_exome.deepvariant.feature_name_exonic_selPct6.txt")
# selFeat_names_iwes1 = pd.read_table("/users/qwu24/data/silvio/Qing_Wu/SFARI/batch_jobs/python_script/machine_learning_script/tpot_out/wes_70487_exome.deepvariant.feature_name_exonic_selPct0462.txt")
selFeat_names_iwes1_lst = selFeat_names_iwes1['selected_features'].tolist() # 1153

# get full variant annotation matrix
VxA_anno_file_exonic = '/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.var_anno_dummy_df.txt'

# region: get SPxG input data file from wes12
GxSP_matrix_file_exonic_wes12 = "/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het_by_sample_matrix_cleaned_2GxSP.txt"
V2G_unconsec_lst_file_exonic_wes12 ='/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.unconsecutive_indices_splited.txt'
V2G_lst_file_exonic_wes12 ='/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.txt'

X_exonic_wes12, y_exonic_wes12, _, _ = generate_input_df_SPARK_wes12(GxSP_matrix_file=GxSP_matrix_file_exonic_wes12,
                                                                                 var2gene_lst_file=V2G_lst_file_exonic_wes12, var2gene_unconsec_lst_file=V2G_unconsec_lst_file_exonic_wes12,
                                                                                 VxA_anno_file=VxA_anno_file_exonic)
print(X_exonic_wes12.shape, len(y_exonic_wes12)) # exonic: (43203, 19118) 43203 


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

    y_filtered = y_filtered.astype(int)

    return X_filtered, y_filtered

def create_model(model_name, rs=0):

    if model_name == "Decision Trees":
        model = make_pipeline(DecisionTreeClassifier(random_state=rs))
    elif model_name == "Random Forest":
        model = make_pipeline(RandomForestClassifier(random_state=rs))
    elif model_name == "XGBoost":
        model = XGBClassifier(use_label_encoder =False, random_state=rs, verbosity=0)
    elif model_name == "L1 Logistic Regression":
        model = make_pipeline(LogisticRegression(penalty="l1", solver="saga", max_iter=1000, n_jobs=-1, random_state=rs))
    elif model_name == "L2 Logistic Regression":
        model = make_pipeline(LogisticRegression(penalty="l2", solver='saga', max_iter=1000, n_jobs=-1, random_state=rs))
    elif model_name == "Multilayer Perceptron":
        model = make_pipeline(MLPClassifier(random_state=rs))
    elif model_name == "RBF Support Vector Classifier":
        model = make_pipeline(svm.SVC(C=5.0, kernel='rbf', probability=True, random_state=rs))
    elif model_name == "Linear Support Vector Classifier":
        model = make_pipeline(LinearSVC(C=0.1, penalty="l2", dual=False, random_state=rs))

    return  model


# analysis
X_cleaned, y_filtered = get_X_by_selected_feature(set_name="selFeat", feature_lst=selFeat_names_wes12_lst, X=X_exonic_wes12, y=y_exonic_wes12)
training_X, testing_X, training_y, testing_y = train_test_split(X_cleaned, y_filtered, random_state=17)


# model = make_pipeline(LinearSVC(C=0.1, penalty="l2", dual=False, random_state=17))
model = make_pipeline(LogisticRegression(penalty="l2", solver='liblinear', random_state=17))
model.fit(training_X, training_y)
prediction = model.predict(testing_X)
# probabilities = model.decision_function(testing_X)
probabilities = model.predict_proba(testing_X)[:,1]

#export to csv file
out_df = pd.DataFrame(list(zip(probabilities, testing_X['is_female'], testing_y)))
out_file = "/users/qwu24/data/silvio/Qing_Wu/SFARI/batch_jobs/python_script/machine_learning_script/wes1_wes2_combined.deepvariant.exonic.L2LR_liblinear.prob_sex_on_test_selected_genes.csv"
out_df.to_csv(out_file, index=False)


# test models
model_name_set = ["Decision Trees", "Random Forest", "XGBoost", "L1 Logistic Regression", "L2 Logistic Regression", "Linear Support Vector Classifier"]
names, fpr_lst, tpr_lst, auc_lst, precision_lst, recall_lst, auprc_lst = list(), list(), list(), list(), list(), list(), list()
for mn in model_name_set:

    print("Model Name: ", mn)

    probabilities_mn, acc_mn, f1_mn, precision_mn, recall_mn, auc_mn = list(), list(), list(), list(), list(), list()

    for i in range(20):

        model = create_model(model_name=mn, rs=i)
        model.fit(training_X, training_y)
        prediction = model.predict(testing_X)

        if hasattr(model, "predict_proba"):
            probabilities = model.predict_proba(testing_X)[:,1]
        else:
            probabilities = model.decision_function(testing_X)

        acc_mn.append(accuracy_score(testing_y, prediction)) 
        f1_mn.append(f1_score(testing_y, prediction))
        precision_mn.append(precision_score(testing_y, prediction))
        recall_mn.append(recall_score(testing_y, prediction))
        auc_mn.append(roc_auc_score(testing_y, probabilities))

        probabilities_mn.append(probabilities)

    # metrics
    auc = np.mean(auc_mn)

    print("Accuracy: %.4f | F1: %.4f | Precision: %.4f | Recall: %.4f | ROC-AUC: %.4f " % (np.mean(acc_mn), np.mean(f1_mn), np.mean(precision_mn), np.mean(recall_mn), auc))

    # for plot
    names.append(mn)
    prob = np.mean(probabilities_mn, axis=0)
    
    fpr, tpr, _ = roc_curve(testing_y, prob)
    fpr_lst.append(fpr)
    tpr_lst.append(tpr)
    auc_lst.append(auc)

    auprc = average_precision_score(testing_y, prob)
    auprc_lst.append(auprc)
    precision, recall, _ = precision_recall_curve(testing_y, prob)
    precision_lst.append(precision)
    recall_lst.append(recall)

# add RBF Support Vector Classifier
model_SVM_rbf = make_pipeline(svm.SVC(C=5.0, kernel='rbf', probability=True, random_state=0))
model_SVM_rbf.fit(training_X, training_y)
prediction = model_SVM_rbf.predict(testing_X)
probabilities = model_SVM_rbf.predict_proba(testing_X)[:,1]
names.append("RBF Support Vector Classifier")
fpr, tpr, _ = roc_curve(testing_y, probabilities)
fpr_lst.append(fpr)
tpr_lst.append(tpr)

# region: plot roc_auc
plt.figure()
color = plt.cm.rainbow(np.linspace(0, 1, len(names)))
for i, c in zip(range(len(names)), color):
    plt.plot(
        fpr_lst[i], tpr_lst[i], label="{} ({:.4f})".format(names[i], auc_lst[i]), color=c, lw=1
    )
plt.plot([0, 1], [0, 1], "k--", lw=2)
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Receiver operating characteristic")
plt.legend(loc="lower right")
plt.savefig("/users/qwu24/data/silvio/Qing_Wu/SFARI/batch_jobs/python_script/machine_learning_script/tpot_out/wes1_wes2_combined.deepvariant.exonic.models.roc_auc.pdf")
plt.close



# region: plot precision-recall
plt.figure()
for i in range(len(names)):
    plt.plot(
        recall_lst[0], precision_lst[0], label="{} ({:.4f})".format(names[i], auprc_lst[i]), color="b", lw=1
    )
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision Recall Tradeoff')
plt.legend(loc="lower right")
plt.savefig("/users/qwu24/data/silvio/Qing_Wu/SFARI/batch_jobs/python_script/machine_learning_script/tpot_out/wes1_wes2_combined.deepvariant.exonic.models.precisio_recall.png")
plt.close


# region: test LinearSVC model parameters
def test_LinearSVC_model_par(p="l2", c=1, rs=0):
    model = make_pipeline(LinearSVC(C=c, penalty=p, random_state=rs))
    return model

for p in ["l1", "l2"]:
    print("Penalty: {}".format(p))
    for c in [1e-2, 1e-1, 0.5, 1., 5., 10., 15., 20., 25.]:
        print("C=%.4f" % (c))

        acc_mn, f1_mn, precision_mn, recall_mn, auc_mn = list(), list(), list(), list(), list(),
        for i in range(20):
            model = test_LinearSVC_model_par(p=p, c=c, rs=i)
            model.fit(training_X, training_y)
            prediction = model.predict(testing_X)

            if hasattr(model, "predict_proba"):
                probabilities = model.predict_proba(testing_X)[:,1]
            else:
                probabilities = model.decision_function(testing_X)

            acc_mn.append(accuracy_score(testing_y, prediction)) 
            f1_mn.append(f1_score(testing_y, prediction))
            precision_mn.append(precision_score(testing_y, prediction))
            recall_mn.append(recall_score(testing_y, prediction))
            auc_mn.append(roc_auc_score(testing_y, probabilities))

        print("Accuracy: %.4f | F1: %.4f | Precision: %.4f | Recall: %.4f | ROC-AUC: %.4f " % (np.mean(acc_mn), np.mean(f1_mn), np.mean(precision_mn), np.mean(recall_mn), np.mean(auc_mn)))

# endregion



# region：test Logistic Regression model parameters
def test_LogisticRegression_model_par(C, l1_ratio, rs=0):
    model = make_pipeline(LogisticRegression(C=C, l1_ratio=l1_ratio, penalty="elasticnet", solver="saga", max_iter=1000, n_jobs=-1, random_state=rs))
    return model

for r in (0.3, 0.4):
    
    r = 1.0

    print("l1-ratio: %.1f" % (r))
    for c in [1e-2, 1e-1, 0.5, 1., 5., 10., 15., 20., 25.]:
        print("C=%.4f" % (c))
        probabilities_mn, acc_mn, f1_mn, precision_mn, recall_mn, auc_mn = list(), list(), list(), list(), list(), list()
        for i in range(10):
            model = test_LogisticRegression_model_par(C=c, l1_ratio=r, rs=i)
            model.fit(training_X, training_y)
            prediction = model.predict(testing_X)

            if hasattr(model, "predict_proba"):
                probabilities = model.predict_proba(testing_X)[:,1]
            else:
                probabilities = model.decision_function(testing_X)

            acc_mn.append(accuracy_score(testing_y, prediction)) 
            f1_mn.append(f1_score(testing_y, prediction))
            precision_mn.append(precision_score(testing_y, prediction))
            recall_mn.append(recall_score(testing_y, prediction))
            auc_mn.append(roc_auc_score(testing_y, probabilities))

            probabilities_mn.append(probabilities)

        # metrics
        auc = np.mean(auc_mn)

        print("Accuracy: %.4f | F1: %.4f | Precision: %.4f | Recall: %.4f | ROC-AUC: %.4f " % (np.mean(acc_mn), np.mean(f1_mn), np.mean(precision_mn), np.mean(recall_mn), auc))

# endregion
'''
# region: test SVC model parameters
def test_SVC_model_par(c, rs=0):
    model = make_pipeline(svm.SVC(C=c, kernel="rbf", gamma="scale", random_state=rs))
    return model

for c in [1e-2, 1e-1, 0.5, 1., 5., 10., 15., 20., 25.]:
    print("C=%.4f" % (c))
    probabilities_mn, acc_mn, f1_mn, precision_mn, recall_mn, auc_mn = list(), list(), list(), list(), list(), list()
    for i in range(1):
        # rs only work when probability=True
        model = test_SVC_model_par(c=c, rs=i)
        model.fit(training_X, training_y)
        prediction = model.predict(testing_X)

        if hasattr(model, "predict_proba"):
            probabilities = model.predict_proba(testing_X)[:,1]
        else:
            probabilities = model.decision_function(testing_X)

        acc_mn.append(accuracy_score(testing_y, prediction)) 
        f1_mn.append(f1_score(testing_y, prediction))
        precision_mn.append(precision_score(testing_y, prediction))
        recall_mn.append(recall_score(testing_y, prediction))
        auc_mn.append(roc_auc_score(testing_y, probabilities))

        probabilities_mn.append(probabilities)

    # metrics
    auc = np.mean(auc_mn)

    print("Accuracy: %.4f | F1: %.4f | Precision: %.4f | Recall: %.4f | ROC-AUC: %.4f " % (np.mean(acc_mn), np.mean(f1_mn), np.mean(precision_mn), np.mean(recall_mn), auc))

# endregion


'''
