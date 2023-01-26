import numpy as np
import pandas as pd
from sklearn import metrics
from matplotlib import pyplot as plt

# helper function
def get_X_by_selected_feature(set_name, feature_lst, X, y, filter_y=True):
    
    # SFARI gene, DDG2P genes and selected gene feature from top 4% select Percentile()
    SFARI_gene_df = pd.read_csv("/users/qwu24/data/silvio/Qing_Wu/SFARI/reference_data/SFARI_gene/SFARI-Gene_genes_01-11-2022release_03-22-2022export.csv")
    SFARI_gene_lst = SFARI_gene_df['gene-symbol'].tolist() # 1031

    DDG2P_gene_df = pd.read_csv("/users/qwu24/data/silvio/Qing_Wu/SFARI/reference_data/DDG2P/DDG2P_7_4_2022.csv")
    DDG2P_gene_lst = DDG2P_gene_df['gene symbol'].tolist() # 2580

    if set_name == "selFeat":
        # select genes by select percentile 4
        selected_columns = [g for g in feature_lst if g in X.columns.tolist()] # 765
        X_filtered = X[selected_columns]
    elif set_name == "SFARI_genes":
        # select SFARI genes column
        selected_columns = [g for g in SFARI_gene_lst if g in X.columns.tolist()] # 1019
        selected_columns.insert(0, 'is_female')
        selected_columns.insert(0, "PGS")
        X_filtered = X[selected_columns]
    elif set_name == "DDG2P_genes":
        # select DDG2P genes column
        selected_columns = [g for g in DDG2P_gene_lst if g in X.columns.tolist()] # 2547
        selected_columns.insert(0, 'is_female')
        selected_columns.insert(0, "PGS")
        X_filtered = X[selected_columns]
    elif set_name == "selFeat_SFARI_genes":
        selFeat_SFARI_gene_lst = list(set(feature_lst + SFARI_gene_lst)) # 1717
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
        selected_columns.insert(0, "PGS")
        X_filtered = X[selected_columns]
    elif set_name == "selFeat_SFARI_DDG2P_genes":
        selFeat_SFARI_DDG2P_gene_lst = list(set(feature_lst + SFARI_gene_lst + DDG2P_gene_lst)) # 3424
        # select genes by select percentile 4, SFARI genes and DDG2P genes
        selected_columns = [g for g in selFeat_SFARI_DDG2P_gene_lst  if g in X.columns.tolist()]  
        X_filtered = X[selected_columns]
    elif set_name == "full_gene_lst":
        X_filtered = X
    else:
        raise ValueError("Invalid set name.")

    if filter_y == True:

        # check if all participants have mutations in listed genes
        sum_row = X_filtered.iloc[:,2::].sum(axis=1) 
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

    else:
        y_filtered = y

    return X_filtered, y_filtered


def evaluate(y_test, y_pred, y_pred_score=None):

    accuracy = metrics.accuracy_score(y_true=y_test, y_pred=y_pred)
    f1 = metrics.f1_score(y_true=y_test, y_pred=y_pred)
    precision = metrics.precision_score(y_true=y_test, y_pred=y_pred)
    recall = metrics.recall_score(y_true=y_test, y_pred=y_pred)

    if y_pred_score is None:
        fpr, tpr, thresholds = metrics.roc_curve(y_test, y_pred, pos_label=1)
    else:
        fpr, tpr, thresholds = metrics.roc_curve(y_true=y_test, y_score=y_pred_score, pos_label=1)
    auc = metrics.auc(fpr, tpr)
    
    score = {}
    score['accuracy'] = accuracy
    score['f1'] = f1
    score['precision'] = precision
    score['sensitivity'] = recall
    score['auc'] = auc

    return score

def plot_roc_auc(fpr_lst, tpr_lst, auc_lst, names, file_name="/users/qwu24/data/silvio/Qing_Wu/SFARI/batch_jobs/python_script/net/train/wes12.deepvariant.sparse-net.roc_auc.pdf"):

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
    plt.savefig(file_name)
    plt.close

    # pass