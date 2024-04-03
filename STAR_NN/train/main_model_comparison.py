import os, sys
from os.path import join, realpath, basename, dirname
current_dir = dirname(realpath('__file__'))
sys.path.insert(0, dirname(current_dir))

import random
import numpy as np
import pandas as pd
import tensorflow as tf

from importlib.machinery import SourceFileLoader
from STAR_NN_pipeline import STAR_NN_Pipeline
from model.basic_DNN_model import Basic_DNN
from model.subnet_merged_model import SubnetMerged_DNN
from model.traditional_ML_models import ML_models
from model.model_utils import plot_roc_auc

print("\nIn main.py:\ncurrent_dir", current_dir, "\n")

# environment setting
os.environ['KMP_DUPLICATE_LIB_OK'] = 'True'

fpr_lst, tpr_lst, auc_lst, names = list(), list(), list(), list()

params_file_full_path = join(current_dir, "model_params_tanh_selFeat.py")
params = SourceFileLoader(fullname=params_file_full_path, path=params_file_full_path).load_module()
print(params.data['id'])


# region: traditional ML models
ML_model = ML_models(data_params=params.data)

fpr_lst, tpr_lst, auc_lst, names = ML_model.run()
print(pd.DataFrame(list(zip(names, auc_lst))))
# endregion

# region: tanh, selFeat
model_name = "Sparse-Net"
pipeline = STAR_NN_Pipeline(data_params=params.data, model_params=params.model)
print(model_name)

fpr, tpr, auc_score = pipeline.run()
fpr_lst.append(fpr)
tpr_lst.append(tpr)
auc_lst.append(auc_score)
names.append(model_name)
# endregion

# region: basic DNN
model_name = "basic_DNN"
basic_dnn_model = Basic_DNN(data_params=params.data)
print(model_name)

fpr, tpr, auc_score = basic_dnn_model.run()
fpr_lst.append(fpr)
tpr_lst.append(tpr)
auc_lst.append(auc_score)
names.append(model_name)
# endregion

# region: DeepDEP
model_name = "subnet_merged_model"
subnet_merged_model = SubnetMerged_DNN(data_params=params.data)
print(model_name)

fpr, tpr, auc_score = subnet_merged_model.run()
fpr_lst.append(fpr)
tpr_lst.append(tpr)
auc_lst.append(auc_score)
names.append(model_name)
# endregion



# roc-auc
plot_roc_auc(fpr_lst, tpr_lst, auc_lst, names, file_name="/PATH/SFARI/batch_jobs/python_script/net/train/wes12.deepvariant.selFeat.roc_auc.model_comparison_i1.avg.pdf")
