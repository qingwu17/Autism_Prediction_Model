import os, sys
from os.path import join, realpath, basename, dirname
current_dir = dirname(realpath('__file__'))
sys.path.insert(0, dirname(current_dir))

import random
import numpy as np
import tensorflow as tf

from importlib.machinery import SourceFileLoader
from STAR_NN_pipeline import STAR_NN_Pipeline
from model.model_utils import plot_roc_auc

print("\nIn main.py:\ncurrent_dir", current_dir, "\n")

# environment setting
os.environ['KMP_DUPLICATE_LIB_OK'] = 'True'

fpr_lst, tpr_lst, auc_lst, names = list(), list(), list(), list()

# region: tanh, selFeat
params_file_full_path = join(current_dir, "model_params_tanh_selFeat.py")
params = SourceFileLoader(fullname=params_file_full_path, path=params_file_full_path).load_module()
pipeline = STAR_NN_Pipeline(data_params=params.data, model_params=params.model)

print(params.data['id'])

fpr, tpr, auc_score = pipeline.run()
fpr_lst.append(fpr)
tpr_lst.append(tpr)
auc_lst.append(auc_score)
names.append(params.data['id'])
# endregion

# region: tanh, SFARI genes
params_file_full_path = join(current_dir, "model_params_tanh_SFARI.py")
params = SourceFileLoader(fullname=params_file_full_path, path=params_file_full_path).load_module()
pipeline = STAR_NN_Pipeline(data_params=params.data, model_params=params.model)

print(params.data['id'])

fpr, tpr, auc_score = pipeline.run()
fpr_lst.append(fpr)
tpr_lst.append(tpr)
auc_lst.append(auc_score)
names.append(params.data['id'])
# endregion

# region: tanh, selFeat_SFARI
params_file_full_path = join(current_dir, "model_params_tanh_selFeat_SFARI.py")
params = SourceFileLoader(fullname=params_file_full_path, path=params_file_full_path).load_module()
pipeline = STAR_NN_Pipeline(data_params=params.data, model_params=params.model)

print(params.data['id'])

fpr, tpr, auc_score = pipeline.run()
fpr_lst.append(fpr)
tpr_lst.append(tpr)
auc_lst.append(auc_score)
names.append(params.data['id'])
# endregion

# region: tanh, full_gene
params_file_full_path = join(current_dir, "model_params_tanh_full_genes.py")
params = SourceFileLoader(fullname=params_file_full_path, path=params_file_full_path).load_module()
pipeline = STAR_NN_Pipeline(data_params=params.data, model_params=params.model)

print(params.data['id'])

fpr, tpr, auc_score = pipeline.run()
fpr_lst.append(fpr)
tpr_lst.append(tpr)
auc_lst.append(auc_score)
names.append(params.data['id'])
# endregion

plot_roc_auc(fpr_lst, tpr_lst, auc_lst, names, file_name="/PATH/SFARI/batch_jobs/python_script/net/train/wes12.deepvariant.roc_auc.sparse-net.pack4_wval.avg.pdf")
