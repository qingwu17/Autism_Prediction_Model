import os, sys
from os.path import join, realpath, basename, dirname
current_dir = dirname(realpath('__file__'))
sys.path.insert(0, dirname(current_dir))

import random
import numpy as np
import tensorflow as tf

from importlib.machinery import SourceFileLoader
from STAR_NN_pipeline import STAR_NN_Pipeline

print("\nIn main.py:\ncurrent_dir", current_dir, "\n")

# environment setting
os.environ['KMP_DUPLICATE_LIB_OK'] = 'True'

# call parameter file
params_file_full_path = join(current_dir, "model_params_tanh_selFeat.py")
params = SourceFileLoader(fullname=params_file_full_path, path=params_file_full_path).load_module()

# ==> to ./pipeline_train_validate.py
pipeline = STAR_NN_Pipeline(data_params=params.data, model_params=params.model)

pipeline.run()