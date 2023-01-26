import os, sys
import datetime

import numpy as np
import pandas as pd
import tensorflow as tf

import keras

from keras import backend as K
from keras import regularizers, activations, initializers, constraints
from keras.models import Model
from keras.layers import Input, Layer, Dense, Activation
from sklearn.base import BaseEstimator

from model.model_utils import evaluate

from os.path import join, realpath, basename, dirname
current_dir = dirname(realpath('__file__'))
sys.path.insert(0, dirname(current_dir))


class STAR_NN(BaseEstimator):

    def __init__(self, model_params):

        self.model_params = model_params['model_params']

        self.num_epoch = model_params['fitting_params']['epoch']
        self.batch_size = model_params['fitting_params']['batch_size']
        self.verbose = model_params['fitting_params']['verbose']

        pid = os.getpid()
        time_stamp = '_{0:%b}-{0:%d}_{0:%H}-{0:%M}-{0:%S}'.format(datetime.datetime.now())
        self.save_filename = join(current_dir, model_params['fitting_params']['save_name'] + str(pid) + time_stamp)


    def fit(self, X_train, y_train, X_val, y_val, X_test, y_test, col_name_w_variant_type, random_seed):
        
        # print("\nIn STAR_NN_model.STAR_NN.fit(), model: ")

        model = get_compiled_model(
            n_features=X_train.shape[1],
            col_name_w_variant_type=col_name_w_variant_type,
            sparse=True,
            sparse_first_layer=True,
            activation=self.model_params['activation'],
            use_bias=self.model_params['use_bias'],
            kernel_initializer=self.model_params['kernel_initializer'],
            w_reg=self.model_params['reg'],            
            )

        # print("self.model=model, model=self.get_compiled_model")
        self.model = model


        early_stopping = tf.keras.callbacks.EarlyStopping(
            monitor='val_loss', 
            verbose=self.verbose,
            patience=50,
            mode='min',
            restore_best_weights=True)


        def scheduled_lr_decay(epoch, lr):
            if epoch < 300:
                return lr
            else:
                decay = 1e-5/30
                return float(lr * 1/(1+decay*epoch))

        lr_decay = tf.keras.callbacks.LearningRateScheduler(scheduled_lr_decay, verbose=0)


        # fit model 
        tf.keras.backend.clear_session()

        tf.random.set_seed(random_seed)


        # print("Fit model ...")
        history = self.model.fit(
            x=X_train, y=y_train, 
            batch_size=self.batch_size, epochs=self.num_epoch, 
            callbacks=[early_stopping, lr_decay],
            validation_data=(X_val, y_val),
            verbose=self.verbose, use_multiprocessing=True)


        # output
        test_results = self.model.evaluate(X_test, y_test, verbose=0)
        for name, value in zip(model.metrics_names, test_results):
            print(name, ': ', value)

        test_prediction = self.model.predict(X_test)
        print("Test prediction:", test_prediction)


        # get model weights
        if export_model_weights == True:
            # variant level weight
            variant_level_weight = model.layers[1].get_weights()[0]
            # variant_level_bias = model.layers[1].get_weights()[1]

            # gene level weight
            gene_weight = model.layers[2].get_weights()[0]
            # gene_level_bias = model.layers[2].get_weights()[1]

            # save variant level weight to output file
            variant_level_weight_df = pd.DataFrame(list(zip(col_name_w_variant_type, variant_level_weight)), columns=['variant', 'weight'])
            variant_level_weight_df_filename = "../_saving_dir/wes12.deepvariant.selFeat1489.variant_level_weight." + str(random_seed) + ".csv"
            variant_level_weight_df.to_csv(variant_level_weight_df_filename, index=False)
            # save gene level weight to output file
            gene_level = [item[0] for item in col_name_w_variant_type]
            gene_level_weight_df = pd.DataFrame(list(zip(gene_level[0:len(gene_level):3], gene_weight[:,0])), columns=['gene', 'weight'])
            gene_level_weight_df_filename = "../_saving_dir/wes12.deepvariant.selFeat1489.gene_level_weight." + str(random_seed) + ".csv"
            gene_level_weight_df.to_csv(gene_level_weight_df_filename, index=False)

        return history


    def predict(self, X):
        predictions = self.model.predict(X)
        return predictions

    def evaluate(self, X_test, y_test):
        results = self.model.evaluate(X_test, y_test, verbose=0)
        return results



# assume the inputs are connected to the layer nodes according to a 3 to 1 pattern
# the first node is connect to the first 3 inputs (PTV, MisAB, MisC), the second node is connect to the second 3 inputs and so on.
class Sparse_mapping(Layer):
    def __init__(self, units, activation=None, use_bias=None, kernel_initializer='glorot_uniform', bias_initializer='zeros', W_regularizer=None, bias_regularizer=None, kernel_constraint=None, bias_constraint=None, **kwargs):
        
        self.units = units
        self.activation = activation
        self.activation_fn = activations.get(activation)
        self.use_bias = use_bias

        self.kernel_initializer = initializers.get(kernel_initializer) 
        self.bias_initializer = initializers.get(bias_initializer)         
        
        self.W_regularizer = W_regularizer
        self.kernel_regularizer = regularizers.get(W_regularizer)
        self.bias_regularizer = regularizers.get(bias_regularizer)

        self.kernel_constraint = constraints.get(kernel_constraint)
        self.bias_constraint = constraints.get(bias_constraint)
        super(Sparse_mapping, self).__init__(**kwargs)

    # the number of weights equals the number of inputs to the layer
    def build(self, input_shape):

        print("\nin Sparse_mapping.build(input_shape):", input_shape)

        # create a trainable weight variable for this layer
        input_dimension = input_shape[1]
        self.kernel_shape = (input_dimension, self.units)
        self.n_inputs_per_node = input_dimension / self.units

        print("\nin compile_models.Diagnal.build()")
        print("input dimension {} self.units {}".format(input_dimension, self.units)) # input dimension 57375 self.units 19125
        print("n_inputs_per_node {}".format(self.n_inputs_per_node)) # n_inputs_per_node 3.0
        print("kernel shape", self.kernel_shape, "\n") # kernel shape (57375, 19125)

        self.kernel = self.add_weight(name = "kernel", shape=(input_dimension,), initializer=self.kernel_initializer, regularizer=self.kernel_regularizer, constraint=self.kernel_constraint, trainable=True)

        if self.use_bias:
            self.bias = self.add_weight(name = "bias", shape=(self.units,), initializer=self.bias_initializer, regularizer=self.bias_regularizer, constraint=self.bias_constraint, trainable=True)
        else:
            self.bias = None

        super(Sparse_mapping, self).build(input_shape)


    def call(self, x, mask=None):

        n_features = x.shape[1]
        print("\nin Sparse_mapping.call(), n_feautures:", n_features) # in Sparse_mapping.call(), n_feautures: 57375

        kernel = K.reshape(self.kernel, (1, n_features))
        print("kernel.shape", kernel.shape) # kernel.shape (1, 57375)

        mult = x*kernel
        mult = K.reshape(mult, (-1, int(self.n_inputs_per_node)))
        mult = K.sum(mult, axis=1)
        output = K.reshape(mult, (-1, self.units))
        print("output shape:", output.shape) # output shape: (None, 19125)

        if self.use_bias:
            output = K.bias_add(output, self.bias)
        if self.activation_fn is not None:
            output = self.activation_fn(output)

        return output


    def get_config(self):
        config = {
            'units':self.units,
            'activate':self.activation,
            'use_bias': self.use_bias
        }
        base_config = super(Sparse_mapping, self).get_config()
        return dict(list(base_config.items()) + list(config.items()))




def get_compiled_model(n_features, col_name_w_variant_type, sparse, sparse_first_layer, activation, use_bias, kernel_initializer, w_reg):
    
    # print("\nIn STAR_NN_model.get_compiled_model()")

    features = col_name_w_variant_type
    genes = col_name_w_variant_type.levels[0]
    n_genes = len(genes)

    METRICS = [
        tf.keras.metrics.BinaryAccuracy(name='accuracy'),
        tf.keras.metrics.Precision(name='precision'),
        tf.keras.metrics.Recall(name='recall'),
        tf.keras.metrics.AUC(name='auc'),
    ]

    # compile model
    inputs = Input(shape=(n_features, ), dtype="float32", name="inputs")

    if sparse_first_layer:
        print("\nSparse_first_layer...\n")
        layer1 = Sparse_mapping(input_shape=(n_features, ), units=n_genes, activation=activation, use_bias=use_bias, kernel_initializer=kernel_initializer, **{})
    else:
        print("\nDense_first_layer...\n")
        layer1 = Dense(input_shape=(n_features, ), units=n_genes, activation=activation, use_bias=use_bias, kernel_initializer=kernel_initializer)    

    layer1_out = layer1(inputs)
    # layer2_out = Dense(units = 128, activation="relu", use_bias=True, bias_initializer=tf.keras.initializers.GlorotUniform())(layer1_out)
    final_out = Dense(units = 1, activation="sigmoid", use_bias=True, bias_initializer=tf.keras.initializers.GlorotUniform())(layer1_out)
    model = Model(inputs=[inputs], outputs=final_out)

    print("final_out shape:", final_out.shape)  # (None,1)
    print("Compiling ... ")
    model.compile(
        optimizer=tf.keras.optimizers.Adam(learning_rate=3*1e-5),
        loss=tf.keras.losses.binary_crossentropy,
        metrics=METRICS)
    
    return model