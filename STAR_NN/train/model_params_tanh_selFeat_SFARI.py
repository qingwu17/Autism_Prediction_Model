# parameters used in the model creation and data extraction

selected_genes = "selFeat_SFARI_genes"
data = {
    'id': 'selFeat_SFARI_genes', 
    'type': 'PTV_MisAB_MisC',
    'params': {
        'variant_type_lst': ["PTV", "MisAB", "MisC"],
        'selected_genes': selected_genes,
        'validation': True
    }
}

model = {
    'id': 'sparse-net',
    'type': 'classification_binary',
    'params': {
        'model_params': {
            'data_params': data,
            'activation': 'tanh',
            'use_bias': True,
            'kernel_initializer': 'lecun_uniform',
            'reg': 0.001,
        },
        'fitting_params': dict(
                            epoch=5000,
                            batch_size=128,
                            verbose=0,
                            save_name='net'
        ),
    }
}