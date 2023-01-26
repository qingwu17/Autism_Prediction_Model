import numpy as np
import pandas as pd

from os.path import join

data_path_wes12 = "/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/pub/WES12/DeepVariant/"
data_path_iWESv1 = "/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/pub/iWES_v1/"

# region: input gene name files, including SFARI gene, DDG2P genes and selected gene feature from top 4% select Percentile()
SFARI_gene_df = pd.read_csv("/users/qwu24/data/silvio/Qing_Wu/SFARI/reference_data/SFARI_gene/SFARI-Gene_genes_01-11-2022release_03-22-2022export.csv")
SFARI_gene_lst = SFARI_gene_df['gene-symbol'].tolist() # 1031

DDG2P_gene_df = pd.read_csv("/users/qwu24/data/silvio/Qing_Wu/SFARI/reference_data/DDG2P/DDG2P_7_4_2022.csv")
DDG2P_gene_lst = DDG2P_gene_df['gene symbol'].tolist() # 2580

selFeat_names_wes12 = pd.read_csv("/users/qwu24/data/silvio/Qing_Wu/SFARI/batch_jobs/python_script/machine_learning_script/tpot_out/wes12.deepvariant.feature_name_union.csv")
selFeat_names_wes12_lst = selFeat_names_wes12['selected_features'].tolist() # 1489  

selFeat_names_exonic_wes12 = pd.read_csv("/users/qwu24/data/silvio/Qing_Wu/SFARI/batch_jobs/python_script/machine_learning_script/tpot_out/wes12.deepvariant.feature_name_exonic.csv")
selFeat_names_exonic_wes12_lst = selFeat_names_exonic_wes12['selected_features'].tolist() # 956

selFeat_names_PTV_MisAB_wes12 = pd.read_csv("/users/qwu24/data/silvio/Qing_Wu/SFARI/batch_jobs/python_script/machine_learning_script/tpot_out/wes12.deepvariant.feature_name_PTV_MisAB.csv")
selFeat_names_PTV_MisAB_wes12_lst = selFeat_names_PTV_MisAB_wes12['selected_features'].tolist() # 669

# endregion

# region: SPxG matrix dataframe with variants to gene index list (consecutive and unconsecutive version)
# PTVs
wes12_GxSP_from_PTVs = 'wes1_wes2_combined.deepvariant.rare1pct_variants_het_by_sample_matrix_cleaned_2GxSP.PTVs.txt'
wes12_PTVs_V2G_unconsec_lst_filename ='wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.unconsecutive_indices_splited.PTVs.txt'
wes12_PTVs_V2G_lst_filename ='wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.PTVs.txt'

# MisAB
wes12_GxSP_from_MisAB = 'wes1_wes2_combined.deepvariant.rare1pct_variants_het_by_sample_matrix_cleaned_2GxSP.MisAB.txt'
wes12_MisAB_V2G_unconsec_lst_filename ='wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.unconsecutive_indices_splited.MisAB.txt'
wes12_MisAB_V2G_lst_filename ='wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.MisAB.txt'

# # MisC(nif)
# wes12_GxSP_from_MisC = 'wes1_wes2_combined.deepvariant.rare1pct_variants_het_by_sample_matrix_cleaned_2GxSP.non_impactful.txt'
# wes12_nif_V2G_unconsec_lst_filename ='wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.unconsecutive_indices_splited.non_impactful.txt'
# wes12_nif_V2G_lst_filename ='wes1_wes2_combined.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.non_impactful.txt'

# wes12_GxSP_matrix_file_lst = (wes12_GxSP_from_PTVs, wes12_GxSP_from_MisAB, wes12_GxSP_from_MisC)
# wes12_V2G_unconsec_file_lst = (wes12_PTVs_V2G_unconsec_lst_filename, wes12_MisAB_V2G_unconsec_lst_filename, wes12_nif_V2G_unconsec_lst_filename)
# wes12_V2G_file_lst = (wes12_PTVs_V2G_lst_filename, wes12_MisAB_V2G_lst_filename, wes12_nif_V2G_lst_filename)

wes12_GxSP_matrix_file_lst = (wes12_GxSP_from_PTVs, wes12_GxSP_from_MisAB)
wes12_V2G_unconsec_file_lst = (wes12_PTVs_V2G_unconsec_lst_filename, wes12_MisAB_V2G_unconsec_lst_filename)
wes12_V2G_file_lst = (wes12_PTVs_V2G_lst_filename, wes12_MisAB_V2G_lst_filename)

wes12_paired_SPxG_file_lst = list(zip(wes12_GxSP_matrix_file_lst, wes12_V2G_unconsec_file_lst, wes12_V2G_file_lst))

# endregion

# region: iWES1: PTVs, MisAB, MisC
# PTVs
iwes1_GxSP_matrix_file_exonic_PTVs = "wes_70487_exome.deepvariant.rare1pct_variants_het_by_sample_matrix_cleaned_2GxSP.PTVs.txt"
iwes1_var2gene_unconsec_lst_file_exonic_PTVs ='wes_70487_exome.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.unconsecutive_indices_splited.PTVs.txt'
iwes1_var2gene_lst_file_exonic_PTVs ='wes_70487_exome.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.PTVs.txt'

# MisAB
iwes1_GxSP_matrix_file_exonic_MisAB = "wes_70487_exome.deepvariant.rare1pct_variants_het_by_sample_matrix_cleaned_2GxSP.MisAB.txt"
iwes1_var2gene_unconsec_lst_file_exonic_MisAB ='wes_70487_exome.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.unconsecutive_indices_splited.MisAB.txt'
iwes1_var2gene_lst_file_exonic_MisAB ='wes_70487_exome.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.MisAB.txt'

# # MisC(nif)
# iwes1_GxSP_matrix_file_exonic_MisC = "wes_70487_exome.deepvariant.rare1pct_variants_het_by_sample_matrix_cleaned_2GxSP.non_impactful.txt"
# iwes1_iwes1_var2gene_unconsec_lst_file_exonic_MisC ='wes_70487_exome.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.unconsecutive_indices_splited.non_impactful.txt'
# var2gene_lst_file_exonic_MisC ='wes_70487_exome.deepvariant.rare1pct_variants_het.idx_var2gene_lst_df.non_impactful.txt'

# iwes1_GxSP_matrix_file_lst = (iwes1_GxSP_matrix_file_exonic_PTVs, iwes1_GxSP_matrix_file_exonic_MisAB, iwes1_GxSP_matrix_file_exonic_MisC)
# iwes1_V2G_unconsec_file_lst = (iwes1_var2gene_unconsec_lst_file_exonic_PTVs, iwes1_var2gene_unconsec_lst_file_exonic_MisAB, iwes1_iwes1_var2gene_unconsec_lst_file_exonic_MisC)
# iwes1_V2G_file_lst = (iwes1_var2gene_lst_file_exonic_PTVs, iwes1_var2gene_lst_file_exonic_MisAB, var2gene_lst_file_exonic_MisC)

iwes1_GxSP_matrix_file_lst = (iwes1_GxSP_matrix_file_exonic_PTVs, iwes1_GxSP_matrix_file_exonic_MisAB)
iwes1_V2G_unconsec_file_lst = (iwes1_var2gene_unconsec_lst_file_exonic_PTVs, iwes1_var2gene_unconsec_lst_file_exonic_MisAB)
iwes1_V2G_file_lst = (iwes1_var2gene_lst_file_exonic_PTVs, iwes1_var2gene_lst_file_exonic_MisAB)

iwes1_paired_SPxG_file_lst = list(zip(iwes1_GxSP_matrix_file_lst, iwes1_V2G_unconsec_file_lst, iwes1_V2G_file_lst))
# endregion


def get_SPxis_case_is_female_header_2D_arr(SP_col_header_file):
    
    # read SP column header file, columns: sp_id, is_case, is_female
    SP_col_header_lst = []
    with open(SP_col_header_file) as tmp_file:
        for line in tmp_file:
            SP_col_header_lst.append(line.rstrip().split("\t"))
    SP_col_header_2D_array = np.array(SP_col_header_lst)
    SP_col_header_df = pd.DataFrame(SP_col_header_2D_array)
    # add column names to SPxV matrix
    col_names = ['SP_ID', 'is_female', 'is_case']
    SP_col_header_df.columns = col_names
    print("Shape of SP_id_header numpy array:", SP_col_header_2D_array.shape) # (43203, 3)
    return SP_col_header_2D_array

def get_GxSP_2D_arr(GxSP_matrix_file):
    # read gene by sp matrix file
    GxSP_matrix_file = GxSP_matrix_file
    GxSP_matrix_lst = []
    with open(GxSP_matrix_file) as tmp_file:
        for line in tmp_file:
            GxSP_matrix_lst.append(line.rstrip().split('\t'))
    GxSP_matrix = np.array(GxSP_matrix_lst)
    print("Shape of GxSP matrix numpy array:", GxSP_matrix.shape) # (19626, 43203)
    return GxSP_matrix

def get_gene_lst_of_GxSP(var2gene_lst_file):
    # read gene list, read variant to gene index file, each line represent the index of variant to one gene
    # the first vector of each index line is the gene name
    V2G_idx_file = var2gene_lst_file
    gene_lst = []
    V2G_idx_lst = []
    with open(V2G_idx_file) as tmp_file:
        for line in tmp_file:
            var2gene_idx_w_gene = line.rstrip().split('\t')
            gene_name = var2gene_idx_w_gene[0]
            gene_lst.append(gene_name)
            # get the idx of variant to the gene of that row
            var2gene_idx = [int(x)-1 for x in var2gene_idx_w_gene[1:]]
            # var2gene_idx = var2gene_idx[1:]
            V2G_idx_lst.append(var2gene_idx)
    print("Length of gene converted by rare exonic variants:", len(gene_lst)) # 18825
    return gene_lst, V2G_idx_lst

def get_SPxPGS_header_df(SPxPGS_header_file):
    # read risk score generated from GWAS, columns: sp_id, risk score (centered 0)
    SPxPGS_header_lst = []
    with open(SPxPGS_header_file) as tmp_file:
        for line in tmp_file:
            SPxPGS_header_lst.append(line.rstrip().split("\t"))
    SPxPGS_header_2D_array = np.array(SPxPGS_header_lst)
    SPxPGS_header_df = pd.DataFrame(SPxPGS_header_2D_array)
    print("Shape of SP_id by PGS pandas dataframe:", SPxPGS_header_df.shape) # (43227, 2)
    return SPxPGS_header_df

def get_exclude_sp_id_idx(idx_sp2exclude_file):
    # read sp idx file 
    idx_sp2exclude_lst = []
    with open(idx_sp2exclude_file) as tmp_file:
        for line in tmp_file:
            idx_sp2exclude_lst.append(int(line.rstrip())-1)
    print(len(idx_sp2exclude_lst))
    return idx_sp2exclude_lst

def read_SPxG_matrix_table(data_path, SPxis_case_is_female_header_filename, GxSP_matrix_filename, var2gene_unconsec_lst_filename, var2gene_lst_filename, idx_sp2exclude_filename, PGS_filename, iWESv1_subset=False):
    
    # GxSP_matrix_filename = paired_SPxG_file_lst[0][0]
    # var2gene_unconsec_lst_filename = paired_SPxG_file_lst[0][1]
    # var2gene_lst_filename = paired_SPxG_file_lst[0][2]

    GxSP_file_full_path = join(data_path, GxSP_matrix_filename)
    GxSP_matrix = get_GxSP_2D_arr(GxSP_file_full_path)

    # transpose VxSP matrix into SPxV matrix
    SPxG_matrix = GxSP_matrix.T
    print("Shape of SPxG matrix numpy array:", SPxG_matrix.shape) # (43203, 19626)

    # paste sp_id, is_case, is_female to SPxV matrix
    SPxis_case_is_female_header_file_full_path = join(data_path, SPxis_case_is_female_header_filename)
    SPxis_case_is_female_header_matrix = get_SPxis_case_is_female_header_2D_arr(SPxis_case_is_female_header_file_full_path) # (43203, 3)
    SPxG_matrix_full = np.hstack((SPxis_case_is_female_header_matrix, SPxG_matrix))

    # change numpy 2D array to pd.DataFrame
    SPxG_matrix_full = pd.DataFrame(SPxG_matrix_full)

    # add column names to SPxV matrix
    var2gene_unconsec_lst_file_full_path = join(data_path, var2gene_unconsec_lst_filename)
    gene_lst_w_unconsecutive_indices_splited, _ = get_gene_lst_of_GxSP(var2gene_lst_file=var2gene_unconsec_lst_file_full_path) 
    col_names = [['SP_ID', 'is_female', 'is_case'], gene_lst_w_unconsecutive_indices_splited]
    SPxG_matrix_full.columns = sum(col_names, [])
    print("Shape of SPxG with header column:", SPxG_matrix_full.shape) # (43203, 19629)

    # combined same gene name columns into one
    new_SPxG_matrix_full = SPxG_matrix_full.groupby(level=0, axis=1).sum()
    # print(new_SPxG_matrix_full) # [43203 rows x  columns]

    # order columns by original gene name list
    var2gene_lst_file_full_path = join(data_path, var2gene_lst_filename)
    gene_lst, _ = get_gene_lst_of_GxSP(var2gene_lst_file=var2gene_lst_file_full_path) # 19117
    gene_lst.insert(0, 'is_female')
    gene_lst.insert(0, 'is_case')
    gene_lst.insert(0, 'SP_ID')

    SPxG_matrix_full = new_SPxG_matrix_full[gene_lst]
    # print(SPxG_matrix_full) # [43203 rows x  columns]

    # exclude selected sp from SPxG_matrix_full df
    idx_sp2exclude_file_full_path = join(data_path, idx_sp2exclude_filename)
    idx_sp2excl = get_exclude_sp_id_idx(idx_sp2exclude_file_full_path)
    SPxG_matrix_full = SPxG_matrix_full.drop(idx_sp2excl)
    print("Shape of SPxG matrix after droping individuals without PGS:", SPxG_matrix_full.shape) # (40871, )

    # add PGS columns
    SPxPGS_header_file_full_path = join(data_path, PGS_filename)
    SPxPGS_header_df = get_SPxPGS_header_df(SPxPGS_header_file_full_path)  
    # paste sp_id_PGS to SPxV matrix
    SPxG_matrix_full.insert(loc=2, column="PGS", value=SPxPGS_header_df[5].astype(float).tolist())
    print("Shape of SPxG matrix with PGS column: ", SPxG_matrix_full.shape) 

    if iWESv1_subset == True:
        idx_sp2exclude_lst_subset = get_exclude_sp_id_idx(idx_sp2exclude_file='/users/qwu24/data/silvio/Qing_Wu/SFARI/hail_out/SPARK/pub/iWES_v1/iWES_v1.deepvariant.rare1pct_variants_het_by_sample_matrix_cleaned.spid2subset_wPGS_27430.txt') # 640
        # reset index to 0-based
        SPxG_matrix_full = SPxG_matrix_full.reset_index(drop=True)
        # exclude selected sp from VxSP matrix
        SPxG_matrix_full = SPxG_matrix_full.drop(idx_sp2exclude_lst_subset)
    else:
        SPxG_matrix_full = SPxG_matrix_full

    # prepare input file
    samples_SPID = SPxG_matrix_full['SP_ID']
    X = SPxG_matrix_full.iloc[:,2::].astype(np.float32)
    y = SPxG_matrix_full.iloc[:,1].astype(np.float32)
    # print(X.dtypes) # float32
    print("Shape of input data wo SP_id and label :", X.shape) 

    # get input data(x), labels, sample_name, genes 
    X.index = samples_SPID
    y.index = samples_SPID

    # col_name_X = X.columns

    return X, y

def filter_SPxG_matrix_by_selected_feature(set_name, feature_lst, X):
    
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
        selected_columns.insert(0, "PGS")
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

    # # check if all participants have mutations in listed genes
    # sum_row = X_filtered.iloc[:,2::].sum(axis=1) 
    # # print((sum_row > 0).value_counts()) # all participants have selected type of mutations in listed genes
    # # if not all participants have mutations in selected features, remove the participants
    # row_filter = sum_row > 0
    # if any(row_filter) > 0:
    #     X_filtered = X_filtered.loc[row_filter,:]
    #     y_filtered = y[row_filter]
    # # print(X_filtered.shape, len(y_filtered))
    # # print(y_filtered.value_counts())

    # # check if all genes have mutations in certain people among the whole dataset
    # sum_col = X_filtered.sum(axis=0)
    # # print((sum_col > 0).value_counts()) # all listed genes get mutated in certain participants
    # # if not all geen have mutations in all participants, remove the gene columns
    # col_filter = sum_col > 0
    # if any(col_filter) > 0:
    #     X_filtered = X_filtered.loc[:,col_filter]
    # # print(X_filtered.shape)
    # # print(len(X_filtered.columns.tolist()))

    # y_filtered = y.astype(int)

    # print(X_filtered.shape)

    return X_filtered

def combine_SPxG_matrix_table(SPxG_matrix_lst, labels, X_col_name_lst, variant_type_lst):
    
    # SPxG_matrix_lst = X_lst
    # labels = y_lst

    # flatten gene nemas in the list to one single list
    col_names_lst_set = [set(list(col_names)) for col_names in X_col_name_lst]
    # make a union and unique list of gene name
    col_names = list(set.union(*col_names_lst_set))
    print("Length of union genes plus PGS and is_female:", len(col_names))

    tmp_df = pd.DataFrame(index=col_names)

    # combine data frame by gene name
    combined_SPxG_df_lst = list()
    for X, y, X_col_name in zip(SPxG_matrix_lst, labels, X_col_name_lst):
        df = pd.DataFrame(X, columns=X_col_name)
        df = df.T.join(tmp_df, how="right")
        df = df.T
        df = df.fillna(0)
        combined_SPxG_df_lst.append(df)

    # outer combined data frame and got variant type as first level of column names
    combined_SPxG_df_0 = pd.concat(combined_SPxG_df_lst, keys=variant_type_lst, join="outer", axis=1)
    # switch first level and second level column index, genes on the first level and variant type on the second level
    combined_SPxG_df = combined_SPxG_df_0.swaplevel(i=0, j=1, axis=1)
    # order the columns by is_female+PGS+gene name order
    col_order = combined_SPxG_df.columns.levels[0].tolist()
    col_order.insert(0, col_order.pop(col_order.index("is_female")))
    col_order.insert(0, col_order.pop(col_order.index("PGS")))
    combined_SPxG_df = combined_SPxG_df.reindex(columns=col_order, level=0)

    # check if all columns in PGS are the same, if so, keep values in one column and set the rest to 0
    if all((y == combined_SPxG_df['PGS']['PTV']).all() for y in combined_SPxG_df['PGS'].to_numpy().T):
        combined_SPxG_df['PGS']['MisAB'].values[:] = 0
        # combined_SPxG_df['PGS']['MisC'].values[:] = 0

    # check if all columns in is_female are the same, if so, keep values in one column and set the rest to 0
    if all((y == combined_SPxG_df['is_female']['PTV']).all() for y in combined_SPxG_df['is_female'].to_numpy().T):
        combined_SPxG_df['is_female']['MisAB'].values[:] = 0
        # combined_SPxG_df['is_female']['MisC'].values[:] = 0
    
    # check if labels are identical from merged SPxG matrix
    if all((y == labels[0]).all() for y in labels):
        combined_labels = labels[0]

    # check if all participants have mutations in listed genes
    sum_row = combined_SPxG_df.iloc[:,4::].sum(axis=1) 
    # print((sum_row > 0).value_counts()) # all participants have selected type of mutations in listed genes
    # if not all participants have mutations in selected features, remove the participants
    row_filter = sum_row > 0
    if any(row_filter) > 0:
        combined_SPxG_df_filtered = combined_SPxG_df.loc[row_filter,:]
        combined_labels_filtered = combined_labels[row_filter]
    # print(combined_SPxG_df_filtered)
    # print(combined_labels_filtered)

    # get SP_ID and gene name list for combined SPxG matrix
    combined_samples_SPID = combined_SPxG_df_filtered.index
    col_name_w_variant_type = combined_SPxG_df_filtered.columns
    
    return combined_SPxG_df_filtered, combined_labels_filtered, combined_samples_SPID, col_name_w_variant_type



class Combined_SPxG_from_variant_type():
    def __init__(self, variant_type_lst, selected_genes, training=True):

        if training:
            paired_SPxG_file_lst=wes12_paired_SPxG_file_lst
            data_path = data_path_wes12
            SPxis_case_is_female_header_filename = 'wes1_wes2_combined.deepvariant.rare1pct_variants_het_by_sample_matrix_cleaned.sp_sex_case_43203.txt'
            idx_sp2exclude_filename='wes1_wes2_combined.deepvariant.rare1pct_variants_het_by_sample_matrix_cleaned.spid2excl_noPGS.txt'
            PGS_filename='wes1_wes2_combined.deepvariant.rare1pct_variants_het_by_sample_matrix_cleaned.spid2PGS_40871.txt'
            iWESv1_subset=False

            X_lst, y_lst, samples_SPID_lst, X_col_name_lst = list(), list(), list(), list()
            for files in paired_SPxG_file_lst:
                
                print("\nInput SPxG matrix file: ", files[0])

                input_SPxG_matrix, labels = read_SPxG_matrix_table(
                    data_path=data_path, 
                    SPxis_case_is_female_header_filename=SPxis_case_is_female_header_filename, 
                    GxSP_matrix_filename=files[0], var2gene_unconsec_lst_filename=files[1], var2gene_lst_filename=files[2], 
                    idx_sp2exclude_filename=idx_sp2exclude_filename,
                    PGS_filename=PGS_filename,
                    iWESv1_subset=False
                )

                print(input_SPxG_matrix.shape, labels.shape)

                if (selected_genes is not None):
                    input_SPxG_matrix = filter_SPxG_matrix_by_selected_feature(set_name=selected_genes, feature_lst=selFeat_names_wes12_lst, X=input_SPxG_matrix)

                print(input_SPxG_matrix.shape, labels.shape)
                X_lst.append(input_SPxG_matrix), y_lst.append(labels), X_col_name_lst.append(input_SPxG_matrix.columns) 

        else:
            paired_SPxG_file_lst=iwes1_paired_SPxG_file_lst
            data_path=data_path_iWESv1
            SPxis_case_is_female_header_filename='wes_70487_exome.deepvariant.rare1pct_variants_het_by_sample_matrix_cleaned.sp_sex_case_70464.txt'
            idx_sp2exclude_filename='iWES_v1.deepvariant.rare1pct_variants_het_by_sample_matrix_cleaned.spid2excl_noPGS.txt'
            PGS_filename='iWES_v1.deepvariant.rare1pct_variants_het_by_sample_matrix_cleaned.spid2PGS_68301.txt'
            iWESv1_subset=True

            X_lst, y_lst, samples_SPID_lst, X_col_name_lst = list(), list(), list(), list()
            for files in paired_SPxG_file_lst:
                
                print("\nInput SPxG matrix file: ", files[0])

                input_SPxG_matrix, labels, samples_SPID = read_SPxG_matrix_table(
                    data_path=data_path, 
                    SPxis_case_is_female_header_filename=SPxis_case_is_female_header_filename, 
                    GxSP_matrix_filename=files[0], var2gene_unconsec_lst_filename=files[1], var2gene_lst_filename=files[2], 
                    idx_sp2exclude_filename=idx_sp2exclude_filename,
                    PGS_filename=PGS_filename,
                    iWESv1_subset=True
                )

                if (selected_genes is not None):
                    input_SPxG_matrix, labels = filter_SPxG_matrix_by_selected_feature(set_name=selected_genes, feature_lst=selFeat_names_wes12_lst, X=input_SPxG_matrix, y=labels)

                print(input_SPxG_matrix.shape, labels.shape, samples_SPID.shape)
                X_lst.append(input_SPxG_matrix), y_lst.append(labels), samples_SPID_lst.append(samples_SPID), X_col_name_lst.append(input_SPxG_matrix.columns) 
        
        combined_X, combined_labels, combined_samples_SPID_lst, combined_X_col_name_lst = combine_SPxG_matrix_table(X_lst, y_lst, X_col_name_lst, variant_type_lst)
        print("Shape of combined_X:", combined_X.shape)

        self.X = combined_X
        self.y = combined_labels
        self.samples_SPID = combined_samples_SPID_lst
        self.X_col_name_lst = combined_X_col_name_lst

    def get_data(self):
        return self.X, self.y, self.samples_SPID, self.X_col_name_lst
