import pandas as pd
from anndata import AnnData
import numpy as np
import os

def Anndata_to_csv(adata,export_dir):
    if os.path.isdir(export_dir):
        pass
    else:
        os.mkdir(export_dir)
    adata.write_csvs(export_dir, skip_data = False)

def csv_to_Anndata(import_dir): # import dir = location of your .csv files
    X = pd.read_csv(import_dir + "/X.csv", header = None)
    obs = pd.read_csv(import_dir + "/obs.csv", header = 0, index_col = 0)
    var = pd.read_csv(import_dir + "/var.csv", header = 0, index_col = 0)
    obsm_df = pd.read_csv(import_dir + "/obsm.csv", header = 0)
    X.columns = var.index

    adata_import = AnnData(X = X, obs = obs, var = var)
    X_pca = list()
    X_umap = list()
    spatial = list()

    for colname in obsm_df.columns:
        if 'X_pca' in colname:
            X_pca.append(obsm_df[colname].to_list())
        elif 'X_umap' in colname:
            X_umap.append(obsm_df[colname].to_list())
        elif 'spatial' in colname:
            spatial.append(obsm_df[colname].to_list())

    if len(X_pca) != 0:
        X_pca = np.swapaxes(np.array(X_pca, dtype = np.float32),0,1)
        adata_import.obsm['X_pca'] = X_pca
    if len(X_umap) != 0:
        X_umap = np.swapaxes(np.array(X_umap, dtype = np.float32),0,1)
        adata_import.obsm['X_umap'] = X_umap
    if len(spatial) != 0:
        spatial = np.swapaxes(np.array(spatial), 0, 1)
        adata_import.obsm['spatial'] = spatial
        
    return adata_import # returns adata object