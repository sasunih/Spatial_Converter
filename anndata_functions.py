import pandas as pd
from anndata import AnnData
import numpy as np
import os
import warnings

def Anndata_to_csv(adata,export_dir):
    if os.path.isdir(export_dir):
        pass
    else:
        os.mkdir(export_dir)
    adata.write_csvs(export_dir, skip_data = False)

    if 'spatial' in adata.uns.keys():
        if 'hires' in list(adata.uns['spatial'].values())[0]['images'].keys():
            uns_spatial = list(adata.uns['spatial'].values())[0]['images']['hires']
            uns_spatial_r = uns_spatial[:,:,0]
            uns_spatial_g = uns_spatial[:,:,1]
            uns_spatial_b = uns_spatial[:,:,2]

            np.savetxt(export_dir + "/uns/spatial_hires_r.csv", uns_spatial_r, delimiter=",", fmt="%d")
            np.savetxt(export_dir + "/uns/spatial_hires_g.csv", uns_spatial_g, delimiter=",", fmt="%d")
            np.savetxt(export_dir + "/uns/spatial_hires_b.csv", uns_spatial_b, delimiter=",", fmt="%d")

        if 'lowres' in list(adata.uns['spatial'].values())[0]['images'].keys():
            uns_spatial = list(adata.uns['spatial'].values())[0]['images']['lowres']
            uns_spatial_r = uns_spatial[:,:,0]
            uns_spatial_g = uns_spatial[:,:,1]
            uns_spatial_b = uns_spatial[:,:,2]

            np.savetxt(export_dir + "/uns/spatial_lowres_r.csv", uns_spatial_r, delimiter=",", fmt="%d")
            np.savetxt(export_dir + "/uns/spatial_lowres_g.csv", uns_spatial_g, delimiter=",", fmt="%d")
            np.savetxt(export_dir + "/uns/spatial_lowres_b.csv", uns_spatial_b, delimiter=",", fmt="%d")
        
        spot_size = pd.DataFrame(list(adata.uns['spatial'].values())[0]['scalefactors'], index = [0])
        spot_size.to_csv(export_dir + "/uns/spot_size.csv", header = True, index = False)
        
def csv_to_Anndata(import_dir): # import dir = location of your .csv files
    X = pd.read_csv(import_dir + "/X.csv", header = None)
    obs = pd.read_csv(import_dir + "/obs.csv", header = 0, index_col = 0)
    var = pd.read_csv(import_dir + "/var.csv", header = 0, index_col = 0)
    if os.path.exists(import_dir + "/obsm.csv"):
        with open(import_dir + "/obsm.csv", 'r') as f:
            first_line = f.readline().strip()
        if first_line != '':
            obsm_df = pd.read_csv(import_dir + "/obsm.csv", header = 0)
        
    X.columns = var.index

    adata_import = AnnData(X = X, obs = obs, var = var)
    X_pca = list()
    X_umap = list()
    spatial = list()

    try:
        obsm_df
    except NameError:
        pass
    else:
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
            
    if os.path.exists(import_dir + "/uns/spot_size.csv"):
        if os.path.exists(import_dir + "/uns/spatial_hires_r.csv"):
            r_array = np.loadtxt(import_dir + '/uns/spatial_hires_r.csv', delimiter = ',')
            g_array = np.loadtxt(import_dir + '/uns/spatial_hires_g.csv', delimiter = ',')
            b_array = np.loadtxt(import_dir + '/uns/spatial_hires_b.csv', delimiter = ',')

            uns_hires_image = np.array([r_array, g_array, b_array]).astype(int)
            uns_hires_image = np.transpose(uns_hires_image, (1, 2, 0))
        else:
            warnings.warn('High resolution image missing')
            uns_hires_image = None

        if os.path.exists(import_dir + "/uns/spatial_lowres_r.csv"):
            r_array = np.loadtxt(import_dir + '/uns/spatial_lowres_r.csv', delimiter = ',')
            g_array = np.loadtxt(import_dir + '/uns/spatial_lowres_g.csv', delimiter = ',')
            b_array = np.loadtxt(import_dir + '/uns/spatial_lowres_b.csv', delimiter = ',')

            uns_lowres_image = np.array([r_array, g_array, b_array]).astype(int)
            uns_lowres_image = np.transpose(uns_lowres_image, (1, 2, 0))
        else:
            warnings.warn('Low resolution image missing')
            uns_lowres_image = None
            
        spot_size = pd.read_csv(import_dir + '/uns/spot_size.csv', header = 0, index_col = False)
        uns_spatial_dict = {'visium':{'images':{'hires':uns_hires_image,'lowres':uns_lowres_image}, 'scalefactors':spot_size.to_dict(orient = 'records')[0]}}
        uns_spatial_dict['visium']['images'] = {k:v for k, v in uns_spatial_dict['visium']['images'].items() if v is not None}
        adata_import.uns['spatial'] = uns_spatial_dict
        
    return adata_import # returns adata object
