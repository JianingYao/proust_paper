import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import os
import sys
from sklearn.metrics.cluster import adjusted_rand_score
import STAGATE


# the location of R (used for the mclust clustering)
os.environ['R_HOME'] = '/Library/Frameworks/R.framework/Resources'
os.environ['R_USER'] = '/Users/jianingyao/opt/anaconda3/envs/STAGATE/lib/python3.9/site-packages/rpy2'

datasets = [
     'Visium_CKp25_rep1',
     'Visium_CKp25_rep2',
     'Visium_CKp25_rep3',
     'Visium_CKp25_rep4']

for dataset in datasets:
    n_clusters = 20
    print("################################ " + str(dataset) + " #####################################")
    file_fold = '/Users/jianingyao/Desktop/Research/Biostatistics_JHU/PhD/Data/proust_datasets/Visium-mouse_DSB_coronal_forebrain_brain/'+ str(dataset)
    adata = sc.read_visium(file_fold, count_file='filtered_feature_bc_matrix.h5', load_images=True)
    adata.var_names_make_unique()
    adata.obsm['spatial'] = adata.obsm['spatial'].astype(int)
    adata

    # Normalization
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    STAGATE.Cal_Spatial_Net(adata, rad_cutoff=300)
    STAGATE.Stats_Spatial_Net(adata)

    import tensorflow as tf
    tf.compat.v1.disable_eager_execution()
    # adata = STAGATE.train_STAGATE(adata, alpha=0.5, pre_resolution=0.2,
    #                               n_epochs=1000, save_attention=True)
    adata = STAGATE.train_STAGATE(adata, alpha=0)

    adata = STAGATE.mclust_R(adata, used_obsm='STAGATE', num_cluster=n_clusters)

    adata.write_h5ad("./Result/adata/" + str(dataset) + "_clusters.h5ad")

    obs_df = adata.obs.dropna()

    adata.obsm['spatial'][:, 1] = -1 * adata.obsm['spatial'][:, 1]

    dpi = 500
    fig, ax = plt.subplots(figsize=(5, 5), dpi=dpi)
    sc.pl.embedding(adata,
                    basis="spatial",
                    color="mclust",
                    size=40,
                    show=False,
                    ax=ax)
    plt.savefig('./Result/plot/' + dataset + '_ref_stagate.png',
                dpi=dpi, bbox_inches='tight')
    plt.close()