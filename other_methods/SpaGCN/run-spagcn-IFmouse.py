import os, csv, re
import pandas as pd
import numpy as np
import scanpy as sc
import math
import SpaGCN as spg
from scipy.sparse import issparse
import random, torch
import warnings

warnings.filterwarnings("ignore")
import matplotlib.colors as clr
import matplotlib.pyplot as plt
from sklearn import metrics
import SpaGCN as spg
# In order to read in image data, we need to install some package. Here we recommend package "opencv"
# inatll opencv in python
# !pip3 install opencv-python
import cv2

# Read original data and save it to h5ad
from scanpy import read_10x_h5
import imageio


datasets = [
     'Visium_CKp25_rep1',
     'Visium_CKp25_rep2',
     'Visium_CKp25_rep3',
     'Visium_CKp25_rep4']


for i in range(len(datasets)):
    dataset = datasets[i]
    print("################################ " + str(dataset) + " #####################################")
    file_fold = '/Users/jianingyao/Desktop/Research/Biostatistics_JHU/PhD/Data/proust_datasets/Visium-mouse_DSB_coronal_forebrain_brain/'+ str(dataset)
    adata = sc.read_visium(file_fold, count_file='filtered_feature_bc_matrix.h5', load_images=True)
    adata.var_names_make_unique()
    adata.obsm['spatial'] = adata.obsm['spatial'].astype(int)
    adata

    # spatial = pd.read_csv(file_fold + "/spatial/tissue_positions_list.csv", sep=",", header=None, na_filter=False, index_col=0)
    # adata.obs["x1"] = spatial[1]
    # adata.obs["x2"] = spatial[2]
    # adata.obs["x3"] = spatial[3]
    # adata.obs["x4"] = spatial[4]
    # adata.obs["x5"] = spatial[5]
    # adata.obs["x_array"] = adata.obs["x2"]
    # adata.obs["y_array"] = adata.obs["x3"]
    # adata.obs["x_pixel"] = adata.obs["x4"]
    # adata.obs["y_pixel"] = adata.obs["x5"]

    # Select captured samples
    # adata = adata[adata.obs["x1"] == 1]
    adata.var_names = [i.upper() for i in list(adata.var_names)]
    adata.var["genename"] = adata.var.index.astype("str")

    x_array = adata.obs["array_row"].tolist()
    y_array = adata.obs["array_col"].tolist()
    x_pixel = adata.obsm["spatial"].astype(int)[:, 0].tolist()
    y_pixel = adata.obsm["spatial"].astype(int)[:, 1].tolist()

    # Calculate adjacent matrix
    s = 1
    b = 49
    adj=spg.calculate_adj_matrix(x=x_pixel, y=y_pixel, histology=False)

    spg.prefilter_genes(adata, min_cells=3)  # avoiding all genes are zeros
    spg.prefilter_specialgenes(adata)
    # Normalize and take log for UMI
    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)

    p = 0.5
    # Find the l value given p
    l = spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)

    n_clusters = 20
    r_seed = t_seed = n_seed = 100
    # Search for suitable resolution
    res = spg.search_res(adata, adj, l, n_clusters, start=0.7, step=0.1, tol=5e-3, lr=0.05, max_epochs=20, r_seed=r_seed,
                     t_seed=t_seed, n_seed=n_seed)

    clf = spg.SpaGCN()
    clf.set_l(l)
    # Set seed
    random.seed(r_seed)
    torch.manual_seed(t_seed)
    np.random.seed(n_seed)
    # Run
    clf.train(adata, adj, init_spa=True, init="louvain", res=res, tol=5e-3, lr=0.05, max_epochs=200, n_clusters=n_clusters)
    y_pred, prob = clf.predict()
    adata.obs["pred"] = y_pred
    adata.obs["pred"] = adata.obs["pred"].astype('category')
    # Do cluster refinement(optional)
    # shape="hexagon" for Visium data, "square" for ST data.
    adj_2d = spg.calculate_adj_matrix(x=x_array, y=y_array, histology=False)
    refined_pred = spg.refine(sample_id=adata.obs.index.tolist(), pred=adata.obs["pred"].tolist(), dis=adj_2d, shape="hexagon")
    adata.obs["refined_pred"] = refined_pred
    adata.obs["refined_pred"] = adata.obs["refined_pred"].astype('category')

    adata.write_h5ad("./Result/adata/" + str(dataset) + "_clusters.h5ad")

    adata.obsm['spatial'][:, 1] = -1 * adata.obsm['spatial'][:, 1]

    dpi = 500
    fig, ax = plt.subplots(figsize=(5, 5), dpi=dpi)
    sc.pl.embedding(adata,
                    basis="spatial",
                    color="refined_pred",
                    size=40,
                    show=False,
                    ax=ax)
    plt.savefig('./Result/plot/' + dataset + '_ref_spagcn.png',
                dpi=dpi, bbox_inches='tight')
    plt.close()

