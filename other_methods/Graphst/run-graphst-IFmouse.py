import os
import torch
import pandas as pd
import scanpy as sc
from sklearn import metrics
import multiprocessing as mp
import matplotlib.pyplot as plt

from GraphST import GraphST
from GraphST.utils import clustering

device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
os.environ['R_HOME'] = '/users/jyao/.conda/envs/proust-310/lib/R'

# DLPFC
clusters = ['20'] * 4
datasets = [
     'Visium_CKp25_rep1',
     'Visium_CKp25_rep2',
     'Visium_CKp25_rep3',
     'Visium_CKp25_rep4']

datasets = [
     'Visium_CKp25_rep4']


for i in range(len(datasets)):
    dataset = datasets[i]
    n_clusters = clusters[i]
    print("################################ IF-" + str(dataset) + " #####################################")
    file_fold = '/dcs04/hicks/data/jianing/datasets/proust_datasets/Visium-mouse_DSB_coronal_forebrain_brain/' + str(dataset)

    adata = sc.read_visium(file_fold, count_file='filtered_feature_bc_matrix.h5', load_images=True)
    adata.var_names_make_unique()
    adata.obsm['spatial'] = adata.obsm['spatial'].astype(int)
    adata

    # define model
    model = GraphST.GraphST(adata, device=device)

    # train model
    adata = model.train()

    # set radius to specify the number of neighbors considered during refinement
    radius = 10

    tool = 'mclust'  # mclust, leiden, and louvain

    # clustering
    if tool == 'mclust':
        clustering(adata, n_clusters, radius=radius, method=tool, refinement=True)  # For DLPFC dataset, we use optional refinement step.
    elif tool in ['leiden', 'louvain']:
        clustering(adata, n_clusters, radius=radius, method=tool, start=0.1, end=2.0, increment=0.01, refinement=True)

    # plotting spatial clustering result
    ax = sc.pl.spatial(adata, img_key="hires", color=["domain"], show=False)
    plt.savefig('../Result/GraphST/' + str(dataset) + '_refPred.png', dpi=600)
    plt.close()

    adata.write_h5ad("/fastscratch/myscratch/jyao/proust_datasets/adata/GraphST/" + str(dataset) + "_clusters.h5ad")



