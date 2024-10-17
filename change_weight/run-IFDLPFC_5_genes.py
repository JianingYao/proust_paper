import pandas as pd
import scanpy as sc
import torch
from sklearn import metrics
import multiprocessing as mp
import pickle
import matplotlib.pyplot as plt
import cv2
import imageio
import numpy as np
from numpy import newaxis
import sys
import os

sys.path.append('/users/jyao/projects/proust')
from proust.Train import *
from proust.cluster import *
from proust.prep import *


os.environ['R_HOME'] = '/users/jyao/.conda/envs/proust-310/lib/R'
# os.environ['R_HOME'] = '/jhpce/shared/jhpce/core/R/4.0.3/lib64/R/'

seed = 1998

dataset = '2720'

n_clusters = 7
print("################################ IF-" + str(dataset) + " #####################################")
# define radius to limit the number of neighbors during refinement
radius = 50

adata = sc.read("/fastscratch/myscratch/jyao/adata_final/proust/" + str(dataset) + "_results.h5ad")
adata = clustering(adata, n_clusters=n_clusters, gene_pcs=5, image_pcs=5, radius=radius, refinement=True, seed=seed)

adata = sc.read(f"/fastscratch/myscratch/jyao/proust_datasets/adata/{dataset}_clusters.h5ad")
file_fold = '/dcs04/hicks/data/jianing/datasets/proust_datasets/Visium-DLPFC/' + str(dataset)
df_meta = pd.read_csv(file_fold + '/metadata_sorted.csv')
df_meta_layer = df_meta['label']
adata.obs['ground_truth'] = df_meta_layer.values
adata = adata[~pd.isnull(adata.obs['ground_truth'])]
ARI = metrics.adjusted_rand_score(adata.obs['cluster_profile'], adata.obs['ground_truth'])

adata.uns['ARI'] = ARI
print('Dataset:', dataset)
print('ARI:', ARI)
ax=sc.pl.embedding(adata,
                basis="spatial",
                color="cluster_profile",
                size=40,
                show=False,
                title=['ARI=%.4f' % ARI])
ax.set_aspect('equal', 'box')
plt.savefig(f'/users/jyao/projects/proust/Result/{dataset}_ARI.png', dpi=600)
plt.close()


