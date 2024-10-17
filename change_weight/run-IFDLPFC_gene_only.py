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
gene_pcs_values = list(range(30, 41))
for gene_pcs in gene_pcs_values:
    adata = sc.read("/fastscratch/myscratch/jyao/adata_final/proust/" + str(dataset) + "_results.h5ad")
    clustering(adata, n_clusters=n_clusters, gene_pcs=gene_pcs, radius=radius, refinement=True, seed=seed)

ari_results = []
for gene_pcs in gene_pcs_values:
    adata = sc.read(f"/fastscratch/myscratch/jyao/adata_final/proust_gene_only/{dataset}_gene{gene_pcs}_clusters.h5ad")
    file_fold = '/dcs04/hicks/data/jianing/datasets/proust_datasets/Visium-DLPFC/' + str(dataset)
    df_meta = pd.read_csv(file_fold + '/metadata_sorted.csv')
    df_meta_layer = df_meta['label']
    adata.obs['ground_truth'] = df_meta_layer.values
    adata = adata[~pd.isnull(adata.obs['ground_truth'])]
    ARI = metrics.adjusted_rand_score(adata.obs['cluster_gene'], adata.obs['ground_truth'])
    print(f"Dataset: {dataset}, PC value: {gene_pcs}, ARI: {ARI}")
    ari_results.append({'PC_value': gene_pcs, 'ARI': ARI})

ari_df = pd.DataFrame(ari_results)

# Save the DataFrame to a CSV file
ari_df.to_csv(f'./ari_results_{dataset}_gene_weight.csv', index=False)
ari_df = pd.read_csv(f'./ari_results_{dataset}_gene_weight.csv')
# Plot the ARI values as a line chart
plt.figure(figsize=(8, 5))
plt.plot(ari_df['PC_value'], ari_df['ARI'], marker='o', linestyle='-', color='b', label='Gene only')
plt.axhline(y=0.4943, color='red', linestyle='--', linewidth=1, label='Proust with top 30 PCs of gene')
plt.xticks(ari_df['PC_value'])
plt.title(f"ARI using transcriptomics only across different PC Value for {dataset}")
plt.xlabel("PC Value")
plt.ylabel("ARI")
plt.legend(loc='best')
plt.tight_layout()
# Save the plot as a PNG image
plt.savefig(f'./plot/ari_vs_pc_{dataset}.png', dpi=300)
plt.close()

