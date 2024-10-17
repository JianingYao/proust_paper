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
rec_img = adata.obsm['rec_img']
rec_img = adata.obsm['rec_img'].reshape(rec_img.shape[0], -1)
rec_gene = adata.obsm['rec_gene']

pca_gene = PCA(n_components=50)
pca_gene.fit(rec_gene)
cumulative_explained_variance_gene = np.cumsum(pca_gene.explained_variance_ratio_)
pca_img = PCA(n_components=50)
pca_img.fit(rec_img)
cumulative_explained_variance_img = np.cumsum(pca_img.explained_variance_ratio_)

plt.figure(figsize=(12, 5))

# Subplot 1: Cumulative Explained Variance for rec_gene
plt.subplot(1, 2, 1)
plt.plot(
    range(1, 51), 
    cumulative_explained_variance_gene[:50], 
    marker='o', 
    markersize=3,  # Smaller dot size
    linestyle='-', 
    color='blue'
)
plt.title('Transcriptomics', fontsize=12)  # Smaller title font size
plt.xlabel('Number of Principal Components', fontsize=10)  # Smaller x-axis label font size
plt.ylabel('Cumulative Explained Variance', fontsize=10)  # Smaller y-axis label font size
plt.grid(True)

# Subplot 2: Cumulative Explained Variance for rec_img
plt.subplot(1, 2, 2)
plt.plot(
    range(1, 51), 
    cumulative_explained_variance_img[:50], 
    marker='o', 
    markersize=3,  # Smaller dot size
    linestyle='-', 
    color='blue'
)
plt.title('Image', fontsize=12)  # Smaller title font size
plt.xlabel('Number of Principal Components', fontsize=10)  # Smaller x-axis label font size
plt.ylabel('Cumulative Explained Variance', fontsize=10)  # Smaller y-axis label font size
plt.grid(True)

# Adjust layout and save the plot
plt.tight_layout()
plt.savefig(f'./plot/{dataset}_cumulative_explained_variance.png', dpi=300)
plt.close()


