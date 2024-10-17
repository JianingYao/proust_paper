import scanpy as sc
import os
import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score
from sklearn.decomposition import PCA

datasets = [
     'Visium_CKp25_rep1',
     'Visium_CKp25_rep2',
     'Visium_CKp25_rep3',
     'Visium_CKp25_rep4']


######################### Proust #########################
file_fold = '/fastscratch/myscratch/jyao/adata_final/proust'
for dataset in datasets:
    file_path = os.path.join(file_fold, f"{dataset}_clusters.h5ad")
    adata = sc.read_h5ad(file_path)
    X = adata.obsm['profile_gene_img']
    labels = adata.obs['cluster_profile']
    score = silhouette_score(X, labels)
    print(f"Proust - {dataset}: Silhouette Score = {score}")



######################### GraphST #########################
file_fold = '/fastscratch/myscratch/jyao/adata_final/GraphST'
for dataset in datasets:
    file_path = os.path.join(file_fold, f"{dataset}_clusters.h5ad")
    adata = sc.read_h5ad(file_path)
    X = adata.obsm['emb_pca']
    labels = adata.obs['domain']
    score = silhouette_score(X, labels)
    print(f"GraphST - {dataset}: Silhouette Score = {score}")



######################### SpaGCN #########################
pc = 50
pca = PCA(n_components=pc)
file_fold = '/fastscratch/myscratch/jyao/adata_final/SpaGCN'
for dataset in datasets:
    file_path = os.path.join(file_fold, f"{dataset}_clusters.h5ad")
    adata = sc.read_h5ad(file_path)
    X = adata.X.toarray()
    X = pca.fit_transform(X)
    labels = adata.obs['refined_pred']
    score = silhouette_score(X, labels)
    print(f"SpaGCN - {dataset}: Silhouette Score = {score}")



######################### STAGATE #########################
file_fold = '/fastscratch/myscratch/jyao/adata_final/STAGATE'
for dataset in datasets:
    file_path = os.path.join(file_fold, f"{dataset}_clusters.h5ad")
    adata = sc.read_h5ad(file_path)
    X = adata.obsm['STAGATE']
    labels = adata.obs['mclust']
    score = silhouette_score(X, labels)
    print(f"STAGATE - {dataset}: Silhouette Score = {score}")


