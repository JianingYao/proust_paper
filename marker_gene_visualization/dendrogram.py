import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import numpy as np
from matplotlib import rcParams

sc.set_figure_params(dpi=500)

dataset = '6432'
adata = sc.read('/fastscratch/myscratch/jyao/adata_final/proust/' + str(dataset) + "_coded.h5ad")

adata.uns['log1p'] = {'base': None}
adata = adata[~pd.isnull(adata.obs['ground_truth'])]
adata = adata[:, adata.var['highly_variable']]


##### heatmap
sc.tl.rank_genes_groups(adata, "cluster_profile", method="wilcoxon", inplace=True)
ax = sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, show_gene_labels=True, figsize=(8, 7), show=False)
plt.savefig('./plot/cluster_gene_heatmap/' + dataset + '_proust_heatmap.png', dpi=600)

sc.tl.rank_genes_groups(adata, "ground_truth", method="wilcoxon", inplace=True)
ax = sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, show_gene_labels=True, figsize=(8, 7), show=False)
plt.savefig('./plot/cluster_gene_heatmap/' + dataset + '_annotation_heatmap.png', dpi=600)









