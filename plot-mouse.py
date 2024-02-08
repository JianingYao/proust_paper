import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
import os
from sklearn import metrics

# Load data
dataset = 'Visium_CKp25_rep3'


markers = ['Cst7', 'Lgals3bp', 'Lpl', 'H2-D1']

domain_order = [str(i) for i in range(1, 21)]

print("################################ IF-" + str(dataset) + " #####################################")
adata = sc.read("../Result-raw/" + str(dataset) + "_clusters.h5ad")
adata.n = adata.X.toarray()[:, ]

gene_expr = pd.DataFrame(adata.n, columns=adata.var.index, index=adata.obs.index).stack().reset_index()
gene_expr = gene_expr.rename(columns={"level_1": "gene", 0: "expression", 'level_0': 'spot'})
cluster_profile_df = pd.DataFrame({"cluster_profile": adata.obs["label"].values.astype('category')}, index=adata.obs.index)

merged_domains = cluster_profile_df.reset_index().rename(columns={"index": "spot"})

merged_df = pd.merge(gene_expr, merged_domains, on='spot')

for marker in markers:
     df_plot = merged_df[merged_df['gene'] == marker]
     highest_median = df_plot.groupby("cluster_profile")["expression"].median().max()
     highest_index = df_plot.groupby("cluster_profile")["expression"].median().idxmax()
     hpc = ['9', '11', '4', '5']
     # color_fine = {g: "orange" if g == highest_index else "lightgrey" for g in domain_order}
     color_fine = {g: "orange" if g in hpc else "lightgrey" for g in domain_order}
     ax = sns.boxplot(x=df_plot["cluster_profile"], y=df_plot["expression"], palette=color_fine)
     ax.set_xticklabels(domain_order)
     ax.set_xlabel("Spatial domains")
     ax.set_ylabel("Expression")
     ax.set_title(marker)
     plt.savefig('../plots/box/' + dataset + '_' + marker + '_cluster-profile_box.png', dpi=600)
     plt.close()


adata.uns['log1p'] = {'base': None}
for gene in markers:
    sc.pl.spatial(adata, img_key="lowres", color=gene, show=False, color_map = 'magma')
    plt.savefig('../plots/heatmapGenes/' + dataset + '_' + gene + '_.png', dpi=600)






