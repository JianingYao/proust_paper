import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
import os
from sklearn import metrics

dataset = 'Visium_CKp25_rep3'
markers = ['Cst7', 'Lgals3bp', 'Lpl', 'H2-D1']

domain_order = [str(i) for i in range(1, 21)]

print("################################ IF-" + str(dataset) + " #####################################")
adata = sc.read("/fastscratch/myscratch/jyao/adata_final/proust/" + str(dataset) + "_clusters.h5ad")

# Create violin plot
adata.n = adata.X.toarray()[:, ]
gene_expr = pd.DataFrame(adata.n, columns=adata.var.index, index=adata.obs.index).stack().reset_index()
gene_expr = gene_expr.rename(columns={"level_1": "gene", 0: "expression", 'level_0': 'spot'})
cluster_profile_df = pd.DataFrame({"cluster_profile": adata.obs["cluster_profile"].values.astype('category')}, index=adata.obs.index)
merged_domains = cluster_profile_df.reset_index().rename(columns={"index": "spot"})
merged_df = pd.merge(gene_expr, merged_domains, on='spot')

# Create a violin plot
for marker in markers:
    # Profile
    df_plot = merged_df[merged_df['gene'] == marker]
    highest_median = df_plot.groupby("cluster_profile")["expression"].median().max()
    highest_index = df_plot.groupby("cluster_profile")["expression"].median().idxmax()
    hpc = ['19', '15', '12', '18']
    # color_fine = {g: "orange" if g == highest_index else "lightgrey" for g in domain_order}
    color_fine = {g: "orange" if g in hpc else "lightgrey" for g in domain_order}
    ax = sns.boxplot(x=df_plot["cluster_profile"], y=df_plot["expression"], palette=color_fine)
    ax.set_xticklabels(domain_order)
    ax.set_xlabel("Spatial domains by Proust")
    ax.set_ylabel("Expression")
    ax.set_title(marker)
    plt.savefig('./plot/' + dataset + '_proust_profile_' + marker + '.png', dpi=600)
    plt.close()



