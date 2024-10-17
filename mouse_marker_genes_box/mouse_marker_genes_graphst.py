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

domain_order = [str(i) for i in range(1, 21) if i != 14]

print("################################ IF-" + str(dataset) + " #####################################")
adata = sc.read("/fastscratch/myscratch/jyao/adata_final/GraphST/" + str(dataset) + "_clusters.h5ad")

# Create violin plot
adata.n = adata.X.toarray()[:, ]
gene_expr = pd.DataFrame(adata.n, columns=adata.var.index, index=adata.obs.index).stack().reset_index()
gene_expr = gene_expr.rename(columns={"level_1": "gene", 0: "expression", 'level_0': 'spot'})
graphst_df = pd.DataFrame({"graphst": adata.obs["domain"].values.astype('category')}, index=adata.obs.index)
graphst_df = graphst_df.reset_index().rename(columns={"index": "spot"})
merged_df = pd.merge(gene_expr, graphst_df, on='spot')

# Create a violin plot
for marker in markers:
    df_plot = merged_df[merged_df['gene'] == marker]
    highest_median = df_plot.groupby("graphst")["expression"].median().max()
    highest_index = df_plot.groupby("graphst")["expression"].median().idxmax()
    hpc = ['19', '18']
    # color_fine = {g: "orange" if g == highest_index else "lightgrey" for g in domain_order}
    color_fine = {g: "orange" if g in hpc else "lightgrey" for g in domain_order}
    ax = sns.boxplot(x=df_plot["graphst"], y=df_plot["expression"], palette=color_fine)
    ax.set_xticklabels(domain_order)
    ax.set_xlabel("Spatial domains by GraphST")
    ax.set_ylabel("Expression")
    ax.set_title(marker)
    plt.savefig('./plot/' + dataset + '_graphst_' + marker + '.png', dpi=600)
    plt.close()




