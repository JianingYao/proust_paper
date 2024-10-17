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
markers = [marker.upper() for marker in markers]

domain_order = [str(i) for i in range(0, 19)]

print("################################ IF-" + str(dataset) + " #####################################")
adata = sc.read("/fastscratch/myscratch/jyao/adata_final/SpaGCN/" + str(dataset) + "_clusters.h5ad")

# Create violin plot
adata.n = adata.X.toarray()[:, ]
gene_expr = pd.DataFrame(adata.n, columns=adata.var.index, index=adata.obs.index).stack().reset_index()
gene_expr = gene_expr.rename(columns={"level_1": "gene", 0: "expression", 'level_0': 'spot'})
spagcn_df = pd.DataFrame({"spagcn": adata.obs["refined_pred"].values.astype('category')}, index=adata.obs.index)
spagcn_df = spagcn_df.reset_index().rename(columns={"index": "spot"})
merged_df = pd.merge(gene_expr, spagcn_df, on='spot')

# Create a violin plot
for marker in markers:
    df_plot = merged_df[merged_df['gene'] == marker]
    df_plot["spagcn"] = df_plot["spagcn"].astype(str) 
    hpc = ['3', '13']
    color_fine = {g: "orange" if g in hpc else "lightgrey" for g in domain_order}
    ax = sns.boxplot(
        x=df_plot["spagcn"], 
        y=df_plot["expression"], 
        order=domain_order,  
        palette=[color_fine.get(g, "lightgrey") for g in domain_order] 
    )
    ax.set_xticklabels(domain_order)
    ax.set_xlabel("Spatial domains by SpaGCN")
    ax.set_ylabel("Expression")
    ax.set_title(marker)
    plt.savefig('./plot/' + dataset + '_spagcn_' + marker + '.png', dpi=600)
    plt.close()






