import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import os

datasets = ['2720', '6432', '6522', '8667']
markers = ['MOBP', 'KRT17', 'PCP4', 'RORB', 'HPCAL1', 'AQP4']

base_domain_order = ['Layer1', 'Layer2', 'Layer3', 'Layer4', 'Layer5', 'Layer6', 'WM']
group_order = ['Manual annotation', 'Proust', 'GraphST']
color_fine = {
    'Manual annotation': (0.4, 0.76, 0.65),
    'Proust': (0.99, 0.55, 0.38),
    'GraphST': (0.55, 0.63, 0.80)
}

for dataset in datasets:
    print(f"################################ IF-{dataset} #####################################")
    adata = sc.read(f"/fastscratch/myscratch/jyao/adata_final/proust/{dataset}_coded.h5ad")
    adata_graphst = sc.read(f"/fastscratch/myscratch/jyao/adata_final/GraphST/{dataset}_coded.h5ad")
    adata.obs['graphst'] = adata_graphst.obs['domain']
    adata = adata[~pd.isnull(adata.obs['ground_truth'])]
    domain_order = [d for d in base_domain_order if d in adata.obs['ground_truth'].unique()]
    # Create gene expression DataFrame
    gene_expr = pd.DataFrame(
        adata.X.toarray(), 
        columns=adata.var.index, 
        index=adata.obs.index
    ).stack().reset_index()
    gene_expr = gene_expr.rename(columns={"level_1": "gene", 0: "expression", 'level_0': 'spot'})
    # Create DataFrames for each annotation type
    proust = pd.DataFrame({"Proust": adata.obs["cluster_profile"].values}, index=adata.obs.index)
    graphst = pd.DataFrame({"GraphST": adata.obs["graphst"].values}, index=adata.obs.index)
    truth = pd.DataFrame({"Manual annotation": adata.obs["ground_truth"].values}, index=adata.obs.index)
    # Merge the domains into a single DataFrame
    merged_domains = proust.merge(truth, left_index=True, right_index=True)
    merged_domains = merged_domains.merge(graphst, left_index=True, right_index=True)
    merged_domains = merged_domains.reset_index().melt(id_vars=['index'], var_name='Group', value_name='label')
    merged_domains = merged_domains.rename(columns={'index': 'spot'})
    # Merge the gene expression data with the domain annotations
    merged_df = pd.merge(gene_expr, merged_domains, on='spot')
    # plot
    for marker in markers:
        ax = sns.violinplot(
            x="label", 
            y="expression", 
            data=merged_df[merged_df["gene"] == marker], 
            hue='Group', 
            order=domain_order, 
            hue_order=group_order, 
            palette=color_fine, 
            inner=None
        )
        ax.set(xlabel="Spatial domains", ylabel="Gene expression", title=marker)
        ax.get_legend().remove()
        plt.tight_layout() 
        plt.savefig(f'./plot/3groups_{dataset}_{marker}_violin.png', dpi=600)
        plt.close()



