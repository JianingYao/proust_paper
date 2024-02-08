import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
import os


# Load data

clusters = ['7']
datasets = ['6432']

markers = ['MOBP', 'KRT17', 'PCP4', 'RORB', 'HPCAL1', 'AQP4']
# markers = ['PCP4']
domain_order = ['Layer1', 'Layer2', 'Layer3', 'Layer4', 'Layer5', 'Layer6', 'WM']
# domain_order = ['Layer1', 'Layer2', 'Layer3', 'Layer4', 'Layer5', 'Layer6']
group_order = ['Manual annotation', 'Proust']
rgb_values = sns.color_palette("Set2", 3)
color_fine = {'Manual annotation': (0.4, 0.7607843137254902, 0.6470588235294118),
              'Proust': (0.9882352941176471, 0.5529411764705883, 0.3843137254901961),
              'GraphST': (0.5529411764705883, 0.6274509803921569, 0.796078431372549)}

for i in range(len(datasets)):
    dataset = datasets[i]
    print("################################ IF-" + str(dataset) + " #####################################")
    file_fold = '/Users/jianingyao/Desktop/Research/Biostatistics_JHU/Stephanie/Data/Visium-DLPFC/' + str(dataset)
    adata = sc.read('../adatas/proust/' + str(dataset) + "_coded.h5ad")
    #######
    # adata_graphst = sc.read('../adatas/graphst/' + str(dataset) + "_coded.h5ad")
    # adata.obs['graphst'] = adata_graphst.obs['label']
    #######
    adata = adata[~pd.isnull(adata.obs['ground_truth'])]
    # adata = adata[adata.obs['label'] != 'WM']


    # Create violin plot
    adata.n = adata.X.toarray()[:, ]
    gene_expr = pd.DataFrame(adata.n, columns=adata.var.index, index=adata.obs.index).stack().reset_index()
    gene_expr = gene_expr.rename(columns={"level_1": "gene", 0: "expression", 'level_0': 'spot'})
    proust = pd.DataFrame({"Proust": adata.obs["label"].values.astype('category')}, index=adata.obs.index)
    #######
    # graphst = pd.DataFrame({"GraphST": adata.obs["graphst"].values.astype('category')}, index=adata.obs.index)
    #######
    truth = pd.DataFrame({"Manual annotation": adata.obs["ground_truth"].values.astype('category')}, index=adata.obs.index)

    merged_domains = proust.merge(truth, left_index=True, right_index=True)
    #######
    # merged_domains = merged_domains.merge(graphst, left_index=True, right_index=True)
    #######
    merged_domains = merged_domains.reset_index().melt(id_vars=['index'], var_name='Group', value_name='label')
    merged_domains = merged_domains.rename(columns={'index': 'spot'})

    merged_df = pd.merge(gene_expr, merged_domains, on='spot')

    # Create a violin plot
    for marker in markers:
        ax = sns.violinplot(x="label", y="expression", data=merged_df[merged_df["gene"] == marker], hue='Group', split=True,
                       order=domain_order, hue_order=group_order, palette=color_fine, inner=None)
        ax.set(xlabel="Spatial domains", ylabel="Gene expression")
        ax.get_legend().remove()
        plt.title(marker)
        plt.savefig('../plots/violin/' + dataset + '_' + marker + '_violin.png', dpi=600)
        plt.close()








