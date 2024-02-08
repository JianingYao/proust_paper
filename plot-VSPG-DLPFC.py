# scanpy tutorial
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
from matplotlib import rcParams

sc.set_figure_params(dpi=500)

markers = ['MOBP', 'MBP', "CCK", "KRT17", "CUX2", "NTNG2", "PCP4", "TRABD2A", "BCL11B", "PVALB", "RORB",
           "ADCYAP1", "ENC1", "HPCAL1", "FABP7", "AQP4", "RELN", "SNAP25"]

# for dataset in datasets:
dataset = '6432'
print("################################ IF-" + str(dataset) + " #####################################")
adata = sc.read('../adatas/proust/' + str(dataset) + "_coded.h5ad")
# adata = adata[adata.obs['label'] != 'WM']

ax = sc.pl.stacked_violin(adata, markers, groupby='label', swap_axes=True, use_raw=False, show=False, figsize=(5, 7))
plt.savefig('../plots/violin/' + dataset + '_stacked_violin.png', dpi=600)
plt.close()



# scanpy tutorial
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import numpy as np
from matplotlib import rcParams

sc.set_figure_params(dpi=500)

clusters = ['7', '7', '7', '7']
datasets = ['2720',  '6432', '6522', '8667']
domain_order = ['Layer1', 'Layer2', 'Layer3', 'Layer4', 'Layer5', 'Layer6', 'WM']


def prefilter_genes(adata, min_counts=None, max_counts=None,min_cells=10, max_cells=None):
    if min_cells is None and min_counts is None and max_cells is None and max_counts is None:
        raise ValueError('Provide one of min_counts, min_genes, max_counts or max_genes.')
    id_tmp = np.asarray([True]*adata.shape[1],dtype=bool)
    id_tmp = np.logical_and(id_tmp, sc.pp.filter_genes(adata.X, min_cells=min_cells)[0]) if min_cells is not None else id_tmp
    id_tmp = np.logical_and(id_tmp, sc.pp.filter_genes(adata.X, max_cells=max_cells)[0]) if max_cells is not None else id_tmp
    id_tmp = np.logical_and(id_tmp, sc.pp.filter_genes(adata.X, min_counts=min_counts)[0]) if min_counts is not None else id_tmp
    id_tmp = np.logical_and(id_tmp, sc.pp.filter_genes(adata.X, max_counts=max_counts)[0]) if max_counts is not None else id_tmp
    adata._inplace_subset_var(id_tmp)


def prefilter_specialgenes(adata,Gene1Pattern="ERCC",Gene2Pattern="MT-"):
    id_tmp1 = np.asarray([not str(name).startswith(Gene1Pattern) for name in adata.var_names], dtype=bool)
    id_tmp2 = np.asarray([not str(name).startswith(Gene2Pattern) for name in adata.var_names], dtype=bool)
    id_tmp = np.logical_and(id_tmp1,id_tmp2)
    adata._inplace_subset_var(id_tmp)


def prep_gene(adata):
    prefilter_genes(adata, min_cells=3)  # avoiding all genes are zeros
    print("gene number is ", len(adata.var))
    prefilter_specialgenes(adata)
    print("gene number is ", len(adata.var))
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, zero_center=False, max_value=10)


# for dataset in datasets:
dataset = '6432'
print("################################ IF-" + str(dataset) + " #####################################")
file_fold = '/Users/jianingyao/Desktop/Research/Biostatistics_JHU/Stephanie/Data/Visium-DLPFC/' + str(dataset)
adata = sc.read('../adatas/proust/' + str(dataset) + "_coded.h5ad")

adata.uns['log1p'] = {'base': None}
adata = adata[~pd.isnull(adata.obs['ground_truth'])]
adata = adata[:, adata.var['highly_variable']]

groupby = 'label'  

unique_groups = adata.obs[groupby].unique()

# Perform pairwise comparisons for each group against all other groups
pairwise_results = []
for group1 in unique_groups:
    for group2 in unique_groups:
        if group1 != group2:
            comparison_key = f"rank_genes_{group1}_vs_{group2}"
            sc.tl.rank_genes_groups(adata, groupby=groupby, groups=[group1], reference=group2, key_added=comparison_key, method="wilcoxon")
            pairwise_results.append(comparison_key)

# Extract the top 5 differentially expressed genes for each comparison and combine the results
de_genes_per_group = {group: set() for group in unique_groups}
n_top_genes = 5
for result_key in pairwise_results:
    group1 = result_key.split("_")[2]
    top_genes = pd.DataFrame(adata.uns[result_key]['names']).iloc[:n_top_genes].values.flatten()
    de_genes_per_group[group1].update(top_genes)

for group in unique_groups:
    de_genes_per_group[group] = list(de_genes_per_group[group])

for group, genes in de_genes_per_group.items():
    print(f"Top 5 differentially expressed genes for {group}: {', '.join(genes)}")

all_de_genes = set()
for genes in de_genes_per_group.values():
    all_de_genes.update(genes)
all_de_genes = list(all_de_genes)

expression_data = adata[:, all_de_genes].X

adata_de = sc.AnnData(X=expression_data, obs=adata.obs, var=adata.var.loc[all_de_genes])
sc.tl.rank_genes_groups(adata_de, groupby='label', method='wilcoxon')
sc.pl.rank_genes_groups_heatmap(adata_de, groupby='label', n_genes=3)
