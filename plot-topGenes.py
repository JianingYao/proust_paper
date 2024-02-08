# scanpy tutorial
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import numpy as np
from matplotlib import rcParams

sc.set_figure_params(dpi=500)

clusters = ['7', '7', '7', '7']
datasets = ['6432']
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
# adata = adata[adata.obs['label'] != 'WM']

adata_org = sc.read_visium(file_fold, count_file='filtered_feature_bc_matrix.h5', load_images=True)
adata_org.var_names_make_unique()
adata_org.obs['ground_truth'] = adata.obs['ground_truth']
adata_org.obs['label'] = adata.obs['label']
prep_gene(adata_org)
adata_org = adata_org[~pd.isnull(adata_org.obs['ground_truth'])]
adata_org = adata_org[:, adata_org.var['highly_variable']]


adata.uns['log1p'] = {'base': None}
adata = adata[~pd.isnull(adata.obs['ground_truth'])]
adata = adata[:, adata.var['highly_variable']]


##### heatmap
sc.tl.rank_genes_groups(adata, "label", method="wilcoxon", inplace=True)
sc.pl.rank_genes_groups(adata, groupby='label', method='wilcoxon', n_genes=200)
plt.show()
ax = sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, show_gene_labels=True, figsize=(8, 7), show=False)
plt.savefig('../plots/heatmapGenes/' + dataset + '_proust_heatmap.png', dpi=600)

sc.tl.rank_genes_groups(adata, "ground_truth", method="wilcoxon", inplace=True)
ax = sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, show_gene_labels=True, figsize=(8, 7), show=False)
plt.savefig('../plots/heatmapGenes/' + dataset + '_annotation_heatmap.png', dpi=600)




sc.tl.dendrogram(adata_org, groupby='label')
sc.tl.rank_genes_groups(adata_org, "label", method="wilcoxon", inplace=True)
n_genes = 5
top_genes = pd.DataFrame(adata_org.uns['rank_genes_groups']['names']).iloc[:n_genes, :].values.flatten()

# gene_expression = pd.DataFrame(adata_org[:, np.unique(top_genes)].X.todense(), columns=np.unique(top_genes))
# for col in gene_expression.columns:
#     gene_expression[col].plot.hist(bins=50)
#     plt.title(col)
#     plt.show()

expression_proust = adata_org[:, top_genes].X.toarray()[:, ]
z_scores = (expression_proust - np.mean(expression_proust, axis=0)) / np.std(expression_proust, axis=0)


# global_min = np.min(z_scores)
# global_max = np.max(z_scores)
# rescaled_z_scores = (z_scores - np.min(z_scores, axis=0)) / (np.max(z_scores, axis=0) - np.min(z_scores, axis=0))
# rescaled_z_scores = rescaled_z_scores * (global_max - global_min) + global_min
# adata_proust = sc.AnnData(X=rescaled_z_scores, obs=adata_org.obs, var=adata_org.var.loc[top_genes])

adata_proust = sc.AnnData(X=z_scores, obs=adata_org.obs, var=adata_org.var.loc[top_genes])

# gene_expression = pd.DataFrame(adata_proust.X)
# for col in gene_expression.columns:
#     gene_expression[col].plot.hist(bins=50)
#     plt.title(col)
#     plt.show()

adata_proust.uns['dendrogram_label'] = adata_org.uns['dendrogram_label']
adata_proust.uns['rank_genes_groups'] = adata_org.uns['rank_genes_groups']
adata_proust.var_names_make_unique()
# sc.pl.heatmap(adata_proust, groupby='label', var_names=top_genes)
sc.pl.rank_genes_groups_heatmap(adata_proust, n_genes=n_genes, show_gene_labels=True, figsize=(8, 7), show=True, cmap="coolwarm", vmin=-1.9, vmax=2, vcenter=0)
plt.show()



sc.tl.dendrogram(adata_org, groupby='ground_truth')
sc.tl.rank_genes_groups(adata_org, "ground_truth", method="wilcoxon", inplace=True)
top_genes = pd.DataFrame(adata_org.uns['rank_genes_groups']['names']).iloc[:n_genes, :].values.flatten()
expression_truth = adata_org[:, top_genes].X.toarray()[:, ]
z_scores = (expression_truth - np.mean(expression_truth, axis=0)) / np.std(expression_truth, axis=0)

# global_min = np.min(z_scores)
# global_max = np.max(z_scores)
# rescaled_z_scores = (z_scores - np.min(z_scores, axis=0)) / (np.max(z_scores, axis=0) - np.min(z_scores, axis=0))
# rescaled_z_scores = rescaled_z_scores * (global_max - global_min) + global_min
# adata_truth = sc.AnnData(X=rescaled_z_scores, obs=adata.obs, var=adata.var.loc[top_genes])

adata_truth = sc.AnnData(X=z_scores, obs=adata_org.obs, var=adata_org.var.loc[top_genes])
adata_truth.uns['dendrogram_ground_truth'] = adata_org.uns['dendrogram_ground_truth']
adata_truth.uns['rank_genes_groups'] = adata_org.uns['rank_genes_groups']
adata_truth.var_names_make_unique()
# sc.pl.heatmap(adata_proust, groupby='label', var_names=top_genes)
sc.pl.rank_genes_groups_heatmap(adata_truth, n_genes=n_genes, show_gene_labels=True, figsize=(8, 7), show=True, cmap="coolwarm", vmin=-1.9, vmax=2)
plt.show()




from sklearn.preprocessing import StandardScaler
gene_feat = adata.X.toarray()[:, ]
scaler = StandardScaler()
feat_new = scaler.fit_transform(gene_feat)
adata.X = feat_new

sc.pl.spatial(adata, img_key="hires", color="KRT17", show=False) ## layer6
plt.savefig('../plots/heatmapGenes/' + dataset + '_layer6_gene.png', dpi=600)
sc.pl.spatial(adata, img_key="hires", color="MOBP", show=False) ##wm
plt.savefig('../plots/heatmapGenes/' + dataset + '_wm_gene.png', dpi=600)
sc.pl.spatial(adata, img_key="hires", color="PCP4", show=False) ## layer5
plt.savefig('../plots/heatmapGenes/' + dataset + '_layer5_gene.png', dpi=600)
sc.pl.spatial(adata, img_key="hires", color="ENC1", show=False) ## layer3
plt.savefig('../plots/heatmapGenes/' + dataset + '_layer3_gene.png', dpi=600)
sc.pl.spatial(adata, img_key="hires", color="NEFM", show=False) ## layer 4
plt.savefig('../plots/heatmapGenes/' + dataset + '_layer4_gene.png', dpi=600)
sc.pl.spatial(adata, img_key="hires", color="SPARC", show=False) ## layer1
plt.savefig('../plots/heatmapGenes/' + dataset + '_layer1_gene.png', dpi=600)
sc.pl.spatial(adata, img_key="hires", color="HPCAL1", show=False) ## layer2
plt.savefig('../plots/heatmapGenes/' + dataset + '_layer2_gene.png', dpi=600)

# plot individual clusters
color_fine = {'Layer1': (0.4980392156862745, 0.788235294117647, 0.4980392156862745),
              'Layer3': (0.7450980392156863, 0.6823529411764706, 0.8313725490196079),
              'WM': (0.9921568627450981, 0.7529411764705882, 0.5254901960784314), 'Layer6': (1.0, 1.0, 0.6),
              'Layer5': (0.2196078431372549, 0.4235294117647059, 0.6901960784313725),
              'Layer2': (0.9411764705882353, 0.00784313725490196, 0.4980392156862745),
              'Layer4': (0.7490196078431373, 0.3568627450980392, 0.09019607843137253)}
sc.pl.spatial(adata, img_key="hires", color="label", groups=['Layer1'], palette=color_fine, show=False)
plt.savefig('../plots/heatmapGenes/' + dataset + '_layer1_cluster.png', dpi=600)
sc.pl.spatial(adata, img_key="hires", color="label", groups=['Layer2'], palette=color_fine, show=False)
plt.savefig('../plots/heatmapGenes/' + dataset + '_layer2_cluster.png', dpi=600)
sc.pl.spatial(adata, img_key="hires", color="label", groups=['Layer3'], palette=color_fine, show=False)
plt.savefig('../plots/heatmapGenes/' + dataset + '_layer3_cluster.png', dpi=600)
sc.pl.spatial(adata, img_key="hires", color="label", groups=['Layer4'], palette=color_fine, show=False)
plt.savefig('../plots/heatmapGenes/' + dataset + '_layer4_cluster.png', dpi=600)
sc.pl.spatial(adata, img_key="hires", color="label", groups=['Layer5'], palette=color_fine, show=False)
plt.savefig('../plots/heatmapGenes/' + dataset + '_layer5_cluster.png', dpi=600)
sc.pl.spatial(adata, img_key="hires", color="label", groups=['Layer6'], palette=color_fine, show=False)
plt.savefig('../plots/heatmapGenes/' + dataset + '_layer6_cluster.png', dpi=600)
sc.pl.spatial(adata, img_key="hires", color="label", groups=['WM'], palette=color_fine, show=False)
plt.savefig('../plots/heatmapGenes/' + dataset + '_WM_cluster.png', dpi=600)




