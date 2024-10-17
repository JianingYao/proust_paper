import scanpy as sc
import matplotlib.pyplot as plt

dataset = 'Visium_CKp25_rep3'

## Proust
print("################################ IF-" + str(dataset) + " #####################################")
adata = sc.read('/fastscratch/myscratch/jyao/adata_final/proust/' + str(dataset) + "_clusters.h5ad")
sc.pp.neighbors(adata, use_rep='profile_gene_img', n_neighbors=10)
sc.tl.umap(adata, min_dist=0.15)
ax = plt.rcParams["figure.figsize"] = (5, 5)
sc.pl.umap(adata, color=["cluster_profile"], title=["Proust"], show=False)
plt.savefig('./plot/' + dataset + "_proust_umap.png", dpi=600)
plt.close()