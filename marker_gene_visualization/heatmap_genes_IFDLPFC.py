import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import numpy as np
from matplotlib import rcParams

sc.set_figure_params(dpi=500)

dataset = '6432'
adata = sc.read('/fastscratch/myscratch/jyao/adata_final/proust/' + str(dataset) + "_coded.h5ad")

adata.obsm['spatial'][:, 1] = -1 * adata.obsm['spatial'][:, 1]

sc.pl.spatial(adata, img_key="hires", color="KRT17", show=False) ## layer6
plt.savefig('./plot/gene_heatmap/' + dataset + '_layer6_gene.png', dpi=600)
sc.pl.spatial(adata, img_key="hires", color="MOBP", show=False) ##wm
plt.savefig('./plot/gene_heatmap/' + dataset + '_wm_gene.png', dpi=600)
sc.pl.spatial(adata, img_key="hires", color="PCP4", show=False) ## layer5
plt.savefig('./plot/gene_heatmap/' + dataset + '_layer5_gene.png', dpi=600)
sc.pl.spatial(adata, img_key="hires", color="ENC1", show=False) ## layer3
plt.savefig('./plot/gene_heatmap/' + dataset + '_layer3_gene.png', dpi=600)
sc.pl.spatial(adata, img_key="hires", color="NEFM", show=False) ## layer 4
plt.savefig('./plot/gene_heatmap/' + dataset + '_layer4_gene.png', dpi=600)
sc.pl.spatial(adata, img_key="hires", color="SPARC", show=False) ## layer1
plt.savefig('./plot/gene_heatmap/' + dataset + '_layer1_gene.png', dpi=600)
sc.pl.spatial(adata, img_key="hires", color="HPCAL1", show=False) ## layer2
plt.savefig('./plot/gene_heatmap/' + dataset + '_layer2_gene.png', dpi=600)





