# scanpy tutorial
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import numpy as np
from matplotlib import rcParams

sc.set_figure_params(dpi=500)

# for dataset in datasets:
dataset = '6432'
adata = sc.read('/fastscratch/myscratch/jyao/adata_final/proust/' + str(dataset) + "_coded.h5ad")

# plot individual clusters
color_fine = {'Layer1': (0.4980392156862745, 0.788235294117647, 0.4980392156862745),
              'Layer3': (0.7450980392156863, 0.6823529411764706, 0.8313725490196079),
              'WM': (0.9921568627450981, 0.7529411764705882, 0.5254901960784314), 'Layer6': (1.0, 1.0, 0.6),
              'Layer5': (0.2196078431372549, 0.4235294117647059, 0.6901960784313725),
              'Layer2': (0.9411764705882353, 0.00784313725490196, 0.4980392156862745),
              'Layer4': (0.7490196078431373, 0.3568627450980392, 0.09019607843137253)}
adata.obsm['spatial'][:, 1] = -1 * adata.obsm['spatial'][:, 1]
sc.pl.spatial(adata, img_key="hires", color="cluster_profile", groups=['Layer1'], palette=color_fine, show=False)
plt.savefig('./plot/cluster_heatmap/' + dataset + '_layer1_cluster.png', dpi=600)
sc.pl.spatial(adata, img_key="hires", color="cluster_profile", groups=['Layer2'], palette=color_fine, show=False)
plt.savefig('./plot/cluster_heatmap/' + dataset + '_layer2_cluster.png', dpi=600)
sc.pl.spatial(adata, img_key="hires", color="cluster_profile", groups=['Layer3'], palette=color_fine, show=False)
plt.savefig('./plot/cluster_heatmap/' + dataset + '_layer3_cluster.png', dpi=600)
sc.pl.spatial(adata, img_key="hires", color="cluster_profile", groups=['Layer4'], palette=color_fine, show=False)
plt.savefig('./plot/cluster_heatmap/' + dataset + '_layer4_cluster.png', dpi=600)
sc.pl.spatial(adata, img_key="hires", color="cluster_profile", groups=['Layer5'], palette=color_fine, show=False)
plt.savefig('./plot/cluster_heatmap/' + dataset + '_layer5_cluster.png', dpi=600)
sc.pl.spatial(adata, img_key="hires", color="cluster_profile", groups=['Layer6'], palette=color_fine, show=False)
plt.savefig('./plot/cluster_heatmap/' + dataset + '_layer6_cluster.png', dpi=600)
sc.pl.spatial(adata, img_key="hires", color="cluster_profile", groups=['WM'], palette=color_fine, show=False)
plt.savefig('./plot/cluster_heatmap/' + dataset + '_WM_cluster.png', dpi=600)