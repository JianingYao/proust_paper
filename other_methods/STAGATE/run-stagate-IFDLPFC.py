import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import os
import sys
from sklearn.metrics.cluster import adjusted_rand_score
import STAGATE


# the location of R (used for the mclust clustering)
os.environ['R_HOME'] = '/Library/Frameworks/R.framework/Resources'
os.environ['R_USER'] = '/Users/jianingyao/opt/anaconda3/envs/STAGATE/lib/python3.9/site-packages/rpy2'


n_clusters = 7
l = ['2720', '6432', '6522', '8667']

for i in range(len(l)):
    section_id = l[i]

    input_dir = '/Users/jianingyao/Desktop/Research/Biostatistics_JHU/PhD/Data/proust_datasets/Visium-DLPFC/' + str(section_id)
    adata = sc.read_visium(path=input_dir, count_file='filtered_feature_bc_matrix.h5')
    adata.var_names_make_unique()
    adata

    #Normalization
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # construct the spatial network
    STAGATE.Cal_Spatial_Net(adata, rad_cutoff=150)
    STAGATE.Stats_Spatial_Net(adata)

    # run STAGATE
    # disable default activate eager execution
    import tensorflow as tf
    tf.compat.v1.disable_eager_execution()
    adata = STAGATE.train_STAGATE(adata, alpha=0)

    adata = STAGATE.mclust_R(adata, used_obsm='STAGATE', num_cluster=n_clusters)

    adata.write_h5ad("./Result/adata/" + str(section_id) + "_clusters.h5ad")

    obs_df = adata.obs.dropna()

    # manual annotation
    df_meta = pd.read_csv(input_dir + '/metadata_sorted.csv')
    df_meta_layer = df_meta['label']
    adata.obs['ground_truth'] = df_meta_layer.values
    # filter out NA nodes
    adata = adata[~pd.isnull(adata.obs['ground_truth'])]
    adata.obsm['spatial'][:, 1] = -1 * adata.obsm['spatial'][:, 1]

    color_fine = {7: (0.4980392156862745, 0.788235294117647, 0.4980392156862745),
                    1: (0.7450980392156863, 0.6823529411764706, 0.8313725490196079),
                    2: (0.9921568627450981, 0.7529411764705882, 0.5254901960784314),
                    3: (1.0, 1.0, 0.6),
                    4: (0.2196078431372549, 0.4235294117647059, 0.6901960784313725),
                    5: (0.9411764705882353, 0.00784313725490196, 0.4980392156862745),
                    6: (0.7490196078431373, 0.3568627450980392, 0.09019607843137253)}

    ARI = adjusted_rand_score(adata.obs['mclust'], adata.obs['ground_truth'])
    adata.uns['ARI'] = ARI
    print('Dataset:', section_id)
    print('ARI:', ARI)
    ax = sc.pl.embedding(adata,
                            basis="spatial",
                            color="mclust",
                            size=40,
                            palette=color_fine,
                            show=False,
                            title=['ARI=%.4f' % ARI])
    ax.set_aspect('equal', 'box')
    plt.savefig('./Result/plot/coded_' + section_id + '_ref_stagate.png', dpi=600)
    plt.close()