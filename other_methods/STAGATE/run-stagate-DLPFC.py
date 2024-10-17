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


clusters = [7, 7, 7, 7, 5, 5, 5, 5, 7, 7, 7, 7]
l = ['151507', '151508', '151509', '151510', '151669', '151670', '151671', '151672', '151673', '151674',
        '151675', '151676']
for i in range(len(l)):
    n_clusters = clusters[i]
    section_id = l[i]

    input_dir = '/Users/jianingyao/Desktop/Research/Biostatistics_JHU/PhD/Data/proust_datasets/spatialDLPFC/' + str(section_id)
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
    df_meta = pd.read_csv(input_dir + '/metadata.tsv', sep='\t')
    df_meta_layer = df_meta['layer_guess']
    adata.obs['ground_truth'] = df_meta_layer.values
    # filter out NA nodes
    adata = adata[~pd.isnull(adata.obs['ground_truth'])]
    adata.obsm['spatial'][:, 1] = -1 * adata.obsm['spatial'][:, 1]

    ARI = adjusted_rand_score(adata.obs['mclust'], adata.obs['ground_truth'])
    adata.uns['ARI'] = ARI
    print('Dataset:', section_id)
    print('ARI:', ARI)
    ax = sc.pl.embedding(adata,
                            basis="spatial",
                            color="mclust",
                            size=40,
                            palette="Accent",
                            show=False,
                            title=['ARI=%.4f' % ARI])
    ax.set_aspect('equal', 'box')
    plt.savefig('./Result/plot/' + section_id + '_ref_stagate.png', dpi=600)
    plt.close()