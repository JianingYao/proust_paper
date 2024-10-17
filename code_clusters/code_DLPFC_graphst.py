import os
import pandas as pd
import scanpy as sc
from sklearn import metrics
import multiprocessing as mp
import matplotlib.pyplot as plt
import seaborn as sns


datasets = ['151507', '151508', '151509', '151510', '151669', '151670', '151671', '151672', '151673', '151674', '151675', '151676']

for i in range(len(datasets)):
    dataset = datasets[i]
    adata = sc.read('/fastscratch/myscratch/jyao/adata_final/GraphST/' + str(dataset) + "_clusters.h5ad")
    file_fold = '/dcs04/hicks/data/jianing/datasets/proust_datasets/spatialDLPFC/' + str(dataset)
    df_meta = pd.read_csv(file_fold + '/metadata.tsv', sep='\t')
    df_meta_layer = df_meta['layer_guess']
    adata.obs['ground_truth'] = df_meta_layer.values
    adata = adata[~pd.isnull(adata.obs['ground_truth'])]
    ARI = metrics.adjusted_rand_score(adata.obs['domain'], adata.obs['ground_truth'])
    adata.obsm['spatial'][:, 1] = -1 * adata.obsm['spatial'][:, 1]
    color_fine = {"2": (0.4980392156862745, 0.788235294117647, 0.4980392156862745),
                "7": (0.7450980392156863, 0.6823529411764706, 0.8313725490196079),
                "3": (0.9921568627450981, 0.7529411764705882, 0.5254901960784314),
                "4": (0.4, 0.4, 0.4),
                "1": (0.2196078431372549, 0.4235294117647059, 0.6901960784313725),
                "6": (0.9411764705882353, 0.00784313725490196, 0.4980392156862745),
                "5": (0.7490196078431373, 0.3568627450980392, 0.09019607843137253)}
    adata.uns['ARI'] = ARI
    print('Dataset:', dataset)
    print('ARI:', ARI)
    ax=sc.pl.embedding(adata,
                    basis="spatial",
                    color="domain",
                    size=40,
                    palette=color_fine,
                    show=False,
                    title=['ARI=%.4f' % ARI])
    ax.set_aspect('equal', 'box')
    plt.savefig('/users/jyao/projects/proust/raw_result_final/GraphST/coded/' + dataset + '_coded_graphst.png', dpi=600)
    plt.close()