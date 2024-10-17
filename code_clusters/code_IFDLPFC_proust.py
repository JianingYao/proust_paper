import os
import pandas as pd
import scanpy as sc
from sklearn import metrics
import multiprocessing as mp
import matplotlib.pyplot as plt
import seaborn as sns

datasets = ['2720', '6432', '6522', '8667']
clusters = [7, 7, 7, 7]

# 2720
dataset = '2720'
adata = sc.read('/fastscratch/myscratch/jyao/adata_final/proust/' + str(dataset) + "_clusters.h5ad")
labels = adata.obs['cluster_profile']
df = labels.to_frame(name='cluster_profile')
df['cluster_profile'] = df['cluster_profile'].replace({'1': 'WM'})
df['cluster_profile'] = df['cluster_profile'].replace({'2': 'Layer6'})
df['cluster_profile'] = df['cluster_profile'].replace({'4': 'Layer5'})
df['cluster_profile'] = df['cluster_profile'].replace({'6': 'Layer4'})
df['cluster_profile'] = df['cluster_profile'].replace({'3': 'Layer3'})
df['cluster_profile'] = df['cluster_profile'].replace({'7': 'Layer2'})
df['cluster_profile'] = df['cluster_profile'].replace({'5': 'Layer1'})
adata.obs['cluster_profile'] = df
adata.obs['cluster_profile'] = adata.obs['cluster_profile'].astype('category')

adata.write_h5ad("/fastscratch/myscratch/jyao/adata_final/proust/" + str(dataset) + "_coded.h5ad")

# 6432
dataset = '6432'
adata = sc.read('/fastscratch/myscratch/jyao/adata_final/proust/' + str(dataset) + "_clusters.h5ad")
labels = adata.obs['cluster_profile']
df = labels.to_frame(name='cluster_profile')
df['cluster_profile'] = df['cluster_profile'].replace({'4': 'Layer1'})
df['cluster_profile'] = df['cluster_profile'].replace({'6': 'Layer2'})
df['cluster_profile'] = df['cluster_profile'].replace({'5': 'Layer3'})
df['cluster_profile'] = df['cluster_profile'].replace({'1': 'Layer4'})
df['cluster_profile'] = df['cluster_profile'].replace({'7': 'Layer5'})
df['cluster_profile'] = df['cluster_profile'].replace({'3': 'Layer6'})
df['cluster_profile'] = df['cluster_profile'].replace({'2': 'WM'})
adata.obs['cluster_profile'] = df
adata.obs['cluster_profile'] = adata.obs['cluster_profile'].astype('category')

adata.write_h5ad("/fastscratch/myscratch/jyao/adata_final/proust/" + str(dataset) + "_coded.h5ad")

# 6522
dataset = '6522'
adata = sc.read('/fastscratch/myscratch/jyao/adata_final/proust/' + str(dataset) + "_clusters.h5ad")
labels = adata.obs['cluster_profile']
df = labels.to_frame(name='cluster_profile')
df['cluster_profile'] = df['cluster_profile'].replace({'3': 'Layer1'})
df['cluster_profile'] = df['cluster_profile'].replace({'6': 'Layer2'})
df['cluster_profile'] = df['cluster_profile'].replace({'4': 'Layer3'})
df['cluster_profile'] = df['cluster_profile'].replace({'2': 'Layer4'})
df['cluster_profile'] = df['cluster_profile'].replace({'7': 'Layer5'})
df['cluster_profile'] = df['cluster_profile'].replace({'1': 'Layer6'})
df['cluster_profile'] = df['cluster_profile'].replace({'5': 'WM'})
adata.obs['cluster_profile'] = df
adata.obs['cluster_profile'] = adata.obs['cluster_profile'].astype('category')

adata.write_h5ad("/fastscratch/myscratch/jyao/adata_final/proust/" + str(dataset) + "_coded.h5ad")

# 8667
dataset = '8667'
adata = sc.read('/fastscratch/myscratch/jyao/adata_final/proust/' + str(dataset) + "_clusters.h5ad")
labels = adata.obs['cluster_profile']
df = labels.to_frame(name='cluster_profile')
df['cluster_profile'] = df['cluster_profile'].replace({'6': 'Layer1'})
df['cluster_profile'] = df['cluster_profile'].replace({'3': 'Layer2'})
df['cluster_profile'] = df['cluster_profile'].replace({'5': 'Layer3'})
df['cluster_profile'] = df['cluster_profile'].replace({'1': 'Layer4'})
df['cluster_profile'] = df['cluster_profile'].replace({'2': 'Layer5'})
df['cluster_profile'] = df['cluster_profile'].replace({'4': 'Layer6'})
df['cluster_profile'] = df['cluster_profile'].replace({'7': 'WM'})
adata.obs['cluster_profile'] = df
adata.obs['cluster_profile'] = adata.obs['cluster_profile'].astype('category')

adata.write_h5ad("/fastscratch/myscratch/jyao/adata_final/proust/" + str(dataset) + "_coded.h5ad")


# add ground truth to adata
for i in range(len(datasets)):
    dataset = datasets[i]
    adata = sc.read("/fastscratch/myscratch/jyao/adata_final/proust/" + str(dataset) + "_coded.h5ad")
    file_fold = '/dcs04/hicks/data/jianing/datasets/proust_datasets/Visium-DLPFC/' + str(dataset)
    df_meta = pd.read_csv(file_fold + '/metadata_sorted.csv')
    df_meta_layer = df_meta['label']
    adata.obs['ground_truth'] = df_meta_layer.values
    adata = adata[~pd.isnull(adata.obs['ground_truth'])]
    ARI = metrics.adjusted_rand_score(adata.obs['cluster_profile'], adata.obs['ground_truth'])
    adata.write_h5ad("/fastscratch/myscratch/jyao/adata_final/proust/" + str(dataset) + "_coded.h5ad")
    # plot with ARI
    # adata.obsm['spatial'][:, 1] = -1 * adata.obsm['spatial'][:, 1]
    rgb_values = sns.color_palette("Accent", len(adata.obs['ground_truth'].unique()))
    color_fine = dict(zip(list(adata.obs['ground_truth'].unique()), rgb_values))
    color_fine = {'Layer1': (0.4980392156862745, 0.788235294117647, 0.4980392156862745), 
                  'Layer2': (0.9411764705882353, 0.00784313725490196, 0.4980392156862745), 
                  'Layer3': (0.7450980392156863, 0.6823529411764706, 0.8313725490196079), 
                  'Layer4': (0.7490196078431373, 0.3568627450980392, 0.09019607843137253),
                  'Layer5': (0.2196078431372549, 0.4235294117647059, 0.6901960784313725),
                  'Layer6': (1.0, 1.0, 0.6),
                  'WM': (0.9921568627450981, 0.7529411764705882, 0.5254901960784314)}
    adata.uns['ARI'] = ARI
    print('Dataset:', dataset)
    print('ARI:', ARI)
    ax=sc.pl.embedding(adata,
                    basis="spatial",
                    color="cluster_profile",
                    size=40,
                    palette=color_fine,
                    show=False,
                    title=['ARI=%.4f' % ARI])
    ax.set_aspect('equal', 'box')
    plt.savefig('/users/jyao/projects/proust/raw_result_final/proust/coded/' + dataset + '_coded_proust.png', dpi=600)
    plt.close()