import pandas as pd
import scanpy as sc
from sklearn import metrics
import multiprocessing as mp
import matplotlib.pyplot as plt
import seaborn as sns

samples = ['2720', '6432', '6522', '8667']

ari_results = {'Dataset': samples, 'PCA': [], 'Proust without CSL': [], 'Proust': []}

def calculate_ari(data_type, folder_name):
    for dataset in samples:
        adata = sc.read(f'/fastscratch/myscratch/jyao/adata_final/{folder_name}/{dataset}_clusters.h5ad')
        file_fold = f'/dcs04/hicks/data/jianing/datasets/proust_datasets/Visium-DLPFC/{dataset}'
        df_meta = pd.read_csv(file_fold + '/metadata_sorted.csv')
        df_meta_layer = df_meta['label']
        adata.obs['ground_truth'] = df_meta_layer.values
        adata = adata[~pd.isnull(adata.obs['ground_truth'])]
        ARI = metrics.adjusted_rand_score(adata.obs['cluster_profile'], adata.obs['ground_truth'])
        ari_results[data_type].append(ARI)
        # print(f"{data_type}, Dataset: {dataset}, ARI: {ARI}")

calculate_ari('PCA', 'proust_pca')
calculate_ari('Proust without CSL', 'proust_noCSL')
calculate_ari('Proust', 'proust')

df_ari = pd.DataFrame(ari_results)
df_ari.to_csv('ari_results.csv', index=False)

# Plot for ARI comparison
plt.figure(figsize=(8, 6))
df_ari_melted = df_ari.melt(id_vars='Dataset', value_vars=['PCA', 'Proust without CSL', 'Proust'],
                            var_name='Method', value_name='ARI')

sns.barplot(data=df_ari_melted, x='Dataset', y='ARI', hue='Method')
plt.xlabel('Dataset')
plt.ylabel('Adjusted Rand Index (ARI)')
plt.legend(title='Method', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# Save the plot
plt.savefig('./plot/ari_comparison_barplot.png')
plt.close()
