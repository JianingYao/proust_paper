import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
from matplotlib import rcParams

sc.set_figure_params(dpi=500)

datasets = ['2720',  '6432', '6522', '8667']

markers = ['MOBP', 'MBP', "CCK", "KRT17", "CUX2", "NTNG2", "PCP4", "TRABD2A", "BCL11B", "PVALB", "RORB",
           "ADCYAP1", "ENC1", "HPCAL1", "FABP7", "AQP4", "RELN", "SNAP25"]

base_cluster_order = ['Layer1', 'Layer2', 'Layer3', 'Layer4', 'Layer5', 'Layer6', 'WM']

### Proust
for i in range(len(datasets)):
    dataset = datasets[i]
    print("################################ IF-" + str(dataset) + " #####################################")
    adata = sc.read("/fastscratch/myscratch/jyao/adata_final/proust/" + str(dataset) + "_coded.h5ad")
    if dataset == '8667':
        adata = adata[adata.obs['cluster_profile'] != 'WM']
    cluster_order = [layer for layer in base_cluster_order if layer in adata.obs['cluster_profile'].unique()]
    adata.obs['cluster_profile'] = pd.Categorical(adata.obs['cluster_profile'], categories=cluster_order, ordered=True)
    # plot
    ax = sc.pl.stacked_violin(
        adata,
        markers,
        groupby='cluster_profile',
        swap_axes=True,
        use_raw=False,
        show=False,
        figsize=(5, 7)
    )
    plt.savefig('./plot/' + dataset + '_proust_stacked_violin.png', dpi=600)
    plt.close()


### GraphST
for i in range(len(datasets)):
    dataset = datasets[i]
    print("################################ IF-" + str(dataset) + " #####################################")
    adata = sc.read("/fastscratch/myscratch/jyao/adata_final/GraphST/" + str(dataset) + "_coded.h5ad")
    if dataset == '6522':
        adata = adata[adata.obs['domain'] != 'Layer2']
    if dataset == '8667':
        adata = adata[adata.obs['domain'] != 'WM']
    cluster_order = [layer for layer in base_cluster_order if layer in adata.obs['domain'].unique()]
    adata.obs['domain'] = pd.Categorical(adata.obs['domain'], categories=cluster_order, ordered=True)
    # plot
    ax = sc.pl.stacked_violin(
        adata,
        markers,
        groupby='domain',
        swap_axes=True,
        use_raw=False,
        show=False,
        figsize=(5, 7)
    )
    plt.savefig('./plot/' + dataset + '_graphst_stacked_violin.png', dpi=600)
    plt.close()