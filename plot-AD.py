import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
from matplotlib import rcParams
import seaborn as sns
import numpy as np

datasets = [
'VIFAD1_V10A27004_D1_Br3880',
'VIFAD3_V10T31036_C1_Br3873'
]

markers = ['MOBP', 'MBP', "SNAP25"]

for dataset in datasets:
    print("################################ IF-" + str(dataset) + " #####################################")
    adata = sc.read("../Result-raw/" + str(dataset) + "_clusters.h5ad")

    adata.obs['label'] = adata.obs['label'].astype(str).astype('category')
    sc.pl.stacked_violin(adata, markers, groupby='label', swap_axes=True, use_raw=False, show=False, figsize=(5, 2))
    plt.savefig('../plots/violin/' + dataset + '_stacked_violin.png', dpi=600)
    plt.close()



