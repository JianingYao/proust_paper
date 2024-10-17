import os, csv, re
import pandas as pd
import numpy as np
import scanpy as sc
import math
from scipy.sparse import issparse
import matplotlib.colors as clr
import matplotlib.pyplot as plt
from sklearn import metrics

datasets = [
     'Visium_CKp25_rep1',
     'Visium_CKp25_rep2',
     'Visium_CKp25_rep3',
     'Visium_CKp25_rep4']

for i in range(len(datasets)):
    dataset = datasets[i]
    print("################################ " + str(dataset) + " #####################################")
    adata = sc.read('/fastscratch/myscratch/jyao/adata_final/GraphST/' + str(dataset) + '_clusters.h5ad')
    adata.obsm['spatial'][:, 1] = -1 * adata.obsm['spatial'][:, 1]
    dpi = 500
    fig, ax = plt.subplots(figsize=(5, 5), dpi=dpi)
    sc.pl.embedding(adata,
                    basis="spatial",
                    color="domain",
                    size=40,
                    show=False,
                    ax=ax)
    plt.savefig('../../../raw_result_final/GraphST/' + dataset + '_ref_graphst.png',
                dpi=dpi, bbox_inches='tight')
    plt.close()


