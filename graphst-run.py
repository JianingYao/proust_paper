import os
import torch
import pandas as pd
import scanpy as sc
from sklearn import metrics
import multiprocessing as mp
import matplotlib.pyplot as plt

from GraphST import GraphST
from GraphST.utils import clustering

device = torch.device("mps")

os.environ["PYTORCH_ENABLE_MPS_FALLBACK"] = '1'
os.environ['R_HOME'] = '/Library/Frameworks/R.framework/Resources'

# IF-DLPFC
clusters = ['7', '7', '7', '7']
datasets = ['2720',  '6432', '6522', '8667']
datasets = ['2720']

# DLPFC
# clusters = [7, 7, 7, 7, 5, 5, 5, 5, 7, 7, 7, 7]
# datasets = ['151507', '151508', '151509', '151510', '151669', '151670', '151671', '151672', '151673', '151674',
#      '151675', '151676']


datasets = [
     'Visium_CKp25_rep1',
     'Visium_CKp25_rep2',
     'Visium_CKp25_rep3',
     'Visium_CKp25_rep4']

datasets = [
'VIFAD1_V10A27004_D1_Br3880',
'VIFAD2_V10A27106_B1_Br3854',
'VIFAD2_V10A27106_C1_Br3873',
'VIFAD2_V10A27106_D1_Br3880',
'VIFAD3_V10T31036_B1_Br3854',
'VIFAD3_V10T31036_C1_Br3873',
'VIFAD3_V10T31036_D1_Br3880']


for i in range(len(datasets)):
    dataset = datasets[i]
    n_clusters = 9
    print("################################ IF-" + str(dataset) + " #####################################")
    # file_fold = '/Users/jianingyao/Desktop/Research/Biostatistics_JHU/Stephanie/Data/Visium-DLPFC/' + str(dataset)
    # file_fold = '/Users/JianingYao/Desktop/Research/Biostatistics_JHU/Stephanie/Code/DeepST-main/Data/1.DLPFC/' + str(dataset)
    # file_fold = '/Users/jianingyao/Desktop/Research/Biostatistics_JHU/Stephanie/Data/Visium-mouse_DSB_coronal_forebrain_brain/'+ str(dataset)
    file_fold = '/Users/jianingyao/Desktop/Research/Biostatistics_JHU/Stephanie/Data/Visium-AD/' + str(dataset)

    adata = sc.read_visium(file_fold, count_file='filtered_feature_bc_matrix.h5', load_images=True)
    adata.var_names_make_unique()
    adata.obsm['spatial'] = adata.obsm['spatial'].astype(int)
    adata

    # define model
    model = GraphST.GraphST(adata, device=device)

    # train model
    adata = model.train()

    # set radius to specify the number of neighbors considered during refinement
    radius = 50

    tool = 'mclust'  # mclust, leiden, and louvain

    # clustering
    # if tool == 'mclust':
    #     clustering(adata, n_clusters, radius=radius, method=tool, refinement=False)  # For DLPFC dataset, we use optional refinement step.
    # elif tool in ['leiden', 'louvain']:
    #     clustering(adata, n_clusters, radius=radius, method=tool, start=0.1, end=2.0, increment=0.01, refinement=False)
    #
    # # plotting spatial clustering result`
    # ax = sc.pl.spatial(adata, img_key="hires", color=["domain"], show=False)
    # plt.savefig('./result/' + str(dataset) + '_rawPred.png', dpi=600)
    # plt.close()

    if tool == 'mclust':
        clustering(adata, n_clusters, radius=radius, method=tool, refinement=True)  # For DLPFC dataset, we use optional refinement step.
    elif tool in ['leiden', 'louvain']:
        clustering(adata, n_clusters, radius=radius, method=tool, start=0.1, end=2.0, increment=0.01, refinement=False)

    # plotting spatial clustering result`
    # ax = sc.pl.spatial(adata, img_key="hires", color=["domain"], show=False)
    # plt.savefig('./GraphST/result/' + str(dataset) + '_refPred.png', dpi=600)
    # plt.close()

    adata.obsm['spatial'][:, 1] = -1 * adata.obsm['spatial'][:, 1]
    dpi = 500
    fig, ax = plt.subplots(figsize=(5, 5), dpi=dpi)
    sc.pl.embedding(adata,
                    basis="spatial",
                    color="domain",
                    size=40,
                    show=False,
                    ax=ax)
    plt.savefig('./result/' + dataset + '_ref_graphst.png',
                dpi=dpi, bbox_inches='tight')
    plt.close()

    adata.write_h5ad("../../adatas/graphst/" + str(dataset) + "_clusters.h5ad")


