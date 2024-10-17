import pandas as pd
import scanpy as sc
import torch
from sklearn import metrics
import multiprocessing as mp
import pickle
import matplotlib.pyplot as plt
import cv2
import imageio
import numpy as np
from numpy import newaxis
import sys
import os

sys.path.append('/users/jyao/projects/proust')
from proust.Train import *
from proust.cluster import *
from proust.prep import *


device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
os.environ['R_HOME'] = '/users/jyao/.conda/envs/proust-310/lib/R'

seed = 50

clusters = ['7'] * 7
samples = [
'VIFAD1_V10A27004_D1_Br3880',
'VIFAD2_V10A27106_B1_Br3854',
'VIFAD2_V10A27106_C1_Br3873',
'VIFAD2_V10A27106_D1_Br3880',
'VIFAD3_V10T31036_B1_Br3854',
'VIFAD3_V10T31036_C1_Br3873',
'VIFAD3_V10T31036_D1_Br3880']


for i in range(len(samples)):
    dataset = samples[i]
    n_clusters = clusters[i]
    print('################################ IF-' + str(dataset) + ' #####################################')
    file_fold = '/dcs04/hicks/data/jianing/datasets/proust_datasets/Visium-AD/' + str(dataset)
    adata = sc.read_visium(file_fold, count_file='filtered_feature_bc_matrix.h5', load_images=True)
    adata.var_names_make_unique()
    ####################### 5 channels #######################
    img_raw = imageio.volread('/dcs04/hicks/data/jianing/datasets/proust_datasets/Visium-AD/' + str(dataset) + '/images/Capture_Areas/' + str(dataset) + '.tif')     # IF images for DAPI, Abeta, pTau, GFAP, MAP2
    image_abeta = cv2.imread('/dcs04/hicks/data/jianing/datasets/proust_datasets/Visium-AD/' + str(dataset) + '/images/VistoSeg/' + str(dataset) + '_Abeta.png', 0)
    image_ptau = cv2.imread('/dcs04/hicks/data/jianing/datasets/proust_datasets/Visium-AD/' + str(dataset) + '/images/VistoSeg/' + str(dataset) + '_pTau.png', 0)
    image_abeta_exp = np.expand_dims(image_abeta, axis=0)
    image_ptau_exp = np.expand_dims(image_ptau, axis=0)
    image = np.concatenate([img_raw[[0, 3, 4], :, :], image_abeta_exp, image_ptau_exp], axis=0)
    ####################### 5 channels #######################
    x_pixel = adata.obsm['spatial'][:, 1].astype(int)
    y_pixel = adata.obsm['spatial'][:, 0].astype(int)
    # Extract img features
    Img_learn(adata, image, device=device)
    # extract_img(image, adata)
    print("Finish extracting image features!")
    model = proust(adata, device=device, random_seed=seed)
    adata = model.train()
    adata.write_h5ad("/fastscratch/myscratch/jyao/proust_datasets/adata/" + str(dataset) + "_results.h5ad")
    adata = sc.read("/fastscratch/myscratch/jyao/proust_datasets/adata/" + str(dataset) + "_results.h5ad")
    # define radius to limit the number of neighbors during refinement
    radius = 10
    adata = clustering(adata, n_clusters=n_clusters, radius=radius, refinement=True, seed=seed)  # if 'refinement' is set as 'True', the clustering result would be improved.

    adata.obsm['spatial'][:, 1] = -1 * adata.obsm['spatial'][:, 1]
    adata.obsm['spatial'] = adata.obsm['spatial'].astype(int)

    dpi = 500
    fig, ax = plt.subplots(figsize=(5, 5), dpi=dpi)
    sc.pl.embedding(adata,
                    basis="spatial",
                    color="cluster_profile",
                    size=40,
                    show=False,
                    ax=ax)
    plt.savefig(f'/users/jyao/projects/proust/Result/{dataset}_refPred.png', dpi=dpi, bbox_inches='tight')
    plt.close()



