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
# os.environ['R_HOME'] = '/jhpce/shared/jhpce/core/R/4.0.3/lib64/R/'

seed = 1998

clusters = ['7', '7', '7', '7']
samples = ['2720',  '6432', '6522', '8667']

for i in range(len(samples)):
    dataset = samples[i]
    n_clusters = clusters[i]
    print("################################ IF-" + str(dataset) + " #####################################")
    file_fold = '/dcs04/hicks/data/jianing/datasets/proust_datasets/Visium-DLPFC/' + str(dataset)
    adata = sc.read_visium(file_fold, count_file='filtered_feature_bc_matrix.h5', load_images=True)
    adata.var_names_make_unique()
    image = imageio.volread("/dcs04/hicks/data/jianing/datasets/proust_datasets/Visium-DLPFC/" + str(dataset) + "/" + str(dataset) + ".tif")
    image = image[1:6, :, :]
    x_pixel = adata.obsm['spatial'][:, 1].astype(int)
    y_pixel = adata.obsm['spatial'][:, 0].astype(int)
    # Extract img features
    Img_learn(adata, image, device=device)
    print("Finish extracting image features!")
    adata.obsm['rec_img'] = adata.obsm['img_feat']
    # Preprocess gene expression
    prep_gene(adata)
    adata = adata[:, adata.var['highly_variable']]
    adata.obsm['rec_gene'] = adata.X.toarray()[:, ]
    construct_interaction(adata)
    # define radius to limit the number of neighbors during refinement
    radius = 50
    adata = clustering(adata, n_clusters=n_clusters, radius=radius, refinement=True, seed=seed)  # if 'refinement' is set as 'True', the clustering result would be improved.

