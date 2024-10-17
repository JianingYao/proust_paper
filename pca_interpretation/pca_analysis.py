import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.decomposition import PCA
import os
import pandas as pd

datasets = ['151509', '151674']

# Define colors for each cluster
color_fine = {
    "2": (0.4980392156862745, 0.788235294117647, 0.4980392156862745),
    "7": (0.7450980392156863, 0.6823529411764706, 0.8313725490196079),
    "3": (0.9921568627450981, 0.7529411764705882, 0.5254901960784314),
    "4": (0.4, 0.4, 0.4),
    "1": (0.2196078431372549, 0.4235294117647059, 0.6901960784313725),
    "6": (0.9411764705882353, 0.00784313725490196, 0.4980392156862745),
    "5": (0.7490196078431373, 0.3568627450980392, 0.09019607843137253)
}

num_pc_pairs = 5

for dataset in datasets:
    print(f"################################ IF-{dataset} #####################################")
    # Read the dataset
    adata = sc.read(f'/fastscratch/myscratch/jyao/adata_final/proust/{dataset}_clusters.h5ad')
    latent_rep = adata.obsm['profile_gene_img']
    
    pca = PCA(n_components=num_pc_pairs + 1)  
    latent_pca = pca.fit_transform(latent_rep)
    adata.obsm['X_pca_custom'] = latent_pca
    
    output_dir = './plot/'
    os.makedirs(output_dir, exist_ok=True)
    pdf_path = f'{output_dir}{dataset}_proust_pca_pairs_img_pca.pdf'
    with PdfPages(pdf_path) as pdf:
        for i in range(1, num_pc_pairs + 1):
            plt.figure(figsize=(5, 5))
            sc.pl.embedding(
                adata,
                basis='X_pca_custom',  
                color="cluster_profile",  
                components=f'{i},{i + 1}',  
                title=f"Proust - {dataset} - PC{i} vs. PC{i + 1}",
                show=False,
                palette=color_fine
            )
            pdf.savefig() 
            plt.close()

    print(f"Saved PCA pairs plots to {pdf_path}")



