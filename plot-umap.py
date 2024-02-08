import scanpy as sc
import matplotlib.pyplot as plt

datasets = ['151509', '151674']

color_fine = {'Layer1': (0.4980392156862745, 0.788235294117647, 0.4980392156862745),
              'Layer3': (0.7450980392156863, 0.6823529411764706, 0.8313725490196079),
              'WM': (0.9921568627450981, 0.7529411764705882, 0.5254901960784314), 'Layer6': (1.0, 1.0, 0.6),
              'Layer5': (0.2196078431372549, 0.4235294117647059, 0.6901960784313725),
              'Layer2': (0.9411764705882353, 0.00784313725490196, 0.4980392156862745),
              'Layer4': (0.7490196078431373, 0.3568627450980392, 0.09019607843137253)}



for dataset in datasets:
    print("################################ IF-" + str(dataset) + " #####################################")
    adata = sc.read("/Users/jianingyao/Desktop/Research/Biostatistics_JHU/Stephanie/Proust/Thesis/Result-raw/Mouse/Visium_CKp25_rep3_clusters.h5ad")
    # adata = adata[adata.obs['label_refined'] != 'WM']
    sc.pp.neighbors(adata, use_rep='profile', n_neighbors=15)
    sc.tl.umap(adata)
    ax = plt.rcParams["figure.figsize"] = (5, 5)
    sc.pl.umap(adata, color=["label"], title=["Proust"], show=False)
    # sc.pl.umap(adata, color=["ground_truth", "label_refined"], title=['Manual annotation', "Proust"], show=False, palette=color_fine)
    plt.savefig('../plots/umap/' + dataset + "_umap.png", dpi=600)
    plt.close()