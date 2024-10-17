import scanpy as sc
import matplotlib.pyplot as plt

datasets = ['151509', '151676']

color_fine = {"2": (0.4980392156862745, 0.788235294117647, 0.4980392156862745),
                "7": (0.7450980392156863, 0.6823529411764706, 0.8313725490196079),
                "3": (0.9921568627450981, 0.7529411764705882, 0.5254901960784314),
                "4": (0.4, 0.4, 0.4),
                "1": (0.2196078431372549, 0.4235294117647059, 0.6901960784313725),
                "6": (0.9411764705882353, 0.00784313725490196, 0.4980392156862745),
                "5": (0.7490196078431373, 0.3568627450980392, 0.09019607843137253)}
## Proust
for dataset in datasets:
    print("################################ IF-" + str(dataset) + " #####################################")
    adata = sc.read('/fastscratch/myscratch/jyao/adata_final/proust/' + str(dataset) + "_clusters.h5ad")
    sc.pp.neighbors(adata, use_rep='profile_gene_img', n_neighbors=10)
    sc.tl.umap(adata)
    ax = plt.rcParams["figure.figsize"] = (5, 5)
    sc.pl.umap(adata, color=["cluster_profile"], title=["Proust"], show=False, palette=color_fine)
    plt.savefig('./plot/' + dataset + "_proust_umap.png", dpi=600)
    plt.close()