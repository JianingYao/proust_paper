import scanpy as sc
import matplotlib.pyplot as plt

dataset ='6432'

color_fine = {'Layer1': (0.4980392156862745, 0.788235294117647, 0.4980392156862745), 
                'Layer2': (0.9411764705882353, 0.00784313725490196, 0.4980392156862745), 
                'Layer3': (0.7450980392156863, 0.6823529411764706, 0.8313725490196079), 
                'Layer4': (0.7490196078431373, 0.3568627450980392, 0.09019607843137253),
                'Layer5': (0.2196078431372549, 0.4235294117647059, 0.6901960784313725),
                'Layer6': (1.0, 1.0, 0.6),
                'WM': (0.9921568627450981, 0.7529411764705882, 0.5254901960784314)}

## Proust
print("################################ IF-" + str(dataset) + " #####################################")
adata = sc.read('/fastscratch/myscratch/jyao/adata_final/proust/' + str(dataset) + "_coded.h5ad")
sc.pp.neighbors(adata, use_rep='profile_gene_img', n_neighbors=10)
sc.tl.umap(adata, min_dist=0.15)
ax = plt.rcParams["figure.figsize"] = (5, 5)
sc.pl.umap(adata, color=["cluster_profile"], title=["Proust"], show=False, palette=color_fine)
plt.savefig('./plot/' + dataset + "_proust_umap.png", dpi=600)
plt.close()