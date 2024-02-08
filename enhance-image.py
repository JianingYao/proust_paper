import imageio
import cv2
import numpy as np
import scanpy as sc
import math

# mouse dataset
datasets = [
     'Visium_CKp25_rep1',
     'Visium_CKp25_rep2',
     'Visium_CKp25_rep3',
     'Visium_CKp25_rep4']

for dataset in datasets:
    image = imageio.volread(str(dataset) + ".tif")
    image = image[2, :]
    flattened = image.flatten()
    mean = np.mean(flattened)
    std = np.std(flattened)
    for i in range(len(flattened)):
        if flattened[i] < (mean + 6*std):
            flattened[i] = 0
    enhanced = flattened.reshape(image.shape)

    # Define block size
    block_size = (10, 10)
    for i in range(0, enhanced.shape[0], block_size[0]):
        for j in range(0, enhanced.shape[1], block_size[1]):
            block_max = np.max(enhanced[i:i + block_size[0], j:j + block_size[1]])
            enhanced[i:i + block_size[0], j:j + block_size[1]] = block_max

    cv2.imwrite(str(dataset) + "_enhanced.png", enhanced)
    


# AD dataset
protein = 'pTau'

datasets = [
'VIFAD1_V10A27004_D1_Br3880',
'VIFAD2_V10A27106_B1_Br3854',
'VIFAD2_V10A27106_C1_Br3873',
'VIFAD2_V10A27106_D1_Br3880',
'VIFAD3_V10T31036_B1_Br3854',
'VIFAD3_V10T31036_C1_Br3873',
'VIFAD3_V10T31036_D1_Br3880']

for dataset in datasets:
    image = cv2.imread(str(dataset) + "/images/VistoSeg/" + str(dataset) + '_' + protein + ".png", 0)
    flattened = image.flatten()
    mean = np.mean(flattened)
    std = np.std(flattened)
    for i in range(len(flattened)):
        if flattened[i] < (mean + 6*std):
            flattened[i] = 0
    enhanced = flattened.reshape(image.shape)

    # Define block size
    block_size = (25, 25)
    for i in range(0, enhanced.shape[0], block_size[0]):
        for j in range(0, enhanced.shape[1], block_size[1]):
            block_max = np.max(enhanced[i:i + block_size[0], j:j + block_size[1]])
            enhanced[i:i + block_size[0], j:j + block_size[1]] = block_max

    cv2.imwrite(str(dataset) + "/images/VistoSeg/" + str(dataset) + '_' + protein + "_enhanced.png", enhanced)