# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 18:13:21 2022

@author: Alicia Mosquera
"""
cellnum = 10


#%%
import json
d = json.load(open("../data/datasets.json","r"))

path = 'C:/Users/34646/Documents/Copenhagen/Clases/TFM/newdata/sperm0006'+str(d['data'][cellnum-1])

if 0 < cellnum < 8:
    d = d['dataset2']
elif 8 < cellnum < 16:
    d = d['dataset1']
   
s = d['s']
thresh = d['thresh']    


#%%
from functions import scale, reshapeToEvenSides

import cv2 as cv
import numpy as np
from skimage import io
import matplotlib.pyplot as plt
from numpy.fft import fft2, ifft2
import matplotlib.ticker as ticker

# Perform segmentation of the whole dataset to run in the script 'neck_motion.py'
I = io.imread(path+'/'+path[55:]+'.tif')
I = reshapeToEvenSides(I)
    
M = np.real(ifft2(scale(fft2(I),s,0,0)))

print(np.quantile(M,0.99)/255)
K = np.zeros(M.shape)
K[M>255*thresh] = 1

io.imsave(path+'/segmentation/'+path[55:]+'_segmentation_head_gaussian.tif', K)


#%%
# EXAMPLE

# Read the manual segmentation
path += '/segmentation/tuning/slice_t64_z1'
planenum = 2 # ['z', 'y', 'x'][planenum]
seg_manual = io.imread(path+'_segmentManual.tiff', as_gray=True)
seg_manual[seg_manual>0] = 1

# Read original image
I = io.imread(path+'.tiff', as_gray=True)
I = reshapeToEvenSides(I)
seg_manual = reshapeToEvenSides(seg_manual)

# Gaussian smoothing through scale-space
M = np.real(ifft2(scale(fft2(I),s,0,0)))

# Thresholding the image
print(np.quantile(M,0.99)/255)
K = np.zeros(M.shape)
K[M>255*thresh] = 1


plt.figure()
plt.suptitle(path[93:])

ax1 = plt.subplot(221)
ax1.imshow(I, cmap='gray', vmin=0, vmax=255)
ax1.set_title('Original image')
ax1.set_xlabel(['x', 'z', 'y'][planenum])
ax1.set_ylabel(['y', 'x', 'z'][planenum])
ax1.xaxis.set_major_locator(ticker.NullLocator())
ax1.yaxis.set_major_locator(ticker.NullLocator())

ax2 = plt.subplot(222)
ax2.imshow(seg_manual, cmap='gray', vmin=0, vmax=1)
ax2.set_title('Manual segmentation')
ax2.set_xlabel(['x', 'z', 'y'][planenum])
ax2.set_ylabel(['y', 'x', 'z'][planenum])
ax2.xaxis.set_major_locator(ticker.NullLocator())
ax2.yaxis.set_major_locator(ticker.NullLocator())

ax3 = plt.subplot(223)
ax3.imshow(M, cmap='gray', vmin=0, vmax=255)
ax3.set_title('Gaussian smoothing')
ax3.set_xlabel(['x', 'z', 'y'][planenum])
ax3.set_ylabel(['y', 'x', 'z'][planenum])
ax3.xaxis.set_major_locator(ticker.NullLocator())
ax3.yaxis.set_major_locator(ticker.NullLocator())

ax4 = plt.subplot(224)
ax4.imshow(K, cmap='gray', vmin=0, vmax=1)
ax4.set_title('Thresholding segmentation')
ax4.set_xlabel(['x', 'z', 'y'][planenum])
ax4.set_ylabel(['y', 'x', 'z'][planenum])
ax4.xaxis.set_major_locator(ticker.NullLocator())
ax4.yaxis.set_major_locator(ticker.NullLocator())

# Save thresholding segmentation
cv.imwrite(path+'_segmentThresh.tiff', K)