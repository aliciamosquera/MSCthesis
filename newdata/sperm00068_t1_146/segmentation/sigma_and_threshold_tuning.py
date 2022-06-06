# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 16:02:45 2022

@author: Alicia Mosquera
"""
import re
import glob
import cv2 as cv
import numpy as np
from skimage import io
import matplotlib.pyplot as plt
from numpy.fft import fft2, ifft2
from mpl_toolkits.axes_grid1 import make_axes_locatable

import sys
sys.path.insert(0, 'C:/Users/34646/Documents/Copenhagen/Clases/TFM/newdata/sperm00068_t1_146/segmentation')
from functions import *

# Gaussian smoothing through scale-space by choosing two parameters: sigma and threshold.
# Both these parameters will be optimised through finding the best precision and
# recall (sensitivity) of the confusion matrix

#%%
# Read the original 4D image (I)
path = 'C:/Users/34646/Documents/Copenhagen/Clases/TFM/newdata/sperm00068_t1_146'
I = io.imread(path+'/'+path[55:]+'.tif')

# Crop the image in any odd dimension to be able to apply the Gaussian scale-space
if np.any(np.array(I.shape)%2 != 0):
    # cropping out the last rows/columns (rather than the beginning ones) is important to
    # ensure index consistency btw the original image and the cropped one
    if I.shape[0]%2 != 0:
        I = I[:-1,:,:,:]
    if I.shape[1]%2 != 0: #because I have images with z8 and no z0, we remove a dimension from the beginning
        I = I[:,1:,:,:] 
    if I.shape[2]%2 != 0: 
        I = I[:,:,:-1,:]
    if I.shape[3]%2 != 0:
        I = I[:,:,:,:-1]

# Obtain the confusion matrices of the three manually segmented slices in the XY plane
files = glob.glob(path+'/segmentation/tuning/slice_t*_z*_segmentManual.tiff')
sig=2 # 3 for first dataset, 2 for the second
thresh=0.420 # 0.455 for first dataset, 0.420 for the second

prec = []; rec = []
fig = plt.figure()
ax = fig.add_subplot()

M = np.real(ifft2(scale(fft2(I),sig,0,0)))

for count,thresh in enumerate(range(400,500,5)): #sig in range(0,11,1), thresh in range(400,500,5)
    thresh = thresh*0.001
    
    # Obtain the smoothed image (M) and the segmented image after thresholding (K)  
    K = np.zeros(M.shape)
    K[M>255*thresh] = 1
    
    cfmatrix_tot = np.zeros((2,2))
    for f in files:
        t = int(re.findall(r'slice_t(\d+)', f)[0])
        z = int(re.findall(r'z(\d+)', f)[0])-1 #-1 because I removed the first dimension from z
        smooth = K[t][z]
        
        manual = cv.imread(f, 0) # flag=0 reads in grayscale mode
        manual[manual>0] = 1  
        manual = reshapeToEvenSides(manual)
        
        cfmatrix = confusionMatrix(manual, smooth)
        cfmatrix_tot += cfmatrix
    
    # Plot the whole precision vs recall
    precision, recall = paramConfusionMatrix(cfmatrix_tot)
    #p = ax.scatter(x=recall, y=precision, c=sig, vmin=0, vmax=10, s=25, cmap='viridis')
    p = ax.scatter(x=recall, y=precision, c=thresh*255, vmin=0.4*255, vmax=0.5*255, s=25, cmap='viridis')
    
    prec = np.append(prec, precision)
    rec = np.append(rec, recall)
    
    # The intensity for the 99% quantile should be close to the threshold
    print(count+1)
    print(np.quantile(M,0.99)/255)
    print(precision+recall)

ax.plot(rec, prec, color='b', linestyle='dotted', label='$\sigma$='+str(int(sig))+' pixel')

divider = make_axes_locatable(ax)
cbax = divider.append_axes("right", size="5%", pad=0.05)
cbar = fig.colorbar(p, cax=cbax, label='Threshold')

ax.set_xlabel('Recall')
ax.set_ylabel('Precision')
ax.legend()