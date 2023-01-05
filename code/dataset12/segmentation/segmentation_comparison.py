# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 10:04:33 2022

@author: Alicia Mosquera
"""
import re
import glob
from skimage import io
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import sys
sys.path.insert(0, 'C:/Users/34646/Documents/Copenhagen/Clases/TFM/newdata/sperm00068_t1_146/segmentation')
from functions import plot_zx, plot_yz, confusionMatrix, paramConfusionMatrix, plotConfusionMatrix

#%%
path = 'C:/Users/34646/Documents/Copenhagen/Clases/TFM/newdata/sperm00068_t1_146'
planenum = 0
plane = ['z', 'y', 'x'][planenum]

files = glob.glob(path+'/segmentation/tuning/slice_t*_'+plane+'*.png') 
for i,f in enumerate(files):
    t = int(re.findall(r'slice_t(\d+)', f)[0])
    ax = int(re.findall(plane+'(\d+)', f)[0])
    
    volume_time = io.imread(path+'/segmentation/'+path[55:]+'_segmentation_head.tif')
    if plane=='z':
        imagej = volume_time[t][ax]
    elif plane=='y':
        imagej = plot_zx(volume_time, t, ax)
    else:
        imagej = plot_yz(volume_time, t, ax) 
        
    orig = io.imread(path+'/segmentation/tuning/slice_t'+str(t)+'_'+plane+str(ax)+'.tiff')
    thresh = io.imread(path+'/segmentation/tuning/slice_t'+str(t)+'_'+plane+str(ax)+'_segmentThresh.tiff')
    manual = io.imread(path+'/segmentation/tuning/slice_t'+str(t)+'_'+plane+str(ax)+'_segmentManual.tiff', as_gray=True)
    
    manual[manual>0] = 1 #so the masks are binary
    imagej[imagej>0] = 1
    
    orig = orig[1:,1:] #so they have the same dimensions as the thresholded image
    manual = manual[1:,1:]
    imagej = imagej[1:,1:]
       
    # Visual comparison
    plt.figure(i)
    plt.suptitle(f[93:-4])
    
    ax1 = plt.subplot(321)
    ax1.imshow(orig, cmap='gray', vmin=0, vmax=255)
    ax1.set_title('Original image')
    ax1.set_xlabel(['x', 'z', 'y'][planenum])
    ax1.set_ylabel(['y', 'x', 'z'][planenum])
    ax1.xaxis.set_major_locator(ticker.NullLocator())
    ax1.yaxis.set_major_locator(ticker.NullLocator())
    
    ax2 = plt.subplot(322)
    ax2.imshow(manual, cmap='gray', vmin=0, vmax=1)
    ax2.set_title('Manual')
    ax2.set_xlabel(['x', 'z', 'y'][planenum])
    ax2.set_ylabel(['y', 'x', 'z'][planenum])
    ax2.xaxis.set_major_locator(ticker.NullLocator())
    ax2.yaxis.set_major_locator(ticker.NullLocator())

    ax3 = plt.subplot(323)
    ax3.imshow(thresh, cmap='gray', vmin=0, vmax=1)
    ax3.set_title('Gaussian-smoothing (s=3)')
    ax3.set_xlabel(['x', 'z', 'y'][planenum])
    ax3.set_ylabel(['y', 'x', 'z'][planenum])
    ax3.xaxis.set_major_locator(ticker.NullLocator())
    ax3.yaxis.set_major_locator(ticker.NullLocator())
    
    ax4 = plt.subplot(324)
    ax4.imshow(imagej, cmap='gray', vmin=0, vmax=1)
    ax4.set_title('ImageJ')
    ax4.set_xlabel(['x', 'z', 'y'][planenum])
    ax4.set_ylabel(['y', 'x', 'z'][planenum])
    ax4.xaxis.set_major_locator(ticker.NullLocator())
    ax4.yaxis.set_major_locator(ticker.NullLocator())
    
    # Construct the confusion matrix for gaussian and imagej against the manual as 'ground truth'
    cfmatrix_thresh = confusionMatrix(manual, thresh)
    prec_thresh, recall_thresh = paramConfusionMatrix(cfmatrix_thresh)
    cfmatrix_imagej  = confusionMatrix(manual, imagej)
    prec_imagej, recall_imagej = paramConfusionMatrix(cfmatrix_imagej)
    
    ax5 = plt.subplot(325)
    plotConfusionMatrix(ax5, cfmatrix_thresh, 'Gaussian-smoothing')
    ax6 = plt.subplot(326)
    plotConfusionMatrix(ax6, cfmatrix_imagej, 'ImageJ')