# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 11:36:01 2022

@author: Alicia Mosquera

Python script to manually choose the label from the connected
components analysis in each time frame. 
"""
cellnum = 10
# Type of segmentation to analyze (Gaussian/ImageJ)
Gaussian = True


#%%
import json
import numpy as np
d = json.load(open("../data/datasets.json","r"))

path = 'C:/Users/34646/Documents/Copenhagen/Clases/TFM/newdata/sperm0006'+str(d['data'][cellnum-1])+'/segmentation'

if 0 < cellnum < 8:
    d = d['dataset2']
elif 8 < cellnum < 16:
    d = d['dataset1']  
time_step = d['data_shape'][2]/d['scanning_freq'] #s 
voxel = np.array(d['voxel']) #um (x,y,z)


#%%
from functions import read

import cc3d #connected components in 3D
from skimage import io
import matplotlib.pyplot as plt

def choose_label(volume_time, time_frame):
    labels = cc3d.connected_components(volume_time[time_frame])
    if np.max(labels)==0:
        chosen_label=0
    else:
        plt.figure()
        ax = plt.axes(projection ='3d')
        ax.set_box_aspect([1,1,0.5])
    
        for l in range(1, np.max(labels)+1): #assuming label=0 is always the background
            k, j, i = np.where(labels==l) #k refers to the z direction, j to the y and i to the x
            j = y-j #The y axis is inverted with respect to how it is seen when printing the image (bc imshow flips the y axis)
            ax.plot3D(i*voxel[0], j*voxel[1], k*voxel[2], 'o', label='Label %i'%l)
    
        ax.legend()
        ax.set_xlabel('x/$\mu m$'); ax.set_ylabel('y/$\mu m$'); ax.set_zlabel('z/$\mu m$')
        ax.text2D(0.05, 0.95, 't = %i' % (time_frame+1), transform=ax.transAxes, size=14) #placement (0,0) would be the bottom left, (0,1) would be the top left
        ax.set_xlim(0,x*voxel[0]); ax.set_ylim(0,y*voxel[1])
        ax.set_title('Connected components in 3D')
        plt.show()
    
        chosen_label = int(input('Which label do you choose?\n'))
    return labels, chosen_label


# Read the segmented cell and head
#t, z, y, x are the same for the head and whole cell segmentations
_, t, z, y, x = read(path+'/'+path[55:-13]+'_segmentation_hea'+['d','d_gaussian'][Gaussian]+'.tif')
volume_time_head = io.imread(path+'/'+path[55:-13]+'_segmentation_hea'+['d','d_gaussian'][Gaussian]+'.tif')
#volume_time_cell = io.imread('segmentation/sperm00068_t1_146_segmentation_cell.tif')


#file_c = open('label_cell.txt','w')
file_h = open(path+'/label_hea'+['d','d_gaussian'][Gaussian]+'.txt','w')

for time_frame in range(0,t):
    # Do connected components of the cell and head segmentations, and choose desired label manually
    #labels_cell,  chosen_label_cell = choose_label(volume_time_cell, time_frame)
    labels_head,  chosen_label_head = choose_label(volume_time_head, time_frame)
    
    # Save the labels for all the timeframes for the main head and cell   
    #file_c.write(str(chosen_label_cell)+'\n')
    file_h.write(str(chosen_label_head)+'\n')
    
#file_c.close()  
file_h.close()