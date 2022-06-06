  # -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 17:03:17 2022

@author: Alicia Mosquera
"""
import os
import glob
import cv2 as cv
import matplotlib.pyplot as plt
from numpy.random import default_rng
rng = default_rng()

import sys
sys.path.insert(0, 'C:/Users/34646/Documents/Copenhagen/Clases/TFM/newdata/sperm00068_t1_146/segmentation')
from functions import read, plot_yz, plot_zx


#%%
# Create single image
def createImg(filename, path):
    path += '/segmentation/tuning/slice_t'
    volume_time, t, z, y, x = read(filename)
    
    # Create the random time and z coordinates, also y and x for the perpendicular slices
    t_chosen = rng.integers(low=0, high=t, size=3)
    z_chosen = rng.integers(low=0, high=z, size=1)[0]
    y_chosen = rng.integers(low=0, high=y, size=1)[0]
    x_chosen = rng.integers(low=0, high=x, size=1)[0]
    
    # XY plane
    img = volume_time[t_chosen[0]][z_chosen]
    cv.imwrite(path+str(t_chosen[0])+'_z'+str(z_chosen)+'.tiff', img)
    plt.figure()
    plt.imshow(img, cmap='gray')
    plt.title('slice_t'+str(t_chosen[0])+'_z'+str(z_chosen))
    plt.xlabel('x'); plt.ylabel('y')
    plt.savefig(path+str(t_chosen[0])+'_z'+str(z_chosen)+'.png')
    
    # YZ plane
    img = plot_yz(volume_time, t_chosen[1], x_chosen)
    cv.imwrite(path+str(t_chosen[1])+'_x'+str(x_chosen)+'.tiff', img)
    plt.figure()
    plt.imshow(img, cmap='gray')
    plt.title('slice_t'+str(t_chosen[1])+'_x'+str(x_chosen))
    plt.xlabel('y'); plt.ylabel('z')
    plt.savefig(path+str(t_chosen[1])+'_x'+str(x_chosen)+'.png')
    
    # ZX plane
    img = plot_zx(volume_time, t_chosen[2], y_chosen)
    cv.imwrite(path+str(t_chosen[2])+'_y'+str(y_chosen)+'.tiff', img)
    plt.figure()
    plt.imshow(img, cmap='gray')
    plt.title('slice_t'+str(t_chosen[2])+'_y'+str(y_chosen))
    plt.xlabel('z'); plt.ylabel('x')
    plt.savefig(path+str(t_chosen[2])+'_y'+str(y_chosen)+'.png')
    
    return 0


# Create set of images without duplicates
def replaceDuplicates(array1, array2):
    global t
    for i in range(len(array1)):
        for j in range(len(array1)):
            if i!=j and array1[i] == array1[j]:
                if array2[i] == array2[j]:
                    value = array1[i]
                    while value == array1[i]:
                        value = rng.integers(low=0, high=t)
                    array1[i] = value
    return array1, array2

def createSet(filename, n_img):
    #volume_time: 4d array of data
    #n_img: number of images to produce for exampler segmentation
    volume_time, t, z, y, x = read(filename)
    
    # Remove any previous files
    # files_tiff = glob.glob('./cases/slice_t*.tiff')
    # files_png = glob.glob('./cases/slice_t*.png')
    # files = files_tiff + files_png
    # for f in files:
    #     os.remove(f)
    
    # Create the time and z coordinates arrays for all images
    t_array = rng.integers(low=0, high=t, size=n_img)
    z_array = rng.integers(low=0, high=z, size=n_img)
    # Replace any image used more than once
    t_array, z_array = replaceDuplicates(t_array, z_array)
    
    for i in range(n_img, path):
        path += '/segmentation/tuning/slice_t'+str(t_array[i])+'_z'+str(z_array[i])
        img = volume_time[t_array[i]][z_array[i]]
       
        cv.imwrite(path+'.tiff', img)
        plt.figure(i)
        plt.imshow(img, cmap='gray')
        plt.title('slice_t'+str(t_array[i])+'_z'+str(z_array[i]))
        plt.savefig(path+'.png')
    
    return 0

#%%
path = 'C:/Users/34646/Documents/Copenhagen/Clases/TFM/newdata/sperm00068_t1_146'
filename = path+'/'+path[55:]+'.tif'

# Show all images created
files = glob.glob(path+'/segmentation/tuning/slice_t*.tiff')
for i,f in enumerate(files):
    img = cv.imread(f, -1)
    plt.figure(i)
    plt.imshow(img, cmap='gray')
    plt.title(f[93:-5])
    plt.axis('off')
plt.show()