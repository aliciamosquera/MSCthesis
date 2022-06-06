# -*- coding: utf-8 -*-
"""
Created on Mon May  9 15:04:43 2022

@author: Alicia Mosquera
"""
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from numpy.random import default_rng
rng = default_rng()

plt.style.use('C:/Users/34646/Documents/Copenhagen/Clases/TFM/thesis.mplstyle')

import sys
sys.path.insert(0, 'C:/Users/34646/Documents/Copenhagen/Clases/TFM/newdata/sperm00068_t1_146/segmentation')
from functions import read

volume_time, t, z_pix, y_pix, x_pix = read('C:/Users/34646/Documents/Copenhagen/Clases/TFM/data/data_9-2.tif') #Debería ser t,z,y,x
d = json.load(open("C:/Users/34646/Downloads/github/datasets.json","r"))
voxel = np.array(d['dataset0']['voxel'])

x, y, z = [x_pix, y_pix, z_pix]*voxel
zx_scale = voxel[2]/voxel[0]
print(f'Slices in Z: {z_pix} --> {z} um')
print('Time steps:', t)
print('Data type:', volume_time.dtype)
print()
print(f"Rows or Height of frame: {y} um")
print(f"Columns or Width of frame: {x} um")


#%%
def plot_zy(fullvolume_time, t, x):
    _, z, y, _ = fullvolume_time.shape
    img = np.zeros((y,z), dtype=np.uint16)
    for zi in range(z):
        for yi in range(y):
            img[yi][zi] = fullvolume_time[t][zi][yi][x]
    return img

def plot_xz(fullvolume_time, t, y):
    _, z, _, x = fullvolume_time.shape
    img = np.zeros((z,x), dtype=np.uint16)
    for zi in range(z):
        for xi in range(x):
            img[zi][xi] = fullvolume_time[t][zi][y][xi]
    return img


#%%
# How the image looks in 3D without resizing
fig, ax = plt.subplots(2, 2, gridspec_kw={'height_ratios': [y, z], 'width_ratios': [x, z]})

ax[0,0].imshow(volume_time[28][1], cmap='gray', extent=[0,x,y,0])
ax[0,0].set_title('$t=87$ s, $z\sim 0~\mu$m') #29*3 = 87 s (is the 29th frame, with 3 s time step, the first frame is already 3 s bc that's the time it took to take the picture)
ax[0,0].set_xlabel('$x/\mu m$'); ax[0,0].set_ylabel('$y/\mu m$')

img = plot_zy(volume_time, 28, 350)
ax[0,1].imshow(img, cmap='gray', extent=[0,z,y,0])
ax[0,1].set_title('$t=87$ s, $x\sim 101~\mu$m')
ax[0,1].set_xlabel('$z/\mu m$'); ax[0,1].set_ylabel('$y/\mu m$')

img = plot_xz(volume_time, 28, 200)
ax[1,0].imshow(img, cmap='gray', extent=[0,x,z,0])
ax[1,0].set_title('$t=87$ s, $y\sim 58~\mu$m')
ax[1,0].set_xlabel('$x/\mu m$'); ax[1,0].set_ylabel('$z/\mu m$')

ax[1,1].xaxis.set_major_locator(ticker.NullLocator())
ax[1,1].yaxis.set_major_locator(ticker.NullLocator())
ax[1,1].spines.clear()


#%%
# Resizing z axis to have the same resolution as x and y without scaling

#Because I don't know how to deal with the .3 pixel yet, I'm going to assume the voxel scale z/x is 10.
#Therefore, there should be a 9 black pixel separation between z-slices.
fullvolume_time = np.zeros((t, int(z_pix*zx_scale), y_pix, x_pix), dtype=np.uint16)
print(fullvolume_time.shape, fullvolume_time.dtype)

for ti in range(t):
    counter = 0
    for zi in range(0, (z_pix*10), 10):
        fullvolume_time[ti][zi] = np.copy(volume_time[ti][counter])
        counter+=1
    
fig, ax = plt.subplots(2, 2, gridspec_kw={'height_ratios': [y_pix, z_pix*10], 'width_ratios': [x, z]})
        
ax[0,0].imshow(fullvolume_time[28][10], cmap='gray', extent=[0,x,y,0])
ax[0,0].set_title('$t=87$ s, $z\sim 0~\mu$m')
ax[0,0].set_xlabel('$x/\mu m$'); ax[0,0].set_ylabel('$y/\mu m$')

img = plot_zy(fullvolume_time, 28, 350)
ax[0,1].imshow(img, cmap='gray', extent=[0,z,y,0], vmin=0, vmax=255)
ax[0,1].set_title('$t=87$ s, $x\sim 101~\mu$m')
ax[0,1].set_xlabel('$z/\mu m$'); ax[0,1].set_ylabel('$y/\mu m$')

img = plot_xz(fullvolume_time, 28, 200)
ax[1,0].imshow(img, cmap='gray', extent=[0,x,z,0])
ax[1,0].set_title('$t=87$ s, $y\sim 58~\mu$m')
ax[1,0].set_xlabel('$x/\mu m$'); ax[1,0].set_ylabel('$z/\mu m$')

ax[1,1].xaxis.set_major_locator(ticker.NullLocator())
ax[1,1].yaxis.set_major_locator(ticker.NullLocator())
ax[1,1].spines.clear()