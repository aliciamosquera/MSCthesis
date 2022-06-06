# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 11:52:49 2022

@author: Alicia Mosquera
"""
cellnum = 10
# Type of segmentation to analyze (Gaussian/ImageJ)
Gaussian = True
neck_exist = False


#%%
import json
import numpy as np
d = json.load(open("C:/Users/34646/Downloads/github/datasets.json","r"))

path = 'C:/Users/34646/Documents/Copenhagen/Clases/TFM/newdata/sperm0006'+str(d['data'][cellnum-1])+'/segmentation'

if 0 < cellnum < 8:
    d = d['dataset2']
elif 8 < cellnum < 16:
    d = d['dataset1']  
time_step = d['data_shape'][2]/d['scanning_freq'] #s 
voxel = np.array(d['voxel']) #um (x,y,z)


#%%
import sys
sys.path.insert(0, 'C:/Users/34646/Documents/Copenhagen/Clases/TFM/newdata/sperm00068_t1_146/segmentation')
from functions import reshapeToEvenSides, read

import cc3d #connected components in 3D
import pandas as pd
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from skimage import io, morphology
#import matplotlib.animation as animation

plt.style.use('C:/Users/34646/Documents/Copenhagen/Clases/TFM/thesis.mplstyle')


def PCA(img_h, origin, file_h, count, eigvec_old):
    k,j,i = np.where(img_h!=0)
    i = i*voxel[0] - origin[0]
    j = j*voxel[1] - origin[1]
    k = k*voxel[2] - origin[2]

    head_full = np.array([i,j,k])
    covMatrix = np.cov(head_full, bias=False) #bias is False to calculate the sample covariance and not the population one
    eigval_ns, eigvec_ns = np.linalg.eig(covMatrix) #ns: non-sorted
    
    # Sort axes from biggest to smallest
    idx_sorted = np.argsort(eigval_ns)[::-1] #sort the eigenvalues in descending order
    eigval_s = eigval_ns[idx_sorted]
    eigvec_s = eigvec_ns[:,idx_sorted]
    
    # Choose the orientation of the axes (PCA only gives the direction)
    if len(eigvec_old)>1: #ensures eigvec_old exists
        for m in range(3):
            orient = np.sign(np.dot(eigvec_old[:,m],eigvec_s[:,m]))
            if orient != 0:
                eigvec_s[:,m] *= orient
    
    # Make sure the coordinate axis is always x,y,z and not x,y,-z
    eigvec_s[2] = np.dot(np.cross(eigvec_s[0], eigvec_s[1]),eigvec_s[2])*eigvec_s[2]
    
    # The vectors shown are the eigenvectors of the covariance matrix scaled by
    # the square root of the corresponding eigenvalue
    if np.all(eigval_s):
        vec_s = eigvec_s/np.sqrt(eigval_s)
    else: #in case any eigval_s = 0
        vec_s = np.copy(eigvec_s)
        if eigval_s[0]>0:
            vec_s[:,0] = vec_s[:,0]/np.sqrt(eigval_s[0])
            if eigval_s[1]>0:
                vec_s[:,1] = vec_s[:,1]/np.sqrt(eigval_s[1])
    
    return i,j,k, eigvec_s, vec_s, eigval_s

def find_neck(time_frame, img_c, img_h, it):
    # I dilate two times and compute the difference
    dilation1 = ndi.binary_dilation(img_h, morphology.ball(radius=2), iterations=it)
    dilation2 = ndi.binary_dilation(img_h, morphology.ball(radius=2), iterations=(it+3))
    band1 = dilation2^dilation1
    
    #? do connected components here to assure that we'll do the mean of only the neck?
    # use a mask of the x projection function??
    
    # Where the dilation ring and the cell segmentation overlap, there is the neck
    neck = band1&img_c
    neck_coord = np.where(neck!=0)

    if len(neck_coord[0])>0:
        neck_coord = np.mean(neck_coord, axis=1)
        neck_coord = neck_coord[::-1] #to specify in x,y,z and not z,y,x
    else: #in case there's no overlapping
        neck_coord = np.array([0,0,0]) 
    neck_coord = neck_coord*voxel

    return neck_coord


#%%
# Read the data of the segmented head and the whole cell
# t, z, y, x are the same for the head and whole cell segmentations
_, t, z_pix, y_pix, x_pix = read(path+'/'+path[55:-13]+'_segmentation_hea'+['d','d_gaussian'][Gaussian]+'.tif')
x,y,z = np.array([x_pix,y_pix,z_pix])*voxel

print(f'Slices in Z: {z_pix} --> {z} um')
print(f'Time frames: {t} --> {t*time_step} s')
print()
print(f"Rows or Height of frame: {y} um")
print(f"Columns or Width of frame: {x} um")

if neck_exist: volume_time_cell = io.imread(path+'/'+path[55:-13]+'_segmentation_cell.tif')
volume_time_head = io.imread(path+'/'+path[55:-13]+'_segmentation_hea'+['d','d_gaussian'][Gaussian]+'.tif')

if neck_exist: volume_time_cell = reshapeToEvenSides(volume_time_cell)
volume_time_head = reshapeToEvenSides(volume_time_head)
t, z_pix, y_pix, x_pix = volume_time_head.shape

# Read the txt with the labels chosen for the connected components
if neck_exist: file_c = pd.read_csv(path+'/label_cell_done.txt', names=['cell'])
file_h = pd.read_csv(path+'/label_hea'+['d','d_gaussian'][Gaussian]+'_done.txt', names=['head'])

# Start the lists to later convert into a dataframe
idx         = np.zeros(t)
t_all       = np.zeros(t)
origin_all  = np.zeros((t,3))
eigvec1_all = np.zeros((t,3))
eigvec2_all = np.zeros((t,3))
eigvec3_all = np.zeros((t,3))
necki_all   = np.zeros((t,3))
neckh_all   = np.zeros((t,3))


#%%
# Find center of mass of the midpiece and perform PCA of the head in all time frames
# Segmented points in blue, the neck points in orange, and their coordinate system 
fig   = plt.figure()
ax    = plt.axes(projection ='3d')
ax.set_box_aspect([1,0.5,0.3])
fig_h = plt.figure()
ax_h  = plt.axes(projection ='3d')

history_d = []
history_x = []; history_y = []; history_z = []
history_angles = []
history_vol = []; history_volz = []
eigvec_old = []; vec_old = []

for count,ti in enumerate(np.arange(0,t)): #enumerate in case the range is not from 0 on
    if neck_exist: img_cell = volume_time_cell[ti]
    img_head = volume_time_head[ti]
    
    # Do connected components of the cell and head segmentations, and use the labels obtained manually
    if neck_exist: labels_cell = cc3d.connected_components(img_cell)
    labels_head = cc3d.connected_components(img_head)
    stats = cc3d.statistics(labels_head)
    
    chosen_label_cell=1
    if neck_exist: chosen_label_cell = file_c['cell'][ti]
    chosen_label_head = file_h['head'][ti]
    
    # We continue only if the wanted cell and head exist in the time frame
    if chosen_label_cell!=0 and chosen_label_head!=0:
        # Remove the parts of the image which do not belong to the chosen label
        img_head[np.where(labels_head!=chosen_label_head)] = 0    
        img_head = img_head.astype('bool')
        if neck_exist: img_cell[np.where(labels_cell!=chosen_label_cell)] = 0
        if neck_exist: img_cell = img_cell.astype('bool')
        
        if len(img_head[img_head>0]) > 5: #necessary to have at least 5 pixels segmented with the chosen label, at the specified time frame 
            # Perform the PCA of the head
            origin = stats.get('centroids')[chosen_label_head][::-1]*voxel #same as origin = np.mean(img_head, axis=1) in PCA function
            i,j,k, eigvec, vec, eigval = PCA(img_head, origin, file_h, count, eigvec_old)    
        
            # Save volume of cell for later histogram
            pixels_on = len(img_head[img_head>0])
            pixels_on_z = 0
            for zi in range(z_pix):
                if np.any((img_head>0)[zi]):
                    pixels_on_z += 1
            history_volz = np.append(history_volz, pixels_on_z)  
            #history_vol = np.append(history_vol, pixels_on*np.prod(voxel)) #um^3
            if np.all(eigval):
                history_vol = np.append(history_vol, 4/3*np.pi*np.prod(eigval))
            
            # Find the center point of the neck
            neck = np.zeros(3)
            if neck_exist: neck = find_neck(ti, img_cell, img_head, 2)
            #tail = find_neck(ti, img_cell, img_head, 4)
            
            # In image coordinates
            if (ti%1 == 0):
                #ax.plot3D(i, j, k, '.')
                #ax.plot3D(*origin, 'c.')
                ax.plot(origin[1], origin[2], '.', color='tab:blue', zdir='x', zs=0, markersize=2)
                ax.plot(origin[0], origin[2], '.', color='tab:blue', zdir='y', zs=45, markersize=2)
                ax.plot(origin[0], origin[1], '.', color='tab:blue', zdir='z', zs=4, markersize=2)
                for l in range(vec.shape[1]):
                    ax.quiver(*origin, *vec[:,l], color=['r','g','k'][l])
                ax.set_xlim(0,130); ax.set_ylim(5,45); ax.set_zlim(4,12)
                ax.set_xlabel('x/$\mu m$'); ax.set_ylabel('y/$\mu m$'); ax.set_zlabel('z/$\mu m$')
                        
            if np.all(neck==0): # if there isn't a neck present
                new_neck = np.zeros(3)
            else: # if there's a neck
                # Calculate the neck point in head coordinates and its distance to the head center
                p = neck-origin
                new_neck = eigvec.T@p #the matrix with the eigenvectors is the rotation matrix
    
                d = np.linalg.norm(p)
                history_d.append(d)
            
                # In head coordinates
                if (ti%1 == 0):
                    #ax_h.plot(*new_neck, '.',color='chocolate')                
                    for l in range(eigvec.shape[1]):
                        new_vec = 3*eigvec.T@(eigvec[:,l]) #the 3 is to make the axes bigger in the representation
                        #ax_h.quiver(0,0,0, *new_vec, color=['r','g','k'][l])
                    history_x.append(new_neck[0]); history_y.append(new_neck[1]); history_z.append(new_neck[2])
                ax_h.plot(history_y, history_z, 'r', zdir='x', zs=-30, markersize=3)
                ax_h.plot(history_x, history_z, 'g', zdir='y', zs=5, markersize=3)
                ax_h.plot(history_x, history_y, 'k', zdir='z', zs=-5, markersize=3)
                #ax_h.plot(history_x, history_y, history_z, color='orange')
                ax_h.set_xlabel('X/$\mu m$'); ax_h.set_ylabel('Y/$\mu m$'); ax_h.set_zlabel('Z/$\mu m$')
                ax_h.set_xlim(-30,15); ax_h.set_ylim(-5,5); ax_h.set_zlim(-5,12)
    
            # Save them to know the orientation of the next eigenvector
            eigvec_old = eigvec
            vec_old = vec


#%%        
            # Saving data for the dataframe
            idx[ti]           = ti
            t_all[ti]         = ti*time_step #so the time is in seconds
            origin_all[ti,:]  = origin[:]
            eigvec1_all[ti,:] = eigvec[:,0]
            eigvec2_all[ti,:] = eigvec[:,1]
            eigvec3_all[ti,:] = eigvec[:,2]
            necki_all[ti,:]   = neck[:]
            neckh_all[ti,:]   = new_neck[:]

# Save the dataframe while removing the steps where there's no segmentation
mask = idx>0
df = pd.DataFrame()
df['idx']             = idx[mask]
df['t']               = t_all[mask] #s
df['head_coord']      = origin_all[mask].tolist() #um
df['head_axis_eigv1'] = eigvec1_all[mask].tolist() #um
df['head_axis_eigv2'] = eigvec2_all[mask].tolist() #um
df['head_axis_eigv3'] = eigvec3_all[mask].tolist() #um
df['neck_coord_i']    = necki_all[mask].tolist() #um
df['neck_coord_h']    = neckh_all[mask].tolist() #um

for i in range(len(origin_all)): origin_all[i]=origin_all[i] + i%2 * np.array([0, 0, voxel[2]])
df['head_coord_plus1'] = origin_all[mask].tolist() #um

df.to_csv(path[:-13]+'/realdata_'+['imagej','gaussian'][Gaussian]+'_param.txt')


#%%
#a = plt.hist(volumee_gaus, range=[0,9], bins=9)
#b = plt.hist(volumee_imgj, range=[0,9], bins=9)

# # Histogram of distances btw neck and head
fig_histd = plt.figure()
ax_histd = fig_histd.add_subplot(111)
counts, bins, _ = ax_histd.hist(history_d, color='tab:blue', alpha=1, label='Experimental data', bins=17)
ax_histd.fill_betweenx(np.array([40,0]), 9.5, 10.9, facecolor='tab:orange', alpha=0.5, label='Expected range')
ax_histd.set_xlabel('Distance/$\mu m$'); ax_histd.set_ylabel('Counts')
ax_histd.set_ylim(0,np.max(counts)+2)
ax_histd.legend()
ax_histd.set_title('Between the head\'s and the midpiece\'s centroids')
#ax_histd.set_title(['ImageJ','Gaussian'][Gaussian]+' segmentation')

# Histogram of volume of the cell in each timeframe
fig_histv = plt.figure()
ax_histv = fig_histv.add_subplot(111)
ax_histv.hist(history_vol)
ax_histv.set_xlabel('Head\'s volume/$\mu m^3$'); ax_histv.set_ylabel('Counts')
ax_histv.set_title(['ImageJ','Gaussian'][Gaussian]+' segmentation')
#width = a[1][1:]-a[1][:-1]
#labels = a[1][:-1]
#ax_histv.bar(labels, a[0], width, label='Dataset 1')
#ax_histv.bar(labels, b[0], width, label='Dataset 2', bottom=a[0]) #bottom for stacked histogram
#ax_histv.legend()

# Histogram of z-slices of the cell in each timeframe
fig_histvz = plt.figure()
ax_histvz = fig_histvz.add_subplot(111)
counts, bins = np.histogram(history_volz, bins=int(np.ptp(history_volz)+1))
bins = np.arange(np.min(history_volz),np.max(history_volz)+2)
#ax_histvz.hist(bins[:-1], bins, weights=counts, align='left')
ax_histvz.set_xlabel('# of z slices containing the cell'); ax_histvz.set_ylabel('Counts')
ax_histvz.set_title(['ImageJ','Gaussian'][Gaussian]+' segmentation')

width = (a[1][1:]-a[1][:-1])/1.3
labels = a[1][:-1]
ax_histvz.bar(labels, a[0], width, alpha = 0.7, label='Gaussian')
ax_histvz.bar(labels, b[0], width, alpha = 0.7, label='ImageJ')#, bottom=a[0]) #bottom for stacked histogram
ax_histvz.legend()