# -*- coding: utf-8 -*-
"""
Created on Mon May  9 22:33:28 2022

@author: Alicia Mosquera
"""
import json
import numpy as np
import pandas as pd
from skimage import io
import matplotlib.pyplot as plt
from sklearn.cluster import SpectralClustering

plt.style.use('C:/Users/34646/Documents/Copenhagen/Clases/TFM/thesis.mplstyle')
d = json.load(open("C:/Users/34646/Downloads/github/datasets.json","r"))

prop = 'Head CoM'
def prop_real(prop):
    if prop=='Head CoM':
        return('Head centre')
    else:
        return prop

d = d['dataset0']
time_step = d['data_shape'][2]/d['scanning_freq'] #s
voxel = np.array(d['voxel']) #um (x,y,z)


#%%
#The annotations are saved as a csv file
#The dataframe is tweaked slightly to show only the parameters we are intereseted in
df = pd.read_csv('png/labels_data_9-2_png.csv')
#df = df.drop(['349', '175'], axis=1)
df['Frame'] = df.file_name.replace('data_9-2_|.png', '', regex=True).astype('int64')

df['x_pix'] = df['x'] - 7 #The png image had a 7 pixel border
df['y_pix'] = df['y'] - 7
df['z_pix'] = 30-df['Frame']%30
df['t_pix'] = ((df['Frame']-df['z_pix'])/30).astype('int64') + 1
df['t'] = df['t_pix']*time_step

df = df.reindex(columns=['file_name', 'label', 'Frame', 'x_pix', 'y_pix', 'z_pix', 't_pix', 't'])

#Sort the df by increasing value of z to better find the mirrored clusters (and for each z, increasing in t)
df = df.sort_values(by=['z_pix', 't'])


#%%
#Print the shape of one image to know the size of all and fix the axis limits easily
#As well as write the x, y and z coordinates in real units instead of pixel
im = io.imread('png/data_9-2_1.png')

ymax_pix = im.shape[0] - 7.   #maximum is 313 in tif
scaley = 313/ymax_pix * voxel[1] #um
xmax_pix = im.shape[1] - 7.   #maximum is 650 in tif
scalex = 650/xmax_pix * voxel[0] #um

#Convert pixel resolution into real measures
df['x'] = df['x_pix']*scalex
df['y'] = df['y_pix']*scaley
df['z'] = df['z_pix']*voxel[2]

tmax = np.max(df['t'])       #maximum is 106 in slices, 318 in s
zmax = np.max(df['z'])       #um
ymax = np.max(df['y'])       #um
xmax = np.max(df['x'])       #um

print('Max t:', tmax)
print('Max z:', zmax)
print('Max y:', ymax)
print('Max x:', xmax)
#print('scale x and y:', scaley, scalex)


#%%
# Spectral clustering
mask = (df['label']==prop)

#Try to find the clusters with scikit learn
X = []
length = len(df.x[mask])
for i in range(length):
    X = np.append(X, (df.t[mask].iloc[i], df.y[mask].iloc[i], df.z[mask].iloc[i]))
X = X.reshape((length, 3))

clustering = SpectralClustering(n_clusters=3,
                                assign_labels='discretize',
                                random_state=0).fit(X)

clust0 = (clustering.labels_==0)
clust1 = (clustering.labels_==1)
clust2 = (clustering.labels_==2)


#%%
fig1 = plt.figure()
ax1 = fig1.add_subplot(projection="3d")
ax1.set_title('3D Clustering (t,y,z)')
ax1.set_xlabel('x/μm'); ax1.set_ylabel('y/μm'); ax1.set_zlabel('z/μm')
#ax1.set_xlim(0,xmax)
ax1.set_ylim(0,ymax); ax1.set_zlim(0,zmax)
ax1.text2D(0.05, 0.90, '%s \n%i points' % (prop_real(prop), len(df.x[mask])), transform=ax1.transAxes, fontsize=12)

ax1.scatter(df.x[mask][clust0], df.y[mask][clust0], df.z[mask][clust0], s=5, color='b')
ax1.scatter(df.x[mask][clust1], df.y[mask][clust1], df.z[mask][clust1], s=5, color='orange')
ax1.scatter(df.x[mask][clust2], df.y[mask][clust2], df.z[mask][clust2], s=5, c='g')


fig2 = plt.figure()
ax2 = fig2.add_subplot(projection="3d")
ax2.set_xlabel('x/μm'); ax2.set_ylabel('y/μm'); ax2.set_zlabel('z/μm')
#ax2.set_xlim(0,xmax)
ax2.set_ylim(0,ymax); ax2.set_zlim(0,zmax)
ax2.text2D(0.05, 0.90, '%s \n%i points' % (prop_real(prop), len(df.x[mask])), transform=ax2.transAxes, fontsize=12)

p = ax2.scatter(df.x[mask], df.y[mask], df.z[mask], c=df.t[mask], s=5, cmap='viridis')

cbax = fig2.add_axes([0.1, 0.2, 0.02, 0.6]) #Adding the colorbar to the left (x start, y start, width, length)
cbar = fig2.colorbar(p, cax=cbax, label='$t/s$')
cbar.ax.yaxis.set_ticks_position('left')
cbar.ax.yaxis.set_label_position('right')
#Adding a colorbar to the bottom
#cbar = fig2.colorbar(p, ax=[ax2], location='bottom', shrink=0.65, label='t')
#cbar.ax.invert_xaxis()

plt.show()


#%%
# Measure the experimental size of the head's cell
def measure_dist(df,label1,label2):
    d = []
    for i in range(np.min(df.Frame),np.max(df.Frame+1)):
        mask1 = (df['Frame']==i) & (df['label']==label1)
        mask2 = (df['Frame']==i) & (df['label']==label2)
        if sum(mask1*np.ones(len(mask1)))==1 and sum(mask2*np.ones(len(mask2)))==1:
            x1 = df.x_pix[mask1]
            x2 = df.x_pix[mask2]
            y1 = df.y_pix[mask1]
            y2 = df.y_pix[mask2]
            d.append(np.linalg.norm(np.array([x1,y1]) - np.array([x2,y2])))
    return(len(d), np.mean(d), np.std(d)/np.sqrt(len(d)))


#%%
# Measure the vesicle's linear and curvilinear speeds
def calc_speed(head_coord, time_frame):
    deltat = np.reshape(np.diff(time_frame).tolist()*3, (3,-1,)).T
    diff = np.diff(head_coord, axis=0)
    velocity = diff/deltat
    speed = np.linalg.norm(velocity, axis=1)
    return velocity, speed

def find_centre_of_frame(df, prop):
    mask1 = (df['label']==prop)
    time_frame = []
    head_coord = []
    for i in range(np.min(df.t_pix), np.max(df.t_pix)+1):
        mask = mask1 & (df['t_pix']==i)
        if sum(mask*np.ones(len(mask)))>0:
            time_frame.append(df.t[mask].iloc[0])
            head_x, head_y, head_z = [df.x[mask], df.y[mask], df.z[mask]]
            if len(head_coord) == 0:
                head_coord = np.array([np.mean(head_x), np.mean(head_y), np.mean(head_z)])
            else:
                head_coord = np.row_stack([head_coord, np.array([np.mean(head_x), np.mean(head_y), np.mean(head_z)])])
    return(np.array(time_frame), head_coord)

time_frame, head_coord = find_centre_of_frame(df, 'Midpiece')
VSL = np.linalg.norm((head_coord[-1]-head_coord[0])/time_frame[-1])
velocity, speed = calc_speed(head_coord, time_frame)
VCL, sVCL = [np.mean(speed), np.std(speed)]
LIN = VSL/VCL
print('VSL=%.4f $\mu m s^{-1}$, VCL=(%.4f \pm %.2f) $\mu m s^{-1}$' %(VSL,VCL,sVCL))
print('LIN=%.2f percent' %(LIN*100))