# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 16:39:27 2022

@author: Alicia Mosquera
"""
#%%
import cc3d
import numpy as np
import pandas as pd
from skimage import io
import matplotlib.pyplot as plt

plt.style.use('C:/Users/34646/Documents/Copenhagen/Clases/TFM/thesis.mplstyle')

def read(filename):
    volume_time = io.imread(filename)
    t,z,y,x = volume_time.shape
    return volume_time,t,z,y,x

volume_time,t,z,y,x = read('sperm00068_t1_146/sperm00068_t1_146.tif')
segmented_head = io.imread('sperm00068_t1_146/segmentation/sperm00068_t1_146_segmentation_head.tif')
labels_h = pd.read_csv('sperm00068_t1_146/segmentation/label_head_gaussian_done.txt', names=['head'])

#%%
# Estimate the noise as the std of homogeneous areas (both background and foreground)
# DONE WITH IMAGEJ
df1 = pd.read_csv('noise_std_00068.csv') #dat0:9-2 (10 Hz), dat1:00068 (40 Hz), dat2:00064 (70 Hz)

mean = df1['Mean'].to_numpy()
std = df1['StdDev'].to_numpy()

mean_b = np.mean(mean[:10])
mean_f = np.mean(mean[10:])
std_b = np.sqrt(np.mean(std[:10]**2))
std_f= np.sqrt(np.mean(std[10:]**2))

print('The background mean is %.2f and the foreground %.2f' % (mean_b, mean_f))
print('SBR is %.2f' % (mean_f/mean_b))
print('The background noise is %.2f and the foreground %.2f' % (std_b, std_f))
print('SNR in foreground is %i' % (mean_f/std_f))

#%%
# Histogram of the whole dataframe
df = pd.read_csv('hist_00068.csv')
df_param = pd.read_csv('hist_00068_param.csv')

hist = df['count'].to_numpy()
bin_start = df['bin start'].to_numpy()
bin_width = df_param['Bin Width'][0]/1000
mean = df_param['Mean'][0]/1000
std = df_param['StdDev'][0]/1000

#%%
fig, ax1 = plt.subplots()
ax1.set_xlabel('Intensity')
ax1.set_ylabel('Counts')
ax1.set_yscale('log')
ax1.step(bin_start, hist, where='post', color='tab:green', alpha=1, label='Dataset 1')
ax1.legend()

#ax2 = ax1.twinx()
#ax2.set_ylabel('Counts')
#ax2.bar(bin_start, hist, width=bin_width, color='tab:blue', label='Counts')
#ax1.tick_params(axis ='y', labelcolor=(0.12156862745098039, 0.4666666666666667, 0.7058823529411765,1)) 
#ax2.tick_params(axis ='y', labelcolor=(0.12156862745098039, 0.4666666666666667, 0.7058823529411765, 0.5)) 

textstr = '$\mu=%.i$, $\sigma=%.1f$' %(mean,std)
props = dict(boxstyle='round', facecolor='white', alpha=0.25)
#ax1.text(0.6, 0.96, textstr, transform=ax1.transAxes, fontsize=14, verticalalignment='top', bbox=props)
#fig.legend(loc="upper right", bbox_to_anchor=(1,1), bbox_transform=ax1.transAxes)

#%%
# Check for the histogram of the cell in each z-slice

# fig = plt.figure()
# ax1 = fig.add_subplot(211)
# ax2 = fig.add_subplot(212)
# for ti in range(56,57):
#     img = volume_time[ti]
#     img_head = segmented_head[ti]
#     labels_head = cc3d.connected_components(img_head)
#     chosen_label_head = labels_h['head'][ti]
    
#     img_head[np.where(labels_head!=chosen_label_head)] = 0    
#     img_head = img_head.astype('bool')
#     masked = img_head * img
    
    
#     # for zi in range(z):
#     #     ax1.imshow(img[zi], cmap='gray')
#     #     ax2.hist(masked[zi].flatten(), range=(1,np.max(img)), label=f'z={zi+1}', density=True)
#     hist_data = np.reshape(masked, (z,x*y)).T
#     ax2.hist(hist_data, range=(1,np.max(img)), label=np.arange(z+1), density=True)#, stacked=True)

#     ax1.set_title(f't={ti+1}')
#     ax1.set_xlabel('x/pixels'); ax1.set_ylabel('y/pixels')
#     ax2.set_xlabel('Intensity of pixel'); ax2.set_ylabel('Probability density')
#     ax2.legend()

#%%
# I have the histogram of the pixel intensities in the cell for each z-slice
# I want to create a slider for all the time frames

import plotly.graph_objects as go
import plotly.io as pio
pio.renderers.default='browser' #'svg' for showing the figure in Spyder

zmin = 3; zmax = 6
show_cell = True # choose between histogram of the cell or the background

# Create figure
fig = go.Figure()

# Add traces, one ti for each slider step and the histogram of all the zi
for ti in range(t):
    img = volume_time[ti]
    img_head = segmented_head[ti]
    labels_head = cc3d.connected_components(img_head)
    chosen_label_head = labels_h['head'][ti]
       
    img_head[np.where(labels_head!=chosen_label_head)] = 0
    img_head = img_head.astype('bool')
    if show_cell==False:
        img_head = ~img_head
    masked = img_head * img
    
    for zi in range(zmin,zmax):
        fig.add_trace(go.Histogram(
            x=masked[zi].flatten(),
            histnorm='probability',
            name='z'+str(zi+1),
            opacity=0.8,
            xbins=dict(
                start=1.0,
                end=np.max(img)),
            visible=ti<1 # Make 1st time step trace visible
        ))

# Create and add slider
steps = []
for ti in range(t):
    step = dict(
        method="update",
        args=[{"visible": [False] * t*(zmax-zmin)},
              {"title": "Pixel intensities in " + ['background','cell'][show_cell] + " at time step " + str(ti)}],  # layout attribute
        label='t={}'.format(ti)
    )
    for zi in range(zmin,zmax):
        step["args"][0]["visible"][(zmax-zmin)*ti+(zi-zmin)] = True  # Toggle i'th trace to "visible"
    steps.append(step)

slider = [dict(
        active=0,
        pad={"t": 25},
        steps=steps
)]

# Update figure        
fig.update_layout(
    xaxis_title_text='Pixel intensity', # xaxis label
    yaxis_title_text='Probability of counts', # yaxis label
    bargap=0.2, # gap between bars of adjacent location coordinates
    bargroupgap=0.1, # gap between bars of the same location coordinates
    sliders=slider,
    barmode='stack'
)

fig.show()
#%%
# I have the histogram of the pixel intensities in the cell for each z-slice
# I want to create a slider for all the time frames

# from matplotlib.widgets import Slider, Button

# Create the Slider
# time_ax = plt.axes([0.25, 0.1, 0.65, 0.03])
# slider = Slider(
#     ax=time_ax,
#     label='Time frame',
#     valmin=1,
#     valmax=t,
#     valinit=1,
#     valstep=1,
#     initcolor='none'
# )

# fig, ax = plt.subplots()
# # Fixing bin edges
# HIST_BINS = np.linspace(50, 200, 20)
# def update(ti):
#     # Update the image's time frame
#     img = volume_time[ti]
#     img_head = segmented_head[ti]
#     labels_head = cc3d.connected_components(img_head)
#     chosen_label_head = labels_h['head'][ti]
    
#     img_head[np.where(labels_head!=chosen_label_head)] = 0    
#     img_head = img_head.astype('bool')
#     masked = img_head * img
#     hist_data = np.reshape(masked, (z,x*y)).T

#     ax.hist(hist_data, HIST_BINS, label=np.arange(z+1), density=True)

#     # Redraw the figure to ensure it updates
#     fig.canvas.draw_idle()

# slider.on_changed(update)
# plt.show()
#%%
# # I have the histogram of the pixel intensities in the cell for each z-slice
# # I want to create an animation for all the time frames

# # There's a problem with trying to animate the histogram of multiple data, I tried to copy https://matplotlib.org/stable/gallery/animation/animated_histogram.html

# import matplotlib.animation as animation

# # Fixing bin edges
# HIST_BINS = np.linspace(50, 200, 20)

# # histogram our data with numpy
# img = volume_time[0]
# img_head = segmented_head[0]
# labels_head = cc3d.connected_components(img_head)
# chosen_label_head = labels_h['head'][0]

# img_head[np.where(labels_head!=chosen_label_head)] = 0    
# img_head = img_head.astype('bool')
# masked = img_head * img
# hist_data = np.reshape(masked, (z,x*y)).T
# #n, _ = np.histogram(hist_data, HIST_BINS)

# def prepare_animation(bar_container):

#     def animate(ti):
#         img = volume_time[ti]
#         img_head = segmented_head[ti]
#         labels_head = cc3d.connected_components(img_head)
#         chosen_label_head = labels_h['head'][ti]
        
#         img_head[np.where(labels_head!=chosen_label_head)] = 0    
#         img_head = img_head.astype('bool')
#         masked = img_head * img
        
#         #hist_data = np.reshape(masked, (z,x*y)).T
#         returned = []
#         for zi in range(z):
#             hist_data = masked[zi].flatten()
#             n, _ = np.histogram(hist_data, HIST_BINS)
#             for count, rect in zip(n, bar_container[zi].patches):
#                 rect.set_height(count)
                
#             returned = np.append(returned, bar_container[zi].patches)
#             return returned
#     return animate

# fig, ax = plt.subplots()
# _, _, bar_container = ax.hist(hist_data, HIST_BINS, label=np.arange(z+1), density=True)
# ax.set_ylim(top=0.5)  # set safe limit to ensure that all data is visible.

# ani = animation.FuncAnimation(fig, prepare_animation(bar_container), t, repeat=False, blit=True)
# plt.show()