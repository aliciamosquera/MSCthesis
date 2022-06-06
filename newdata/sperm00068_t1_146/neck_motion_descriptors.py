# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 00:30:20 2022

@author: Alicia Mosquera
"""
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.style.use('C:/Users/34646/Documents/Copenhagen/Clases/TFM/thesis.mplstyle')


#%%
cellnum = 10
# Type of segmentation to analyze (Gaussian/ImageJ)
Gaussian = True

d = json.load(open("C:/Users/34646/Downloads/github/datasets.json","r"))
path = 'C:/Users/34646/Documents/Copenhagen/Clases/TFM/newdata/sperm0006'+str(d['data'][cellnum-1])

def convert_to_array(x):
    return np.array(x.str.strip('[]').str.split(',').to_list(), dtype=float)

df = pd.read_csv(path+'/realdata_'+['imagej','gaussian'][Gaussian]+'_param.txt') #simulation_param.txt
#df = df.iloc[::2] # to resample to half the scanning frequency
idx          = df['idx'].to_numpy()
print('All consecutive time frames?', ~np.any(np.diff(idx, axis=0)-1))
time_frame   = df['t'].to_numpy()
head_coord   = convert_to_array(df['head_coord'])
eigvec1      = convert_to_array(df['head_axis_eigv1'])
eigvec2      = convert_to_array(df['head_axis_eigv2'])
eigvec3      = convert_to_array(df['head_axis_eigv3'])
neck_coord_i = convert_to_array(df['neck_coord_i']) #Neck coordinates wrt the image axes
neck_coord_h = convert_to_array(df['neck_coord_h']) #Neck coordinates wrt the head axes
n_points     = len(time_frame)

head_coord_plus1 = convert_to_array(df['head_coord_plus1']) 


#----------- HEAD'S SPEED (VSL, VCL) -----------------
def calc_speed(head_coord, time_frame):
    deltat = np.reshape(np.diff(time_frame).tolist()*3, (3,-1,)).T
    diff = np.diff(head_coord, axis=0)
    velocity = diff/deltat
    speed = np.linalg.norm(velocity, axis=1)
    return velocity, speed

VSL = np.linalg.norm((head_coord[-1]-head_coord[0])/time_frame[-1]-time_frame[0])
velocity, speed = calc_speed(head_coord, time_frame)
VCL, sVCL = [np.mean(speed), np.std(speed)]
LIN = VSL/VCL
print('VSL=%.4f $\mu m$, VCL=%.4f \pm %.2f $\mu m$' %(VSL,VCL,sVCL))
print('LIN=%.2f percent' %(LIN*100))

#VSLt.append(VSL)
#VCLt.append(VCL)
#LINt.append(LIN*100)

fig_speed, ax_speed = plt.subplots()
ax_speed.set_xlabel('Time/s')
ax_speed.set_ylabel('Speed/$\mu m\ s^{-1}$')
ax_speed.plot(time_frame[1:], speed, c='tab:blue', alpha=0.7, label='data')
ax_speed.plot(time_frame[1:], savgol_filter(speed, window_length=21, polyorder=2, deriv=0), c='k', linestyle='-', label='smoothed data')
#ax_speed.hlines(np.mean(speed), time_frame[1], time_frame[-1], linestyles='--', color='tab:red', label='VCL (mean)')

_,speed_plus1 = calc_speed(head_coord_plus1, time_frame)
#ax_speed.plot(time_frame[1:], speed_plus1, 'tab:orange', alpha=0.7, label='disp. centroid')
#ax_speed.plot(time_frame[1:], savgol_filter(speed_plus1, window_length=21, polyorder=2, deriv=0), c='k', linestyle='--', label='smoothed disp. centroid')
#ax_speed.hlines(np.mean(speed_plus1), time_frame[1], time_frame[-1], linestyles='-', color='g', label='mean')
#ax_speed.set_title('Centroid displaced 1 pixel in $z$')
ax_speed.legend()

#divider_speed = make_axes_locatable(ax_speed)
#cbax_speed = divider_speed.append_axes("right", size="5%", pad=0.05)
#cbar_speed = fig_speed.colorbar(p, cax=cbax_speed, ticks=np.arange(2, 20+1, 2), label='Cell')

#----------- MSD for time lags -----------------
def calc_msd(r):
    N = r.shape[0]
    MSD = []
    std = []
    for n in range(1,N): #n = 1,..., N-1
        suma = []
        for i in range(N-n): #i = 0,...,N-n-1
            suma = np.append(suma, np.linalg.norm(r[i+n] - r[i])**2)
        mean = np.sum(suma)/(N-n)
        variance = np.mean((suma-mean)**2)
        std = np.append(std, np.sqrt(variance))
        MSD = np.append(MSD, mean)
    return MSD, std

MSD_headcoord, std = calc_msd(head_coord)
fig_msd, ax_msd = plt.subplots()
ax_msd.set_xlabel('Time lag/s')
ax_msd.set_ylabel('MSD/$\mu m^2$')
ax_msd.set_yscale('log')
ax_msd.set_xscale('log')
ax_msd.set_title('Mean squared displacement for different cells')# and $\pm\sigma$ interval')

#ax_msd.fill_between(time_frame[1:], MSD_headcoord-std, MSD_headcoord+std, color=c[cellnum-1], alpha=.2, linewidth=0)
p = ax_msd.scatter(x=time_frame[1:], y=MSD_headcoord, c=np.ones(len(MSD_headcoord))*cellnum, vmin=1-0.5, vmax=20+0.5, s=15, cmap='tab20c_r')

a = np.arange(0,7**3,0.5)
ax_msd.plot(a, a, linestyle='dotted', color='k', label='slope=1')
ax_msd.plot(a, a**2, linestyle='dashed', color='k', label='slope=2')
ax_msd.legend()

divider = make_axes_locatable(ax_msd)
cbax = divider.append_axes("right", size="5%", pad=0.05)
cbar = fig_msd.colorbar(p, cax=cbax, ticks=np.arange(2, 20+1, 2), label='Cell')

#%%
print()
def R_x(alpha):
    return np.array([[1,0,0], [0,np.cos(alpha),-np.sin(alpha)], [0,np.sin(alpha),np.cos(alpha)]])
def R_y(beta):
    return np.array([[np.cos(beta),0,np.sin(beta)], [0,1,0], [-np.sin(beta),0,np.cos(beta)]])
def R_z(gamma):
    return np.array([[np.cos(gamma),-np.sin(gamma),0], [np.sin(gamma),np.cos(gamma),0], [0,0,1]])

#----------- Tait-Bryan ANGLES ------------------------
def tait_angles(eigvec1, eigvec2, eigvec3): #rotation zyx
    alpha, beta, gamma = [], [], []
    for i in range(1,len(eigvec1)):
        Ri = np.column_stack([eigvec1[i-1],eigvec2[i-1],eigvec3[i-1]]) #the eigvec matrix is the rotation matrix
        Rf = np.column_stack([eigvec1[i],eigvec2[i],eigvec3[i]])
        #To test the code
        #Ri = R_z(np.pi/4)@R_y(0)@R_x(0)
        #Rf = R_z(np.pi/2)@R_y(0)@R_x(0)
        
        R = np.matmul(Rf,Ri.T)
        if round(R[2,0]**2,4) != 1: #otherwise the rounding off error throws the if off
            beta1 = np.arctan2(-R[2,0], np.sqrt(1-R[2,0]**2))
            beta2 = np.pi-beta1
            
            alpha1 = np.arctan2(R[2,1]/np.cos(beta1), R[2,2]/np.cos(beta1))
            alpha2 = np.arctan2(R[2,1]/np.cos(beta2), R[2,2]/np.cos(beta2))
            
            gamma1 = np.arctan2(R[1,0]/np.cos(beta1), R[0,0]/np.cos(beta1))
            gamma2 = np.arctan2(R[1,0]/np.cos(beta2), R[0,0]/np.cos(beta2))
        else:
            gamma1, gamma2 = 0,0 #Gimbal lock, could be anything so we set it to zero and calculate alpha
            if R[2,0] == -1:
                beta1, beta2 = [np.pi/2]*2
                alpha1, alpha2 = [gamma1 + np.arctan2(R[0,1], R[0,2])]*2
            else:
                beta1, beta2 = [-np.pi/2]*2
                alpha1, alpha2 = [gamma1 + np.arctan2(-R[0,1], -R[0,2])]*2
 
        alpha = np.append(alpha, alpha1)
        beta = np.append(beta, beta1)
        gamma = np.append(gamma, gamma1)
    return alpha,beta,gamma

alpha,beta,gamma = tait_angles(eigvec1, eigvec2, eigvec3)
#alphat = np.concatenate((alphat, alpha), axis=0)
#betat = np.concatenate((betat, beta), axis=0)
#gammat = np.concatenate((gammat, gamma), axis=0)
#%%

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_box_aspect([1,0.65,0.5])
#ax.set_title('Cell '+str(cellnum))
ax.set_xlabel('Angle/degrees'); ax.set_zlabel('Counts')
ax.tick_params(axis ='y', labelcolor=(1,1,1))
ax.set_xlim(xmin=-180,xmax=180)
nbins = 30
for y,lab,z in zip(np.array([alpha, beta, gamma])*180/np.pi, ['Roll (X)','Pitch (Y)','Yaw (Z)'], [0,1,2]):
    ys,bins = np.histogram(y, bins=nbins, range=(-180,180))
    xs = (bins[:-1] + bins[1:])/2
    ax.bar(xs, ys, zs=z, zdir='y', width=10, color=['r','g','k'][z] ,alpha=0.8, label=lab)
yax = ax.get_yaxis()
yax = yax.set_visible(False)
ax.legend()

# fig, ax = plt.subplots()
# ax.set_title('Cell '+str(cellnum))
# ax.set_xlabel('Angle/degrees'); ax.set_ylabel('Counts')
# nbins = 30
# for y,lab in zip(np.array([alpha, beta, gamma])*180/np.pi, ['Yaw (z)','Pitch (y)','Roll (x)']):
#     ys,bins = np.histogram(y, bins=nbins, range=(-180,180))
#     xs = (bins[:-1] + bins[1:])/2
#     ax.plot(xs, ys, label=lab)
# ax.legend()

#----------- ANGULAR DISPLACEMENT ---------------------
# def calc_angdisp(head_coord):
#     angdisp = []
#     for i in range(1,len(head_coord)-1):
#         vi = (head_coord[i]-head_coord[i-1])
#         vf = (head_coord[i+1]-head_coord[i])
#         dotprod = np.dot(vf/np.linalg.norm(vf), vi/np.linalg.norm(vi))
#         if -1<=dotprod<=1:
#             angdisp = np.append(angdisp, np.arccos(dotprod))
#         else:
#             print('The dot product is out of bounds for arccos')
#     return angdisp

# angdisp = calc_angdisp(head_coord)
# MAD = np.mean(angdisp)
# print('MAD_red=%.4f rad' %(MAD))

#----------- ANGULAR DISPLACEMENT ---------------------
# def calc_angdisp(eigvec):
#     angdisp = []
#     for i in range(len(eigvec)-1):
#         dotprod = np.dot(eigvec[i+1], eigvec[i]) #no need to divide by the norms bc they are 1
#         if -1<=dotprod<=1:
#             angdisp = np.append(angdisp, np.arccos(dotprod))
#         else:
#             print('The dot product is out of bounds for arccos')
#     return angdisp

# angdisp_red = calc_angdisp(eigvec1)

#----------- HEAD'S ANGLE OF ROTATION -----------------
angle = []
for k in range(1,n_points):
    Ri = np.column_stack([eigvec1[k-1],eigvec2[k-1],eigvec3[k-1]])
    Rf = np.column_stack([eigvec1[k],eigvec2[k],eigvec3[k]])
    R = Ri@Rf.T
    anglei = np.arccos((np.trace(R)-1)/2)
    if np.isnan(anglei):
        angle = np.append(angle, np.pi)#problems with 180º, gives nan
    else:
        angle = np.append(angle, anglei) 
    
angle_mean = np.mean(angle)

print('Angle mean: ', angle_mean)
fig_angle = plt.figure()
ax_angle = fig_angle.add_subplot()
ax_angle.plot(time_frame[1:], angle, 'orange', label='data')
ax_angle.hlines(angle_mean, time_frame[1], time_frame[-1], label='mean')
ax_angle.set_xlabel('Time/s')
ax_angle.set_ylabel('Angle/rad')
ax_angle.legend()

#angle = np.linspace(0,2*np.pi,n_points+2)[1:-1] #bc I get a nan and wrong values
#------------- RADIUS OF CURVE ----------------
r_x = []; r_y = []
for k in range(n_points-1):
    r_x = np.append(r_x, head_coord[k][0]/np.cos(angle[k]))
    r_y = np.append(r_y, head_coord[k][1]/np.sin(angle[k]))
    #print(head_coord[k][0], head_coord[k][1])
#print(np.cos(angle))

#----------- SPEED OF ROTATION -----------------


#----------- NOISE OF NECK MOVEMENT ------------


#%%
# Save results
motion_descrip = pd.DataFrame()
motion_descrip['speed'] = speed
motion_descrip['velocity'] = velocity.tolist()
motion_descrip['angle'] = angle
motion_descrip['r_x'] = r_x
motion_descrip['r_y'] = r_y

#%%
plt.scatter(VCLt_dat1, VSLt_dat1, c='tab:blue', marker='s', s=50, label='Dataset 1')
plt.scatter(VCLt_dat2, VSLt_dat2, c='tab:green', marker='d', s=50, label='Dataset 2')
plt.scatter(VCLt_dat2_resample, VSLt_dat2_resample, c='tab:orange', marker='*', s=50, label='Dataset 2 resampled')

plt.plot(VCLt_dat1, np.mean(VSLt_dat1/VCLt_dat1)*VCLt_dat1, c='tab:blue')
plt.plot(VCLt_dat2, np.mean(VSLt_dat2/VCLt_dat2)*VCLt_dat2, c='tab:green')
plt.plot(VCLt_dat2_resample, np.mean(VSLt_dat2_resample/VCLt_dat2_resample)*VCLt_dat2_resample, c='tab:orange')

plt.legend()
plt.xlabel('VCL/$\mu m\ s^{-1}$'); plt.ylabel('VSL/$\mu m\ s^{-1}$')
plt.ylim(0,5)