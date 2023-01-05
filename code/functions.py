# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 16:25:44 2022

@author: Alicia Mosquera
"""
import numpy as np
import seaborn as sns
from skimage import io
from numpy import newaxis
from numpy.matlib import repmat
import matplotlib.pyplot as plt

def read(filename):
    volume_time = io.imread(filename)
    t,z,y,x = volume_time.shape
    return volume_time,t,z,y,x

#%%
def scale(I,s,dr,dc):
    """ 
    Gaussian Scale-Space using the Fourier Domain.
    
    Parameters
    ----------
    I : array of size (rows,cols)
        The Fourier transform of a matrix or an image.
    s : float
        The standard deviation of the Gaussian (must be larger than 0).
    dr : float
        The derivative order in the direction of the rows.
    dc : float
        The derivative order in the direction of the columns.

    Returns
    -------
    Is : array of size (rows,cols)
        The Fourier transform of the gaussian scaled matrix or image.

    -------
    This is an implementation of the Gaussian scale space on matrices.
    The convolution is implemented in the Fourier domain and for that
    reason the number of rows and columns of the matrix must be powers
    of 2.
    Fractional valued dr and dc are possible, but be warned the result
    will probably be complex.
    The complexity of this algorithm is O(n) where n is the total number
    of elements of I.
    
    To calculate an image of scale (variance) 2^2 use,
    I2 = real(ifft2(scale(fft2(I),2,0,0)));
    To derive an image once in the direction of rows at scale 1^2 do,
    I2 = real(ifft2(scale(fft2(I),1,1,0)));
    
    Copyright: Jon Sporring, January 1, 1996
    """
    # Expand dimensions in case I is a 1D horizontal array, so it has to return rows=1
    rows,cols = I[newaxis,:].shape[-2:]
    
    if s < 0:
      raise ValueError('s must be larger than or equal to 0')
      
    elif s == 0:
        if (dr == 0) and (dc == 0):
            Is = I
        else:
            Is = np.zeros([rows,cols])
            
    elif s == float('inf'):
        Is = np.zeros([rows,cols])
        Is[0,0] = I[0,0]

    else:
        G = np.zeros([rows,cols])
        if ((rows != 1) and (rows%2 != 0)) or ((cols != 1) and (cols%2 != 0)):
            raise ValueError('The matrix must have side lengths that are either 1 or even.')
        else:
            
            # Calculate the Fourier transform of a gaussian fct.
            if (rows > 1) and (cols > 1):
                mat1 = repmat((np.arange(rows/2+1)[newaxis,:].T/rows)**2, 1, cols//2+1)
                mat2 = repmat((np.arange(cols/2+1)/cols)**2, rows//2+1, 1)
                G[:rows//2+1, :cols//2+1] = np.exp(-(mat1+mat2)*(s*2*np.pi)**2/2)
                G[rows//2:, :cols//2+1] = np.flipud(G[1:rows//2+1, 0:cols//2+1])
                G[:rows//2+1, cols//2:] = np.fliplr(G[0:rows//2+1, 1:cols//2+1])
                G[rows//2:, cols//2:] = np.fliplr(np.flipud(G[1:rows//2+1, 1:cols//2+1]))
            else:
                val = np.max([rows,cols])
                G[:val/2+1] = np.exp(-(np.arange(val/2+1)[newaxis,:].T/val)**2*(s*2*np.pi)**2/2)
                if rows > 1:
                    G[val//2:] = np.flipud(G[1:val//2+1])
                else:
                    G[val//2:] = np.fliplr(G[1:val//2+1])
                
            # Calculate the Differentiation matrix
            if (rows > 1) and (cols > 1):
                x = np.hstack((np.arange(rows/2), np.arange(-rows/2,0,1)))/rows
                y = np.hstack((np.arange(cols/2), np.arange(-cols/2,0,1)))/cols
                DG = (x**dr)[newaxis,:].conj().T*(y**dc)*(2j*np.pi)**(dr+dc)
            else:
                if rows > 1:
                    x = np.hstack((np.arange(rows/2), np.arange(-rows/2,0,1)))[newaxis,:]/rows
                    DG = (2j*np.pi*x.T)**dr
                else:
                    y = np.hstack((np.arange(cols/2), np.arange(-cols/2,0,1)))/cols
                    DG = (2j*np.pi*y)**dc
                    
            Is = I*G*DG
    return Is

#%%
def confusionMatrix(true, pred):
    if true.shape != pred.shape:
        raise ValueError('both images must have the same dimensions')
    true = true.flatten()
    pred = pred.flatten()
    
    k = len(np.unique(true))
    if k == 1:
        raise ValueError('there must be more than one outcome in the true image')
    else:
        cfmatrix = np.zeros((k,k))
        for i in range(len(true)):
            cfmatrix[int(true[i]), int(pred[i])] += 1

    return cfmatrix

def paramConfusionMatrix(cfmatrix):
    TN = cfmatrix[0,0]
    FP = cfmatrix[0,1]
    FN = cfmatrix[1,0]
    TP = cfmatrix[1,1]
    
    if (TP+FP)==0 or (TP+FN)==0:
        precision = 0
        recall = 0
    else:
        precision = TP/(TP+FP)
        recall = TP/(TP+FN)
    
    return precision, recall

def plotConfusionMatrix(ax, cf_matrix, case):
    group_names = ['True Neg','False Pos','False Neg','True Pos']
    group_counts = ["{0:0.0f}".format(value) for value in cf_matrix.flatten()]
    group_percentages = ["{0:.2%}".format(value) for value in cf_matrix.flatten()/np.sum(cf_matrix)]
    labels = [f"{v1}\n{v2}\n{v3}" for v1, v2, v3 in zip(group_names,group_counts,group_percentages)]
    labels = np.array(labels).reshape(2,2)

    sns.heatmap(cf_matrix, annot=labels, fmt='', cmap='Blues')

    ax.set_xlabel('\n'+case+' segmentation')
    ax.set_ylabel('Manual segmentation');

    # Ticket labels - List must be in alphabetical order
    ax.xaxis.set_ticklabels(['Background','Cell'])
    ax.yaxis.set_ticklabels(['Background','Cell'])

    # Display the visualization of the Confusion Matrix.
    plt.show()
    return 0

#%%
def reshapeToEvenSides(I): 
    if np.any(np.array(I.shape)%2 != 0):
        # cropping out the last rows/columns is important, rather than the beginning ones
        if I.shape[0]%2 != 0:
            I = I[:-1]
        if (len(I.shape)>1) and (I.shape[1]%2 != 0):
            I = I[:,:-1]
        if (len(I.shape)>2) and (I.shape[2]%2 != 0):
            I = I[:,:,:-1]
        if (len(I.shape)>3) and (I.shape[3]%2 != 0):
            I = I[:,:,:,:-1]
    return I

def plot_yz(fullvolume_time, t, x):
    _, z, y, _ = fullvolume_time.shape
    img = np.zeros((z,y), dtype=np.uint16)
    for zi in range(z):
        for yi in range(y):
            img[zi][yi] = fullvolume_time[t][zi][yi][x]
    return img

def plot_zx(fullvolume_time, t, y):
    _, z, _, x = fullvolume_time.shape
    img = np.zeros((x,z), dtype=np.uint16)
    for zi in range(z):
        for xi in range(x):
            img[xi][zi] = fullvolume_time[t][zi][y][xi]
    return img