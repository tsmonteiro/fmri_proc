#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 14:28:03 2019

@author: u0101486
"""

# Bootstrap model order selection for ICA


import nibabel as nib
import numpy as np
import nilearn.signal as sgn


from matplotlib import pyplot as plt


TR = 2.5

funcData = nib.load('/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/RS010/proc_data_native.nii')
funcData = funcData.get_fdata()

nX, nY, nZ, nT = funcData.shape


funcData = np.reshape(funcData, (nX*nY*nZ, nT))

meanFunc = np.mean(funcData, axis=1)
funcData  = funcData[meanFunc>200,:]

funcData = sgn.clean(np.transpose(funcData), detrend=True, standardize=True,
                     confounds=None, low_pass=None, high_pass=None, t_r=TR, ensure_finite=False)
funcData = np.transpose(funcData)

#%%


from sklearn.decomposition import PCA
import scipy.stats as sst


pca = PCA(svd_solver='auto', n_components=150)
prinComps = pca.fit_transform(funcData)

pcFit = pca.fit(funcData)

def entropy(X):
    #value,counts = np.unique(labels, return_counts=True)
    h = np.histogram(X, bins=25, normed=True)
    #sst.entropy(h[0], base=None)
    return np.rank(X)#sst.entropy(h[0], base=None)



compEntropy = []

for c in range(150):
    smoothComp = prinComps[:,c]
    smoothCompA = np.where((smoothComp)<.5, 0, smoothComp)
    e1 = entropy(smoothCompA)
    smoothCompB = np.where((smoothComp)<.5, smoothComp, 0)
    e2 = entropy(smoothCompB)
    
    if e2 > 0:
        compEntropy.append(e1 / e2)

    
plt.plot(compEntropy)
