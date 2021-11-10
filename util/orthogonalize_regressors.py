#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 12:05:42 2020

Creates and orthogonalizes regressors prior to 3dREMLfit

@author: u0101486
"""



import numpy as np

import os

from sklearn.linear_model import TheilSenRegressor
from scipy.stats import zscore
import scipy.signal as sig

import argparse


def whiten(X, method='zca'):
    """
    Whitens the input matrix X using specified whitening method.
    Inputs:
        X:      Input data matrix with data examples along the first dimension
        method: Whitening method. Must be one of 'zca', 'zca_cor', 'pca',
                'pca_cor', or 'cholesky'.
    """
    X = X.reshape((-1, np.prod(X.shape[1:])))
    X_centered = X - np.mean(X, axis=0)
    Sigma = np.dot(X_centered.T, X_centered) / X_centered.shape[0]
    W = None
    
    if method in ['zca', 'pca', 'cholesky']:
        U, Lambda, _ = np.linalg.svd(Sigma)
        if method == 'zca':
            W = np.dot(U, np.dot(np.diag(1.0 / np.sqrt(Lambda + 1e-5)), U.T))
        elif method =='pca':
            W = np.dot(np.diag(1.0 / np.sqrt(Lambda + 1e-5)), U.T)
        elif method == 'cholesky':
            W = np.linalg.cholesky(np.dot(U, np.dot(np.diag(1.0 / (Lambda + 1e-5)), U.T))).T
    elif method in ['zca_cor', 'pca_cor']:
        V_sqrt = np.diag(np.std(X, axis=0))
        P = np.dot(np.dot(np.linalg.inv(V_sqrt), Sigma), np.linalg.inv(V_sqrt))
        G, Theta, _ = np.linalg.svd(P)
        if method == 'zca_cor':
            W = np.dot(np.dot(G, np.dot(np.diag(1.0 / np.sqrt(Theta + 1e-5)), G.T)), np.linalg.inv(V_sqrt))
        elif method == 'pca_cor':
            W = np.dot(np.dot(np.diag(1.0/np.sqrt(Theta + 1e-5)), G.T), np.linalg.inv(V_sqrt))
    else:
        raise Exception('Whitening method not found.')

    return np.dot(X_centered, W.T)


parser = argparse.ArgumentParser(description='Run CanICA estimation [with some melodic structure]')

# Required options
reqoptions = parser.add_argument_group('Required arguments')

reqoptions.add_argument('-regressors',  dest="regFile", required=True, help='Directory where images are to be saved' )
reqoptions.add_argument('-signal',  dest="sigFile", required=False, default=None, help='EPI file' )
reqoptions.add_argument('-out',     dest="outFile", required=True, help='EPI file' )

args = parser.parse_args()

regFile = args.regFile
sigFile = args.sigFile
outFile = args.outFile
#%%

#import matplotlib.pyplot as plt

#regFile = '/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/RS001/regressors_FAC_CC_RP24.txt'
#sigFile = '/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/RS001/non_noise_ics.txt'



regressors = np.loadtxt(regFile)

nr = regressors.shape[1]


#for r in range(nr):
    #regressors[:,r] = sig.detrend( regressors[:,r] )
    #regressors[:,r] = zscore( regressors[:,r] )
    
from sklearn.linear_model import LinearRegression

#print(regressors.shape)
if not sigFile == None:
    sigg = np.loadtxt( sigFile )
    
    ns = sigg.shape[1]
        
    #for r in range(ns):
        #sigg[:,r] = sig.detrend( sigg[:,r] )
        #sigg[:,r] = zscore( sigg[:,r] )
    
    #print(sig.shape)
    
    # Removing shared variance
    print('======================================================')
    print('Removing shared variance')
    print('======================================================')
    
    for r in range(nr):
        #r = 45
        y = np.copy(regressors[:,r])
        X = np.copy(sigg)
        
        
        #reg =  SGDRegressor(loss='huber',penalty='elasticnet',alpha=0.00001).fit(X,y)
        reg =  LinearRegression().fit(X,y)
        print('{:.01f} - {:.02f}'.format(reg.score(X,y),np.std(y-reg.predict(X))))
    
        # Residuals
        #plt.close('all')
        #plt.plot(regressors[:,r])
        #plt.plot(y-reg.predict(X))
        #break
        regressors[:,r] = y-reg.predict(X)
        
    
        
        #plt.plot(regressors[:,r], color=(0,0,0))
        #break

np.savetxt(outFile, regressors, fmt='%0.5f')
#%% Orthogonalize

#import scipy.stats as ssp
#nr = regressors.shape[1]
#covMat  = np.zeros((nr,nr))
#wcovMat = np.zeros((nr,nr))

#print( idx )
#wRegressors = whiten(regressors, method='cholesky')

#nRegressors = []
#
#for i in range(nr):
#    print(np.std(wRegressors[:,i], axis=0))
#    #wRegressors[:,i] = wRegressors[:,i] - np.mean(wRegressors[:,i])
#    if np.std(wRegressors[:,i], axis=0) > 0.1:
#        nRegressors.append( wRegressors[:,i] )
#    for j in range(nr):
#        r,p = ssp.pearsonr(regressors[:,i], regressors[:,j])
#        covMat[i,j] = r
#       covMat[j,i] = r
#        
#        r,p = ssp.pearsonr(wRegressors[:,i], wRegressors[:,j])
#        wcovMat[i,j] = r
#       wcovMat[j,i] = r
#
#plt.subplot(121)
#plt.imshow(covMat, vmax=1, vmin=-1)
#plt.colorbar()

#plt.subplot(122)
#plt.imshow(wcovMat, vmax=1, vmin=-1)
#plt.colorbar()

#nRegressors = np.transpose(nRegressors)
#np.savetxt(outFile, np.transpose(pRegressors), fmt='%0.5f')
#np.savetxt(outFile, nRegressors, fmt='%0.5f')