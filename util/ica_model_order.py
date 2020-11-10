#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 22 15:36:33 2020

@author: u0101486
"""

import numpy as np

import argparse

parser = argparse.ArgumentParser(description='Estimate the number of ICA components to be estimated by MELODIC')

# Required options
reqoptions = parser.add_argument_group('Required arguments')

reqoptions.add_argument('-eig_file', dest="eigFile", required=True, help='The PCA eigenvectors estmiated by AFNI\'s 3dpc' )
reqoptions.add_argument('-outfile', dest="outFile", required=False, help='DEPRECATED' )


args = parser.parse_args()


eigFile = args.eigFile


#%%
#eigFile = '/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/tmp/B1_31/pc_var_eig.1D'
# The eig file is composed of 4 columns
# #Num.  --Eigenvalue--  -Var.Fraction-  -Cumul.Fract.-
eigs = np.loadtxt( eigFile )


val = []

nComps = eigs.shape[0]

modelOrder  = 0
sumTotal    = 0

lowChangeCount = 0
for i in range( nComps ):
    if i > 9:
        m0 = np.mean( eigs[i-10:i-1,1]  )
        #val = 100 * eigs[i,1]/m0
        #val = 100 * (eigs[i,3] - np.median( eigs[i-8:i-1,3], axis=0  ))
        val = 100 * (eigs[i,3] - eigs[i-1,3]  )
        modelOrder = i
        #print('{} - {:.5f} / {:.2f}%'.format(i, val, eigs[i,3]*100))
        # New components are explaining little variance (delta < 0.3%)
        # - OR -
        # 80% of total variance explained
        if val < 0.3:
            # Sometimes there is some "noise", so we want to be stably below this threshold
            lowChangeCount += 1
        

        if (lowChangeCount >= 3 and (eigs[i,3]*100) > 50) or (eigs[i,3]*100) > 80:
            break

if modelOrder == 0:
    modelOrder  = int( np.ceil(nComps/2) )

print( int(np.round(modelOrder)) )
#np.savetxt(args.outFile, int(modelOrder), fmt='%d')
#
#print('Extracting {} IC components responsible for {:.01f}% of total variance'.format(
#        int(modelOrder), eigs[modelOrder-1,3]*100 ))
