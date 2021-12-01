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

reqoptions.add_argument('-eig_file', dest="eigFile", required=True, help='The PCA eigenvectors estimated by AFNI\'s 3dpc' )
reqoptions.add_argument('-outfile', dest="outFile", required=False, help='DEPRECATED' )
reqoptions.add_argument('-motion', dest="motion", required=True, default="high", help='Motion file' )


args = parser.parse_args()


eigFile = args.eigFile


qualityGrade = 1


rp = np.loadtxt( args.motion )

drp     = np.concatenate((np.zeros((1,6), dtype=float), np.diff( rp, axis=0 )) )
drp_rad = np.concatenate((np.zeros((1,6), dtype=float), np.diff( np.deg2rad(rp), axis=0 )) )

fd_rot = abs(drp_rad[:,0]*50) + abs(drp_rad[:,1]*50) + abs(drp_rad[:,2]*50)
fd_mot = abs(drp[:,3]) + abs(drp[:,4]) + abs(drp[:,5])


fd = fd_rot + fd_mot

if np.mean(fd) > 0.25 and np.max(rp) > 3:
    qualityGrade = 3
elif np.mean(fd) > 0.25 or np.max(rp) > 3:
        qualityGrade = 2

# Somewhat by experience, higher motion subjects tend to
#if qualityGrade == 1:
#    thr = 0.175
#   ev  = 90
#elif qualityGrade == 2:
#    thr = 0.125
#    ev  = 90
#else:
thr = 0.15
ev  = 75


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
    if i > 30:

        val = 100 * (eigs[i,3] - eigs[i-1,3]  )
        modelOrder = i
        #print('{} - {:.5f} / {:.2f}%'.format(i, val, eigs[i,3]*100))
        # New components are explaining little variance (delta < 0.2%)
        # - OR -
        # 80% of total variance explained

        # High-motion/lower quality data sets likely require more stringent
        # parameters here.
        if val < thr:
            # Sometimes there is some "noise", so we want to be stably below this threshold
            lowChangeCount += 1


        if (lowChangeCount >= 3 and i > 80) or (eigs[i,3]*100) > ev or i > 150:
            break

if modelOrder == 0:
    modelOrder  = int( np.ceil(nComps/2) )

print( int(np.round(modelOrder)) )
#np.savetxt(args.outFile, int(modelOrder), fmt='%d')
#
#print('Extracting {} IC components responsible for {:.01f}% of total variance'.format(
#        int(modelOrder), eigs[modelOrder-1,3]*100 ))
