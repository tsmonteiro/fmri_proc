#!/usr/bin/env python


import os
import argparse
import sys

import nibabel as nib
from builtins import str

import matplotlib.pyplot as plt
import numpy as np

import nipype.algorithms.confounds as npalg
import nilearn.plotting as nlp
import nilearn.image as nimg
import nilearn.signal as sgn


from nilearn import datasets
from nilearn.input_data import NiftiMapsMasker
from nilearn.connectome import ConnectivityMeasure



parser = argparse.ArgumentParser(description='Save QA check Plots')

# Required options                    
reqoptions = parser.add_argument_group('Required arguments')

reqoptions.add_argument('-o', '-out', dest="outDir", required=True, help='Directory where images are to be saved' )
reqoptions.add_argument('-a', '-in', dest="inDir", required=True, help='Dir where EPI + masks are stored [MNI SPACE]' )
reqoptions.add_argument('-m', '-mot', dest="motest", required=True, help='Motion estimate file' )
reqoptions.add_argument('-f', '-func', dest="funcImg", required=True, help='Functional Data' )



args = parser.parse_args()


outDir = args.outDir
inDir = args.inDir
mFile = args.motest
rp = np.loadtxt(inDir + '/' + mFile)
funcImgFile = inDir + '/' + args.funcImg


funcImg = nib.load(funcImgFile)
data = np.array(funcImg.get_data())

X,Y,Z,N = data.shape
nVox = X*Y*Z
# Reshape to voxels x time
data = np.reshape(data, (X*Y*Z, N))

dvar =  np.sqrt( np.sum( np.power( np.diff( data, axis=1 ), 2 ), axis=0 ) )
dvar = 100 * dvar / np.median(dvar)
dvar = np.insert( dvar, 0, 100 )


fd1 =  np.abs(np.diff(np.deg2rad( rp[:,0] ))*50) + np.abs(np.diff(np.deg2rad( rp[:,1] ))*50) + np.abs(np.diff(np.deg2rad( rp[:,2] ))*50) 
fd2 =  np.abs(np.diff(rp[:,3] )) + np.abs(np.diff(rp[:,4] )) + np.abs(np.diff(rp[:,5] ))

fd1 = np.insert(fd1, 0, 0)
fd2 = np.insert(fd2, 0, 0)


arp = np.max( np.sum( np.abs(rp[:,3:5]), axis=1) )

fd = fd1 + fd2


mot = np.concatenate(( [np.mean(fd,axis=0)],[arp],[np.std(fd,axis=0)],[np.std(dvar,axis=0)]))

print(mot)

# Save framewise displacement and DVARs numbers
# Can be used later for group QA summary
np.savetxt(outDir+'/motion_summary.txt',mot, delimiter=',', fmt='%.3f')




