#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 21 08:47:56 2020

@author: u0101486
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 11:07:00 2019

@author: Thiago Monteiro
"""
import os.path
import os


import numpy as np
import nibabel as nib


from nilearn.input_data import NiftiMapsMasker
from scipy import signal


import argparse

parser = argparse.ArgumentParser(description='Run CanICA estimation [with some melodic structure]')

# Required options
reqoptions = parser.add_argument_group('Required arguments')

reqoptions.add_argument('-o', '-out', dest="outDir", required=True, help='Directory where images are to be saved' )
reqoptions.add_argument('-i', '-in',  dest="inFile", required=True, help='EPI file' )
reqoptions.add_argument('-icmap',  dest="icMap", required=True, help='EPI Mask [func space]' )
reqoptions.add_argument('-nic',  dest="nIC", required=True, help='EPI Mask [func space]' )


args = parser.parse_args()


subDir        = args.outDir #baseDir + '/RS011'

nComps        = int(args.nIC)
compFile      = args.icMap #subDir + '/nat_mask.nii'
fName         = args.inFile #subDir + '/tmp_melodic.nii'



compNii = nib.load(compFile)
masker = NiftiMapsMasker(compNii, standardize=True, detrend=False)
timeSeries = masker.fit_transform(fName, confounds=None)
#%%

if not os.path.exists(subDir + '/report'):
    os.mkdir(subDir + '/report')

for ic in range(nComps):
    np.savetxt(subDir + '/report/t' + str(ic+1) + '.txt', timeSeries[:,ic],
           fmt='%0.6f', newline=' ')
    
    freqs, psd = signal.welch(timeSeries[:,ic])
    
    np.savetxt(subDir + '/report/f' + str(ic+1) + '.txt', psd,
           fmt='%0.6f', newline=' ')

np.savetxt(subDir + '/melodic_mix', timeSeries,
       fmt='%0.9f', newline=os.linesep, encoding='ascii')