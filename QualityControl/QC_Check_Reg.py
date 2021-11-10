#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 11:59:11 2020

@author: u0101486
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 12:26:49 2019

@author: u0101486
"""

# Aggregate QC measures

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

import configparser

config = configparser.ConfigParser()
config.read('/home/luna.kuleuven.be/u0101486/workspace/fmri_proc/params.ini')

#PATH to qclib
sys.path.append(config['PATHS']['QCLIB_PATH'])


import qclib.group_plots as gp
import qclib.common_funcs as cf

import nibabel as nib
from nilearn import plotting
import nilearn.image as nimg

#TODO Make it easier to switch between projects

#project = 'CRUNCH'
project = 'RepImpact'
#project = 'CAI_China'

baseDir       = '/home/luna.kuleuven.be/u0101486/workspace/data/' + project + '/tmp/'
outDir      = '/home/luna.kuleuven.be/u0101486/workspace/data/' + project + '/Quality_Control/Reg/'


if not os.path.isdir(outDir):
    os.mkdir(outDir)


i = 0
for sub in sorted(os.listdir(baseDir)):
    #print(sub)
    subDir        = baseDir + '/' + sub
    t1File        = subDir + '/anat/anat_proc.nii'
    meanFuncFile  = subDir + '/anat/mean_func_data_nds.nii'
        
    
    if os.path.isfile(t1File) and os.path.isfile(meanFuncFile):
        plt.close('all')
        fig = plt.figure(figsize=(16,8), dpi=200, facecolor='k', edgecolor='k')
        print('Reading ' + sub)
        anatImage = nib.load(t1File)
        anatData = anatImage.get_fdata()
        anatData = anatData / np.quantile(anatData, .9)
        
        anatImage = nib.Nifti1Image(anatData, anatImage.affine, header=anatImage.header)
        
        
        funcImage = nib.load(meanFuncFile)
        funcData = np.array( funcImage.get_fdata(), dtype=np.float64)
        funcData = np.array( funcData / np.quantile(funcData, .99) )
        
        funcImage = nib.Nifti1Image(funcData, affine=funcImage.affine, header=funcImage.header)
        
        
        
        sp1=plt.subplot(311)
        display = plotting.plot_anat(anatImage, title=sub, draw_cross=False, vmin=0.1, vmax=0.95, axes=sp1,
                                     display_mode='z', cut_coords=np.linspace(0,80,8), black_bg=True)
        #plotting.plot_epi(funcImage, draw_cross=False, vmin=.200, vmax=.900, figure=fig,  cmap='gray')
        display.add_overlay(funcImage, threshold=0.35, cmap='gist_rainbow', alpha=0.95)
        display.add_contours(anatImage, levels=np.quantile(anatData, [.7,.9]), colors='k', alpha=0.8, linewidths=0.8)
        
        sp2=plt.subplot(312)
        display = plotting.plot_anat(anatImage, title=sub, draw_cross=False, vmin=0.1, vmax=0.95, axes=sp2,
                                     display_mode='x', cut_coords=np.linspace(-60,60,8), black_bg=True)
        #plotting.plot_epi(funcImage, draw_cross=False, vmin=.200, vmax=.900, figure=fig,  cmap='gray')
        display.add_overlay(funcImage, threshold=0.35, cmap='gist_rainbow', alpha=0.95)
        display.add_contours(anatImage, levels=np.quantile(anatData, [.7,.9]), colors='k', alpha=0.8, linewidths=0.8)
        
        sp3=plt.subplot(313)
        display = plotting.plot_anat(anatImage, title=sub, draw_cross=False, vmin=0.1, vmax=0.95, axes=sp3,
                                     display_mode='y', cut_coords=np.linspace(-80,80,8), black_bg=True)
        #plotting.plot_epi(funcImage, draw_cross=False, vmin=.200, vmax=.900, figure=fig,  cmap='gray')
        display.add_overlay(funcImage, threshold=0.35, cmap='gist_rainbow', alpha=0.95)
        display.add_contours(anatImage, levels=np.quantile(anatData, [.7,.9]), colors='k', alpha=0.8, linewidths=0.8)
        
        plt.savefig(outDir + sub + '.png')
        plt.close('all')
        i += 1

        
        
        
        
        
        #imY = axY.imshow(np.transpose(anatImage[:,curY,:]), cmap='gray',  alpha=None,
        #               vmin=0.1, vmax=1.35,  aspect='auto' )

            
       #ovY = axY.imshow(maskedIc, cmap='hot',  interpolation=None, alpha=None,  aspect='auto', vmin=vmin, vmax=vmax)
       
    
    # END for sub sorted(....)    
    
