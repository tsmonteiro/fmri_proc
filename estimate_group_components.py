#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 16:59:58 2020

@author: u0101486
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 11:07:00 2019

@author: Thiago Monteiro
"""
import os.path
import os


from nilearn.decomposition import DictLearning
from nilearn.input_data import NiftiMapsMasker




baseDir = '/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/fix_sample_n/'
#baseDir='/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/GroupICA/'
files = []
#for t in range(3):
#    for i in range(75):
#        subDir = '/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/tmp/N{}_{:02d}/'.format(t,i)
#        subDir += 'proc_data_mni.nii'
#        
#       if os.path.exists(subDir):
#            files.append(subDir)
for file in os.listdir(baseDir):
    if file.endswith('.nii') and file.startswith('N'):
        files.append(baseDir + file)

maskFile = baseDir + 'mask.nii'

#%%

canica = DictLearning(n_components=35, smoothing_fwhm=None,
                n_jobs=12, alpha=10, n_epochs=1, t_r=1.6,
                low_pass=0.15, high_pass=0.01, detrend=True,
                verbose=10, standardize=True)




canica.fit(files, confounds=None)
# Retrieve the independent components in brain space. Directly
# accesible through attribute `components_img_`. Note that this
# attribute is implemented from version 0.4.1. For older versions,
# see note section above for details.
compNii = canica.components_img_
compNii.to_filename(baseDir + '/Group_Comps.nii')

#%%

masker = NiftiMapsMasker(compNii, standardize=True, detrend=True)
timeSeries = masker.fit_transform(files[0], confounds=None)
