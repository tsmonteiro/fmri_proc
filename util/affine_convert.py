#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 14:18:52 2020

@author: u0101486
"""



import sys
import numpy as np

import nitransforms.linear as nitl

sys.path.append('/home/luna.kuleuven.be/u0101486/workspace/fmri_proc/ext/nitransforms/')


inAffine = '/home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/tmp/sub001/anat/anat2func.mat'
ref='/home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/tmp/sub001/anat/anat_proc.nii'
moving='/home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/tmp/sub001/FIX/mean_func.nii.gz'

aff = np.loadtxt(inAffine)
nAff= nitl.Affine(aff, reference=ref)
nAff.to_filename('/home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/tmp/sub001/anat/anat2func_fsl.mat',
                 fmt='fsl', moving=moving)

