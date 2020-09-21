#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 14:01:14 2020

@author: u0101486
"""

from nilearn import surface
import nibabel as nib
import numpy as np
from nilearn.image import resample_img
np.set_printoptions(precision=3, suppress=True)
fmri_img=  '/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/RS001/proc_data_native.nii'
t1_img=  '/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/RS001/t1_frn.nii'
fmri_a_img =  '/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/RS001/00_anat_proc_data_native.nii'
func2anat = '/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/RS001/afunc2anat.mat'
nat_mask=  '/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/RS001/nat_mask.nii'
surf  = '/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/RS001/surf/s15.mesh.thickness.resampled.t1_f.gii'

# TODO pra applicar essa transformacao, melhor usar ANTs em cada volume tendo como referencia
# o espaco da imagem funcional [verificar se funciona, talvez padding seja necessario]
nii2 = nib.load(fmri_img)
print(nii2.affine)
print('')
f2a = np.loadtxt(func2anat)
print(f2a)
print('')
print(np.dot(nii2.affine, f2a))
print('')

#for i in range(3):
#    f2a[i,i] = nii2.affine[i,i]


nii = nib.Nifti1Image(nii2.get_data(), np.dot(nii2.affine, f2a))


#for i in range(3):
#    nii.affine[i,i] = nii2.affine[i,i]

print(nii.header)

nii.to_filename(fmri_a_img)

#%%

texture = surface.vol_to_surf(fmri_img, surf, kind='ball', radius=2.5)
#%%
from nilearn import plotting


delt_score = texture[:,2]-texture[:,1]

thickSurf = '/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/RS001/surf/lh.central.t1_f.gii'

#plotting.plot_surf(thickSurf, output_file='/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/RS001/surf/lh.central.t1_f.png')


plotting.plot_surf_stat_map(surf, delt_score, hemi='left',  colorbar=True,output_file='/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/RS001/surf/lh.central.t1_f.png',
                            vmax=60)