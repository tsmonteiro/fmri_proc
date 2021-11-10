# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 11:07:00 2019

@author: Thiago Monteiro
"""
import os.path
import os

from nilearn.decomposition import CanICA
from nilearn.decomposition import DictLearning

import numpy as np
import nibabel as nib


from nilearn.input_data import NiftiMapsMasker
from scipy import signal
from scipy import stats

import argparse

parser = argparse.ArgumentParser(description='Run CanICA estimation [with some melodic structure]')

# Required options
reqoptions = parser.add_argument_group('Required arguments')

reqoptions.add_argument('-o', '-out', dest="outDir", required=True, help='Directory where images are to be saved' )
reqoptions.add_argument('-i', '-in',  dest="inFile", required=True, help='EPI file' )
reqoptions.add_argument('-m', '-mask',  dest="mask", required=True, help='EPI Mask [func space]' )
reqoptions.add_argument('-n', '-nIC',  dest="nIC", required=True, help='Number of estimated ICs' )
reqoptions.add_argument('-decomp',  dest="decomp", required=False, default='canica', help='Number of estimated ICs' )
reqoptions.add_argument('-dict_init',  dest="dict_init", required=False, default=None, help='Number of estimated ICs' )

args = parser.parse_args()


subDir        = args.outDir #baseDir + '/RS011'

nComps        = int(args.nIC)
maskFile      = args.mask #subDir + '/nat_mask.nii'
fName         = args.inFile #subDir + '/tmp_melodic.nii'


#if args.decomp == 'dictlearning':
#    if args.dict_init == None:
#        canica = DictLearning(n_components=nComps, smoothing_fwhm=None,
#                    mask=maskFile, n_jobs=12, alpha=3, n_epochs=2,
#                    verbose=10, standardize=True)
#    else:
#        canica = DictLearning(n_components=nComps, smoothing_fwhm=None,
#                    mask=maskFile, n_jobs=12, alpha=4, n_epochs=1,
#                    verbose=10, standardize=True, dict_init=args.dict_init)
#
#else:
canica = CanICA(n_components=nComps, smoothing_fwhm=None,
                mask=maskFile, do_cca=True, n_jobs=None, detrend=False,
               threshold=None, verbose=10, standardize=True)


canica.fit(fName, confounds=None)
# Retrieve the independent components in brain space. Directly
# accesible through attribute `components_img_`. Note that this
# attribute is implemented from version 0.4.1. For older versions,
# see note section above for details.
compNii = canica.components_img_
compNii.to_filename(subDir + '/CanICA_Comps.nii')


# TODO Run melodic statistica estimation on this (see last section here https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/MELODIC#Stats)
# TODO Manually clean the image and check the results (possibly compare to the melodic output)

compData = compNii.get_fdata()
maskData = nib.load(maskFile).get_fdata()
nX, nY, nZ, nIC = compData.shape
nVox = nX*nY*nZ
for ic in range(nIC):
    comp = compData[:,:,:,ic]
    comp = np.reshape( comp, (nVox) )


    comp = comp - stats.mode(comp)[0][0]


    comp = (comp - np.mean(comp))/(np.std(comp))



    compData[:,:,:,ic] = np.reshape( comp, (nX, nY, nZ))


zcompNii = nib.Nifti1Image(compData, compNii.header.get_base_affine())
zcompNii.to_filename(subDir + '/CanICA_Comps_Z.nii')
#%%
if args.decomp == 'canica':
    icVar = canica.variance_ - np.min(canica.variance_)
    icVar = icVar/np.sum(icVar)

    np.savetxt(subDir + '/eigenvalues_percent', np.cumsum(icVar),
           fmt='%0.6f', newline=' ')



#masker = NiftiMapsMasker(compNii, standardize=True, detrend=False)
#timeSeries = masker.fit_transform(fName, confounds=None)
##%%
#
#os.mkdir(subDir + '/report')
#
#for ic in range(nComps):
#    np.savetxt(subDir + '/report/t' + str(ic+1) + '.txt', timeSeries[:,ic],
#           fmt='%0.6f', newline=' ')
#    
#    freqs, psd = signal.welch(timeSeries[:,ic])
#    
#    np.savetxt(subDir + '/report/f' + str(ic+1) + '.txt', psd,
#           fmt='%0.6f', newline=' ')

# Write the melodic_mix, so to speak
#np.savetxt(subDir + '/melodic_mix', timeSeries,
#       fmt='%0.9f', newline=os.linesep, encoding='ascii')
