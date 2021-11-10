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

from PIL import Image, ImageDraw





# ++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# END OF function definitions
#
# ++++++++++++++++++++++++++++++++++++++++++++++


parser = argparse.ArgumentParser(description='Save QA check Plots')

# Required options                    
reqoptions = parser.add_argument_group('Required arguments')

reqoptions.add_argument('-o', '-out', dest="outDir", required=True, help='Directory where images are to be saved' )
reqoptions.add_argument('-a', '-in', dest="inDir", required=True, help='Dir where EPI + masks are stored [MNI SPACE]' )

reqoptions.add_argument('-x', '-outf', dest="outFname", required=False, default='greyplot.png', help='Name of the output plot' )


reqoptions.add_argument('-k', '-im1', dest="im1", required=False,  help='Image 1 in normalisation check' )
reqoptions.add_argument('-l', '-im2', dest="im2", required=False,  help='Image 2 in normalisation check' )



reqoptions.add_argument('-d', '-dpi', dest="dpi", required=False, default=120, help='Saved figure DPI' )


args = parser.parse_args()


outDir = args.outDir
inDir = args.inDir



# PNG resolution of the saved file
figDpi=int(args.dpi)

outFile = outDir + '/' + args.outFname

outGFile = outDir + '/' + args.outFname + '.gif'



# Font size and weight for ALL plots
plt.rcParams.update({'font.size': 20, 'font.weight':'bold'} )




# ++++++++++++++++++++++++++++++++++++++++++++++++
#
# END OF PARAMETER/SETUP
#
# ++++++++++++++++++++++++++++++++++++++++++++++++





baseImg = nib.load(args.im1)
data    = np.array(baseImg.get_data())

if len(data.shape) > 3:
	data = np.mean( data, 3)

baseImg = nimg.new_img_like(baseImg, data)


overlayImg = nib.load(args.im2)
data    = np.array(overlayImg.get_data())

if len(data.shape) > 3:
	data = np.mean( data, 3)

overlayImg = nimg.new_img_like(overlayImg, data)



zCuts=(-34,58)
xCuts=(-40,40)
yCuts=(-80,40)





fig = plt.figure(figsize=(30,6), dpi=figDpi, facecolor='w', edgecolor='k')

# Display the different planes
cuts = np.linspace(zCuts[0], zCuts[1], 20)

display1 = nlp.plot_stat_map(overlayImg, bg_img=baseImg, cmap='jet', figure=fig, cut_coords=cuts, display_mode='z', draw_cross=False, threshold=20, vmax=125,  \
		axes=(0,0,1,.33), symmetric_cbar=False, alpha=0.5 )


display1.add_contours(overlayImg, filled=False, levels=[20, 70], colors='k')


cuts = np.linspace(xCuts[0], xCuts[1], 20)

display2 = nlp.plot_stat_map(overlayImg, bg_img=baseImg, cmap='jet', figure=fig, cut_coords=cuts, display_mode='x', draw_cross=False, threshold=20, vmax=125, \
		axes=(0,.33,1,.33), symmetric_cbar=False, alpha=0.5)

display2.add_contours(overlayImg, filled=False, levels=[20, 70], colors='k')

cuts = np.linspace(yCuts[0], yCuts[1], 20)

display3 = nlp.plot_stat_map(overlayImg, bg_img=baseImg, cmap='jet', figure=fig, cut_coords=cuts, display_mode='y', draw_cross=False, threshold=20, vmax=125, \
		axes=(0,.66,1,.33), symmetric_cbar=False, alpha=0.5)

display3.add_contours(overlayImg, filled=False, levels=[20, 70], colors='k')
plt.savefig(outFile)

	
