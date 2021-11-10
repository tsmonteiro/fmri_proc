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

import configparser

config = configparser.ConfigParser()
config.read('params.ini')



sys.path.append(config['PATHS']['qclib_path'])

import qclib.motion_handler as mh
import qclib as qc





def plot_brain_overlay(data, bgImg, rangeVal, figDpi, slAxis='z', thr=0, is_stat=False):
	fig = plt.figure(figsize=(16,12), dpi=figDpi, facecolor='w', edgecolor='k')

	# z = Axial
	# y = coronal
	# x = sagital


	if slAxis == 'z':
		zSlices  = np.linspace(-60, 80, 60)

	if slAxis == 'x':
		zSlices  = np.linspace(-60, 60, 60)
	
	if slAxis == 'y':
		zSlices  = np.linspace(-100, 60, 60)

	nRow = 7
	rowPlace = np.linspace(0, 1, nRow + 1)

	step= int( np.floor( 60/nRow ) )
	t0 = 0
	tf = step

	for i in range(nRow):
		if is_stat == True:
			nlp.plot_stat_map(data, display_mode=slAxis, figure=fig, draw_cross=False,  vmax=.9, cmap='cold_hot', \
					cut_coords=zSlices[t0:tf], axes=(0,rowPlace[i],1,1/nRow), colorbar=True, \
					bg_img=bgImg, threshold=thr)
		else:
			nlp.plot_epi(data, display_mode=slAxis, figure=fig, draw_cross=False, vmin=0, vmax=rangeVal, cmap='gnuplot2', cut_coords=zSlices[t0:tf], axes=(0,rowPlace[i],1,1/nRow), colorbar=True)
		t0 = t0 + step
		tf = tf + step

def vcorrcoef(X,y):
	# Vectorized, faster why to compute correlation for large matrices
	Xm = np.reshape(np.mean(X,axis=1),(X.shape[0],1))
	ym = np.mean(y)
	r_num = np.sum((X-Xm)*(y-ym),axis=1)
	r_den = np.sqrt(np.sum((X-Xm)**2,axis=1)*np.sum((y-ym)**2))
	r = r_num/r_den
	return r

# ++++++++++++++++++++++++++++++++++++++++++++++
#
# END OF function definitions
#
# ++++++++++++++++++++++++++++++++++++++++++++++


parser = argparse.ArgumentParser(description='Save QA check Plots')

# Required options   
# TODO Adjust grouping of arguments                 
reqoptions = parser.add_argument_group('Required arguments')

reqoptions.add_argument('-o', '-out', dest="outDir", required=True, help='Directory where images are to be saved' )
reqoptions.add_argument('-a', '-in', dest="inDir", required=True, help='Dir where EPI + masks are stored [MNI SPACE]' )

reqoptions.add_argument('-f', '-fname', dest="fname", required=True,  help='EPI image name ' )
reqoptions.add_argument('-x', '-outf', dest="outFname", required=True,  help='Name of the output plot' )

reqoptions.add_argument('-b', '-bg', dest="bg", required=True, help='Background Image for plots (e.g. MNI or mean func file path)' )

reqoptions.add_argument('-c', '-plane', dest="plane", required=False, default='axial', help='Plane to slice [Axial (Z)], Sagital (x), Coronal (y)' )

reqoptions.add_argument('-t', '-type', dest="type", required=False, default='std', help='Type of plot [STD], GlobalCorr, MotionCorr' )

reqoptions.add_argument('-e', '-mpe', dest="mpe", required=False, default='motion_estimate.par', help='Motion Estimate File' )
reqoptions.add_argument('-p', '-prog', dest="prog", required=False, default='AFNI',  help='Software which generated the motion estimate: [AFNI], FSL or SPM' )

reqoptions.add_argument('-r', '-range', dest="range", required=False, default='80%', help='Range Plot' )
reqoptions.add_argument('-l', '-thr', dest="thr", required=False, default=0, help='(ABsolute) Threshold for plotting [0.3]' )

reqoptions.add_argument('-d', '-dpi', dest="dpi", required=False, default=120, help='Saved figure DPI' )
reqoptions.add_argument('-s', '-smooth', dest="smooth", required=False, default=0, help='Visualisation smoothness [FWHM, mm]' )

reqoptions.add_argument('-n', '-save_nii', dest="save_nii", required=False, default=0, help='Save the 3d volume' )


args = parser.parse_args()


outDir = args.outDir
inDir = args.inDir
outFile = outDir + '/' + args.outFname
rpFile = inDir + '/' + args.mpe


plane = args.plane.lower()
if plane == 'axial':
	plane ='z'

if plane == 'coronal':
	plane ='y'

if plane == 'sagital':
	plane ='x'



# PNG resolution of the saved file
figDpi=int(args.dpi)

saveNii = int(args.save_nii)==1

rangeVal=args.range
usePctl=False
if rangeVal[-1] == '%':
	usePctl=True
	rangeVal=float(rangeVal[0:-1])


smooth = float(args.smooth)

# Font size and weight for ALL plots
plt.rcParams.update({'font.size': 20, 'font.weight':'bold'} )


funcImgFile = inDir + '/' + args.fname 

thr = float(args.thr)

bgImgPath = args.bg
bgImg = nimg.load_img(args.bg)


plotType = args.type.lower()

# ++++++++++++++++++++++++++++++++++++++++++++++++
#
# END OF PARAMETER/SETUP
#
# ++++++++++++++++++++++++++++++++++++++++++++++++


print('\n\n ===============================================================================================\n\n')
print('Starting Temporal SD/Correlation QC')
print('\n\n ===============================================================================================\n\n')




# =============================================
#
#  3. Temporal Standard Deviation
#
# =============================================

if plotType == 'std':
	funcImg = nib.load(funcImgFile)
	data = np.array(funcImg.get_data())
	data = np.std(data, axis=3)


	X,Y,Z = data.shape
	if usePctl == True:
		rangeVal=np.percentile( np.reshape(data, (X*Y*Z,1)), rangeVal)

	data = nimg.new_img_like(funcImg, data)

	plot_brain_overlay(data, None, rangeVal, figDpi, slAxis=plane)
	
	plt.savefig(outFile)
	if saveNii == True:
		nib.save(data, outDir + '/tSD.nii')



# =============================================
#
#  4. Voxelwise global signal correlation
#
# =============================================

if plotType == 'globalcorr':
	funcImg = nib.load(funcImgFile)
	data = np.array(funcImg.get_data())

	globalSignal = data

	X,Y,Z,N = data.shape
	nVox = X*Y*Z

	data = np.reshape(data, (X*Y*Z, N))
	dataMask = np.std( data, axis=1 )
	np.where(np.abs(data)<thr, 0, data)


	globalSignal = np.mean( np.reshape(globalSignal, (X*Y*Z, N)), axis=0 )
	corrMap = np.zeros((nVox, 1))



	corrMap[dataMask>0,0] = vcorrcoef( data[dataMask>0, :], globalSignal )
		
	corrMap = np.reshape( corrMap, (X,Y,Z) )




	data = nimg.new_img_like(funcImg, corrMap)



	# For visualization purposes, slightly smooths image [5mm FWHM Gaussian kernel]
	#data = nimg.smooth_img(data, 5)
	if smooth > 0:
		data = nimg.smooth_img(data, smooth)
	plot_brain_overlay(data, bgImg, rangeVal, figDpi, slAxis=plane, thr=thr, is_stat=True)

	plt.savefig(outFile)
	if saveNii == True:
		nib.save(data, outDir + '/Global_Corr.nii')



if plotType == 'motioncorr':

	rp = mh.convert_motion(np.loadtxt(rpFile), prog=args.prog)
	funcImg = nib.load(funcImgFile)
	data = np.array(funcImg.get_data())



	X,Y,Z,N = data.shape
	nVox = X*Y*Z

	data = np.reshape(data, (X*Y*Z, N))
	dataMask = np.std( data, axis=1 )

	np.where(np.abs(data)<thr, 0, data)


	corrMapRp = np.zeros((nVox, 1))

	
	for p in range(6):
		re = vcorrcoef( data[dataMask>0.5, :], rp[:,p] )
		tmp = corrMapRp[dataMask>0.5,0]
		comp = [ x>y for (x,y) in zip(abs(re), abs(tmp))]
		tmp = np.where( comp, re, tmp)
		corrMapRp[dataMask>0.5,0] = tmp
		#corrMapRp[dataMask>0.5,0] = vcorrcoef( data[dataMask>0.5, :], rp[:,p] )
	

	corrMapRp = np.reshape( corrMapRp, (X,Y,Z) )


	data = nimg.new_img_like(funcImg, corrMapRp)
	if smooth > 0:
		data = nimg.smooth_img(data, smooth)
	plot_brain_overlay(data, bgImg, rangeVal, figDpi, slAxis=plane, thr=thr, is_stat=True)

	plt.savefig(outFile)

	if saveNii == True:
		nib.save(data, outDir + '/Motion_Corr.nii')

