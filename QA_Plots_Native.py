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

reqoptions.add_argument('-e', '-mpe', dest="plotMP", required=False, default=True, help='Plot motion estimate (rotation and shifts)' )
reqoptions.add_argument('-f', '-fd', dest="plotFD", required=False, default=True, help='Plot Framewise Displacement + Grayplot' )
reqoptions.add_argument('-s', '-sd', dest="plotSD", required=False, default=True, help='Plot Temporal Standard Deviation' )
reqoptions.add_argument('-g', '-gs', dest="plotGS", required=False, default=True, help='Plot Voxelwise Global Signal & Max Motion correlation ' )
reqoptions.add_argument('-i', '-ic', dest="plotIC", required=False, default=True, help='Plot ICA AROMA classification overview ' )
reqoptions.add_argument('-n', '-nm', dest="plotNM", required=False, default=True, help='Plot Normalisation check ' )
reqoptions.add_argument('-x', '-fc', dest="plotFC", required=False, default=True, help='Plot Connectivity Matrix [MSDL Atlas]' )


reqoptions.add_argument('-p', '-fname', dest="fname", required=False, default='proc_data_MNI2mm.nii', help='Name of functional EPI file in MNI space' )
reqoptions.add_argument('-c', '-csf_name', dest="csfname", required=False, default='csf_mask_grp.nii', help='Name of CSF Mask file in MNI space' )
reqoptions.add_argument('-w', '-wm_name', dest="wmname", required=False, default='wm_mask_grp.nii', help='Name of WM Mask file in MNI space' )
reqoptions.add_argument('-z', '-gm_name', dest="gmname", required=False, default='gm_mask_grp.nii', help='Name of GM Mask file in MNI space' )

reqoptions.add_argument('-k', '-im1', dest="image1", required=False, default='', help='Registration QA Image 1' )
reqoptions.add_argument('-l', '-im2', dest="image2", required=False, default='', help='Registration QA Image 2' )


#'/mnt/hgfs/ssd_tmp/RepImpact/Template/Group_T1_MNI_1_5mm.nii'
reqoptions.add_argument('-t', '-tr', dest="tr", required=True, help='TR' )
reqoptions.add_argument('-b', '-bg', dest="bg", required=True, help='Background Image for plots (e.g. MNI or mean func file path)' )

reqoptions.add_argument('-d', '-dpi', dest="dpi", required=False, default=120, help='Saved figure DPI' )


args = parser.parse_args()


outDir = args.outDir
inDir = args.inDir

normIm1 = args.image1
normIm2 = args.image2


TR=args.tr

# PNG resolution of the saved file
figDpi=args.dpi


# What checks should be performed?

print('\nThe following QA plots will be generated and save to [' + outDir + ']:')

if args.plotMP.lower() == 'true' or args.plotMP == '1':
	print('\tMotion Estimates')
	plotMP = 1 
else:
	plotMP = 0 

if args.plotFD.lower() == 'true' or args.plotFD == '1':
	print('\tFramewise Displacement + Gray plot')
	plotFD = 1 
else:
	plotFD = 0 

if args.plotSD.lower() == 'true' or args.plotSD == '1':
	print('\tTemporal Standard deviation')
	plotSD = 1 
else:
	plotSD = 0 


if args.plotGS.lower() == 'true' or args.plotGS == '1':
	print('\tVoxelwise correlation to global signal and motion')
	plotGS = 1 
else:
	plotGS = 0 

if args.plotIC.lower() == 'true' or args.plotIC == '1':
	print('\tICA AROMA overview')
	plotIC = 1 
else:
	plotIC = 0 


if args.plotNM.lower() == 'true' or args.plotNM == '1':
	print('\tNormalisation resutls')
	plotNM = 1 
else:
	plotNM = 0 


if args.plotFC.lower() == 'true' or args.plotFC == '1':
	print('\tFunctional Connectivity [MSDL, between 0.01 and 0.1 Hz]')
	plotFC = 1 
else:
	plotFC = 0 




# Font size and weight for ALL plots
plt.rcParams.update({'font.size': 20, 'font.weight':'bold'} )


print('\n\n')
print('!!!ATTENTION!!!')
print('Make sure that the files below are correct!')

funcImgFile = inDir + '/' + args.fname 
csfFile = inDir + '/' + args.csfname 
gmFile = inDir + '/' + args.gmname 
wmFile = inDir + '/' + args.wmname 

print('EPI functional data: ' + funcImgFile)
print('GM Mask: ' + gmFile)
print('WM Mask: ' + wmFile)
print('CSF Mask: ' + csfFile)

print('\n\n ===============================================================================================\n\n')
print('Starting QA')
print('\n\n ===============================================================================================\n\n')

# Standard output filenames. Can be changed to whatever names

outRPFile  = outDir + '/01_motion_estimates.png'
outFDFile  = outDir + '/02_framewise_displacement.png'
outA2MFile = outDir + '/03_' + normIm1 + '2' + normIm2 + '.png'
outF2AFile = outDir + '/04_func2anat.png'
outSDFile  = outDir + '/05_func_tSD.png'
outGSFile  = outDir + '/06_GlobalSignal_Corr.png'
outCRFile  = outDir + '/07_Motion_Corr.png'
outICFile  = outDir + '/08_Motion_IC_Overview.png'
outNICFile = outDir + '/09_Motion_IC_NoMot_Overview.png'
outFCFile  = outDir + '/10_Connectivity_Matrix.png'

mniTemplate = nimg.load_img(args.bg)


#TODO
# Plot normalisation check anatT12MNI
# FNC <-> Distance relationship
# FNC histogram plot

# ++++++++++++++++++++++++++++++++++++++++++++++++
#
# END OF PARAMETER/SETUP
#
# ++++++++++++++++++++++++++++++++++++++++++++++++







# =============================================
#
#  1. Motion Estimates
#
# =============================================
rpFile = inDir + '/motion_estimate.par'

rp = np.loadtxt(rpFile)

N = rp.shape
N = N[0]

if plotMP == 1:
	print('Plotting motion estimates')
	# Roll, pitch yaw, dS, dL, dP
	fig = plt.figure(figsize=(16,12), dpi=figDpi, facecolor='w', edgecolor='k')
	plt.subplot(211)

	time = []
	[time.append(x+1) for x in range(N)]

	plt.plot(time,  rp[:,0], label='Roll', color=[.6,.6,1])
	plt.plot(time,  rp[:,1], label='Pitch', color=[.0,.0,1])
	plt.plot(time,  rp[:,2], label='Yaw', color=[0,0,.5])

	# plt.xlabel('TIME [s]')
	plt.ylabel('Rotation [degs]')
	plt.title('ROTATION')


	plt.axis((0, N,-1,1))
	plt.yticks([-1, -.5, 0, 0.5, 1])

	plt.legend()

	ax = plt.gca()

	for axis in ['left', 'bottom']:
		ax.spines[axis].set_linewidth(3)

	for axis in ['top', 'right']:
		ax.spines[axis].set_linewidth(0)


	plt.subplot(212)

	plt.plot(time,  rp[:,3], label='dS', color=[1,.6,.6])
	plt.plot(time,  rp[:,4], label='dL', color=[1,.0, 0])
	plt.plot(time,  rp[:,5], label='dP', color=[.5,0,0])

	plt.xlabel('TIME [Volume]')
	plt.ylabel('Shift [mm]')
	plt.title('TRANSLATION')

	plt.axis((0, N,-1,1))

	ax = plt.gca()

	for axis in ['left', 'bottom']:
		ax.spines[axis].set_linewidth(3)

	for axis in ['top', 'right']:
		ax.spines[axis].set_linewidth(0)


	plt.legend()
	fig.tight_layout()


	plt.yticks([-1.5, -1, -.5, 0 ,0.5,1, 1.5])

	plt.savefig(outRPFile)




# =============================================
#
#  2. Framewise displacement + Greyplot
#
# =============================================

if plotFD == 1:
	print('PLotting framewise displacement')
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

	fd = fd1 + fd2


	# Save framewise displacement and DVARs numbers
	# Can be used later for group QA summary
	np.savetxt(outDir+'/framewise_displacement.txt', fd, delimiter=',', fmt='%.3f')
	np.savetxt(outDir+'/dvars.txt', dvar, delimiter=',', fmt='%.3f')



	fig = plt.figure(figsize=(16,12), dpi=figDpi, facecolor='w', edgecolor='k')
	grid = plt.GridSpec(3, 1)



	ax1 = plt.subplot(grid[0, 0])

	time = []
	[time.append(x+1) for x in range(N)]

	ax1.plot(time,  fd, color=[0,0,0])


	ax1.set_ylabel('FD [mm]', color='black')
	ax1.tick_params(axis='y', labelcolor='black')
	ax1.plot([0, N], [0.3, 0.3], linestyle='dotted', color=[0,0,0])
	ax1.axis((0, N,0, 2))



	ax2 = ax1.twinx()
	ax2.set_ylabel('DVAR [a.u.]', color='blue')
	ax2.tick_params(axis='y', labelcolor='blue')
	ax2.plot(time, dvar, color=[0,0,1])

	ax2.plot([0, N], [150, 150], linestyle='dotted', color=[0,0,1])
	ax2.axis((0, N,0, 220))


	for axis in ['left', 'bottom', 'right']:
		ax1.spines[axis].set_linewidth(3)
		ax2.spines[axis].set_linewidth(3)

	for axis in ['top']:
		ax1.spines[axis].set_linewidth(0)
		ax2.spines[axis].set_linewidth(0)





	plt.subplot(grid[1:, 0])



	muData = np.mean(data, axis=1)

	for v in range(nVox):
		if muData[v] != 0:
			data[v,:] = data[v,:] / muData[v]
		else:
			data[v,:] = 0


	csfImg = nib.load(csfFile)
	csfMask = np.array(csfImg.get_data())
	csfMask = np.reshape(csfMask, X*Y*Z)


	wmImg = nib.load(wmFile)
	wmMask = np.array(wmImg.get_data())
	wmMask = np.reshape(wmMask, X*Y*Z)

	gmImg = nib.load(gmFile)
	gmMask = np.array(gmImg.get_data())
	gmMask = np.reshape(gmMask, X*Y*Z)

	nVoxGm = sum((gmMask>0.5) & (muData>0))
	nVoxWm = sum((wmMask>0.5) & (muData>0))
	nVoxCsf = sum((csfMask>0.5) & (muData>0))

	nTotal = nVoxGm + nVoxWm + nVoxCsf

	gmplot = np.reshape(data[(gmMask>0.5) & (muData>0),:], (nVoxGm,N))
	wmplot = np.reshape(data[(wmMask>0.5) & (muData>0),:], (nVoxWm,N))
	csfplot = np.reshape(data[(csfMask>0.5) & (muData>0),:], (nVoxCsf,N))


	allTissuePlot = np.concatenate((gmplot,wmplot,csfplot))



	plt.imshow( allTissuePlot, interpolation='gaussian', extent=[0, N, 0, nVoxGm+nVoxWm+nVoxCsf], aspect='auto',  cmap='YlGnBu_r', vmin=0.9, vmax=1.1 )
	plt.plot([0, N],[nTotal-nVoxGm, nTotal-nVoxGm], color='w', linewidth=4)
	plt.plot([0, N],[nTotal-(nVoxGm+nVoxWm), nTotal-(nVoxGm+nVoxWm)], color='w', linewidth=4)

	plt.plot([0, N],[nTotal-nVoxGm, nTotal-nVoxGm], color='w', linewidth=4)
	plt.plot([0, N],[nTotal-(nVoxGm+nVoxWm), nTotal-(nVoxGm+nVoxWm)], color='w', linewidth=4)


	labels = []
	labels.append('CSF')
	labels.append('WM')
	labels.append('GM')
	plt.yticks([nVoxCsf/2, nVoxWm/2 +nVoxCsf, nVoxGm/2 + (nVoxCsf+nVoxWm) ])
	ax = plt.gca()



	ax.set_yticklabels(labels)


	#plt.colorbar()
	plt.savefig(outFDFile)



# =============================================
#
#  3. Temporal Standard Deviation
#
# =============================================

if plotSD == 1:
	print('Plotting tSD')
	funcImg = nib.load(funcImgFile)
	data = np.array(funcImg.get_data())
	data = np.std(data, axis=3)


	data = nimg.new_img_like(funcImg, data)

	#print( data.shape )
	fig = plt.figure(figsize=(16,12), dpi=figDpi, facecolor='w', edgecolor='k')

	zSlices  = np.linspace(-60, 80, 60)

	nRow = 7
	rowPlace = np.linspace(0, 1, nRow + 1)





	step= int( np.floor( 60/nRow ) )
	t0 = 0
	tf = step

	for i in range(nRow):
		nlp.plot_epi(data, display_mode='z', figure=fig, draw_cross=False, vmin=0, vmax=50, cmap='gnuplot2', cut_coords=zSlices[t0:tf], axes=(0,rowPlace[i],1,1/nRow), colorbar=True)
		t0 = t0 + step
		tf = tf + step


	plt.savefig(outSDFile)



# =============================================
#
#  4. Voxelwise correlation with GS and MP
#
# =============================================

if plotGS == 1:
	print('Plotting GS and MP correlation')
	funcImg = nib.load(funcImgFile)
	data = np.array(funcImg.get_data())

	globalSignal = data

	X,Y,Z,N = data.shape
	nVox = X*Y*Z

	data = np.reshape(data, (X*Y*Z, N))
	dataMask = np.std( data, axis=1 )



	globalSignal = np.mean( np.reshape(globalSignal, (X*Y*Z, N)), axis=0 )
	corrMap = np.zeros((nVox, 1))
	corrMapRp = np.zeros((nVox, 1))




	for v in range(nVox):
		if dataMask[v] > 0.5:
			r = np.corrcoef( data[v,:], globalSignal )
			rMot1 = np.corrcoef( np.transpose( data[v,:] ), rp[:,0] )
			rMot2 = np.corrcoef( np.transpose( data[v,:] ), rp[:,1] )
			rMot3 = np.corrcoef( np.transpose( data[v,:] ), rp[:,2] )
			rMot4 = np.corrcoef( np.transpose( data[v,:] ), rp[:,3] )
			rMot5 = np.corrcoef( np.transpose( data[v,:] ), rp[:,4] )
			rMot6 = np.corrcoef( np.transpose( data[v,:] ), rp[:,5] )

			rMot = [ rMot1[0,1], rMot2[0,1], rMot3[0,1], rMot4[0,1], rMot5[0,1], rMot6[0,1] ]
			rMotAbs = np.abs(rMot)
			maxIdx = np.argmax(rMotAbs)

			corrMap[v] = r[0,1]
			corrMapRp[v] = rMot[maxIdx]
		else:
			corrMap[v] = 0
			corrMapRp[v] = 0

		#if np.mod(v, 25000) == 0:
		#	print(str(v) + ' / ' + str(nVox) )
		
	corrMap = np.reshape( corrMap, (X,Y,Z) )
	corrMapRp = np.reshape( corrMapRp, (X,Y,Z) )


	for p in range(2):
		if p == 0:
			data = nimg.new_img_like(funcImg, corrMap)
		else:
			data = nimg.new_img_like(funcImg, corrMapRp)


		# For visualization purposes, slightly smooths image [5mm FWHM Gaussian kernel]
		data = nimg.smooth_img(data, 5)

		fig = plt.figure(figsize=(16,12), dpi=figDpi, facecolor='w', edgecolor='k')

		zSlices  = np.linspace(-60, 80, 60)

		nRow = 7
		rowPlace = np.linspace(0, 1, nRow + 1)



		step= int( np.floor( 60/nRow ) )
		t0 = 0
		tf = step

		for i in range(nRow):
			nlp.plot_stat_map(data, display_mode='z', figure=fig, draw_cross=False,  vmax=.8, cmap='cold_hot', \
					cut_coords=zSlices[t0:tf], axes=(0,rowPlace[i],1,1/nRow), colorbar=True, \
					bg_img=mniTemplate, threshold=0.15)
			t0 = t0 + step
			tf = tf + step

		if p == 0:
			plt.savefig(outGSFile)
		else:
			plt.savefig(outCRFile)


# =============================================
#
#  5. Plot Independent Components
#
# =============================================

if plotIC == 1:
	baseIcDir = inDir + '/ICA_AROMA'
	mniIcFile = baseIcDir + '/melodic_IC_thr_MNI2mm.nii.gz'
	motionIcs = baseIcDir + '/classified_motion_ICs.txt'
	



	motionIcs = np.loadtxt(motionIcs, delimiter=',')


	funcImg = nib.load(mniIcFile)
	data = np.array(funcImg.get_data())

	X,Y,Z,N=data.shape


	if len(motionIcs) > 40:
		fig = plt.figure(figsize=(16,60), dpi=figDpi, facecolor='w', edgecolor='k')
	else:
		fig = plt.figure(figsize=(16,36), dpi=figDpi, facecolor='w', edgecolor='k')


	nRow = len(motionIcs)
	rowPlace = np.linspace(0, 1, nRow + 1)

	grid = plt.GridSpec(nRow, 3, figure=fig, hspace=0)


	gridI = -1

	for ic in range(len(motionIcs)):


		icMap = data[:,:,:, int(motionIcs[ic]-1)]
		maxCoord = np.argmax(icMap)
		icMap = nimg.new_img_like(funcImg, icMap)
	
	
		#maxCoord
		mx,my,mz = np.unravel_index(maxCoord, icMap.shape)

		mx,my,mz= nimg.coord_transform(mx,my,mz, icMap.affine)

		icSignal = baseIcDir + '/melodic.ica/report/t' + str(int(motionIcs[ic])) + '.txt'
		icSignal = np.loadtxt(icSignal)

		icPwr = baseIcDir + '/melodic.ica/report/f' + str(int(motionIcs[ic])) + '.txt'
		icPwr = np.loadtxt(icPwr)

		nlp.plot_stat_map(icMap, display_mode='ortho', figure=fig, draw_cross=False,  cut_coords=(mx,my,mz), axes=(0,rowPlace[ic],.33,1/nRow), colorbar=False, \
					threshold=0.5, bg_img=mniTemplate, annotate=False)

		ax1 = plt.subplot(grid[gridI, 1])




		minSig   = np.min(icSignal)
		maxSig   = np.max(icSignal)

		icStr = 'IC ' + str(int(motionIcs[ic]))

		plt.text(-10, (minSig+maxSig)/2, icStr, fontsize=12, rotation='vertical', horizontalalignment='left', va='center', color=[1, 0, 0])

		tenSecI  = 0
		hundSecI = 0
		for tp in range(len(icSignal)):
			tenSec = tenSecI*TR
			hundSec = hundSecI*TR

			if tenSec >= 10:
				tenSecI = 0
				plt.plot([tp,tp], [minSig, maxSig], color=[.66, .66, .66])
			else:
				tenSecI += 1


			if hundSec >= 100:
				hundSecI = 0
				plt.plot([tp,tp], [minSig, maxSig], color=[0, 0, 1])
			else:
				hundSecI += 1


		ax1.plot(icSignal, color='k', linewidth=1)
		ax1.axes.get_xaxis().set_visible(False)
		ax1.axes.get_yaxis().set_visible(False)

		ax1 = plt.subplot(grid[gridI, 2])

		maxPwr = np.max(icPwr)
		NP = len(icPwr)

		NPi = np.linspace(0, (1/TR)/2, NP)
		
		pl001 = 0
		pl01  = 0

		li = 0
		for l in range(NP):
			if NPi[l] >= 0.01 and pl001 == 0:
				ax1.plot([l, l], [0, maxPwr], color=[0, 0, 1])
				pl001 = 1
			
			if NPi[l] >= 0.1 and pl01 == 0:
				ax1.plot([l, l], [0, maxPwr], color=[0, 0, 1])
				pl01 = 1


		ax1.plot(icPwr, color='k')
	

		ax1.axes.get_xaxis().set_visible(False)
		ax1.axes.get_yaxis().set_visible(False)

		#if ic == (nRow-1):
		#	plt.title('IC Power')

		gridI -= 1

	plt.subplots_adjust(0, 0, 1, 1, 0, 0)

	plt.savefig(outICFile)
	


	nRow = N-len(motionIcs)
	rowPlace = np.linspace(0, 1, nRow + 1)


	if nRow > 40:
		fig = plt.figure(figsize=(12,60), dpi=figDpi, facecolor='w', edgecolor='k')
	else:
		fig = plt.figure(figsize=(12,36), dpi=figDpi, facecolor='w', edgecolor='k')

	grid = plt.GridSpec(nRow, 3, figure=fig, hspace=0)




	gridI = -1
	idx=0
	for ic in range(N):

		if not np.any( (ic+1) == motionIcs.astype(int) ):
			#print(int(motionIcs[ic]-1) )
			icMap = data[:,:,:, ic]
			maxCoord = np.argmax(icMap)
			icMap = nimg.new_img_like(funcImg, icMap)
	
	
			#maxCoord
			mx,my,mz = np.unravel_index(maxCoord, icMap.shape)

			mx,my,mz= nimg.coord_transform(mx,my,mz, icMap.affine)

			icSignal = baseIcDir + '/melodic.ica/report/t' + str(int(ic+1)) + '.txt'
			icSignal = np.loadtxt(icSignal)

			icPwr = baseIcDir + '/melodic.ica/report/f' + str(int(ic+1)) + '.txt'
			icPwr = np.loadtxt(icPwr)

			nlp.plot_stat_map(icMap, display_mode='ortho', figure=fig, draw_cross=False,  cut_coords=(mx,my,mz), axes=(0,rowPlace[idx],.33,1/nRow), colorbar=False, \
						threshold=0.5, bg_img=mniTemplate, annotate=False)

			idx+=1
			ax1 = plt.subplot(grid[gridI, 1])




			minSig   = np.min(icSignal)
			maxSig   = np.max(icSignal)

			icStr = 'IC ' + str(int(ic))

			plt.text(-10, (minSig+maxSig)/2, icStr, fontsize=12, rotation='vertical', horizontalalignment='left', va='center', color=[1, 0, 0])

			tenSecI  = 0
			hundSecI = 0
			for tp in range(len(icSignal)):
				tenSec = tenSecI*TR
				hundSec = hundSecI*TR

				if tenSec >= 10:
					tenSecI = 0
					plt.plot([tp,tp], [minSig, maxSig], color=[.8, .8, .8])
				else:
					tenSecI += 1


				if hundSec >= 100:
					hundSecI = 0
					plt.plot([tp,tp], [minSig, maxSig], color=[0, 0, 1])
				else:
					hundSecI += 1


			ax1.plot(icSignal, color='k', linewidth=1)
			ax1.axes.get_xaxis().set_visible(False)
			ax1.axes.get_yaxis().set_visible(False)

			ax1 = plt.subplot(grid[gridI, 2])

			maxPwr = np.max(icPwr)
			NP = len(icPwr)

			NPi = np.linspace(0, (1/TR)/2, NP)
		
			pl001 = 0
			pl01  = 0

			li = 0
			for l in range(NP):
				if NPi[l] >= 0.01 and pl001 == 0:
					ax1.plot([l, l], [0, maxPwr], color=[0, 0, 1])
					pl001 = 1
			
				if NPi[l] >= 0.1 and pl01 == 0:
					ax1.plot([l, l], [0, maxPwr], color=[0, 0, 1])
					pl01 = 1


			ax1.plot(icPwr, color='k')
	

			ax1.axes.get_xaxis().set_visible(False)
			ax1.axes.get_yaxis().set_visible(False)

			#if ic == (nRow-1):
			#	plt.title('IC Power')

			gridI -= 1

	plt.subplots_adjust(0, 0, 1, 1, 0, 0)

	plt.savefig(outNICFile)



# =============================================
#
#  6. Plot Connectivity Matrix
#
# =============================================

if plotFC == 1:
	atlas = datasets.fetch_atlas_msdl()
	atlasFname = atlas['maps']
	labels = atlas['labels']

	masker = NiftiMapsMasker(maps_img=atlasFname, standardize=True)
	signal = masker.fit_transform(funcImgFile)

	
	signal = sgn.clean(signal, low_pass=0.1, high_pass=0.01, t_r=1.6)

	corrMeas = ConnectivityMeasure(kind='correlation')
	corrMat = corrMeas.fit_transform([signal])[0]



	fig = plt.figure(figsize=(20,20), dpi=figDpi, facecolor='w', edgecolor='k')

	nlp.plot_matrix(corrMat, vmin=-.8, vmax=.8, labels=labels, grid=True, auto_fit=True, colorbar=True, figure=fig, cmap='cold_white_hot')

	plt.savefig(outFCFile)

# =============================================
#
#  7. Normalisation Results
#
# =============================================


if plotNM == 1:

	
	afunc2anat = nimg.load_img( inDir + '/' + normIm1 )
	data2 =  np.array(afunc2anat.get_data())
	anatFile = inDir + '/' + normIm2

	niimg = nimg.new_img_like(afunc2anat, data2)

	vmax = np.amax(data2)

	vmin = vmax * 0
	vmax = vmax * 0.66

	fig = plt.figure(figsize=(30,6), dpi=figDpi, facecolor='w', edgecolor='k')

	cuts = np.linspace(-34, 58, 20)

	display1 = nlp.plot_epi(niimg, cmap='gray', figure=fig, cut_coords=cuts, display_mode='z', draw_cross=False, vmax=vmax, vmin=vmin, \
			axes=(0,0,1,.33))
	display1.add_edges(anatFile)

	cuts = np.linspace(-40, 40, 20)

	display2 = nlp.plot_epi(niimg, cmap='gray', figure=fig, cut_coords=cuts, display_mode='x', draw_cross=False, vmax=vmax, vmin=vmin, \
			axes=(0,.33,1,.33))
	display2.add_edges(anatFile)


	cuts = np.linspace(-80, 40, 20)

	display3 = nlp.plot_epi(niimg, cmap='gray', figure=fig, cut_coords=cuts, display_mode='y', draw_cross=False, vmax=vmax, vmin=vmin, \
			axes=(0,.66,1,.33))
	display3.add_edges(anatFile)

	plt.savefig(outA2MFile)




