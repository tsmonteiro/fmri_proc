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
from nilearn.input_data import NiftiLabelsMasker
from nilearn.input_data import NiftiMapsMasker
from nilearn.connectome import ConnectivityMeasure


from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
from sklearn.linear_model import LinearRegression

from scipy.interpolate import griddata
import numpy.ma as ma


def least_sq(X,Y):
	# Building the model
	X_mean = np.mean(X)
	Y_mean = np.mean(Y)

	num = 0
	den = 0
	for i in range(len(X)):
	    num += (X[i] - X_mean)*(Y[i] - Y_mean)
	    den += (X[i] - X_mean)**2
	m = num / den
	c = Y_mean - m*X_mean
	return  m*X + c,m,c





def load_nuis(filepath):
	nuis = np.loadtxt(filepath)

	print(filepath)
	print(nuis.shape)
	print(len(nuis.shape))
	if len(nuis.shape)>1:
		d1,d2= nuis.shape
		if d2 > d1:
			nuis = np.transpose(nuis)
	else:
		tmpNuis = np.zeros((len(nuis),1))
		tmpNuis[:,0]=nuis
		nuis=tmpNuis



	return nuis


def calc_dist_matrix(coords):
	#coords=atlas.region_coords
	nNodes = len(coords)

	distMat = np.zeros((nNodes,nNodes))

	for n1 in range(nNodes):
		pt1 = np.asarray(coords[n1])
		for n2 in range(nNodes):
			pt2 = np.asarray(coords[n2])
			distMat[n1,n2] = np.linalg.norm(pt1-pt2)

	return distMat



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
reqoptions.add_argument('-x', '-outf', dest="outFname", required=False, default='fc_mat', help='Name of the output plot' )


reqoptions.add_argument('-f', '-fname', dest="fname", required=False, default='proc_data_MNI2mm.nii', help='Name of functional EPI file in MNI space' )

reqoptions.add_argument('-g', '-cov1', dest="cov1", required=False, default=None, help='File containing nuisance covariates' )
reqoptions.add_argument('-y', '-cov2', dest="cov2", required=False, default=None, help='File containing nuisance covariates' )
reqoptions.add_argument('-j', '-cov3', dest="cov3", required=False, default=None, help='File containing nuisance covariates' )
reqoptions.add_argument('-k', '-cov4', dest="cov4", required=False, default=None, help='File containing nuisance covariates' )

reqoptions.add_argument('-t', '-tr', dest="tr", required=True,  help='TR' )

reqoptions.add_argument('-u', '-high', dest="high", required=False, default=None, help='Highpass cutoff' )
reqoptions.add_argument('-l', '-low', dest="low", required=False, default=None, help='Low-pass cutoff' )

reqoptions.add_argument('-d', '-dpi', dest="dpi", required=False, default=120, help='Saved figure DPI' )


args = parser.parse_args()


outDir = args.outDir
inDir = args.inDir

outFile = outDir + '/' + args.outFname
outFileMsdl = outDir + '/' + args.outFname + '_msdl'
outFileCorr = outDir + '/' + args.outFname + '_corr'
outFileHist = outDir + '/' + args.outFname + '_hist' 

funcImgFile = inDir + '/' + args.fname 

tr = float(args.tr)

if tr > 10:
	tr = tr / 1000

if not args.low == None:
	lpCutoff=float(args.low)
else:
	lpCutoff = None

if not args.high == None:
	hpCutoff=float(args.high)
else:
	hpCutoff = None

# PNG resolution of the saved file
figDpi=int(args.dpi)


nuis = None

if args.cov1 != None:
	nuis = load_nuis(inDir + '/' + args.cov1)

if args.cov2 != None:
	nuis = np.concatenate((nuis, load_nuis(inDir + '/' + args.cov2)), axis=1)

if args.cov3 != None:
	nuis = np.concatenate((nuis, load_nuis(inDir + '/' + args.cov3)), axis=1)

if args.cov4 != None:
	nuis = np.concatenate((nuis, load_nuis(inDir + '/' + args.cov4)), axis=1)

# Font size and weight for ALL plots
plt.rcParams.update({'font.size': 20, 'font.weight':'bold'} )




#atlas = datasets.fetch_atlas_msdl()
atlas = datasets.fetch_atlas_schaefer_2018()
atlasFname = atlas.maps #atlas['maps']
labels = atlas.labels #atlas['labels']


atlasMsdl = datasets.fetch_atlas_msdl()
msdlFname = atlasMsdl['maps']
msdlLabels = atlasMsdl['labels']

masker = NiftiLabelsMasker(labels_img=atlasFname, standardize=False)
maskerMsdl = NiftiMapsMasker(maps_img=msdlFname, standardize=False)
signal = masker.fit_transform(funcImgFile)
signalMsdl = maskerMsdl.fit_transform(funcImgFile)


signal = sgn.clean(signal, low_pass=lpCutoff, high_pass=hpCutoff, t_r=tr, confounds=nuis)
signalMsdl = sgn.clean(signalMsdl, low_pass=lpCutoff, high_pass=hpCutoff, t_r=tr, confounds=nuis)

corrMeas = ConnectivityMeasure(kind='correlation', discard_diagonal=True)
corrMat = corrMeas.fit_transform([signal])[0]
corrMatMsdl = corrMeas.fit_transform([signalMsdl])[0]

n1,n2=corrMat.shape

for n in range(n1):
	corrMat[n,n] =  0

fig = plt.figure(figsize=(20,20), dpi=figDpi, facecolor='w', edgecolor='k')

netLen = np.zeros((14,1))

for l in atlas.labels:
	if 'Vis' in str(l):
		i = 0 

	if 'SomMot' in str(l):
		i = 1

	if 'DorsAttn' in str(l):
		i = 2

	if 'SalVentAttn' in str(l):
		i = 3


	if 'Limbic' in str(l):
		i = 4

	if 'Cont' in str(l):
		i = 5

	if 'Default' in str(l):
		i = 6

	if 'RH' in str(l):
		i += 7
	
	netLen[i] += 1

prev = 0
xticks = []

# Notation to create a list with 400 empty string elements
labs = [''] * 400 
netNames = list(['Left Visual', 'Left Somatomotor', 'Left Dorsal Attention', 'Left Sal/Vent Attention', 'Left Limbic', 'Left Control', 'Left DMN', 
'Right Visual', 'Right Somatomotor', 'Right Dorsal Attention', 'Right Sal/Vent Attention', 'Right Limbic', 'Right Control', 'Right DMN'])

for n in range(len(netLen)):
	xtick = netLen[n]/2. + prev    
	xticks.append(xtick)

	prev += netLen[n]
	labs[int(xtick)] = netNames[n]



nlp.plot_matrix(corrMat, vmin=-.8, vmax=.8, labels=labs, grid=False, auto_fit=True, colorbar=True, figure=fig, cmap='cold_white_hot')

prev = 0
for n in range(len(netLen)):
	x0 = prev
	x = netLen[n] + prev

	plt.plot([x0,x0], [x0,x], color='black', linewidth=3)	
	plt.plot([x,x], [x0,x], color='black', linewidth=3)
	plt.plot([x0,x], [x0,x0], color='black', linewidth=3)
	plt.plot([x0,x], [x,x], color='black', linewidth=3)

	prev += netLen[n]

plt.savefig(outFile)



# Saving MSDL
n1,n2=corrMatMsdl.shape

for n in range(n1):
	corrMatMsdl[n,n] =  0

fig = plt.figure(figsize=(20,20), dpi=figDpi, facecolor='w', edgecolor='k')
nlp.plot_matrix(corrMatMsdl, vmin=-.8, vmax=.8, labels=msdlLabels, grid=True, auto_fit=True, colorbar=True, figure=fig, cmap='cold_white_hot')
plt.savefig(outFileMsdl)


np.savetxt(outDir + '/' + args.outFname + '_correlation_matrix.txt', corrMat, fmt='%.5f'  )
np.savetxt(outDir + '/' + args.outFname + '_correlation_matrix_msdl.txt', corrMatMsdl, fmt='%.5f'  )


nNodes = len(atlasMsdl.region_coords)
distMat = calc_dist_matrix(atlasMsdl.region_coords)

distVec = np.reshape(distMat, (nNodes*nNodes, 1))
corrVec = np.reshape(corrMatMsdl, (nNodes*nNodes, 1))

corrVec = corrVec[distVec>0]
distVec = distVec[distVec>0]


fig = plt.figure(figsize=(16,12), dpi=figDpi, facecolor='w', edgecolor='k')

ax = fig.add_subplot(111)

x = np.array(distVec)
y = np.array(corrVec)



y_pred,m,c = least_sq(x,y)

r = np.corrcoef(x,y)
plt.text(0.05, -0.030, 'r = %.3f' % r[0,1], bbox=dict(facecolor='white', alpha=0.5) )

plt.plot(x, y_pred, '-', linewidth=3, color=(0.2,0.2,0.7,0.8))



ax.axvline(c='grey', lw=1)
ax.axhline(c='grey', lw=1)


ax.scatter(x, y, s=10, c='k' )
ax.scatter(np.mean(x), np.mean(y), c='b', s=250, marker='x')

plt.ylabel('FC [correlation]')
plt.xlabel('ROI-to-ROI Distance [mm]')

plt.savefig(outFileCorr)



# Plot histogram
fig = plt.figure(figsize=(16,12), dpi=figDpi, facecolor='w', edgecolor='k')
#plt.hist(corrVec, bins=30, density=True, histtype='bar')

# the histogram of the data
n, bins, patches = plt.hist(corrVec, bins=50, density=True, histtype='bar', linestyle='-',edgecolor=(.2,.2,.5),facecolor=(.6,.6,1),range=(-1,1))

mu = np.mean(corrVec)
sigma = np.std(corrVec)

# add a 'best fit' line
y = ((1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-0.5 * (1 / sigma * (bins - mu))**2))
bin1 = mu-sigma
bin2 = mu-2*sigma
bin3 = mu-3*sigma
y1 = ((1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-0.5 * (1 / sigma * (bin1 - mu))**2))
y2 = ((1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-0.5 * (1 / sigma * (bin2 - mu))**2))
y3= ((1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-0.5 * (1 / sigma * (bin3 - mu))**2))
plt.plot(bins, y, linestyle='-', color=(0,0,0.4),linewidth=3)

plt.plot([mu, mu],[0, np.amax(y)], color=(0,0,0.4), linestyle='--', linewidth=3 )
plt.plot([mu-sigma, mu-sigma],[0, y1], color=(0,0,0.4), linestyle='--', linewidth=3 )
plt.plot([mu+sigma, mu+sigma],[0, y1], color=(0,0,0.4), linestyle='--', linewidth=3 )

plt.plot([mu-2*sigma, mu-2*sigma],[0, y2], color=(0,0,0.4), linestyle='--', linewidth=3 )
plt.plot([mu+2*sigma, mu+2*sigma],[0, y2], color=(0,0,0.4), linestyle='--', linewidth=3 )

plt.plot([mu-3*sigma, mu-3*sigma],[0, y3], color=(0,0,0.4), linestyle='--', linewidth=3 )
plt.plot([mu+3*sigma, mu+3*sigma],[0, y3], color=(0,0,0.4), linestyle='--', linewidth=3 )



plt.savefig(outFileHist)




np.savetxt(outDir + '/' + args.outFname + '_correlation_vector.txt', corrVec, fmt='%.3f'  )
np.savetxt(outDir + '/' + args.outFname + '_distance_vector.txt', distVec, fmt='%.3f'  )
np.savetxt(outDir + '/' + args.outFname + '_correlation_matrix.txt', corrMat, fmt='%.3f'  )
np.savetxt(outDir + '/' + args.outFname + '_correlation_matrix_msdl.txt', corrMatMsdl, fmt='%.3f'  )
