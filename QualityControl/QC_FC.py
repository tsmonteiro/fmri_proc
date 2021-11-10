#!/usr/bin/env python


import os
import argparse
import sys

import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable

import nibabel as nib
from builtins import str

import matplotlib.pyplot as plt
import numpy as np

import nipype.algorithms.confounds as npalg
import nilearn.plotting as nlp
import nilearn.image as nimg
import nilearn.signal as sgn

import nibabel as nib
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


def plot_lg_matrix(matrix, cmap='RdBu_r', ax=None, colorbar=True, ylabel=True, xlabel=True, vmin=-5, vmax=5, cbarlabel='Z Score'):

    if not ax == None:
        plt.sca(ax)

    colors1 = plt.cm.BuPu_r(np.linspace(0., 1, 256))
    colors2 = plt.cm.afmhot_r(np.linspace(0, 1, 256))

    # combine them and build a new colormap
    colors = np.vstack((colors1, colors2))
    mymap  = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    im = plt.imshow( matrix, vmin=vmin, vmax=vmax, cmap=mymap )


    netNames = ['L Visual', 'L Somatomotor', 'L Dorsal Attention', 'L Ventral Attention',
                'L Limbic', 'L Control', 'L DMN', 'L Frontoparietal',
                'R Visual', 'R Somatomotor', 'R Dorsal Attention', 'R Ventral Attention',
                'R Limbic', 'R Control', 'R DMN', 'R Frontoparietal',
                'L Non-Cortical Regions', 'R Non-Cortical Regions']

    netL1    = [0,23,58,84,107,119,147,193,199,222,257,283,311,323,356,389,399,426,452]
    prev     = 0
    xticks   = []
    for d in netL1:
        gray   = 0.5
        gAlpha = 0.8

        xticks.append(prev + (d-prev)/2)

        plt.plot([0,452],[d,d], color=(gray,gray,gray,gAlpha), lw=1)
        plt.plot([0,452],[prev,prev], color=(gray,gray,gray,gAlpha), lw=1)
        plt.plot([d,d], [0,452], color=(gray,gray,gray,gAlpha), lw=1)
        plt.plot([prev,prev],[0,452], color=(gray,gray,gray,gAlpha), lw=1)
        prev = d


    prev = 0

    for d in netL1:
        plt.plot([prev,d],[d,d], 'k', lw=3)
        plt.plot([prev,d],[prev,prev], 'k', lw=3)
        plt.plot([d,d],[prev,d], 'k', lw=3)
        plt.plot([prev,prev],[prev,d], 'k', lw=3)
        prev = d

    plt.plot([0,452],[199,199], 'k', lw=3)
    plt.plot([199,199],[0,452], 'k', lw=3)

    plt.plot([0,452],[399,399], 'k', lw=3)
    plt.plot([399,399],[0,452], 'k', lw=3)
    if xlabel == True:
        plt.xticks(xticks[1:], netNames, fontsize=14, rotation=90 )
    else:
        plt.xticks([] )


    if ylabel == True:
        plt.yticks(xticks[1:], netNames, fontsize=14 )
    else:
        plt.yticks([] )

    ax = plt.gca()

    ax.xaxis.set_ticks_position('bottom')


    for axis in ['left', 'bottom']:
        ax.spines[axis].set_linewidth(3)

    for axis in ['top', 'right', 'left', 'bottom']:
        ax.spines[axis].set_linewidth(0)

    plt.xlim([-1,453])
    plt.ylim([453,-1])

    if colorbar == True:
        #cbar = plt.colorbar()
        #cbar.ax.set_ylabel(cbarlabel)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.15)

        cbar = plt.colorbar(im, cax=cax)
        cbar.ax.set_ylabel(cbarlabel)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# END OF function definitions
#
# ++++++++++++++++++++++++++++++++++++++++++++++


parser = argparse.ArgumentParser(description='Save QA check Plots')

# Required options
reqoptions = parser.add_argument_group('Required arguments')

reqoptions.add_argument( '-out', dest="outDir", required=True, help='Directory where images are to be saved' )
reqoptions.add_argument( '-in', dest="inDir", required=True, help='Dir where EPI + masks are stored [MNI SPACE]' )
reqoptions.add_argument( '-outf', dest="outFname", required=False, default='fc_mat', help='Name of the output plot' )


reqoptions.add_argument( '-fname', dest="fname", required=False, default='proc_data_MNI2mm.nii', help='Name of functional EPI file in MNI space' )

reqoptions.add_argument( '-TR', dest="TR", required=False, default=1.0, help='Name of functional EPI file in MNI space' )


plotoptions = parser.add_argument_group('Plot arguments')
# Plot options
plotoptions.add_argument('-vmin', dest="vmin", required=False, default=-0.8, help='Highpass cutoff' )
plotoptions.add_argument('-vmax', dest="vmax", required=False, default=0.8, help='Low-pass cutoff' )
plotoptions.add_argument('-d', '-dpi', dest="dpi", required=False, default=120, help='Saved figure DPI' )


plotoptions.add_argument('-atlas', dest="atlas", required=True, help='Saved figure DPI' )
plotoptions.add_argument('-atlas_name', dest="atlasName", required=True, help='Saved figure DPI' )
plotoptions.add_argument('-dist', dest="distMat", required=False, default=None, help='Saved figure DPI' )
plotoptions.add_argument('-labels', dest="labels", required=False, default=None, help='Saved figure DPI' )


args = parser.parse_args()

labels = args.labels
outDir = args.outDir
inDir = args.inDir
#TODO Add labels
outFileFC  = outDir + '/' + args.outFname + '_' + args.atlasName

outFileFCMat  = outDir + '/' + args.outFname + '_' + args.atlasName + '.txt'
outFileFCMat2  = outDir + '/' + args.outFname + '_' + args.atlasName + '_node.txt'
outFileCorr = outDir + '/' + args.outFname + '_corr'
outFileHist = outDir + '/' + args.outFname + '_hist'
outFileHistD = outDir + '/' + args.outFname + '_hist_dist'

funcImgFile = inDir + '/' + args.fname

TR = float(args.TR)
# PNG resolution of the saved file
figDpi=int(args.dpi)

vmin = float(args.vmin)
vmax = float(args.vmax)

# Font size and weight for ALL plots
plt.rcParams.update({'font.size': 20, 'font.weight':'bold'} )

# low_pass=0.1, smoothing_fwhm=5
masker   = NiftiLabelsMasker(labels_img=args.atlas, standardize=False, resampling_target=None,
                             t_r=TR, high_pass=0.01, low_pass=0.1, smoothing_fwhm=5)
signal   = masker.fit_transform(funcImgFile)

nScans   = signal.shape[0]



corrMeas = ConnectivityMeasure(kind='correlation')

corrMat  = corrMeas.fit_transform([signal])[0]



n1,n2=corrMat.shape

for n in range(n1):
	corrMat[n,n] =  0


atlas=nib.load(args.atlas)
atlas=atlas.get_fdata()
x,y,z=atlas.shape
atlas[np.where(np.isnan(atlas))]=0
uNodes=np.unique(atlas.reshape((x*y*z)))
uNodes = uNodes[1:]
np.savetxt(outFileFCMat2, uNodes, fmt='%d')


print('Writing ' + outFileFCMat)
np.savetxt(outFileFCMat, corrMat, fmt='%.8f')

fig = plt.figure(figsize=(20,20), dpi=figDpi, facecolor='w', edgecolor='k')

plot_lg_matrix(corrMat, cmap='RdBu_r', ax=None, colorbar=True, ylabel=True, xlabel=True, vmin=vmin, vmax=vmax, cbarlabel='r')

#if labels == None:
#    nlp.plot_matrix(corrMat, vmin=vmin, vmax=vmax, grid=False, auto_fit=True, colorbar=True, figure=fig, cmap='jet')
#else:
#    nlp.plot_matrix(corrMat, vmin=vmin, vmax=vmax, grid=False, auto_fit=True, colorbar=True, labels=labels, figure=fig, cmap='jet')

plt.savefig(outFileFC)

np.savetxt(outFileFCMat, corrMat, fmt='%.5f'  )


plt.close('all')

fig = plt.figure(figsize=(20,20), dpi=figDpi, facecolor='w', edgecolor='k')

fcs = []
for i in range(corrMat.shape[0]):
    for j in range(corrMat.shape[1]):
        if j > i:
            fcs.append( corrMat[i,j] )
import seaborn as sns
sns.distplot(fcs, norm_hist=True, bins=100, color = 'darkblue',
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 4, 'color':'black'})

plt.plot([0,0],[0, 3], lw=3, color=(0,0,0))
#plt.hist(fcs, bins=np.linspace(-1,1,200))

ax = plt.gca()

for axis in ['left', 'bottom']:
    ax.spines[axis].set_linewidth(3)

for axis in ['top', 'right']:
    ax.spines[axis].set_linewidth(0)
plt.xlim([-1,1])
plt.savefig(outFileHist)
plt.close('all')

#%%


#import pandas as pd
#if not args.distMat == None:
#
#    dmat = np.loadtxt(args.distMat, delimiter=',')
#    dvec = []
#
#    for i in range(corrMat.shape[0]):
#        for j in range(corrMat.shape[1]):
#            if j > i:
#
#                dvec.append( dmat[i,j] )
#
#
#
#    df = pd.DataFrame(np.concatenate((np.array(dvec).reshape([-1,1]),
#                                      np.array(fcs).reshape([-1,1])), axis=1),
#                        columns=['Distance','FC'])
#
#    fig = plt.figure(figsize=(20,10), dpi=figDpi, facecolor='w', edgecolor='k')
#    #cmap = sns.cubehelix_palette(as_cmap=True, dark=0, light=1, reverse=True)
#    sns.jointplot(x='Distance', y='FC', data=df, kind='kde', n_levels=60, shade=True, cmap='magma')
#    #g.plot_joint(plt.scatter, c=(0.8,0.8,0.8,0.1), s=2, linewidth=1, marker="+")
#
#    plt.savefig(outFileHistD)
#    plt.close('all')
