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

#PATH to qclib
sys.path.append(config['PATHS']['qclib_path'])

import qclib.motion_handler as mh
import qclib as qc



parser = argparse.ArgumentParser(description='Save QA check Plots')

# Required options
reqoptions = parser.add_argument_group('Required arguments')

reqoptions.add_argument('-o', '-out', dest="outDir", required=True, help='Directory where images are to be saved' )
reqoptions.add_argument('-a', '-in', dest="inDir", required=True, help='Dir where EPI + masks are stored [MNI SPACE]' )

reqoptions.add_argument('-e', '-mpe', dest="mpe", required=True, default='motion_estimate.par', help='Motion Estimate File' )

reqoptions.add_argument('-x', '-outf', dest="outFname", required=False, default='greyplot.png', help='Name of the output plot' )
reqoptions.add_argument('-f', '-fname', dest="fname", required=False, default='proc_data_MNI2mm.nii', help='Name of functional EPI file in MNI space' )

reqoptions.add_argument('-c', '-csf_name', dest="csfname", required=False, default='csf_mask_grp.nii', help='Name of CSF Mask file in MNI space' )
reqoptions.add_argument('-w', '-wm_name', dest="wmname", required=False, default='wm_mask_grp.nii', help='Name of WM Mask file in MNI space' )
reqoptions.add_argument('-z', '-gm_name', dest="gmname", required=False, default='gm_mask_grp.nii', help='Name of GM Mask file in MNI space' )

reqoptions.add_argument('-r', '-range', dest="range", required=False, default=3, help='Range of greyplot' )
reqoptions.add_argument('-n', '-norm', dest="norm", required=False, default='psc', help='Data normalisation [psc], zscore, none' )


reqoptions.add_argument('-t', '-tr', dest="tr", required=True, help='TR' )

reqoptions.add_argument('-d', '-dpi', dest="dpi", required=False, default=120, help='Saved figure DPI' )
reqoptions.add_argument('-p', '-prog', dest="prog", required=False, default='AFNI',  help='Software which generated the motion estimate: [AFNI], FSL or SPM' )


args = parser.parse_args()


outDir = args.outDir
inDir = args.inDir

outFile = outDir + '/' + args.outFname
rpFile = inDir + '/' + args.mpe

normType = args.norm.lower()

tr=float(args.tr)
if tr > 10:
	tr = tr / 1000

# PNG resolution of the saved file
figDpi=int(args.dpi)

rangeVal=args.range
usePctl=False
if rangeVal[-1] == '%':
	usePctl=True
	rangeVal=float(rangeVal[0:-1])
else:
	rangeVal=float(args.range)
	if normType == 'psc':
		rangeVal=rangeVal/100.



# Font size and weight for ALL plots
plt.rcParams.update({'font.size': 20, 'font.weight':'bold'} )


funcImgFile = inDir + '/' + args.fname
csfFile = args.csfname
gmFile = args.gmname
wmFile = args.wmname

print('EPI functional data: ' + funcImgFile)
print('GM Mask: ' + gmFile)
print('WM Mask: ' + wmFile)
print('CSF Mask: ' + csfFile)

# ++++++++++++++++++++++++++++++++++++++++++++++++
#
# END OF PARAMETER/SETUP
#
# ++++++++++++++++++++++++++++++++++++++++++++++++


print('\n\n ===============================================================================================\n\n')
print('Starting QC : Greyplot')
print('\n\n ===============================================================================================\n\n')

rp = mh.convert_motion(np.loadtxt(rpFile), prog=args.prog)


funcImg = nib.load(funcImgFile)
data = np.array(funcImg.get_fdata())

X,Y,Z,N = data.shape
nVox = X*Y*Z
# Reshape to voxels x time
data = np.reshape(data, (X*Y*Z, N))

dvar =  np.sqrt( np.sum( np.power( np.diff( data, axis=1 ), 2 ), axis=0 ) )
dvar = 100 * dvar / np.median(dvar)
dvar = np.insert( dvar, 0, 100 )


fd = mh.framewise_displacement(rp)
totalShift = mh.total_translation(rp)


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
#TODO Change here to whatever the user defines as the FD threshold
ax1.plot([0, N], [0.5, 0.5], linestyle='dotted', color=[0,0,0])
ax1.axis((0, N,0, 2))



ax2 = ax1.twinx()
ax2.set_ylabel('Total Shift [mm]', color='blue')
ax2.tick_params(axis='y', labelcolor='blue')

ax2.plot(time, totalShift, color=[0,0,1])


ax2.axis((0, N,0, 2))


for axis in ['left', 'bottom', 'right']:
	ax1.spines[axis].set_linewidth(3)
	ax2.spines[axis].set_linewidth(3)

ax2.spines['right'].set_color((0,0,1))


for axis in ['top']:
	ax1.spines[axis].set_linewidth(0)
	ax2.spines[axis].set_linewidth(0)





plt.subplot(grid[1:, 0])



muData = np.mean(data, axis=1)
sdData = np.std(data, axis=1)
#for v in range(nVox):
#	data[v,:] = data[v,:] / muData[v]

if normType == 'psc':
	data = data / np.transpose(np.tile(muData, (N,1)))
	data = data - 1

if normType == 'zscore':
	data = (data - np.transpose(np.tile(muData, (N,1))))/np.transpose(np.tile(sdData, (N,1)))

if normType == 'none':
	data = (data - np.transpose(np.tile(muData, (N,1))))


csfImg = nib.load(csfFile)
csfMask = np.array(csfImg.get_fdata())
csfMask = np.reshape(csfMask, X*Y*Z)


wmImg   = nib.load(wmFile)
wmMask  = np.array(wmImg.get_fdata())
wmMask  = np.reshape(wmMask, X*Y*Z)

gmImg   = nib.load(gmFile)
gmMask  = np.array(gmImg.get_fdata())
gmMask  = np.reshape(gmMask, X*Y*Z)

nVoxGm  = sum((gmMask>0.6) & (muData>0))
nVoxWm  = sum((wmMask>0.6) & (muData>0))
nVoxCsf = sum((csfMask>0.6) & (muData>0))

nTotal  = nVoxGm + nVoxWm + nVoxCsf

gmplot  = np.reshape(data[(gmMask>0.6) & (muData>0),:], (nVoxGm,N))
wmplot  = np.reshape(data[(wmMask>0.3) & (muData>0),:], (nVoxWm,N))
csfplot = np.reshape(data[(csfMask>0.3) & (muData>0),:], (nVoxCsf,N))

sdOrd   = np.std( gmplot, axis=1 )
sortIdx = np.argsort((-sdOrd))
gmplot  = gmplot[sortIdx,:]

sdOrd   = np.std( wmplot, axis=1 )
sortIdx = np.argsort((-sdOrd))
wmplot  = wmplot[sortIdx,:]

sdOrd   = np.std( csfplot, axis=1 )
sortIdx = np.argsort((-sdOrd))
csfplot = gmplot[sortIdx,:]



allTissuePlot = np.concatenate((gmplot,wmplot,csfplot))

if usePctl == True:
	rangeVal=np.percentile(allTissuePlot, rangeVal)

import matplotlib.colors as mcolors
colors1 = plt.cm.bone_r(np.linspace(0., 1, 256))
colors2 = plt.cm.inferno(np.linspace(0, 1, 256))

# combine them and build a new colormap
colors = np.vstack((colors1, colors2))
mymap  = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

plt.imshow( allTissuePlot, interpolation='bilinear', extent=[0, N, 0, nVoxGm+nVoxWm+nVoxCsf], aspect='auto',  cmap=mymap, vmin=-rangeVal, vmax=rangeVal )
plt.plot([0, N],[nTotal-nVoxGm, nTotal-nVoxGm], color='k', linewidth=4)
plt.plot([0, N],[nTotal-(nVoxGm+nVoxWm), nTotal-(nVoxGm+nVoxWm)], color='w', linewidth=4)


labels = []
labels.append('CSF')
labels.append('WM')
labels.append('GM')
plt.yticks([nVoxCsf/2, nVoxWm/2 +nVoxCsf, nVoxGm/2 + (nVoxCsf+nVoxWm) ])
ax = plt.gca()



ax.set_yticklabels(labels)
#plt.colorbar()
plt.savefig(outFile)
