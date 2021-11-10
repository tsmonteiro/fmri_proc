#!/usr/bin/env python


import os
import subprocess
import argparse
import sys

import nibabel as nib
from builtins import str

import matplotlib.pyplot as plt
import numpy as np

import nipype.algorithms.confounds as npalg

import nilearn.image as nimg
import nilearn.signal as sgn

import matplotlib.colors as mcolors
from nilearn import datasets
from nilearn.image import index_img


from PIL import Image, ImageDraw, ImageFont

def plot_brain(bgImg, overlayImg, ax, z, vmax, vmin, fig, axMat=None):
    if np.max(bgImg) > 1:
        bgImg = bgImg / np.nanmax(bgImg)

    bgImg[np.where( np.isnan(bgImg) )] = 0
    overlayImg[np.where( np.isnan(overlayImg) )] = 0
    #roy_big_bl
    colors1 = plt.cm.BuPu(np.linspace(0., 1, 256))
    colors2 = plt.cm.inferno(np.linspace(0, 1, 256))

    mymap2  = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors1)
    mymap1  = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors2[120:,:])

    ax.imshow(bgImg[:,:,z], cmap='gray',  aspect=None, interpolation=None, alpha=None,
               vmin=0.1, vmax=1.35)
    ax.set_title('{}'.format(z+1), color=(0.8,0.8,0.3), fontdict={'fontsize':14}, pad=0)
    plt.xticks([])
    plt.yticks([])
    if overlayImg is None:
        dummy=1+1
    else:
        sliceIm = overlayImg[:,:,z]


        sliceIm[sliceIm<vmin] = np.nan

        p1 = ax.imshow(sliceIm, cmap=mymap1,  aspect=None, interpolation=None, alpha=None, vmin=vmin, vmax=vmax)

        sliceIm = -overlayImg[:,:,z]

        sliceIm[sliceIm<vmin] = np.nan

        p2 = ax.imshow(sliceIm, cmap=mymap2,  aspect=None, interpolation=None, alpha=None, vmin=vmin, vmax=vmax)

        if not fig == None:

            nsx = int(np.floor(axMat.shape[0]/2))

            cb  = fig.colorbar(p1, ax=axMat[-1,(nsx):(nsx+nsx-1)], orientation='horizontal', fraction=0.2, ticks=[vmin,vmax], shrink=1 )
            cb2 = fig.colorbar(p2, ax=axMat[-1,0:(nsx-1)], orientation='horizontal', fraction=0.2, ticks=[vmin,vmax], shrink=1 )

            cb.ax.set_xticklabels(labels=('{:.01f}'.format(vmin), '{:.01f}'.format(vmax)), fontsize=16, weight='bold',color='w')
            cb2.ax.set_xticklabels(labels=('{:.01f}'.format(-vmax), '{:.01f}'.format(-vmin)), fontsize=16, weight='bold',color='w')

parser = argparse.ArgumentParser(description='Save QA check Plots')

# Required options
reqoptions = parser.add_argument_group('Required arguments')

reqoptions.add_argument('-o', '-out', dest="outDir", required=True, help='Directory where images are to be saved' )
reqoptions.add_argument('-i', '-in', dest="inDir", required=True, help='Dir where EPI + masks are stored [MNI SPACE]' )

reqoptions.add_argument('-a', '-im1', dest="im1", required=True, help='First Image' )
reqoptions.add_argument('-b', '-im2', dest="im2", required=True, help='Second Image' )

reqoptions.add_argument('-x', '-msg1', dest="msg1", required=False, default=None, help='Title for image 1' )
reqoptions.add_argument('-y', '-msg2', dest="msg2", required=False, default=None, help='Title for iamge 2' )

reqoptions.add_argument('-d', '-dpi', dest="dpi", required=False, default=120, help='Saved figure DPI' )

reqoptions.add_argument('-t', '-type', dest="type", required=False, default='std', help='Type of comparison: [std], mean' )

reqoptions.add_argument('-ord', dest="voxOrder", required=False, default='std', help='How to order voxels' )


args = parser.parse_args()


outDir = args.outDir
inDir = args.inDir

im1 = args.im1
im2 = args.im2
outType=args.type

msg1 = args.msg1
msg2 = args.msg2

voxOrder = args.voxOrder


fileDirs = str.split( __file__, '/')
numDirs = len(fileDirs)

fileDir=''
for d in range(numDirs-1):
	fileDir = os.path.join(fileDir, fileDirs[d])


fontFile = os.path.join('/', fileDir, 'fonts/DAYPBL__.ttf')




# PNG resolution of the saved file
figDpi=int(args.dpi)



# Font size and weight for ALL plots
plt.rcParams.update({'font.size': 20, 'font.weight':'bold'} )


im1Name = im1.split('.')
im2Name = im2.split('.')

outSD  = outDir + '/' + im1Name[0] + '_' + im2Name[0] + '_sd_diff.png'

outSD1  = outDir + '/' + im1Name[0] + '.png'
outSD2  = outDir + '/' + im2Name[0] + '.png'
outSDD  = outDir + '/' + im1Name[0] + '_to_' + im2Name[0] + '.gif'


outG1  = outDir + '/' + im1Name[0] + '_GP.png'
outG2  = outDir + '/' + im2Name[0] + '_GP.png'
outGDD  = outDir + '/' + im1Name[0] + '_to_' + im2Name[0] + '_GP.gif'


if msg1 == None:
	msg1 = im1Name[0]


if msg2 == None:
	msg2 = im2Name[0]

nii1 = nib.load(inDir +'/' + im1)
nii2 = nib.load(inDir +'/' + im2)

data1 = np.array(nii1.get_fdata())
data2 = np.array(nii2.get_fdata())


bgImg = np.mean(data1, axis=3)
bgImg = nimg.new_img_like(nii1,bgImg)

x,y,z = bgImg.shape




nSteps = 50
nRow = 7


origX = float(subprocess.check_output(['3dinfo', '-oi', inDir +'/' + im1]))
origY = float(subprocess.check_output(['3dinfo', '-oj', inDir +'/' + im1]))
origZ = float(subprocess.check_output(['3dinfo', '-ok', inDir +'/' + im1]))

nZ = float(subprocess.check_output(['3dinfo', '-nk', inDir +'/' + im1]))

minPoint = float(subprocess.check_output(['3dinfo', '-Iextent', inDir +'/' + im1]))
maxPoint = float(subprocess.check_output(['3dinfo', '-Sextent', inDir +'/' + im1]))

minX = float(subprocess.check_output(['3dinfo', '-Rextent', inDir +'/' + im1]))
maxX = float(subprocess.check_output(['3dinfo', '-Sextent', inDir +'/' + im1]))


minY = float(subprocess.check_output(['3dinfo', '-Aextent', inDir +'/' + im1]))
maxY = float(subprocess.check_output(['3dinfo', '-Pextent', inDir +'/' + im1]))



oz = (maxPoint-minPoint)/4.
oy = (maxX-minX)/2.
ox = (maxY-minY)/2.



minPoint -= oz
maxPoint -= oz

maxX -= ox
minX -= ox
maxY -= oy
minY -= oy


zSlices  = np.linspace(minPoint, maxPoint, nSteps)



rowPlace = np.linspace(0, 1, nRow + 1)
step= int( np.floor( nSteps/nRow ) )

if outType == 'std':
	sd2 = np.std(data2, axis=3)
	sd2 = nimg.new_img_like(nii2,sd2)
elif outType == 'mean':
	sd2 = np.mean(data2, axis=3)
	sd2 = nimg.new_img_like(nii2,sd2)

t0 = 0
tf = step

if outType == 'std':
	sd1 = np.nanstd(data1, axis=3)
	sd1 = nimg.new_img_like(nii1,sd1)
elif outType == 'mean':
	sd1 = np.nanmean(data1, axis=3)
	sd1 = nimg.new_img_like(nii1,sd1)

sdd = sd1.get_fdata() - sd2.get_fdata()
sdd = nimg.new_img_like(nii1,sdd)

xPos = minX * 0.35
yPos = maxY * 0.2
vmaxSd1 = np.percentile(sd1.get_fdata(), 99)
vmaxSd2 = np.percentile(sd2.get_fdata(), 99)

vminSd1 = np.percentile(sd1.get_fdata(), 90)
vminSd2 = np.percentile(sd2.get_fdata(), 90)

#np.max([np.max(sd1.get_fdata()),np.max(sd2.get_fdata())])
vmaxSd = np.nanmax([vmaxSd1, vmaxSd2])
vminSd = np.nanmax([vminSd1, vminSd2])



vmaxSdd = np.percentile(abs(sdd.get_fdata()), 98)
vminSdd = np.percentile(abs(sdd.get_fdata()), 90)

colors1 = plt.cm.BuPu_r(np.linspace(0., 1, 256))
colors2 = plt.cm.hot(np.linspace(0, 1, 256))

# combine them and build a new colormap
colors = np.vstack((colors1, colors2))
mymap  = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

nCols = int(np.ceil(np.sqrt( nZ )))
nRows = int(nCols)




fig, axMat = plt.subplots(nCols, nRows, figsize=(16,12), dpi=figDpi, facecolor='k', edgecolor='k', sharex=True, sharey=True)
bgImg3d    = bgImg.get_fdata()

currentZ   = 0
sddData = sdd.get_fdata()
for i in range(nCols):
    for j in range(nRows):

    	if currentZ < nZ and outType == 'std':
            if currentZ == (nZ-1):
                plot_brain(bgImg3d, sddData, axMat[i,j], currentZ, vmaxSdd, vminSdd, fig, axMat)
            else:
                plot_brain(bgImg3d, sddData, axMat[i,j], currentZ, vmaxSdd, vminSdd, None)
            currentZ += 1



plt.savefig(outSD, facecolor=fig.get_facecolor(), edgecolor='none', transparent=True)
plt.close('all')


fig, axMat = plt.subplots(nCols, nRows, figsize=(16,12), dpi=figDpi, facecolor='k', edgecolor='k', sharex=True, sharey=True)
currentZ   = 0

sd1 = sd1.get_fdata()
#TODO Add dots for telling what is what
for i in range(nCols):
    for j in range(nRows):
        if currentZ < nZ and outType == 'std':
            if currentZ == (nZ-1):
                plot_brain(bgImg3d, sd1, axMat[i,j], currentZ, vmaxSd, vminSd, fig, axMat)
            else:
                plot_brain(bgImg3d, sd1, axMat[i,j], currentZ, vmaxSd, vminSd, None)
            currentZ += 1

        if currentZ < nZ and outType == 'mean':
            plot_brain(sd1, None, axMat[i,j], currentZ, vmaxSd, vminSd, None)
            currentZ += 1

plt.savefig(outSD1, facecolor=fig.get_facecolor(), edgecolor='none', transparent=True)
plt.close('all')

fig, axMat = plt.subplots(nCols, nRows, figsize=(16,12), dpi=figDpi, facecolor='k', edgecolor='k', sharex=True, sharey=True)
currentZ   = 0

sd2 = sd2.get_fdata()

for i in range(nCols):
    for j in range(nRows):
        if currentZ < nZ and outType == 'std':
            if currentZ == (nZ-1):
                plot_brain(bgImg3d, sd2, axMat[i,j], currentZ, vmaxSd, vminSd, fig, axMat)
            else:
                plot_brain(bgImg3d, sd2, axMat[i,j], currentZ, vmaxSd, vminSd, None)
            currentZ += 1

        if currentZ < nZ and outType == 'mean':
            plot_brain(sd2, None, axMat[i,j], currentZ, vmaxSd, vminSd, None)
            currentZ += 1

plt.savefig(outSD2, facecolor=fig.get_facecolor(), edgecolor='none', transparent=True)
plt.close('all')

im1 = Image.open(outSD1)
im2 = Image.open(outSD2)
w, h = im1.size

fnt = ImageFont.truetype(fontFile, 25)


draw = ImageDraw.Draw(im1)
draw.text( ((w*0.85)/2,0), msg1, fill=(255,255,255), font=fnt)
del draw

draw = ImageDraw.Draw(im2)
draw.text( ((w*0.85)/2,0), msg2, fill=(255,255,255), font=fnt )
del draw

ims = []

ims.append(im1)
ims.append(im2)


ims[0].save(outSDD, format='GIF',
   append_images=ims[1:], save_all=True, duration=1000, loop=0)




X,Y,Z,N = data1.shape
nVox = X*Y*Z


cX = (X - 1) // 2
cY = (Y - 1) // 2
cZ = (Z - 1) // 2

cent = np.zeros((3,1))
cent[0] = cX
cent[1] = cY
cent[2] = cZ

cent = np.transpose(cent)

# Sort voxels by temporal standard deviation
if voxOrder == 'std':
    # Reshape to voxels x time
    data1 = np.reshape(data1, (X*Y*Z, N))
    data2 = np.reshape(data2, (X*Y*Z, N))

    sdOrd  = np.std( data1, axis=1 )

    sortIdx = np.argsort((-sdOrd))

    data1 = data1[sortIdx,:]
    data2 = data2[sortIdx,:]

# Sort voxels by slice
if voxOrder == 'z':
    # Reshape to voxels x time
    data1 = np.transpose( data1, (2,1,0,3) )
    data2 = np.transpose( data2, (2,1,0,3) )

    # Order F --> Last index changes the slowest
    data1 = np.reshape(data1, (X*Y*Z, N))
    data2 = np.reshape(data2, (X*Y*Z, N))


mu    = np.nanmean(data1,axis=1)
mu2   = np.nanmean(data2,axis=1)

thr   = np.percentile(mu, 75)

data1 = data1 - mu[:, np.newaxis]
data2 = data2 - mu2[:, np.newaxis]

muData = np.nanmean( np.nanmean(data1,axis=1) )
sdData = np.nanmean( np.nanstd(data1,axis=1) )
vmin   = muData - 1*sdData
vmax   = muData + 1*sdData

data1 = data1[abs(mu) > thr, :]
data2 = data2[abs(mu) > thr, :]

nVox,un = data1.shape

import matplotlib.colors as mcolors
#colors1 = plt.cm.bone_r(np.linspace(0., 1, 256))
#colors2 = plt.cm.inferno(np.linspace(0, 1, 256))

colors1 = plt.cm.BuPu_r(np.linspace(0., 1, 256))
colors2 = plt.cm.BuGn(np.linspace(0, 1, 256))

# combine them and build a new colormap
colors = np.vstack((colors1, colors2))
mymap  = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

fig = plt.figure(figsize=(16,12), dpi=figDpi, facecolor='w', edgecolor='k')

im = plt.imshow( data1, interpolation='bilinear', extent=[0, N, 0, nVox], aspect='auto',  cmap=mymap, vmin=vmin, vmax=vmax )
if voxOrder == 'z':
    plt.yticks([0, nVox], labels=('Z = 0', 'Z = Max'))
if voxOrder == 'std':
    plt.yticks([0, nVox], labels=('Low tSD', 'High tSD'))
fig.colorbar(im, orientation='vertical')

plt.savefig(outG1)
plt.close('all')

fig = plt.figure(figsize=(16,12), dpi=figDpi, facecolor='w', edgecolor='k')

im=plt.imshow( data2, interpolation='bilinear', extent=[0, N, 0, nVox], aspect='auto',  cmap=mymap, vmin=vmin, vmax=vmax )
if voxOrder == 'z':
    plt.yticks([0, nVox], labels=('Z = 0', 'Z = Max'))
if voxOrder == 'std':
    plt.yticks([0, nVox], labels=('Low tSD', 'High tSD'))
fig.colorbar(im, orientation='vertical')

plt.savefig(outG2)
plt.close('all')


im1 = Image.open(outG1)
im2 = Image.open(outG2)


w, h = im1.size

fnt = ImageFont.truetype(fontFile, 40)


draw = ImageDraw.Draw(im1)
draw.text( ((w*0.8)/2,10), msg1, fill=(0,0,0), font=fnt)
del draw

draw = ImageDraw.Draw(im2)
draw.text( ((w*0.8)/2,10), msg2, fill=(0,0,0), font=fnt )
del draw

ims = []

ims.append(im1)
ims.append(im2)


ims[0].save(outGDD, format='GIF',
   append_images=ims[1:], save_all=True, duration=1000, loop=0)

#os.remove(outG1)
#os.remove(outG2)
