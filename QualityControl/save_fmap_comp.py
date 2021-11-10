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
import nilearn.plotting as nlp
import nilearn.image as nimg
import nilearn.signal as sgn


from nilearn import datasets
from nilearn.image import index_img

import nilearn.masking as nbmsk


from PIL import Image, ImageDraw, ImageFont

parser = argparse.ArgumentParser(description='Save QA check Plots')

# Required options                    
reqoptions = parser.add_argument_group('Required arguments')

reqoptions.add_argument('-o', '-out', dest="outDir", required=True, help='Directory where images are to be saved' )
reqoptions.add_argument('-i', '-in', dest="inDir", required=True, help='Dir where EPI + masks are stored [MNI SPACE]' )

reqoptions.add_argument('-a', '-im1', dest="im1", required=True, help='First Image' )
reqoptions.add_argument('-b', '-im2', dest="im2", required=True, help='Second Image' )
reqoptions.add_argument('-c', '-bg_im', dest="bg_im", required=True, help='Anatomical Image' )

reqoptions.add_argument('-x', '-msg1', dest="msg1", required=False, default=None, help='Title for image 1' )
reqoptions.add_argument('-y', '-msg2', dest="msg2", required=False, default=None, help='Title for image 2' )

reqoptions.add_argument('-d', '-dpi', dest="dpi", required=False, default=120, help='Saved figure DPI' )



args = parser.parse_args()


outDir = args.outDir
inDir = args.inDir

im1 = args.im1
im2 = args.im2
bgImg = args.bg_im



msg1 = args.msg1
msg2 = args.msg2


	


# PNG resolution of the saved file
figDpi=int(args.dpi)



# Font size and weight for ALL plots
plt.rcParams.update({'font.size': 20, 'font.weight':'bold'} )


im1Name = im1.split('.')
im2Name = im2.split('.')

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

bgNii = nib.load(inDir +'/' + bgImg)


data1 = np.array(nii1.get_data())
data2 = np.array(nii2.get_data())
dataBg = np.array(bgNii.get_data())

#bgImg = np.mean(data1, axis=3)
bgImg = nimg.new_img_like(bgNii,dataBg)

x,y,z = bgImg.shape


fig = plt.figure(figsize=(18,14), dpi=figDpi, facecolor='w', edgecolor='k')

nSteps = 50
nRow = 7

origX = float(subprocess.check_output(['3dinfo', '-oi', inDir +'/' + im1]))
origY = float(subprocess.check_output(['3dinfo', '-oj', inDir +'/' + im1]))
origZ = float(subprocess.check_output(['3dinfo', '-ok', inDir +'/' + im1]))



minPoint = float(subprocess.check_output(['3dinfo', '-Iextent', inDir +'/' + im1]))# + origZ
maxPoint = float(subprocess.check_output(['3dinfo', '-Sextent', inDir +'/' + im1]))# + origZ

minX = float(subprocess.check_output(['3dinfo', '-Rextent', inDir +'/' + im1]))# + origX
maxX = float(subprocess.check_output(['3dinfo', '-Sextent', inDir +'/' + im1]))#+ origX


minY = float(subprocess.check_output(['3dinfo', '-Aextent', inDir +'/' + im1]))# + origY
maxY = float(subprocess.check_output(['3dinfo', '-Pextent', inDir +'/' + im1]))# + origY


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


 
t0 = 0
tf = step
if len(data1.shape)>3:
	sd1 = np.mean(data1, axis=3)
else:
	sd1 = data1
sd1 = nimg.new_img_like(nii1,sd1)
xPos = minX * 0.375
yPos = maxY * 0.2



for i in range(nRow):

	display = nlp.plot_epi(sd1, display_mode='z', figure=fig, draw_cross=False, vmin=300, vmax=2200, cmap='gray', cut_coords=zSlices[t0:tf], axes=(0,rowPlace[i],1,1/nRow), colorbar=True, black_bg=True)

	display.add_edges(bgNii, color=(.1,.1,.75))
	
	nSl = tf-t0
	
	for j in range(nSl):
		display.add_markers([(xPos,yPos,zSlices[t0+j])], marker_color='r', marker_size=75)


	t0 = t0 + step
	tf = tf + step


plt.savefig(outSD1)

t0 = 0
tf = step

fig = plt.figure(figsize=(18,14), dpi=figDpi, facecolor='w', edgecolor='k')

if len(data2.shape)>3:
	sd2 = np.mean(data2, axis=3)
else:
	sd2 = data2
#sd2 = np.mean(data2, axis=3)
sd2 = nimg.new_img_like(nii2,sd2)

for i in range(nRow):

	display = nlp.plot_epi(sd2, display_mode='z', figure=fig, draw_cross=False, vmin=300, vmax=2200, cmap='gray', cut_coords=zSlices[t0:tf], axes=(0,rowPlace[i],1,1/nRow), colorbar=True)
	display.add_edges(bgNii, color=(.1,.1,.75))

	nSl = tf-t0
	for j in range(nSl):
		display.add_markers([(xPos,yPos,zSlices[t0+j])], marker_color='w', marker_size=75)


	t0 = t0 + step
	tf = tf + step



plt.savefig(outSD2)


im1 = Image.open(outSD1)
im2 = Image.open(outSD2)

w, h = im1.size

fnt = ImageFont.truetype('/home/fsluser/Documents/rs_proc/QC_funcs/fonts/DAYPBL__.ttf', 25)


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

os.remove(outSD1)
os.remove(outSD2)


