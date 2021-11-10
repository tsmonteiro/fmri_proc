#!/usr/bin/env python


import os
import argparse
import sys


from builtins import str

import matplotlib.pyplot as plt



import numpy as np


from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
from sklearn.linear_model import LinearRegression


def confidence_ellipse(x,y,ax, n_std=3.0,facecolor='none',**kwargs):
	mean_x = np.mean(x)
	mean_y = np.mean(y)

	sd_x   = np.std(x)
	sd_y   = np.std(y)

	r = np.corrcoef(x,y)
	r = r[0,1]
	

	hp = 3.14159/2
	rad_angle=np.arccos(r)
	angle = np.rad2deg(rad_angle)

	ell_radius_x = mean_x + n_std*sd_x/1
	ell_radius_y = mean_y + n_std*sd_y/1

	ellipse = Ellipse((mean_x, mean_y),width=ell_radius_x , height=ell_radius_y, angle=-rad_angle, facecolor='none', **kwargs)


	cos_angle = np.cos(np.radians(rad_angle))
	sin_angle = np.sin(np.radians(rad_angle))

	xc = x - mean_x
	yc = y - mean_y

	xct = xc * cos_angle - yc * sin_angle
	yct = xc * sin_angle + yc * cos_angle 

	rad_cc = (xct**2/(ellipse.width/2.)**2) + (yct**2/(ellipse.height/2.)**2)



	return ax.add_patch(ellipse), rad_cc

#TODO Aggregate plots per measure [FD, DVAR & so on]


parser = argparse.ArgumentParser(description='Save QA check Plots')

# Required options                    
reqoptions = parser.add_argument_group('Required arguments')

reqoptions.add_argument('-o', '-out', dest="outDir", required=True, help='Directory where images are to be saved' )
reqoptions.add_argument('-a', '-in', dest="inDir", required=True, help='Dir where ALL motion summaries are located' )


reqoptions.add_argument('-d', '-dpi', dest="dpi", required=False, default=120, help='Saved figure DPI' )


args = parser.parse_args()


outDir = args.outDir
inDir = args.inDir

# PNG resolution of the saved file
figDpi=int(args.dpi)


# Font size and weight for ALL plots
plt.rcParams.update({'font.size': 20, 'font.weight':'bold'} )

# ++++++++++++++++++++++++++++++++++++++++++++++++
#
# END OF PARAMETER/SETUP
#
# ++++++++++++++++++++++++++++++++++++++++++++++++


#TODO Hardcoded to get age colors, remove later
ageFile = '/mnt/hgfs/ssd_tmp/ages.txt'
ages = np.loadtxt(ageFile, delimiter=',')

ages = ages[0:105,0]
normAges = ages / np.max(ages)


idx=0
for a in ages:
	if a < 35:
		normAges[idx]=1
	elif a < 65:
		normAges[idx]=2
	else:
		normAges[idx]=3
	idx+=1


meanFdAll = []
sdDvarAll = []
maxDistAll= []
sdFdAll   = []

ageColors = []

for fname in os.listdir(inDir):
	if fname.endswith(".txt"):
		motAgg = np.loadtxt(inDir + '/' + fname)
		meanFd = motAgg[0]
		maxDist = motAgg[1]
		sdFd = motAgg[2]
		sdDvar = motAgg[3]

		meanFdAll = np.concatenate( (meanFdAll, [meanFd]) )
		sdDvarAll = np.concatenate( (sdDvarAll, [sdDvar]) )
		maxDistAll = np.concatenate( (maxDistAll, [maxDist]) )
		sdFdAll = np.concatenate( (sdFdAll, [sdFd]) )

		


#TODO subplots for different combination
#TODO linewstyle for different ellipses


fig = plt.figure(figsize=(16,12), dpi=figDpi, facecolor='w', edgecolor='k')
ax = fig.add_subplot(111)

x = meanFdAll
y = sdFdAll


#m,b = np.polyfit(x,y,1)

fit = np.polyfit(x,y,1)
fit_fn = np.poly1d(fit)

xs = np.linspace(0, .6, len(x))


r = np.corrcoef(x,y)
plt.text(0.05, 60, 'r = %.3f' % r[0,1] )
plt.plot(xs, fit_fn(xs), '-', linewidth=3, color=(.5,0,0,0.8))


ax.scatter(x, y, s=100, c=normAges, cmap='gnuplot', vmin=.8, vmax=3.2)



#ax.axvline(c='grey', lw=1)
#ax.axhline(c='grey', lw=1)

ptch, rad_cc1 = confidence_ellipse(x, y, ax, edgecolor=(.5,0.5,1,1),n_std=1.0, linewidth=2, linestyle=':')
ptch, rad_cc2 = confidence_ellipse(x, y, ax, edgecolor=(.1,0.1,1,1),n_std=2.0, linewidth=3, linestyle='--')
ptch, rad_cc3 = confidence_ellipse(x, y, ax, edgecolor=(0,0,0.5,1),n_std=3.0, linewidth=4, linestyle='-')


idx=-1
for r in rad_cc1:
	idx+=1
	if r > 1:
		plt.text(x[idx], y[idx], str(idx+1), {'size':12} )


ax.scatter(np.mean(x), np.mean(y), c='k', s=250, marker='x')

plt.xlabel('Mean FD')
plt.ylabel('FD STD')

plt.savefig(outDir + '/fd_fd_group.png' )

