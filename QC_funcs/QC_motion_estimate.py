#!/usr/bin/env python


import os
import argparse
import sys

import nibabel as nib
from builtins import str

import matplotlib.pyplot as plt
import numpy as np


import nilearn.plotting as nlp
import nilearn.image as nimg

import motion_handler as mh




# ++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# END OF function definitions
#
# ++++++++++++++++++++++++++++++++++++++++++++++


parser = argparse.ArgumentParser(description='Save QA check Plots')

# Required options                    
reqoptions = parser.add_argument_group('Required arguments')

reqoptions.add_argument('-o', '-out', dest="outDir", required=True, help='Directory where images are to be saved' )
reqoptions.add_argument('-a', '-in', dest="inDir", required=True, help='Dir where files are stored' )


reqoptions.add_argument('-e', '-mpe', dest="mpe", required=True, default='motion_estimate.par', help='Motion Estimate File' )
reqoptions.add_argument('-x', '-outf', dest="outFname", required=True, default='motion_estimate.png', help='Name of the output plot' )

reqoptions.add_argument('-t', '-rngT', dest="rngT", required=False, default=1.5, help='Range of translation plot [+-1.5]' )
reqoptions.add_argument('-r', '-rngR', dest="rngR", required=False, default=1.5, help='Range of rotation plot [+-1.5]' )

reqoptions.add_argument('-k', '-tr', dest="tr", required=False, default=0, help='TR ' )


reqoptions.add_argument('-p', '-prog', dest="prog", required=False, default='AFNI',  help='Software which generated the motion estimate: [AFNI], FSL or SPM' )
reqoptions.add_argument('-d', '-dpi', dest="dpi", required=False, default=120, help='Saved figure DPI' )


args = parser.parse_args()


outDir = args.outDir
inDir = args.inDir
outFile = outDir + '/' + args.outFname
rpFile = inDir + '/' + args.mpe

tr = float(args.tr)
if tr > 10:
	tr = tr / 1000

# PNG resolution of the saved file
figDpi=int(args.dpi)
rngT = float(args.rngT)
rngR = float(args.rngR)

# Font size and weight for ALL plots
plt.rcParams.update({'font.size': 20, 'font.weight':'bold'} )


# ++++++++++++++++++++++++++++++++++++++++++++++++
#
# END OF PARAMETER/SETUP
#
# ++++++++++++++++++++++++++++++++++++++++++++++++


print('\n\n ===============================================================================================\n\n')
print('Starting QC : Motion Estimate Plot')
print('\n\n ===============================================================================================\n\n')


rp = mh.convert_motion(np.loadtxt(rpFile), prog=args.prog)

N = rp.shape
N = N[0]


# Roll, pitch yaw, dS, dL, dP
fig = plt.figure(figsize=(16,12), dpi=figDpi, facecolor='w', edgecolor='k')
plt.subplot(211)

time = []
[time.append(x+1) for x in range(N)]


plt.plot([0, N], [0, 0], color=[0.25,0.25,0.25], linewidth=2, linestyle=':')

if tr > 0:
	t = 60/tr # One-minute mark
	while t < (N*tr):
		plt.plot([t, t], [-rngR, rngR], color=[0.75,0.75,0.75], linewidth=1, linestyle='--')
		t += 60/tr

plt.plot(time,  rp[:,0], label='Roll', color=[.6,.6,1], linewidth=3)
plt.plot(time,  rp[:,1], label='Pitch', color=[.0,.0,1], linewidth=3)
plt.plot(time,  rp[:,2], label='Yaw', color=[0,0,.5], linewidth=3)



# plt.xlabel('TIME [s]')
plt.ylabel('Rotation [degs]')
#plt.title('ROTATION')


plt.axis((0, N, -rngR, rngR))
ticks = np.linspace(-rngR, rngR, 5)

plt.yticks(ticks)

plt.legend()

ax = plt.gca()

for axis in ['left', 'bottom' ]:
	ax.spines[axis].set_linewidth(3)

for axis in ['top', 'right']:
	ax.spines[axis].set_linewidth(0)


plt.subplot(212)
plt.plot([0, N], [0, 0], color=[0.25,0.25,0.25], linewidth=2, linestyle=':')

if tr > 0:
	t = 60/tr # One-minute mark
	while t < (N*tr):
		plt.plot([t, t], [-rngT, rngT], color=[0.75,0.75,0.75], linewidth=1, linestyle='--')
		t += 60/tr


plt.plot(time,  rp[:,3], label='dS', color=[1,.6,.6], linewidth=3)
plt.plot(time,  rp[:,4], label='dL', color=[1,.0, 0], linewidth=3)
plt.plot(time,  rp[:,5], label='dP', color=[.5,0,0], linewidth=3)

plt.xlabel('TIME [Volume]')
plt.ylabel('Shift [mm]')
#plt.title('TRANSLATION')

plt.axis((0, N, -rngT, rngT))
ticks = np.linspace(-rngT, rngT, 5)

ax = plt.gca()

for axis in ['left', 'bottom']:
	ax.spines[axis].set_linewidth(3)

for axis in ['top', 'right']:
	ax.spines[axis].set_linewidth(0)


plt.legend()
fig.tight_layout()




plt.savefig(outFile)


