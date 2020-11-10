#!/usr/bin/env python

# Import required modules
from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from builtins import str
import os
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Running aCompCor')

# Required options
reqoptions = parser.add_argument_group('Required arguments')
reqoptions.add_argument('-d', '-dir', dest="dir", required=True, help='Directory name' )
reqoptions.add_argument('-r', '-rp', dest="rp", required=True, help='Motion estimates file [rotation first]' )
reqoptions.add_argument('-t', '-fdt', dest="fdt", required=False, help='Threshold for motion spikes [Default 0.5mm]' )

args = parser.parse_args()

sessDir = args.dir
rpFile = sessDir + '/' + args.rp
#icFile = sessDir + '/ic_tcs.txt'
outFile_12 = sessDir + '/' + 'motion_regressors_12.txt'
outFile_24 = sessDir + '/' + 'motion_regressors_24.txt'
outFile_1224 = sessDir + '/' + 'motion_regressors_12_24.txt'
outSpikeFile = sessDir + '/' + 'temporal_mask_fd.txt'

fdt = args.fdt



if not fdt:
	fdt = 0.5
else:
	fdt = float(fdt)


rp = np.loadtxt( rpFile )
#ics = np.loadtxt( icFile )
drp = np.concatenate((np.zeros((1,6), dtype=float), np.diff( rp, axis=0 )) )
drp_rad = np.concatenate((np.zeros((1,6), dtype=float), np.diff( np.deg2rad(rp), axis=0 )) )

rp2 = np.power(rp, 2)
drp2 = np.power(drp, 2)

mot_12 = np.concatenate( (rp,drp), axis=1)
mot_24 = np.concatenate( (rp2,drp2), axis=1  )
mot_12_24 = np.concatenate( (rp,drp,rp2,drp2), axis=1  )

np.savetxt( outFile_12, mot_12, fmt='%.5f' )
np.savetxt( outFile_24, mot_24, fmt='%.5f' )
np.savetxt( outFile_1224, mot_12_24, fmt='%.5f' )





fd_rot = abs(drp_rad[:,0]*50) + abs(drp_rad[:,1]*50) + abs(drp_rad[:,2]*50)
fd_mot = abs(drp[:,3]) + abs(drp[:,4]) + abs(drp[:,5])


fd = fd_rot + fd_mot



N = fd.shape
N = N[0]

censor_list = np.ones( (N, 1) )

i=-1
for displace in fd:
	i+=1
	if displace > fdt:
		censor_list[i] = 0



np.savetxt( outSpikeFile, censor_list, fmt='%d' )
