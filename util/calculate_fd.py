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
reqoptions.add_argument('-mdir', dest="outdir", required=True, help='Threshold for motion spikes [Default 0.5mm]' )

args = parser.parse_args()

outdir    = args.outdir
rpFile    = outdir + '/motion_estimate.par'
outFileFd = outdir + '/maximum_disp.1d_delt'
outFileD  = outdir + '/maximum_disp.1d'


rp = np.loadtxt( rpFile )

drp = np.concatenate((np.zeros((1,6), dtype=float), np.diff( rp, axis=0 )) )
drp_rad = np.concatenate((np.zeros((1,6), dtype=float), np.diff( np.deg2rad(rp), axis=0 )) )

fd_rot = abs(drp_rad[:,0]*50) + abs(drp_rad[:,1]*50) + abs(drp_rad[:,2]*50)
fd_mot = abs(drp[:,3]) + abs(drp[:,4]) + abs(drp[:,5])


fd = fd_rot + fd_mot

maxDisp = np.sum( rp[:,[3,4,5]], axis=1 )

np.savetxt( outFileFd, fd, fmt='%.5f' )
np.savetxt( outFileD, maxDisp, fmt='%.5f' )
