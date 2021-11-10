#!/usr/bin/env python
'''
Deals with motion coming from different software packages so they can be used uniformly in the QC scripts

'''

import numpy as np

def convert_motion(mpe, prog='AFNI'):

	prog = prog.lower()
	cmpe = mpe
	if prog == 'afni':
		# Nothing to do, that is what I expect
		cmpe = mpe

	if prog == 'spm':
		spmMpe = mpe
		cmpe[:,0] = np.rad2deg(mpe[:,3])
		cmpe[:,1] = np.rad2deg(mpe[:,4])
		cmpe[:,2] = np.rad2deg(mpe[:,5])
		cmpe[:,3] = mpe[:,0]
		cmpe[:,4] = mpe[:,1]
		cmpe[:,5] = mpe[:,2]

	return cmpe



def framewise_displacement(rp):
	fd1 =  np.abs(np.diff(np.deg2rad( rp[:,0] ))*50) + np.abs(np.diff(np.deg2rad( rp[:,1] ))*50) + np.abs(np.diff(np.deg2rad( rp[:,2] ))*50) 
	fd2 =  np.abs(np.diff(rp[:,3] )) + np.abs(np.diff(rp[:,4] )) + np.abs(np.diff(rp[:,5] ))

	fd = fd1 + fd2

	fd = np.insert(fd, 0, 0)
	


	return fd


def total_translation(rp):
	return np.sqrt( np.power(rp[:,3],2) + np.power(rp[:,4],2) + np.power(rp[:,5],2) )


def total_rotation(rp, rot_type='degree'):
	if rot_type == 'radian':
		rp = np.deg2rad(rp)
	

	return np.sqrt( np.power(rp[:,0],2) + np.power(rp[:,1],2) + np.power(rp[:,2],2) )


