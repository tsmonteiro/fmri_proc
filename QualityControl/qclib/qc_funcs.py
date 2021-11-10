import nibabel as nib
from builtins import str

import matplotlib.pyplot as plt
import numpy as np


import nilearn.plotting as nlp
import nilearn.image as nimg

import motion_handler as mh



def motion_estimate(rpFile, outFname, range_trans=1.5, range_rot=1.5, tr=0, prog='AFNI', dpi=150):
	'''motion_estimate(rpFile, outFname, range_trans=1.5, range_rot=1.5, tr=0, prog='AFNI', dpi=150):
	
	   Plots translation and rotation estimates (in separate subplots)
	   

  	   Parameters
	   --------------------------------

	   [REQUIRED]
	   rpFile: string
		   Motion estimate files. File is expected to be in the Volumes x Parameters (rows x columns).
	   outFname: string
		     Output file for the graphics

	   [OPTIONAL]
	   range_trans: int
			range of the plot for translation parameters
	   range_rot: int
			range of the plot for rotation parameters
	   tr: float,
		Repetition time. It is only used to plot faint vertical gray lines every 60 seconds
	   prog: {'AFNI', 'SPM', 'FSL'}
		Program used to generate the motion file . 
	   dpi: int
		Output figure DPI

	'''
	# Font size and weight for ALL plots
	plt.rcParams.update({'font.size': 20, 'font.weight':'bold'} )

	rp = mh.convert_motion(np.loadtxt(rpFile), prog=prog)

	N = rp.shape
	N = N[0]


	# Roll, pitch yaw, dS, dL, dP
	fig = plt.figure(figsize=(16,12), dpi=dpi, facecolor='w', edgecolor='k')
	plt.subplot(211)

	time = []
	[time.append(x+1) for x in range(N)]


	plt.plot([0, N], [0, 0], color=[0.25,0.25,0.25], linewidth=2, linestyle=':')

	if tr > 0:
		t = 60/tr # One-minute mark
		while t < (N*tr):
			plt.plot([t, t], [-range_rot, range_rot], color=[0.75,0.75,0.75], linewidth=1, linestyle='--')
			t += 60/tr

	plt.plot(time,  rp[:,0], label='Roll', color=[.6,.6,1], linewidth=3)
	plt.plot(time,  rp[:,1], label='Pitch', color=[.0,.0,1], linewidth=3)
	plt.plot(time,  rp[:,2], label='Yaw', color=[0,0,.5], linewidth=3)


	plt.ylabel('Rotation [degs]')


	plt.axis((0, N, -range_rot, range_rot))
	ticks = np.linspace(-range_rot, range_rot, 5)

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
			plt.plot([t, t], [-range_trans, range_trans], color=[0.75,0.75,0.75], linewidth=1, linestyle='--')
			t += 60/tr


	plt.plot(time,  rp[:,3], label='dS', color=[1,.6,.6], linewidth=3)
	plt.plot(time,  rp[:,4], label='dL', color=[1,.0, 0], linewidth=3)
	plt.plot(time,  rp[:,5], label='dP', color=[.5,0,0], linewidth=3)

	plt.xlabel('TIME [Volume]')
	plt.ylabel('Shift [mm]')


	plt.axis((0, N, -range_trans, range_trans))
	ticks = np.linspace(-range_trans, range_trans, 5)

	ax = plt.gca()

	for axis in ['left', 'bottom']:
		ax.spines[axis].set_linewidth(3)

	for axis in ['top', 'right']:
		ax.spines[axis].set_linewidth(0)


	plt.legend()
	fig.tight_layout()

	plt.savefig(outFname)


