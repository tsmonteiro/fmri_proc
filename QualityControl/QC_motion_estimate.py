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


import configparser

config = configparser.ConfigParser()
config.read('params.ini')

#PATH to qclib
sys.path.append(config['QCLIB_PATH'])

import motion_handler as mh
import qclib as qc



parser = argparse.ArgumentParser(description='Save motion estimate plots')

# Required options                    
reqoptions = parser.add_argument_group('Function arguments')

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


qc.motion_estimate(rpFile, outFname, range_trans=rngT, range_rot=rngR, tr=tr, prog=args.prog, dpi=figDpi)


