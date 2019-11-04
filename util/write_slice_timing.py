#!/usr/bin/env python


import os
import argparse
import sys

from builtins import str
import numpy as np


parser = argparse.ArgumentParser(description='Save QA check Plots')

# Required options                    
reqoptions = parser.add_argument_group('Required arguments')

reqoptions.add_argument('-t', '-tr', dest="TR", required=True, help='Time to repeat' )
reqoptions.add_argument('-n', '-nsl', dest="NSL", required=True, help='Number of acquisition slices' )
reqoptions.add_argument('-a', '-acq', dest="acq", required=True, help='Acquisition type: [1: sequantial], 2: interleaved [NOT IMPLEMENTED YET]' )
reqoptions.add_argument('-o', '-out', dest="outFile", required=True, help='Output file' )


args = parser.parse_args()


TR 	= float(args.TR)
outFile = args.outFile
nSlices = int(args.NSL)
acq = int(args.acq)


if TR > 10:
	# Likely, this was given in ms
	print('Assuming TR was given in milliseconds. If this is not the case, the script needs to be changed.')
	TR = TR / 1000

if acq == 1:
	timing = np.linspace(0, TR, nSlices+1)
	timing = timing[0:nSlices]	

np.savetxt(outFile, timing, fmt='%.3f')


