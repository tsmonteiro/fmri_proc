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
reqoptions.add_argument('-a', '-acq', dest="acq", required=True, help='Acquisition type: [1: sequantial], 2: interleaved [NOT IMPLEMENTED YET], 3: parallel' )
reqoptions.add_argument('-outms', dest="outFileMs", required=True, help='Output file' )
reqoptions.add_argument('-outs', dest="outFileS", required=True, help='Output file' )


args = parser.parse_args()


TR 	= float(args.TR)
outFileMs = args.outFileMs
outFileS = args.outFileS
nSlices = int(args.NSL)
acq = int(args.acq)


if TR > 10:
    # Likely, this was given in ms
    print('Assuming TR was given in milliseconds. If this is not the case, the script needs to be changed.')
    TR = TR / 1000

if acq == 1:
    timing = np.linspace(0, TR, nSlices+1)
    timing = timing[0:nSlices]
    np.savetxt(outFileS, np.transpose(timing), fmt='%.3f')
    np.savetxt(outFileMs, np.transpose(timing*1000), fmt='%.0f')

# Ascending, parallel MB2, Philips
if acq == 20:
    dt = TR/(nSlices/2)

    # dt = (1/TR) * 2

    timing1 = np.linspace(0, TR-dt, int((nSlices)/2))
    timing2 = np.linspace(TR-dt, 0, int((nSlices)/2))

    # Philips separated slices in packages and interleaves them
    # The assumption is that the packages will be acquired in parallel
    # Practically, this means two sequential slices at a time
    timing  = ( np.concatenate((timing1, timing1 )) )
    # This happens when there are is an odd number of slices
    if len(timing) < nSlices:
        timing2 = np.linspace(0, TR-dt, int((nSlices)/2)+1)

        timing = np.concatenate((timing2, timing1 ))

    #

    np.savetxt(outFileS, np.transpose(timing), fmt='%.9f')
    timing = (timing*(TR*1000)).astype(int) # As int ms
    np.savetxt(outFileMs, np.transpose(timing), fmt='%d')

if acq == 2:

    dt = 1/TR

    timing1 = np.linspace(0, TR-dt, int((nSlices)/2))
    timing2 = np.linspace(dt, TR, int((nSlices)/2))

    timing = np.concatenate((timing1, timing2 ))
    if len(timing) < nSlices:
        timing2 = np.linspace(dt, TR, int((nSlices)/2)+1)

        timing = np.concatenate((timing1, timing2 ))
    #timing = (timing*(TR*1000)).astype(int) # As int ms

    np.savetxt(outFile, np.transpose(timing), fmt='%.3f')

if acq == 29:

    dt = 1/TR

    timing1 = np.linspace(TR-dt, dt, int((nSlices)/2))
    timing2 = np.linspace(TR, 0, int((nSlices)/2))

    timing = np.concatenate((timing2, timing1 ))
    if len(timing) < nSlices:
        timing2 = np.linspace(dt, TR, int((nSlices)/2)+1)

        timing = np.concatenate((timing2, timing1 ))
    #timing = (timing*(TR*1000)).astype(int) # As int ms

    np.savetxt(outFile, np.transpose(timing), fmt='%.5f')

if acq == 33:
    # Assuming here that it was ascending, 3 slices at a times

    timing = np.linspace(0, TR- (TR/(nSlices/3)), int((nSlices/3))) #np.linspace(0, TR-(1/TR), int((nSlices/3)))
    timing = np.concatenate((timing, timing,timing))


    if len(timing) < nSlices:
        timing2 = np.linspace(0, TR- (TR/(nSlices/3)), int((nSlices/3))+1) #np.linspace(0, TR-(1/TR), int((nSlices/3)))
        timing = np.linspace(0, TR- (TR/(nSlices/3)), int((nSlices/3)))
        timing = np.concatenate((timing, timing,timing2))

    #timing = (timing*(TR*1))#.astype(int) # As int ms
    #np.savetxt(outFile, timing.reshape([-1,1]), fmt='%.5f', delimiter=',')

    np.savetxt(outFileS, np.transpose(timing), fmt='%.9f')
    timing = (timing*(TR*1000)).astype(int) # As int ms
    np.savetxt(outFileMs, np.transpose(timing), fmt='%d')


if acq == 32:
    # Assuming here that it was ascending, 3 slices at a times
    dt = (1/TR) * 2

    #timing1 = np.linspace(0, TR-dt, int((nSlices)/2))
    #timing2 = np.linspace(dt, TR, int((nSlices)/2))
    timing1 = np.linspace(0, TR-dt, int((nSlices)/2))
    timing2 = np.linspace(0, TR-dt, int((nSlices)/2))

    #timing = np.concatenate((timing1, timing2 ))
    #timing = np.linspace(0, TR- (TR/(nSlices/2) * 2), int((nSlices/2))) #np.linspace(0, TR-(1/TR), int((nSlices/3)))
    timing = np.concatenate((timing1, timing2))
    timing = (TR- (TR/(nSlices/2) * 2)) - (timing*(TR*1)) # As int ms

    np.savetxt(outFile, timing.reshape([-1,1]), fmt='%.5f', delimiter=',')
