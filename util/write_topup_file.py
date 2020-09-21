#!/usr/bin/env python


import os
import argparse
import sys

from builtins import str
import numpy as np


parser = argparse.ArgumentParser(description='Save QA check Plots')

# Required options
reqoptions = parser.add_argument_group('Required arguments')

reqoptions.add_argument('-o', '-out', dest="outFile", required=True, help='Output file' )


args = parser.parse_args()

outFile = args.outFile

topup = np.zeros((10,4))
for l in range(5):
    topup[l,0] = 0
    topup[l,1] = 1
    topup[l,2] = 0
    topup[l,3] = 0.0267

    topup[l+5,0] = 0
    topup[l+5,1] = -1
    topup[l+5,2] = 0
    topup[l+5,3] = 0.0267


np.savetxt(outFile, topup, fmt='%.4f')
