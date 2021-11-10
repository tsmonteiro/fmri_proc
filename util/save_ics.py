#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 12:05:42 2020

Creates 2 text files, noise and not for ICs

@author: u0101486
"""


import os

import numpy as np


import argparse

parser = argparse.ArgumentParser(description='Run CanICA estimation [with some melodic structure]')

# Required options
reqoptions = parser.add_argument_group('Required arguments')

reqoptions.add_argument('-out',  dest="outDir", required=True, help='Directory where images are to be saved' )
reqoptions.add_argument('-class',  dest="icClass", required=True, help='Directory where images are to be saved' )


args    = parser.parse_args()

outDir  = args.outDir
icClass = args.icClass

#%%

#outDir = '/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/tmp/N1_07/'
#icClass = '/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/tmp/B1_02/FIX/fix4melview_classifier_thr60.txt'
#icClass = '/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/tmp/N1_07/FIX/hand_labels_noise.txt'

#outDir = '/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/tmp/B1_07/'
#icClass = '/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/tmp/B1_07/FIX/hand_labels_noise.txt'

fid   = open(icClass, 'r')
lines = fid.readlines()

if len(lines[-1]) == 1:
    noiseComps = lines[-2]
else:
    noiseComps = lines[-1]
noiseComps = noiseComps[1:-2]
noiseComps = np.array( noiseComps.split(',') ).astype(int)

#print(noiseComps)


##%%


ics  = []
nics = []

noiseIcs = ''

for i in range(1000):
    tcFile = outDir + '/FIX/filtered_func_data.ica/report/t' + str(i) + '.txt'

    if os.path.isfile(tcFile):
        icTc = np.loadtxt(tcFile)

        if len( np.intersect1d( i, noiseComps ) ) > 0:
            noiseIcs= noiseIcs +  str(i) + ','
            nics.append(icTc)
        else:
            ics.append(icTc)

#print(noiseIcs[0:-1])

np.savetxt(outDir + '/noise_ics.txt', np.transpose(np.array(nics)))
np.savetxt(outDir + '/non_noise_ics.txt', np.transpose(np.array(ics)))
