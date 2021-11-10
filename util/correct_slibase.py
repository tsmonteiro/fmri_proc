#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 15:00:01 2020

@author: u0101486
"""

# Hack-ish correction of belt disconnection....

import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Run CanICA estimation [with some melodic structure]')

# Required options
reqoptions = parser.add_argument_group('Required arguments')

reqoptions.add_argument('-slibase', dest="sliFile", required=True, help='Philips physiological recording' )
reqoptions.add_argument('-sliout', dest="sliOut", required=True, help='Philips physiological recording' )

args = parser.parse_args()



sliFile = args.sliFile#'/home/luna.kuleuven.be/u0101486/workspace/data/ConnectEx/tmp/A003/physiological_regressors.slibase.1D'
sliOFile = args.sliOut#'/home/luna.kuleuven.be/u0101486/workspace/data/ConnectEx/tmp/A003/physiological_regressors_c.slibase.1D'

sliOrvt = args.sliOut + '.rvt' #'/home/luna.kuleuven.be/u0101486/workspace/data/ConnectEx/tmp/A003/physiological_regressors_c.slibase.1D'
sliOresp = args.sliOut + '.resp' #'/home/luna.kuleuven.be/u0101486/workspace/data/ConnectEx/tmp/A003/physiological_regressors_c.slibase.1D'
sliOcard = args.sliOut + '.card' #'/home/luna.kuleuven.be/u0101486/workspace/data/ConnectEx/tmp/A003/physiological_regressors_c.slibase.1D'

data = np.loadtxt(sliFile)

r,c = data.shape



skipCols = 13

rvtCols = np.sort(np.concatenate( (np.arange(0,c,skipCols), np.arange(1,c,skipCols),
                                   np.arange(2,c,skipCols), np.arange(3,c,skipCols),
                                   np.arange(4,c,skipCols) ) ) ) 

respCols = np.sort(np.concatenate( (np.arange(5,c,skipCols), np.arange(6,c,skipCols),
                                   np.arange(7,c,skipCols), np.arange(8,c,skipCols) ) ) ) 

cardCols = np.sort(np.concatenate( (np.arange(9,c,skipCols), np.arange(10,c,skipCols),
                                   np.arange(11,c,skipCols), np.arange(12,c,skipCols) ) ) ) 


for k in np.arange(c):
    if np.var(data[:,k]) == 0:
        data[:,k] = np.random.rand(r)
        
        
np.savetxt(sliOFile, data)

np.savetxt(sliOrvt, data[:,list(rvtCols)])
np.savetxt(sliOresp, np.transpose( data[:,list(respCols)]))
np.savetxt(sliOcard, np.transpose( data[:,list(cardCols)]))
        
