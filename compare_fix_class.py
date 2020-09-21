#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 11:18:19 2020

@author: u0101486
"""

import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Compare classifier with hand marked artifacts')

#%% Required options                    
reqoptions = parser.add_argument_group('Required arguments')

reqoptions.add_argument('-t', '-hand', dest="hand_labels", required=True, help='Labels created by hand' )
reqoptions.add_argument('-l', '-thr', dest="thr", required=True, help='Labels created by hand' )
reqoptions.add_argument('-f', '-fix', dest="fix_labels", required=True, help='Labels created by FIX' )
reqoptions.add_argument('-o', '-out', dest="outFile", required=True, help='Save results' )

args = parser.parse_args()


handLabels = args.hand_labels
fixLabels  = args.fix_labels
outFile  = args.outFile

#%%
thr = int(args.thr)
#for thr in [1,5,10,20,30,40,50]:
#handLabels = '/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/RS031/FIX/hand_labels_noise.txt'
#fixLabels = '/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/RS031/FIX/fix4melview_classifier_thr' + str(thr) + '.txt'


with open(handLabels, 'r') as file:
    handLabels = file.readline()

handLabels = str.split( handLabels[1:-2], ',' )

lineCount = 0
nComps = 0
noiseComps = []
signalComps = []
with open(fixLabels, 'r') as file:
    line = file.readline()
    
    while line:
       
       isNoise = 0
       if lineCount > 0:
           entry = str.split(line, ',')
           
           if str.strip(entry[2]) == 'True':
               isNoise = 1
           
       if lineCount > 0 and line[0] != '[':
           nComps += 1
           if isNoise == 1:
               noiseComps.append(lineCount)
           else:
               signalComps.append(lineCount)
                   
       #print("Line {}{}: {}".format(lineCount, isNoise, line.strip()))
       lineCount +=1
       line = file.readline()

handNoise = []
handSignal = []

for c in range(nComps):
    comp = c + 1
    
    isHandLabel = False
    for l in handLabels:
        if comp == int(l):
            isHandLabel = True
    
    if isHandLabel == True:
        handNoise.append(comp)
    else:
        handSignal.append(comp)
    
    
    

tpNoise = 0
fpNoise = 0
tnNoise = 0
fnNoise = 0

tpSignal = 0
fpSignal = 0
tnSignal = 0
fnSignal = 0

for c in noiseComps:
    p = [c == h for h in handNoise]
    if np.any(p):
        tpNoise += 1
    else:
        fpNoise += 1

for c in handNoise:
    p = [c == h for h in noiseComps]
    if not np.any(p):
        fnNoise += 1


for c in handNoise:
    p = [c == h for h in noiseComps]
    if not np.any(p):
        fnNoise += 1
    
for c in signalComps:
    p = [c == h for h in handSignal]
    if np.any(p):
        tpSignal += 1
    else:
        fpSignal += 1

print('Threshold ' + str(thr))
print('\tTrue Positive Rate Noise: {:.01f}%'.format(100*tpNoise/(tpNoise + fnNoise) )  )
print('\tSignal Components Detected: {:.01f}%'.format(100*tpSignal/len(handSignal) )  )
        
        
np.savetxt(outFile, [100*tpNoise/(tpNoise + fnNoise), 100*tpSignal/len(handSignal) ], '%0.2f')

