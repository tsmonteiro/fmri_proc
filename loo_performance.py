#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 11:18:19 2020

@author: u0101486
"""

import numpy as np
import argparse
import os
import matplotlib.pyplot as plt
parser = argparse.ArgumentParser(description='Compare classifier with hand marked artifacts')

#%% Required options                    
reqoptions = parser.add_argument_group('Required arguments')

reqoptions.add_argument('-c', '-classDir', dest="classDir", required=True, help='Directory with LOO results' )
reqoptions.add_argument('-o', '-out', dest="outFile", required=True, help='Save results' )

args = parser.parse_args()


classDir = args.classDir
outFile  = args.outFile



#%%
resultPrefix = 'classifier_'
idThr = []

tpr = []
sigPos = []

#classDir = '/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/fix_classifier_rs/'
#classDir = '/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/fix_classifier_belgium/'
#classDir = '/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/fix_classifier_norway/'
classDir = '/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/fix_classifier_germany_b/'
#classDir = '/home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/fix_classifier_rs/'
#classDir = '/home/luna.kuleuven.be/u0101486/workspace/data/RSPET/fix_classifier/'
for fname in os.listdir(classDir):
    print(fname)
    if fname.startswith(resultPrefix) and fname.endswith(".txt"): 
        thr = str.split(fname, '_')
        thr = str.split(thr[2],'.')
        thr = int(thr[0])
        idThr.append(thr)
        
        classAcc = np.loadtxt(classDir + fname)
        tpr.append(classAcc[0])
        sigPos.append(classAcc[1])
        

#classifier_0_1.txt
#%%
plt.close('all')
uniqueThrs = np.unique(idThr)
yticks = [60, 100]
for thr in uniqueThrs:
    idx = np.where( np.array(idThr) == thr )[0]
    idx = np.array( np.reshape( idx, [-1,1]))
    idx = idx.astype(int)
    tpr = np.array(tpr)
    sigPos = np.array(sigPos)
    
    h1 = plt.scatter(thr, np.mean( tpr[ idx ] ), 80, color=(.4,.4,1), edgecolor=(0,0,0))
    h11 = plt.scatter(thr, np.quantile( tpr[ idx ], .2 ), 80, color=(.8,.8,1), edgecolor=(0,0,0))
    h12 =plt.scatter(thr, np.quantile( tpr[ idx ], .05 ), 80, color=(1,1,1), edgecolor=(0,0,0))
    
    #plt.vlines(thr, np.quantile(tpr, .05), np.quantile(tpr, .95),  color=(.4,.4,1), zorder=0)
    #plt.hlines(np.quantile(tpr, 95), thr-.5, thr+.5, color=(.4,.4,1))
    #plt.hlines(np.quantile(tpr, 5), thr-.5, thr+.5, color=(.4,.4,1))
    
    h2 = plt.scatter(thr, np.mean( sigPos[ idx ] ), 80, color=(1,.6,.6), edgecolor=(0,0,0))
    
    plt.plot([0, thr], [np.mean( tpr[ idx ] ), np.mean( tpr[ idx ] )], color=(.8,.8,1), zorder=0)
    plt.plot([0, thr], [np.mean( sigPos[ idx ] ), np.mean( sigPos[ idx ] )], color=(1,.9,.9), zorder=0)
    
    yticks.append(np.round(np.median( tpr[ idx ] )))
    yticks.append(np.round(np.median( sigPos[ idx ] )))

yticks = np.unique(np.sort(yticks))
plt.xticks(uniqueThrs, labels=uniqueThrs, fontsize=12 )

plt.yticks(yticks, labels=yticks, fontsize=12 )
plt.xlabel('ICA-FIX Threshold [a.u.]', fontsize=14)
plt.ylabel('Performance [%]', fontsize=14)

plt.legend([h1, h11, h12, h2], ['% of True Noise Identified',
           '% of True Noise Identified [80% percentile]','% of True Noise Identified [95% percentile]',
           '% True Signal Identifies'],
           fontsize=11)