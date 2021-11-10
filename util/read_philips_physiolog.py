#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 14:37:47 2020

@author: u0101486
"""

import numpy as np
import argparse

def hex2dec(s):
    return int(s, 16)



parser = argparse.ArgumentParser(description='Run CanICA estimation [with some melodic structure]')

# Required options
reqoptions = parser.add_argument_group('Required arguments')

reqoptions.add_argument('-inDir', dest="inDir", required=True, help='Philips physiological recording' )
reqoptions.add_argument('-TR', dest="TR", required=True, default=2, help='Philips physiological recording' )
reqoptions.add_argument('-vols', dest="nVols", required=True, default=2, help='Philips physiological recording' )

reqoptions.add_argument('-fs', dest="fs", required=False, default=500, help='Philips physiological recording' )

reqoptions.add_argument('-lp', dest="lp", required=False, default=-1, help='Philips physiological recording' )
reqoptions.add_argument('-despike', dest="despike", required=False, default=0, help='Philips physiological recording' )


args = parser.parse_args()


inDir        = args.inDir

lp      = float(args.lp)
despike = int(args.despike)

#inDir      = '/home/luna.kuleuven.be/u0101486/workspace/data/ConnectEx/tmp/A001/'

physioFile = inDir + '/physio.log'
outPhys    = inDir + '/physio_edit.log'
outResp    = inDir + '/physio.resp'
outCard    = inDir + '/physio.card'

fs    = float(args.fs) #500
TR    = float(args.TR) #1
nVols = int(args.nVols) #600
data  = np.loadtxt(physioFile)




ppu  = data[:,4]
resp = data[:,5]

gradX = data[:,6]
gradY = data[:,7]
gradZ = data[:,8]

allGrad   = np.abs( gradX + gradY + gradZ )
scanStart = int(np.argmax(allGrad) + fs*3)

dataOut   = data[scanStart:,:]
np.savetxt( fname=outPhys, X=dataOut.astype(int), fmt='%d' )


markersDec = data[:,9]
markers    = []


for mrk in markersDec:
    markers.append(hex2dec(str(int(mrk))))

markers = np.array( markers )

t0 =np.where(markers == 16)[-1][-1] # Scan start trigger
tf = len(markers) #np.where(markers == 32)[-1][-1] # Scan stop  trigger



agx = abs(gradX[t0:tf])
agy = abs(gradY[t0:tf])
gy = (gradY[t0:tf])
agz = abs(gradZ[t0:tf])

appu  = ppu[t0:tf]
aresp = resp[t0:tf]


# Is this consistent?
t = 3500

scanTriggers = []
numScans     = 0


#matplotlib.pyplot.rcParams['agg.path.chunksize'] = 20000
while t < len(agx):
   if agx[t] > 0 or agy[t] > 0 or agx[t] > 0:
       t2 = t-1;
       while 1 == 1 and t2 < len(agx)-1:
           if (agx[t2] + agy[t2] + agx[t2]) < 10:
               break
           t2 = t2 + 1


       numScans = numScans + 1
       scanTriggers.append(t2)
       
       
       
       t = int(t2 + fs * TR -67)
       #t = t + 499 - 11
   else:
       t = t + 1



dmt = np.mean(np.diff(scanTriggers))


t = len(gradX)-1
while t > 0 and abs(gradX[t]) < 10 and abs(gradY[t]) < 10 and abs(gradZ[t]) < 10:
    t = t - 1

sf  = t #int(scanTriggers[-1] + (TR*fs) - 1)
#sf  = int(tf) #int(scanTriggers[-1] + (TR*fs) - 1)
s0  = int(tf-nVols*fs) #int(sf - nVols*fs )


print('=================================================')
print('There are ' + str(numScans) + ' scans')
print('Average difference between triggers = {:.03f} seconds'.format(dmt/fs))
print((tf-t0)/fs)
#print(len(resp))
#print(s0)
#print(sf)
#print(len(ppu[s0:sf]))
#print(len(resp[s0:sf]))
print('=================================================')

#print(scanTriggers)
#import scipy.signal.medfilt
from scipy.signal import medfilt
from nilearn.signal import clean
if despike > 0:
    ppu[s0:sf]  = medfilt(ppu[s0:sf], kernel_size=3)
    resp[s0:sf] = medfilt(resp[s0:sf], kernel_size=3)
    
if lp > 0:
    ppu[s0:sf] = clean(ppu[s0:sf], detrend=False, standardize=False, t_r=(1/fs),
                       low_pass=lp)
    
    resp[s0:sf] = clean(resp[s0:sf], detrend=False, standardize=False, t_r=(1/fs),
                       low_pass=lp)

np.savetxt( fname=outCard, X=ppu[s0:sf] )
np.savetxt( fname=outResp, X=resp[s0:sf] )

#np.savetxt( fname=outCard, X=ppu[s0+16000:s0+32000] )
#np.savetxt( fname=outResp, X=resp[s0+16000:s0+32000] )

