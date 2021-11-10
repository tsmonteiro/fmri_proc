#!/usr/bin/env python

# Import required modules
from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from builtins import str
from sklearn.decomposition import PCA
import os
import argparse
import numpy as np
import scipy.signal as sgn

parser = argparse.ArgumentParser(description='Running aCompCor')

# Required options
reqoptions = parser.add_argument_group('Required arguments')
reqoptions.add_argument('-compcor', dest="compcorFname", required=True, help='Directory name' )
reqoptions.add_argument('-fix_ics', dest="fixIcs", required=True, help='Directory name' )
reqoptions.add_argument('-mot12', dest="mot12", required=True, help='Directory name' )
reqoptions.add_argument('-mot24', dest="mot24", required=True, help='Directory name' )
reqoptions.add_argument('-expVar', dest="expVar", required=False, default=0.5, help='Repetition Time' )
reqoptions.add_argument('-svar', dest="svar", required=False, default=1, help='Repetition Time' )
reqoptions.add_argument('-out', dest="outFname", required=True, help='Repetition Time' )

args = parser.parse_args()


compCorFile = args.compcorFname
icsFile     = args.fixIcs

expVar = float(args.expVar)
compCor = np.loadtxt( compCorFile )
ics = np.loadtxt( icsFile )
mot12 = np.loadtxt( args.mot12 )
mot24 = np.loadtxt( args.mot24 )
svar = float(args.svar)/100

if compCor.ndim == 1:
    compCor = np.reshape(compCor, [-1,1])

if ics.ndim == 1:
    ics = np.reshape(ics, [-1,1])


print(compCor.shape[1])
print(ics.shape[1])
print(mot12.shape[1])
print(mot24.shape[1])


nuis = np.concatenate( (compCor, ics, mot12, mot24), axis=1)

nNuis = nuis.shape[1]
for n in range(nNuis):
    nuis[:,n] =  sgn.detrend( (nuis[:,n] - np.mean(nuis[:,n])) / np.std(nuis[:,n]) )



pcaDecomp = PCA(n_components=None, copy=True, whiten=False, svd_solver='auto', tol=0.0, iterated_power='auto', random_state=None)

pcaRes = pcaDecomp.fit(nuis)
cumVar = pcaRes.explained_variance_ratio_

retainedComps = 0
sumVar = 0
idx=1

prev=0
for v in cumVar:
    sumVar += v
    print( "{}: Cumulative explained variance {:.01f}%".format(idx, sumVar*100) )
    idx += 1
    chg = sumVar-prev
    if sumVar < expVar and chg >= svar:
        retainedComps += 1
        prev = sumVar
    else:
        break

print("Retaining {} components out of {}.".format(retainedComps, len(cumVar)) )

comps = pcaRes.transform(nuis)
print(comps.shape)
comps = comps[:,0:retainedComps]

np.savetxt( args.outFname, comps, fmt='%.5f' )
