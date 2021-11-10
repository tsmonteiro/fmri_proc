#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 12:37:24 2020

@author: u0101486
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 22 15:36:33 2020

@author: u0101486
"""

import numpy as np

import argparse

from sklearn.decomposition import PCA

parser = argparse.ArgumentParser(description='Run CanICA estimation [with some melodic structure]')

# Required options
reqoptions = parser.add_argument_group('Required arguments')

reqoptions.add_argument('-nuis_file', dest="nuisFile", required=True, help='Directory where images are to be saved' )
reqoptions.add_argument('-varexp', dest="varExp", required=True, help='Directory where images are to be saved' )
reqoptions.add_argument('-outfile', dest="outFile", required=True, help='Directory where images are to be saved' )


args = parser.parse_args()


nuisFile= args.nuisFile
varExp  = float(args.varExp)


##%%
#eigFile = '/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/tmp/N1_01/pc_var_eig.1D'
# The eig file is composed of 4 columns
# #Num.  --Eigenvalue--  -Var.Fraction-  -Cumul.Fract.-
nuis = np.loadtxt( nuisFile )

nNuis = nuis.shape[1]

for i in range( nNuis ):
    if np.nanstd(nuis[:,i]) > 0:
        nuis[:,i] = (nuis[:,i] - np.nanmean(nuis[:,i]))/ np.nanstd(nuis[:,i])
    else:
        nuis[:,i] = nuis[:,i] * 0


pca_obj = PCA(n_components=None)
princomp = pca_obj.fit_transform(nuis)


val = []

nComps = pca_obj.explained_variance_ratio_.shape[0]


modelOrder  = 0
sumTotal    = 0
for i in range( nComps ):
        print('{:.05f} >? {:.1f}'.format(sumTotal, varExp))
        if  sumTotal > varExp:

            break
        modelOrder += 1
        sumTotal   += (100*pca_obj.explained_variance_ratio_[i])

if modelOrder == 0:
    modelOrder  = int( np.ceil(nComps/2) )

print('Extracting {} out of {} components responsible for {:.01f}% of total variance'.format(
        int(modelOrder), int(nNuis), sumTotal ))

np.savetxt(args.outFile, princomp[:,0:modelOrder], fmt='%.5f')
#
