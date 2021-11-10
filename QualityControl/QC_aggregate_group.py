#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 09:39:23 2020

@author: u0101486
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 12:26:49 2019

@author: u0101486
"""

# Aggregate QC measures

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import spearmanr    

import configparser

config = configparser.ConfigParser()
config.read('/home/luna.kuleuven.be/u0101486/workspace/fmri_proc/params.ini')

#PATH to qclib
sys.path.append(config['PATHS']['QCLIB_PATH'])
sys.path.append('/home/luna.kuleuven.be/u0101486/workspace/fmri_proc/ext/nitransforms/nitransforms')


import qclib.common_funcs as gpc
import shutil

from warnings import filterwarnings
filterwarnings('ignore')

import pandas as pd
import statsmodels.formula.api as smf
import multiprocessing as mp
import time
from functools import partial

    


# ======================================================================
# ======================================================================
    
    
#project = 'CRUNCH'
#project = 'CAI_China'

project  = 'RepImpact'


baseDir  = '/home/luna.kuleuven.be/u0101486/workspace/data/' + project + '/tmp/'
qcDir    = '/home/luna.kuleuven.be/u0101486/workspace/data/' + project + '/Quality_Control/FinalProc_20200714/Files/'

distMat  = '/home/luna.kuleuven.be/u0101486/workspace/fmri_proc/atlas/CAT12/lg400_cobra_distance.txt'
distMat  = np.loadtxt(distMat, delimiter=',')

distMat  = distMat[0:400, 0:400]

nNodes   = distMat.shape[0]
distVec  = np.reshape(distMat, (1,nNodes*nNodes))

meanDt   = np.mean(distVec)
stdDt    = np.std(distVec)

distSort = np.argsort(distVec)

distThreshold = 0.1

i = 0
j = 1

while j < distVec.shape[1]:
    dist = np.sqrt( (distVec[0,distSort[0,i]]-distVec[0,distSort[0,j]])**2)
    
    
    if dist < distThreshold or dist == 0:
        distVec[0,distSort[0,j]] = -100
    else:
        i = j
        
    j += 1


#print( len(distVec[np.where(distVec > -100)]) )

#plt.plot(np.sort(distVec[np.where(distVec > -100)]), '.')

##%%
import nibabel as nib

atlas=nib.load('/home/luna.kuleuven.be/u0101486/workspace/fmri_proc/atlas/CAT12/lg400_cobra.nii')
atlas=atlas.get_fdata()
x,y,z=atlas.shape
atlas[np.where(np.isnan(atlas))]=0
uNodes=np.unique(atlas.reshape((x*y*z))).reshape([1,-1]).astype(int)




#statModels = ['NONE', 'SFIX', 'SFIX_D', 'SRP24WM1CSF1', 'SRP24CC', 'SRP9']
procBase   = ['']
statModels = ['FAC_DiC_RP6', 'FAC_DiC_RP24', 'FAC', 'FAC_CC_RP24',
              'FAC_RP24', 'FAC_WM_CSF', 'FAD_CC_RP24', 'FAD_DiC_RP24', 'RP24_WM_CSF']

statModels = ['FAC_CC_RP6', 'FAD_CC_RP6', 'FAC_DiC_RP6', 'FAD_DiC_RP6', 
              'FAC_CC_RP24', 'FAD_CC_RP24', 'FAC_DiC_RP24', 'FAD_DiC_RP24',
              'FAC_WM_CSF_RP6', 'FAD_WM_CSF_RP24', 'RP24_WM_CSF', 'RP24_CC']

statModels= ['FAD_DiC_RP6', 'FAD_DiC_RP24', 'FAC_WM_CSF_RP6']

#statModels = [ 'FAC_DiC_RP24'         ]

#if os.path.isdir(qcDir) == True:
#    shutil.rmtree(qcDir)
#os.mkdir(qcDir)



  
qcFcResults     = []
bic             = [] 
aic             = [] 
subsMdl         = []
models          = []

for proc in procBase:
    for mdl in statModels:
        
        fcMats    = []
        meanMov   = []
        stdMov    = []
        movFcR    = []
        
        center    = []
        timepoint = []
        subIds    = []
        
        incSubs   = []
        excSubs   = []
        
        currentModel =  proc +  mdl
        print( currentModel )
      
        
        for sub in sorted(os.listdir(baseDir)):
            subDir  = baseDir + '/' + sub + '/' + '/QA_' + currentModel + '_AGG'
          
            fdFile  = baseDir + '/' + sub + '/maximum_disp.1d_delt'
            
            fcFile    = subDir + '/FC_' + mdl + '_local_global_cobra.txt'
            nodeFile  = subDir + '/FC_' + mdl + '_local_global_cobra_node.txt'
            
            if os.path.isfile(fcFile): #os.path.isfile(fdFile)
                
                fd = np.loadtxt(fdFile)
                fc = np.loadtxt(fcFile)
                
                nodes = np.loadtxt(nodeFile)
                
                fcN = np.zeros((600,600)) * np.nan
                
                idx1 = 0
                
                for n1 in nodes[1:]:
                    idx2 = 0;    
                    for n2 in nodes[1:]:
                        fcN[int(n1-1), int(n2-1)] = fc[idx1,idx2]
                        idx2 += 1
                        
                    idx1 += 1
                        
                
                                
                print('Processing: ' + sub)
                print('\tMeanFD   =  {:.03f}'.format(np.mean(fd)))
                print('\tMaxFD    =  {:.03f}'.format(np.max(fd)))
                pctScrub = 100*np.count_nonzero( fd>0.4) / len(fd)
                print('\tpctScrub =  {:.03f}%'.format(pctScrub))
                #%
                
                if np.mean(fd) > 0.4 or np.max(fd) > 5 or pctScrub > 33:
                    excId = 0
                    if sub[0] == 'B':
                        excId = 100 + int(sub[3:])
                    if sub[0] == 'N':
                        excId = 500 + int(sub[3:])
                        
                    if sub[1] == '1':
                        excId += 100
                    
                    if sub[1] == '2':
                        excId += 200
                        
                    if sub[1] == '3':
                        excId += 300
                        
                        
                    excSubs.append( excId )
                    continue
                
                excId = 0
                if sub[0] == 'B':
                    excId = 100 + int(sub[3:])
                if sub[0] == 'N':
                    excId = 500 + int(sub[3:])
                    
                if sub[1] == '1':
                    excId += 100
                
                if sub[1] == '2':
                    excId += 200
                    
                if sub[1] == '3':
                    excId += 300
                    
                    
                incSubs.append( excId )
                
                
                if sub[0] == 'B':
                    center.append(1)
                    subIds.append(100 + int(sub[3:]))
                if sub[0] == 'N':
                    center.append(10)
                    subIds.append(200 + int(sub[3:]))
                    
                    
                if sub[1] == '1':
                    timepoint.append(1)
                    
                if sub[1] == '2':
                    timepoint.append(2)
                    
                if sub[1] == '3':
                    timepoint.append(3)
                
                #tmp = fcN[uNodes,:]
                #tmp = tmp[0,:,uNodes]
                
                # No cerebellum included
                tmp = fcN[0:400,:]
                tmp = tmp[:,0:400]
                
                fcMats.append(np.squeeze(tmp))
                
                meanMov.append(np.mean(fd) )
                stdMov.append( np.std(fd) )
        # END for sub sorted(....)    
        
        
        
        #plt.imshow(np.nanmean(np.array(fcMats),axis=0), vmin=-0.4, vmax=0.4, cmap='jet')
        
        #%
        nSubs   = len(subIds)
        subsMdl.append(nSubs)
        
        
        if nSubs == 0:
            print( 'NO SUBJECTS IN:' + currentModel)
            continue
        else:
            models.append(currentModel)
            print(nSubs)
        
        fcMatsA = np.array(fcMats)
        

        
    
        fcMatsA   = np.array(fcMats)
        nSubs     = fcMatsA.shape[0]
        
        distVecU, uniqueIdx = np.unique(distVec, return_index=True)
        nUniqueNodes        = len(distVecU)-2
           
        
        qcFc = []
        
        fcMatList = []
        nodeStr = []
        clusStr = []
        modul   = []
        
        import bct as bct
        for s in range(nSubs):
            sfc = fcMatsA[s,:,:]
            sfc[np.isnan(sfc)]=0
            nodeStr.append(  (bct.strengths_und( sfc )) )
            clusStr.append( np.mean(bct.clustering_coef_wu(sfc)) )
            modul.append( np.mean(bct.eigenvector_centrality_und(sfc) ) )
            fcVec    = np.reshape(fcMatsA[s,:,:], (1,nNodes*nNodes))
            distVecU = np.array(distVec[0,uniqueIdx[2:]]).reshape(1, nUniqueNodes)
                
            
            r     = spearmanr(fcVec[0,uniqueIdx[2:]], distVecU, axis=1)
            r     = np.arctanh( r.correlation )
            
            qcFc.append(r)
            fcMatList.append(fcVec)
            
            
        
        fcFd    = []
        pfcFd    = []
        meanMovA = np.reshape( np.array(meanMov), [-1,1])
        
        C1T1 = np.where( (np.array(center) + np.array(timepoint)) == 2 )
        C1T2 = np.where( (np.array(center) + np.array(timepoint)) == 3 )
        C1T3 = np.where( (np.array(center) + np.array(timepoint)) == 4 )
        C2T1 = np.where( (np.array(center) + np.array(timepoint)) == 11 )
        C2T2 = np.where( (np.array(center) + np.array(timepoint)) == 12 )
        C2T3 = np.where( (np.array(center) + np.array(timepoint)) == 13 )
        
        from scipy.stats.mstats import pearsonr
        

          
        
        centerFd    = []
        timepointFd = []

        fcFdC1T1    = []
        pfcFdC1T1    = []
        
        fcFdC1T2    = []
        pfcFdC1T2    = []
        
        fcFdC1T3    = []
        pfcFdC1T3    = []
        
        fcFdC2T1    = []
        pfcFdC2T1    = []
        
        fcFdC2T2    = []
        pfcFdC2T2    = []
        
        fcFdC2T3    = []
        pfcFdC2T3    = []
        
        
        
        movC1T1 = np.transpose(meanMovA[C1T1,0])
        movC1T2 = np.transpose(meanMovA[C1T2,0])
        movC1T3 = np.transpose(meanMovA[C1T3,0])
        movC2T1 = np.transpose(meanMovA[C2T1,0])
        movC2T2 = np.transpose(meanMovA[C2T2,0])
        movC2T3 = np.transpose(meanMovA[C2T3,0])
        for n1 in range(nNodes):
            
            for n2 in range(nNodes):
                if n2 > n1:
                    
                    

                    fcVec    = np.reshape(fcMatsA[C1T1,n1,n2], [-1, 1])
                    nas = ~np.isnan(fcVec)
                    r = pearsonr(fcVec[nas],movC1T1[nas])
                    fcFdC1T1.append(r[0])
                    pfcFdC1T1.append(r[1].data)
                    
                    fcVec    = np.reshape(fcMatsA[C1T2,n1,n2], [-1, 1])
                    nas = ~np.isnan(fcVec)
                    r = pearsonr(fcVec[nas],movC1T2[nas])
                    fcFdC1T2.append(r[0])
                    pfcFdC1T2.append(r[1].data)
                    
                    fcVec    = np.reshape(fcMatsA[C1T3,n1,n2], [-1, 1])
                    nas = ~np.isnan(fcVec)
                    r = pearsonr(fcVec[nas],movC1T3[nas])
                    fcFdC1T3.append(r[0])
                    pfcFdC1T3.append(r[1].data)
                     
           
            
                    fcVec    = np.reshape(fcMatsA[C2T1,n1,n2], [-1, 1])
                    nas = ~np.isnan(fcVec)
                    r = pearsonr(fcVec[nas],movC2T1[nas])
                    fcFdC2T1.append(r[0])
                    pfcFdC2T1.append(r[1].data)
                    
                    fcVec    = np.reshape(fcMatsA[C2T2,n1,n2], [-1, 1])
                    nas = ~np.isnan(fcVec)
                    r = pearsonr(fcVec[nas],movC2T2[nas])
                    fcFdC2T2.append(r[0])
                    pfcFdC2T2.append(r[1].data)
                    
                    fcVec    = np.reshape(fcMatsA[C2T3,n1,n2], [-1, 1])
                    nas = ~np.isnan(fcVec)
                    r = pearsonr(fcVec[nas],movC2T3[nas])
                    fcFdC2T3.append(r[0])
                    pfcFdC2T3.append(r[1].data)
                    
        N = (nNodes*nNodes)/2 - nNodes
        fcFdC1T1        = np.array(np.mean(fcFdC1T1, axis=0))
        pfcFdC1T1_05    = 100*np.array(np.sum(np.where(np.array(pfcFdC1T1) < 0.05, 1, 0), axis=0))/N
        pfcFdC1T1_001   = 100*np.array(np.sum(np.where(np.array(pfcFdC1T1) < 0.001, 1, 0), axis=0))/N
        

        fcFdC1T2        = np.array(np.mean(fcFdC1T2, axis=0))
        pfcFdC1T2_05    = 100*np.array(np.sum(np.where(np.array(pfcFdC1T2) < 0.05, 1, 0), axis=0))/N
        pfcFdC1T2_001   = 100*np.array(np.sum(np.where(np.array(pfcFdC1T2) < 0.001, 1, 0), axis=0))/N
        
        
        fcFdC1T3        = np.array(np.mean(fcFdC1T3, axis=0))
        pfcFdC1T3_05    = 100*np.array(np.sum(np.where(np.array(pfcFdC1T3) < 0.05, 1, 0), axis=0))/N
        pfcFdC1T3_001   = 100*np.array(np.sum(np.where(np.array(pfcFdC1T3) < 0.001, 1, 0), axis=0))/N
        
        fcFdC2T1        = np.array(np.mean(fcFdC2T1, axis=0))
        pfcFdC2T1_05    = 100*np.array(np.sum(np.where(np.array(pfcFdC2T1) < 0.05, 1, 0), axis=0))/N
        pfcFdC2T1_001   = 100*np.array(np.sum(np.where(np.array(pfcFdC2T1) < 0.001, 1, 0), axis=0))/N
        

        fcFdC2T2        = np.array(np.mean(fcFdC2T2, axis=0))
        pfcFdC2T2_05    = 100*np.array(np.sum(np.where(np.array(pfcFdC2T2) < 0.05, 1, 0), axis=0))/N
        pfcFdC2T2_001   = 100*np.array(np.sum(np.where(np.array(pfcFdC2T2) < 0.001, 1, 0), axis=0))/N
        
        
        fcFdC2T3        = np.array(np.mean(fcFdC2T3, axis=0))
        pfcFdC2T3_05    = 100*np.array(np.sum(np.where(np.array(pfcFdC2T3) < 0.05, 1, 0), axis=0))/N
        pfcFdC2T3_001   = 100*np.array(np.sum(np.where(np.array(pfcFdC2T3) < 0.001, 1, 0), axis=0))/N
            
        ##%%
        nodeStr = np.array(nodeStr)
        
        colors = ['blue', 'green', 'black']
        ''' 
        for t in [1,2,3]:
            plt.fill_between(np.arange(399), np.squeeze(np.nanmean(nodeStr[np.where(np.array(timepoint)==t),:],axis=1))+
                     np.squeeze(np.nanstd(nodeStr[np.where(np.array(timepoint)==t),:],axis=1))/np.sqrt(55),
                     np.squeeze(np.nanmean(nodeStr[np.where(np.array(timepoint)==t),:],axis=1))-
                             np.squeeze(np.nanstd(nodeStr[np.where(np.array(timepoint)==t),:],axis=1))/np.sqrt(55),
                             facecolor=colors[t-1], # The fill color
                             color=colors[t-1],       # The outline color
                             alpha=0.2)
            
            
            plt.plot(np.squeeze(np.nanmean(nodeStr[np.where(np.array(timepoint)==t),:],axis=1)), color=colors[t-1])
        '''
        
        #%
        outFile = qcDir + '/' + currentModel + '.bz2'
        
        
        
        pd.DataFrame.from_records( {'fcMatsA':    fcMatList, 
                                    'Center':     center,
                                    'Timepoint':  timepoint,
                                    'Subject':    subIds,
                                    'nNodes':     nNodes} ).to_pickle(outFile)
        
        outFile = qcDir + '/' + currentModel + '_FC_FD.bz2'
    
        pd.DataFrame.from_records( {'fcFd_C1T1':    [fcFdC1T1], 
                                    'fcFd_C1T2':    [fcFdC1T2],
                                    'fcFd_C1T3':    [fcFdC1T3], 
                                    'fcFd_C2T1':    [fcFdC2T1], 
                                    'fcFd_C2T2':    [fcFdC2T2], 
                                    'fcFd_C2T3':    [fcFdC2T3], 
                                    'pfcFd_C1T1_05':    [pfcFdC1T1_05], 
                                    'pfcFd_C1T2_05':    [pfcFdC1T2_05],
                                    'pfcFd_C1T3_05':    [pfcFdC1T3_05], 
                                    'pfcFd_C2T1_05':    [pfcFdC2T1_05], 
                                    'pfcFd_C2T2_05':    [pfcFdC2T2_05], 
                                    'pfcFd_C2T3_05':    [pfcFdC2T3_05], 
                                    'pfcFd_C1T1_001':   [pfcFdC1T1_001], 
                                    'pfcFd_C1T2_001':    [pfcFdC1T2_001],
                                    'pfcFd_C1T3_001':    [pfcFdC1T3_001], 
                                    'pfcFd_C2T1_001':    [pfcFdC2T1_001], 
                                    'pfcFd_C2T2_001':    [pfcFdC2T2_001], 
                                    'pfcFd_C2T3_001':    [pfcFdC2T3_001]} 
                                     ).to_pickle(outFile)
    
    
    
        dataTable = pd.DataFrame.from_records( {'QC_FC_' + currentModel: qcFc, 
                                                'Strength_' + currentModel: np.mean(np.array(nodeStr), axis=1), 
                                                'Cluster_' + currentModel: np.array(clusStr),
                                                'EigCen_' + currentModel: np.array(modul),
                                                'Center':np.asarray(center),
                                                'Timepoint':np.asarray(timepoint),
                                                'Mean_Fd':np.asarray(meanMov),
                                                'Std_Fd':np.asarray(stdMov),
                                                'Subject':np.asarray(subIds) } )
        
        outFile = qcDir + '/' + currentModel + '_DM.bz2'
        dataTable.to_pickle(outFile)
        
        md        = smf.mixedlm('QC_FC_' + currentModel + ' ~ Center*Timepoint', dataTable, groups=dataTable["Subject"], missing='drop')
        mdf       = md.fit(reml=True)
        
        md2        = smf.mixedlm('Strength_' + currentModel + ' ~ Center*Timepoint', dataTable, groups=dataTable["Subject"], missing='drop')
        mdf2       = md2.fit(reml=True)
        
        md3        = smf.mixedlm('Cluster_' + currentModel + ' ~ Center*Timepoint', dataTable, groups=dataTable["Subject"], missing='drop')
        mdf3       = md3.fit(reml=True)
        
        #md4        = smf.mixedlm('FD_' + currentModel + ' ~ Center*Timepoint', dataTable, groups=dataTable["Subject"], missing='drop')
        #mdf4       = md4.fit(reml=True)
        
        
        outFile = qcDir + '/' + currentModel + '_LMME.bz2'
        pd.DataFrame({'Model':[mdf], 'Model_Strength':[mdf2],
                      'Model_Cluster':[mdf3]}).to_pickle(outFile)
        
        qcFcResults.append(mdf)
        
        mdf       = md.fit(reml=False)
        bic.append(mdf.bic)
        aic.append(mdf.aic)
        #break
        #END for proc,mdl in zip(proc...)
        
#%%
        
try:
     clear = lambda: os.system('clear')
     clear()
except:
    pass
        


statModels = ['FAC_CC_RP6', 'FAD_CC_RP6', 'FAC_DiC_RP6', 'FAD_DiC_RP6', 
              'FAC_CC_RP24', 'FAD_CC_RP24', 'FAC_DiC_RP24', 'FAD_DiC_RP24',
              'FAC_WM_CSF_RP6', 'FAD_WM_CSF_RP24', 'RP24_WM_CSF', 'RP24_CC']


# Change to number of sign. edges...
idx = 0

ct = []
tp = []
offs = [-1.5,-1.25,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75, 1, 1.25]
x = [0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48]

colors = [ (1,.4,.4), (.8,.2,.2), (.6,0,0), 
           (.4,.4,1), (.2,.2,.8), (0,0,.6)]

fig = plt.figure(figsize=(50,12), dpi=72, facecolor='w', edgecolor='k')
qcDir = '/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/Quality_Control/FinalProc_20200714/Files/'



xticks = []
selModels = []
for mdl in statModels:
    inFile = qcDir + '/' + mdl + '_LMME.bz2'
    if os.path.isfile(inFile):
        #df = pd.read_pickle(inFile)
        #print(df["Model"][0].summary())
        #print(df["Model_Strength"][0].summary())
        #print(df["Model_Cluster"][0].summary())
        print('========================================')
        #print('========================================')
        #print('========================================')
        
    inFile = qcDir + '/' + mdl + '_FC_FD.bz2'
    
    if os.path.isfile(inFile):
        fcFd = []
        fcFd_ = []
        selModels.append(mdl)
        df = pd.read_pickle(inFile)
        
        fcFd_.append(df["pfcFd_C1T1_05"].to_numpy()[0])
        fcFd_.append(df["pfcFd_C1T2_05"].to_numpy()[0])
        fcFd_.append(df["pfcFd_C1T3_05"].to_numpy()[0])
        fcFd_.append(df["pfcFd_C2T1_05"].to_numpy()[0])
        fcFd_.append(df["pfcFd_C2T2_05"].to_numpy()[0])
        fcFd_.append(df["pfcFd_C2T3_05"].to_numpy()[0])
        
        fcFd.append(float(fcFd_[0]))
        fcFd.append(float(fcFd_[1]))
        fcFd.append(float(fcFd_[2]))
        fcFd.append(float(fcFd_[3]))
        fcFd.append(float(fcFd_[4]))
        fcFd.append(float(fcFd_[5]))
        
        x_ = []
        
        for k in range(len(fcFd)):
            x_.append( x[idx] + offs[k])
      
        plt.bar(x_, fcFd, width=0.25, color=colors, edgecolor='k')
        
        xticks.append( x_[2] + (x_[3]-x_[2])/2)
        
        
        idx += 1
            


plt.xticks(ticks=xticks, labels=selModels, fontsize=26, rotation=0)
plt.yticks(fontsize=18)
plt.grid(b=True, color=(.6,.6,.6), alpha=0.3, linestyle='--', linewidth=1)
plt.ylabel("Significant Edges [% Total Edges]", fontsize=24, fontweight="bold")
plt.xlabel("Nuisance Regression Model", fontsize=24, fontweight="bold")

ax = plt.gca()


for axis in ['left', 'bottom']:
    ax.spines[axis].set_linewidth(3)
   
for axis in ['top', 'right']:
    ax.spines[axis].set_linewidth(0)

plt.title("p$_{uncorrected}$ < 0.05", fontsize=28, fontweight="bold")
outFile = "/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/Quality_Control/FinalProc_20200714/01_NumSigEdges_05.png"

plt.savefig(outFile)
plt.close("all")
#%%
   
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable


def plot_lg_matrix(matrix, cmap='RdBu_r', ax=None, colorbar=True, ylabel=True, xlabel=True, vmin=-5, vmax=5, cbarlabel='Z Score', title=None):
    
    if not ax == None:
        plt.sca(ax)

    colors1 = plt.cm.BuPu_r(np.linspace(0., 1, 256))
    colors2 = plt.cm.afmhot_r(np.linspace(0, 1, 256))

    # combine them and build a new colormap
    colors = np.vstack((colors1, colors2))
    mymap  = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
    
    im = plt.imshow( matrix, vmin=vmin, vmax=vmax, cmap=mymap )

    
    netNames = ['L Visual', 'L Somatomotor', 'L Dorsal Attention', 'L Ventral Attention',
                'L Limbic', 'L Control', 'L DMN', 'L Frontoparietal', 
                'R Visual', 'R Somatomotor', 'R Dorsal Attention', 'R Ventral Attention',
                'R Limbic', 'R Control', 'R DMN', 'R Frontoparietal', 
                'L Non-Cortical Regions', 'R Non-Cortical Regions']

    netL1    = [0,23,58,84,107,119,147,193,199,222,257,283,311,323,356,389,399,426,452]
    prev     = 0
    xticks   = []
    for d in netL1:
        gray   = 0.5
        gAlpha = 0.8
        
        xticks.append(prev + (d-prev)/2)
        
        plt.plot([0,452],[d,d], color=(gray,gray,gray,gAlpha), lw=1)
        plt.plot([0,452],[prev,prev], color=(gray,gray,gray,gAlpha), lw=1)
        plt.plot([d,d], [0,452], color=(gray,gray,gray,gAlpha), lw=1)
        plt.plot([prev,prev],[0,452], color=(gray,gray,gray,gAlpha), lw=1)
        prev = d
    
    
    prev = 0
        
    for d in netL1:    
        plt.plot([prev,d],[d,d], 'k', lw=3)
        plt.plot([prev,d],[prev,prev], 'k', lw=3)
        plt.plot([d,d],[prev,d], 'k', lw=3)
        plt.plot([prev,prev],[prev,d], 'k', lw=3)
        prev = d
    
    plt.plot([0,452],[199,199], 'k', lw=3)
    plt.plot([199,199],[0,452], 'k', lw=3)
    
    plt.plot([0,452],[399,399], 'k', lw=3)
    plt.plot([399,399],[0,452], 'k', lw=3)
    if xlabel == True:
        plt.xticks(xticks[1:], netNames, fontsize=14, rotation=90 )
    else:
        plt.xticks([] )
        
        
    if ylabel == True:
        plt.yticks(xticks[1:], netNames, fontsize=14 )
    else:
        plt.yticks([] )
    
    ax = plt.gca()
    
    ax.xaxis.set_ticks_position('bottom')
    
    
    for axis in ['left', 'bottom']:
        ax.spines[axis].set_linewidth(3)
   
    for axis in ['top', 'right', 'left', 'bottom']:
        ax.spines[axis].set_linewidth(0)

    plt.xlim([-1,400])
    plt.ylim([400,-1])    
    
    if not title == None:
        plt.title(title)
    
    if colorbar == True:
        #cbar = plt.colorbar()
        #cbar.ax.set_ylabel(cbarlabel)
        
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.15)

        cbar = plt.colorbar(im, cax=cax)
        cbar.ax.set_ylabel(cbarlabel)


statModels= ['FAD_DiC_RP6', 'FAD_DiC_RP24', 'FAC_WM_CSF_RP6']


for mdl in statModels:
    fig = plt.figure(figsize=(15,15), dpi=120, facecolor='w', edgecolor='k')
    inFile = qcDir + '/' + mdl + '.bz2'
    df = pd.read_pickle(inFile)
    
    fcMatsA = df["fcMatsA"].to_numpy()
    #fcMatsA = np.array(fcMats)
    
    nSubs = len(df["Subject"])
    fcSub = []
    for s in range(nSubs):
        fcSub.append( fcMatsA[s])
    fcMatsA  = np.reshape(fcSub, (nSubs,nNodes,nNodes))
    #plt.imshow( np.squeeze(np.nanmean(fcMatsA,axis=0)), vmin=-.4, vmax=0.4 )
    corrMat =  np.squeeze(np.nanmean(fcMatsA,axis=0))
    plot_lg_matrix(corrMat, cmap='RdBu_r', ax=None, colorbar=True, 
                   ylabel=True, xlabel=True, 
                   vmin=-0.5, vmax=0.5, cbarlabel='r', title=mdl)
    
    
    outFile = "/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/Quality_Control/FinalProc_20200714/02_FC_" + mdl + ".png"

    plt.savefig(outFile)
    plt.close("all")
#%%
    
import seaborn as sb
    
idx =0


colors = [ (0,0,0), (.4,.4,1), (1,1,.2) ]
idx=0
fig = plt.figure(figsize=(20,12), dpi=120, facecolor='w', edgecolor='k')        
for mdl in statModels:
    
    inFile = qcDir + '/' + mdl + '.bz2'
    df = pd.read_pickle(inFile)
    
    fcMatsA = df["fcMatsA"].to_numpy()
    #fcMatsA = np.array(fcMats)
    
    nSubs = len(df["Subject"])
    fcSub = []
    for s in range(nSubs):
        fcSub.append( fcMatsA[s])
    fcMatsA  = np.reshape(fcSub, (nSubs*nNodes*nNodes, 1))
    
    mu = np.nanmean(fcMatsA[np.where(fcMatsA != 0 )], axis=0)
    mud = np.nanmedian(fcMatsA[np.where(fcMatsA != 0 )], axis=0)
    #sd = np.nanstd(fcMatsA[np.where(fcMatsA != 0 )], axis=0)
    
    pct = np.nanpercentile(fcMatsA[np.where(fcMatsA != 0 )], q=[1,99], axis=0)
    
    #h0a, = plt.plot([mu, mu],[0, 2.2], color=colors[0], label='Mean')
    #h0b, = plt.plot([mud, mud],[0, 2.2],linestyle='--', color='r', label='Median')
    
    
    #h1,= plt.plot([pct[0], pct[0]],[0, 2.2], color=(.7, .7, .7), linestyle='-.', label='1% Pctile')
    #plt.plot([pct[1], pct[1]],[0, 2.2], color=(.7, .7, .7), linestyle='-.')
    
    #pct = np.nanpercentile(fcMatsA[np.where(fcMatsA != 0 )], q=[5,95], axis=0)
    
    #h2,= plt.plot([pct[0], pct[0]],[0, 2.2], color=(.45,.45,.45), linestyle=':', label='5% Pctile')
    #plt.plot([pct[1], pct[1]],[0, 2.2], color=(.45,.45,.45), linestyle=':', label='5% Pctile')
    
    sb.kdeplot( fcMatsA[np.where(fcMatsA != 0 )],  label=mdl,
                         color=colors[idx], shade=True)

    idx += 1
    #plt.xticks(ticks=xticks, labels=selModels, fontsize=26, rotation=0)
    plt.yticks(fontsize=18)
    plt.grid(b=True, color=(.6,.6,.6), alpha=0.3, linestyle='--', linewidth=1)
    plt.ylabel("Distribution [Normalized]", fontsize=24, fontweight="bold")
    plt.xlabel("Correlation [Pearson's R]", fontsize=24, fontweight="bold")

    plt.xlim([-1,1])
    ax = plt.gca()
    
    
    for axis in ['left', 'bottom']:
        ax.spines[axis].set_linewidth(3)
       
    for axis in ['top', 'right']:
        ax.spines[axis].set_linewidth(0)
    
    #plt.title("p$_{uncorrected}$ < 0.05", fontsize=28, fontweight="bold")
    
    outFile = "/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/Quality_Control/FinalProc_20200714/02_FC_" + mdl + "_Dist.png"
    
    #idx += 1
    
    #plt.legend(handles=[h0a, h0b, h1,h2], fontsize=24)

plt.legend(fontsize=24)
outFile = "/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/Quality_Control/FinalProc_20200714/03_FC_Dist.png"

plt.savefig(outFile)
    #plt.close("all")
    
    
#%%
#md        = smf.mixedlm('QC_FC_' + currentModel + ' ~ Center*Timepoint', dataTable, groups=dataTable["Subject"], missing='drop')
#mdf       = md.fit(reml=True)
#
#md2        = smf.mixedlm('Strength_' + currentModel + ' ~ Center*Timepoint', dataTable, groups=dataTable["Subject"], missing='drop')
#mdf2       = md2.fit(reml=True)
#md3        = smf.mixedlm('Cluster_' + currentModel + ' ~ Center*Timepoint', dataTable, groups=dataTable["Subject"], missing='drop')
#mdf3       = md3.fit(reml=True)
try:
     clear = lambda: os.system('clear')
     clear()
except:
    pass
vc_formula={'Subject': '0 + C(Subject)'}
for mdl in statModels:
    inFile = qcDir + '/' + mdl + '_DM.bz2'
    df = pd.read_pickle(inFile)
    md        = smf.mixedlm('Strength_' + mdl + ' ~ 1 + Center+Timepoint', df, groups=df["Subject"], missing='drop',
                            vc_formula=vc_formula)
    mdf       = md.fit(reml=True, method='powell')
    print('++++++++++++++++++++++++++++++++++++++++++++++++')
    print(mdl)
    print( mdf.summary() )
    print('++++++++++++++++++++++++++++++++++++++++++++++++')

    
    
#%%
# This is the only one which seems to be correctly working in the sense of
# not altering any effects other than the Center
# Implementation comes from  https://github.com/brentp/combat.py

# The other implementations, for some reason, change the signal from 
# Timepoint effect and sometimes create a signigicant center effect
import combat

mdl = 'FAC_WM_CSF_RP6'
inFile = qcDir + '/' + mdl + '_DM.bz2'
df = pd.read_pickle(inFile)



df2 = np.transpose(combat.combat(data=np.transpose(df), batch=np.transpose(df['Center']), model=None))
meas = 'Strength_'

vc_formula={'Subject': '0 + C(Subject)'}

md        = smf.mixedlm(meas + mdl + ' ~ 1 + Center + Timepoint', df, groups=df["Subject"], missing='drop')
mdf       = md.fit(reml=True, method='cg')

md2        = smf.mixedlm(meas + mdl + ' ~ 1 + Timepoint + Center', df2, groups=df["Subject"], 
                         missing='drop')
mdf2       = md2.fit(reml=True, method='powell')




try:
     clear = lambda: os.system('clear')
     clear()
except:
    pass
    
print("Non-Harmonized")
print(mdf.summary())


print(".")
print(".")
print("Harmonized ")
print(mdf2.summary())


