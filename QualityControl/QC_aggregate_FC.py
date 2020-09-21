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

import configparser

from scipy.stats import pearsonr, spearmanr
from statsmodels.stats.multitest import multipletests
import scipy.stats as sst
import h5py

import pandas as pd
import pingouin as pg

from scipy.stats import ttest_ind

from sklearn.metrics import mutual_info_score


config = configparser.ConfigParser()
config.read('/home/luna.kuleuven.be/u0101486/workspace/fmri_proc/params.ini')

#PATH to qclib
sys.path.append(config['PATHS']['QCLIB_PATH'])
sys.path.append('/home/luna.kuleuven.be/u0101486/workspace/fmri_proc/ext/nitransforms/nitransforms')


import qclib.group_plots as gp

import shutil

def calc_MI(x, y, bins):
    c_xy = np.histogram2d(x, y, bins)[0]
    mi = mutual_info_score(None, None, contingency=c_xy)
    return mi


def plot_mean_results(mats=None, odir=None, atlas='aal', motion=None, mdl='NONE'):
    mats = np.array(mats)
    
    fname = str.upper(atlas) + '_' + mdl + '_MeanFC.png'
    gp.plot_fc_mat(np.mean(mats,axis=0), outFile=odir+'/' + fname, atlas=atlas, figDpi=100)
    plt.close('all')
    
    nSub = mats.shape[0]
    
    hist = []
    for s in range(nSub):
        hist.append(np.histogram(mats[s,:,:], bins=30) ) 
    
    fig = plt.figure(figsize=(20,20), dpi=100, facecolor='w', edgecolor='k')

    for h in hist:
        plt.plot( h[1][1:], h[0] )


    fname = str.upper(atlas) + '_' + mdl + '_Histogram.png'
    plt.savefig(odir + '/' + fname)
    plt.close('all')
    
    
    if not motion == None:
        motion = np.array(motion)
        nNodes = mats.shape[1]
        motionMat = np.zeros((nNodes,nNodes))
        pvals = []
        for n1 in range(nNodes):
            for n2 in range(n1+1,nNodes):
                fc = mats[:,n1,n2]
                afc = pearsonr(motion, fc)
                motionMat[n1,n2] = afc[0]
                pvals.append(afc[1])
                
        pcorr = multipletests(pvals, 0.05, 'fdr_bh')
        pcorr = pcorr[1]            
           
        idx=0
        for n1 in range(nNodes):
            for n2 in range(n1+1,nNodes):
                if pcorr[idx] < 0.05:
                    motionMat[n2,n1] = motionMat[n2,n1]
                else:
                    motionMat[n2,n1] = 0
                idx += 1
                
        fname = 'Motion_' + str.upper(atlas) + '_' + mdl + '.png'
        gp.plot_fc_mat(motionMat, outFile=odir+'/' + fname, atlas=atlas, figDpi=100)
        plt.close('all')


def plot_group_results(mats=None, odir=None, opref=None, atlas='aal', groups=None,
                       save_group_mats=False, mdl='NONE'):
    mats = np.array(mats)
    nNodes = mats.shape[1]
    
    
    groupEffect = np.zeros( (nNodes, nNodes) )
    pvals = []

    groups = np.array(groups)
    groupsLbl = np.unique( groups )

    for n1 in range(nNodes):
        for n2 in range(n1+1,nNodes):
            fc = mats[:,n1,n2]
            
            tstat, p = ttest_ind( fc[ np.where(groups == groupsLbl[0])  ], fc[ np.where(groups == groupsLbl[-1])  ], equal_var=False )
            
            groupEffect[n1,n2] = np.mean(fc[ np.where(groups == groupsLbl[0])  ]) - np.mean(fc[ np.where(groups == groupsLbl[-1])  ])
            pvals.append(p)
                
    pcorr = multipletests(pvals, 0.05, 'holm')
    pcorr = pcorr[1]
        
       
    idx=0
    for n1 in range(nNodes):
        for n2 in range(n1+1,nNodes):
            if pcorr[idx] < 0.05:
                groupEffect[n2,n1] = groupEffect[n1,n2]
            else:
                groupEffect[n2,n1] = 0
            idx += 1
                                
    fname = odir + '/' + str.upper(atlas) + '_' + opref + '_' + mdl + '.png'
    gp.plot_fc_mat(groupEffect, outFile=fname, atlas=atlas, vmin=-0.4, vmax=0.4, figDpi=100)
    plt.close('all')
    
    if save_group_mats == True:
        matG1 = np.mean( mats[np.where(groups == groupsLbl[0])], axis=0 )
        matG2 = np.mean( mats[np.where(groups == groupsLbl[-1])], axis=0 )
        

        ofileG1 = odir + '/' + str.upper(atlas) + '_' + opref + '_' + mdl + '_G1.png'
        ofileG2 = odir + '/' + str.upper(atlas) + '_' + opref + '_' + mdl + '_G2.png'
        gp.plot_fc_mat(matG1, outFile=ofileG1, atlas=atlas, figDpi=100)
        plt.close('all')
        
        gp.plot_fc_mat(matG2, outFile=ofileG2, atlas=atlas, figDpi=100)
        plt.close('all')


def plot_regression_results(mats=None, ofile=None, atlas='aal', regressor=None):
    print("TODO. Sorry.")


#END OF FUNCTION DEFINITIONS
# ======================================================================
# ======================================================================
    
    
#TODO Make it easier to switch between projects

#project = 'CRUNCH'
#project = 'RepImpact'
project = 'CAI_China'

baseDir       = '/home/luna.kuleuven.be/u0101486/workspace/data/' + project + '/tmp/'
qcDir = '/home/luna.kuleuven.be/u0101486/workspace/data/' + project + '/Quality_Control/FC/'
ages = np.loadtxt('/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/ages.txt', delimiter=',')[:,0]






statModels = ['NONE', 'SFIX_CC', 'SFIX_D', 'SRP24WM1CSF1', 'SRP24CC', 'SRP9', 'SFIX']
#statModels = ['NONE', 'SFIX', 'SFIX_D', 'SRP24WM1CSF1', 'SRP24CC', 'SRP9']
#statModels = ['NONE']
#statModels = ['SFIX_D']

if os.path.isdir(qcDir) == True:
    shutil.rmtree(qcDir)
os.mkdir(qcDir)

if project == 'RepImpact' and not os.path.isdir(qcDir + '/TimePoint'):
    os.mkdir(qcDir + '/TimePoint')

if project == 'RepImpact' and not os.path.isdir(qcDir + '/Center'):
    os.mkdir(qcDir + '/Center')

   


import bct
for mdl in statModels:
    
    aalMats = []
    aalScMats = []
    histAal = []
    histLg4 = []
    lg4Mats = []
    
    strength = []
    
    incAges = []
    incSubs = []
    incMotion = []
    notProc = []
    incSIDs = []
    incGroups = []
    incCenters = []

    
    fcQc = qcDir + '/Mats_' + mdl + '/'
    if os.path.isdir(fcQc):
        shutil.rmtree(fcQc)
    
    os.mkdir(fcQc)
    
    
    for sub in sorted(os.listdir(baseDir)):
        print(sub)
        subDir  = baseDir + '/' + sub
        fcAal   = subDir + '/QA_' + mdl + '/05_FC_' + mdl + '_correlation_matrix_aal.txt' 
        fcLg4   = subDir + '/QA_' + mdl + '/05_FC_' + mdl + '_correlation_matrix_lg400.txt' 
        motFile = subDir + '/maximum_disp.1d_delt'
        
        scMatFile = '/media/u0101486/Seagate Backup Plus Drive/Crunch_SC_AAL2_dirty/R' + sub[2:] + '_muw__connectome_tractogram10MSIFT2_FOD_AVE99_WM_norm_zerodia_sym_aal2.mat'
        
        if os.path.isfile(fcAal) and os.path.isfile(fcLg4): # and os.path.isfile(scMatFile):
            mot     = np.mean(np.loadtxt(motFile) )
            m = np.loadtxt(motFile)
            nSpikes = len(m[np.where(m>0.5)])
            
            if mot > .4 or nSpikes > 50:
                print(sub + ': Mean = {:02f}, Max = {:02f}, (N = {})'.format(mot, np.max(m), nSpikes) )
                shutil.copyfile(fcAal, fcQc + '/AAL_' + sub + '_HM.txt')
                shutil.copyfile(fcLg4, fcQc + '/LocalGlobal_' + sub + '_HM.txt')
                continue
            
            shutil.copyfile(fcAal, fcQc + '/AAL_' + sub + '.txt')
            shutil.copyfile(fcLg4, fcQc + '/LocalGlobal_' + sub + '.txt')
            
            if project == 'CRUNCH':
                i = int(sub[2:])-1
                incAges.append(ages[i])
                f = h5py.File(scMatFile)
                for k, v in f.items():
                    scMat = np.array(v)
                aalScMats.append( (scMat) )
                
                if ages[i] < 30:
                    incGroups.append(0)
                    
                if ages[i] < 62 and ages[i] >= 30:
                    incGroups.append(1)
                    
                if ages[i] >= 62:
                    incGroups.append(2)
                    
            
            
            if project == 'RepImpact':
                sid = int(sub[3:])-1
                
                incGroups.append( int(sub[1]) )
                if sub[0] == 'B':
                    incCenters.append(0)
                    incSIDs.append( sid + 100 )
                if sub[0] == 'N':
                    incCenters.append(1)
                    incSIDs.append( sid + 200 )
                
            if project == 'CAI_China':
                sid = int(sub[3:])
                
                if sub[3] == '0':
                    incGroups.append(0)
                    incSIDs.append( sid + 100 )
                if sub[3] == '1':
                    incGroups.append(1)
                    incSIDs.append( sid + 200 )
            
            aal = np.loadtxt(fcAal)
            lg4 = np.loadtxt(fcLg4)
            incMotion.append( mot )
            
            #randNet = bct.randmio_und_signed(aal, 1)
            #Cr = bct.clustering_coef_wu(randNet[0])
            #C  = bct.clustering_coef_wu(aal)
            #Lr = ( bct.distance_wei(randNet[0])[0] )
            #L = ( bct.distance_wei(aal)[0] )
            #smallWorld = np.nanmean( (C/Cr) / (L/Lr) )
            #aal = bct.threshold_proportional(aal, 0.25)
            
            
            
            nodeStr = bct.strengths_und(aal)
            nodeStr = nodeStr[94:]
            
            #nodeStr = np.sort(nodeStr[0])[::-1] + np.sort(abs(nodeStr[1]))[::-1]
            #nodeStr = bct.eigenvector_centrality_und(lg4)
            #nodeStr = bct.clustering_coef_wu(lg4)
            #histLg4.append(np.sum(nodeStr))
            
            strength.append( np.mean(nodeStr) )
            
            #nodeStr = bct.strengths_und_sign(aal)
            #nodeStr = np.sort(nodeStr[0])[::-1] + np.sort(abs(nodeStr[1]))[::-1]
            #nodeStr = bct.eigenvector_centrality_und(aal)
            #histAal.append(np.histogram(nodeStr, bins=30))
            
            #histLg4.append(np.histogram(lg4, bins=30))
            #histAal.append(np.histogram(aal, bins=30))
            
            #histLg4.append(nodeStr)
            
            aalMats.append( np.arctanh(aal) )
            lg4Mats.append( np.arctanh(lg4) )

      
        else:
            notProc.append(sub)
        
           
    plt.plot(strength)
    # END for sub sorted(....)    
    #%%
    print( "Plotting group comparisons AAL" )
    plot_group_results(mats=aalMats, odir=qcDir, opref='GroupEffect', mdl=mdl,
                       atlas='aal', groups=incGroups, save_group_mats=True)
    
    print( "Plotting group comparisons LocalGlobal" )
    plot_group_results(mats=lg4Mats, odir=qcDir, opref='GroupEffect', mdl=mdl,
                       atlas='localGlobal', groups=incGroups, save_group_mats=True)
    
    
    print( "Plotting mean FC matrices" )
    plot_mean_results(mats=aalMats, odir=qcDir, atlas='aal', motion=incMotion, mdl=mdl)
    plot_mean_results(mats=lg4Mats, odir=qcDir, atlas='localGlobal', motion=incMotion, mdl=mdl)
    '''
    #%%
    import nilearn.plotting as nlp
    from pingouin import partial_corr
    import pandas as pd
    
    caiDuration = [120,36,10,2,4,120,72,24,7,12,7,4,60,8,36,7.5,9,6,192,192,120,60,12]
    kfas        = [52,15,29,65,60,35,12,50,62,55,17,35,57,22,40,60,62,17,25,65,57,57,70]
    aofas       = [73,58,63,79,66,57,65,89,87,89,29,42,80,78,69,79,85,37,39,72,76,95,79]
    age         = [31,40,31,31,27,29,29,24,25,30,13,26,38,30,33,31,31,41,38,37,25,22,43]
    bmi         = [21.25,21.97,26.12,18.69,20.98,28.34,29.69,20.68,20.99,22.35 ,26.95 ,19.38 ,27.14 ,28.74 ,26.53 ,20.03 ,23.03 ,20.83 ,20.90 ,30.02 ,22.49 ,25.71, 24.22 ]
    vasr        = [0.8,0.9,3.5,0,0.6,0,1.7,0.7,0,0,1.1,0,0,3.4,2.6,1.4,0,0.8,0,1.9,1,0,4.5]

    durMat = np.zeros((120,120))
    pvals = []
    lg4Mats = np.array(lg4Mats)
    aalMats = np.array(aalMats)
    for n1 in range(120):
        print(n1)
        for n2 in range(n1+1,120):
            
            conn = aalMats[26:,n1,n2]
            

            
            #df = pd.DataFrame({'FC':conn, 'Duration':caiDuration, 'Age':age, 'BMI':bmi, 'AOFAS':aofas, 'KFAS':kfas, 'VASR':vasr})
            #res = partial_corr(data=df,x='FC',y='KFAS',covar=['Age'])
            
            afc = pearsonr(aofas, conn)
            durMat[n1,n2] = afc[0] #res["r"][0]
            #pvals.append(res["p-val"][0])
            pvals.append(afc[1])
            
            
                
    
    corrP = multipletests(pvals, 0.05, 'fdr_bh')
    corrP = corrP[1]
    
    idx = 0
    for n1 in range(120):
        for n2 in range(n1+1,120):
            if corrP[idx] < 0.05:
                print('{} x {} :  r = {:.03f} (p = {:.03f})'.format(n1,n2,durMat[n1,n2],corrP[idx]))
                durMat[n2,n1] = durMat[n1,n2]
            idx += 1
            
    gp.plot_fc_mat(durMat, outFile=None, figDpi=70, atlas='aal', vmin=-0.75, vmax=0.75, cmap='jet',
                title=None)
    
            
    #%%        
            
    
    import nilearn.plotting as nlp
    pidx=0
    for n1 in range(6):
        for n2 in range(5):
            ax = plt.subplot(6,5,pidx+1)
            nlp.plot_matrix( lg4Mats[pidx], vmax=0.6, vmin=-0.6, axes=ax )
            plt.title(pidx+1)
            pidx += 1
    
    #%%
    strength = np.array(strength)
    incGroups = np.array(incGroups)
    
    mu1 = np.mean(strength[np.where(incGroups==0)])
    sd1 = np.std(strength[np.where(incGroups==0)])
    mu2 = np.mean(strength[np.where(incGroups==1)])
    sd2 = np.std(strength[np.where(incGroups==1)])
    plt.scatter([1,2], [mu1, mu2]  )
    plt.plot( [1,1], [mu1-sd1, mu1+sd1], color='k' ) 
    plt.plot( [2,2], [mu2-sd2, mu2+sd2], color='k' ) 
    
    
    
    #%%
    
    
    meanAal = np.mean( aalMats, axis=0 )
    meanSAal = np.mean( aalScMats, axis=0 )
    meanSAal = meanSAal / np.max(meanSAal)
    aalMats = np.array(aalMats)
    lg4Mats = np.array(lg4Mats)
    
    
    if project == 'RepImpact':
        incGroups = np.array(incGroups)
        incCenters = np.array(incCenters)
        incSIDs = np.array(incSIDs)
        
        #%%
        import pandas as pd
        import pingouin as pg
        
        centerEffectLg4    = np.zeros( (lg4Mats.shape[1],lg4Mats.shape[1]) )
        timePointEffectLg4 = np.zeros( (lg4Mats.shape[1],lg4Mats.shape[1]) )
        pvalsCenterLg4 = []
        pvalsTimePointLg4 = []
        
        
        for n1 in range(lg4Mats.shape[1]):
            for n2 in range(n1+1,lg4Mats.shape[2]):
                fc = lg4Mats[:,n1,n2]
                df = pd.DataFrame({'Center':incCenters,
                                   'TimePoint':incGroups,
                                   'SIDs':incSIDs,
                                   'FC':fc})
        
            
                # Compute the two-way mixed-design ANOVA
                aov = pg.mixed_anova(dv='FC', within='TimePoint', between='Center', subject='SIDs', data=df)
                
                centerEffectLg4[n1,n2] = np.mean(fc[np.where(incCenters==0)]) - np.mean(fc[np.where(incCenters==1)])
                timePointEffectLg4[n1,n2] = aov['F'][1]
                
                pvalsCenterLg4.append(aov['p-unc'][0])
                pvalsTimePointLg4.append(aov['p-unc'][0])
                
                # Pretty printing of ANOVA summary
                #pg.print_table(aov)
                
        corrPvalsCenterLg4 = multipletests(pvalsCenterLg4, 0.05, 'fdr_bh')
        corrPvalsCenterLg4 = corrPvalsCenterLg4[1]
        
        corrPvalsTimePointLg4 = multipletests(pvalsTimePointLg4, 0.05, 'fdr_bh')
        corrPvalsTimePointLg4 = corrPvalsTimePointLg4[1]
        
        idx=0
        for n1 in range(lg4Mats.shape[1]):
            for n2 in range(n1+1,lg4Mats.shape[2]):
                if corrPvalsTimePointLg4[idx] < 0.05:
                    timePointEffectLg4[n2,n1] = timePointEffectLg4[n2,n1]
                    
                if corrPvalsCenterLg4[idx] < 0.05:
                    centerEffectLg4[n2,n1] = centerEffectLg4[n2,n1]
                
        
        #%%
        
        meanAalT1 = np.mean( aalMats[np.where(incGroups==1)], axis=0 )
        meanAalT2 = np.mean( aalMats[np.where(incGroups==2)], axis=0 )
        meanAalT3 = np.mean( aalMats[np.where(incGroups==3)], axis=0 )
        
        meanLg4T1 = np.mean( lg4Mats[np.where(incGroups==1)], axis=0 )
        meanLg4T2 = np.mean( lg4Mats[np.where(incGroups==2)], axis=0 )
        meanLg4T3 = np.mean( lg4Mats[np.where(incGroups==3)], axis=0 )
        
        meanLg4B = np.mean( lg4Mats[np.where(incCenters==0)], axis=0 )
        meanLg4N = np.mean( lg4Mats[np.where(incCenters==1)], axis=0 )

        meanAalB = np.mean( aalMats[np.where(incCenters==0)], axis=0 )
        meanAalN = np.mean( aalMats[np.where(incCenters==1)], axis=0 )
        
        
        gp.plot_fc_mat(meanAalB, outFile=outFileAalC1, atlas='aal', figDpi=120)
        plt.close('all')
        gp.plot_fc_mat(meanAalN, outFile=outFileAalC2, atlas='aal', figDpi=120)
        plt.close('all')
        
        
        gp.plot_fc_mat(meanLg4B, outFile=outFileLg4C1, atlas='localGlobal', figDpi=120)
        plt.close('all')
        gp.plot_fc_mat(meanLg4N, outFile=outFileLg4C2, atlas='localGlobal', figDpi=120)
        plt.close('all')
        
        gp.plot_fc_mat(centerEffectLg4, outFile=outFileLg4CE, atlas='localGlobal', figDpi=120)
        plt.close('all')
        
        gp.plot_fc_mat(meanAalT1, outFile=outFileAalT1, atlas='aal', figDpi=120)
        plt.close('all')
        
        gp.plot_fc_mat(meanAalT2, outFile=outFileAalT2, atlas='aal', figDpi=120)
        plt.close('all')
        
        gp.plot_fc_mat(meanAalT3, outFile=outFileAalT3, atlas='aal', figDpi=120)
        plt.close('all')
        
        
        gp.plot_fc_mat(meanLg4T1, outFile=outFileLg4T1, atlas='localGlobal', figDpi=120)
        plt.close('all')
        
        gp.plot_fc_mat(meanLg4T2, outFile=outFileLg4T2, atlas='localGlobal', figDpi=120)
        plt.close('all')
        
        gp.plot_fc_mat(meanLg4T3, outFile=outFileLg4T3, atlas='localGlobal', figDpi=120)
        plt.close('all')
    
        gp.plot_fc_mat(timePointEffectLg4, outFile=outFileLg4TE, atlas='localGlobal', vmin=0.5, vmax=5, cmap='hot', figDpi=120)
        plt.close('all')
    
    
    sdAal   = np.std( aalMats, axis=0 )
    
    meanLg4 = np.mean( lg4Mats, axis=0 )
    sdLg4   = np.std( lg4Mats, axis=0 )
    
    if project == 'CRUNCH':
        # Corr with age
        alg4 = np.array(lg4Mats)
        
        nNodes4 = len(meanLg4)
        ageMat = np.zeros((nNodes4,nNodes4))
        pvals = []
        for n1 in range(nNodes4):
            for n2 in range(n1+1,nNodes4):
                fc = alg4[:,n1,n2]
                afc = pearsonr(incAges, fc)
                #lr = linregress(incAges, fc)
                #afc = np.corrcoef(incAges, fc)[0][1]
                ageMat[n1,n2] = afc[0]
                pvals.append(afc[1])
                
        
    
        # Corr with age AAL
        aAal = np.array(aalMats)
        
        nNodesA = len(meanAal)
        ageMatA = np.zeros((nNodesA,nNodesA))
        pvalsA = []
        for n1 in range(nNodesA):
            for n2 in range(n1+1,nNodesA):
                fc = aAal[:,n1,n2]
                afc = pearsonr(incAges, fc)
                #lr = linregress(incAges, fc)
                #afc = np.corrcoef(incAges, fc)[0][1]
                ageMatA[n1,n2] = afc[0]
                pvalsA.append(afc[1])
                
        # Corr with age AAL
        sAal = np.array(aalScMats)
        
        nNodesA = len(meanAal)
        ageMatS = np.zeros((nNodesA,nNodesA))
        pvalsS = []
        for n1 in range(nNodesA):
            for n2 in range(n1+1,nNodesA):
                fc = sAal[:,n1,n2]
                afc = pearsonr(incAges, fc)
                #lr = linregress(incAges, fc)
                #afc = np.corrcoef(incAges, fc)[0][1]
                ageMatS[n1,n2] = afc[0]
                pvalsS.append(afc[1])
                
        
        corrPA = multipletests(pvalsA, 0.05, 'fdr_bh')
        corrPA = corrPA[1]
        
        corrPS = multipletests(pvalsS, 0.05, 'fdr_bh')
        corrPS = corrPS[1]
        
        corrP = multipletests(pvals, 0.05, 'fdr_bh')
        corrP = corrP[1]
        
        idx = 0
        for n1 in range(nNodesA):
            for n2 in range(n1+1,nNodesA):
                if corrPA[idx] < 0.05:
                    ageMatA[n2,n1] = ageMatA[n1,n2]
                idx += 1
        
        idx = 0
        for n1 in range(nNodesA):
            for n2 in range(n1+1,nNodesA):
                if corrPS[idx] < 0.05:
                    ageMatS[n2,n1] = ageMatS[n1,n2]
                idx += 1        
                
        idx = 0
        for n1 in range(nNodes4):
            for n2 in range(n1+1,nNodes4):
                if corrP[idx] < 0.05:
                    ageMat[n2,n1] = ageMat[n1,n2]
                idx += 1
        
        
        
        gp.plot_fc_mat(ageMat, outFile=outFileLg4A, atlas='localGlobal', figDpi=72)
        plt.close('all')
        
        gp.plot_fc_mat(ageMatA, outFile=outFileAalA, atlas='aal', figDpi=72)
        plt.close('all')
        
        gp.plot_fc_mat(ageMatS, outFile=outFileAalAS, atlas='aal', figDpi=72)
        plt.close('all')
    
     # Corr with Motion
    alg4 = np.array(lg4Mats)
    
    nNodes4 = len(meanLg4)
    ageMat = np.zeros((nNodes4,nNodes4))
    pvals = []
    for n1 in range(nNodes4):
        for n2 in range(n1+1,nNodes4):
            fc = alg4[:,n1,n2]
            afc = pearsonr(incMotion, fc)
            #lr = linregress(incAges, fc)
            #afc = np.corrcoef(incAges, fc)[0][1]
            ageMat[n1,n2] = afc[0]
            pvals.append(afc[1])
            
    

    # Corr with age AAL
    aAal = np.array(aalMats)
    
    nNodesA = len(meanAal)
    ageMatA = np.zeros((nNodesA,nNodesA))
    pvalsA = []
    for n1 in range(nNodesA):
        for n2 in range(n1+1,nNodesA):
            fc = aAal[:,n1,n2]
            afc = pearsonr(incMotion, fc)
            #lr = linregress(incAges, fc)
            #afc = np.corrcoef(incAges, fc)[0][1]
            ageMatA[n1,n2] = afc[0]
            pvalsA.append(afc[1])
                
        
    corrPA = multipletests(pvalsA, 0.05, 'fdr_bh')
    corrPA = corrPA[1]
    
    corrP = multipletests(pvals, 0.05, 'fdr_bh')
    corrP = corrP[1]
    
    idx = 0
    for n1 in range(nNodesA):
        for n2 in range(n1+1,nNodesA):
            if corrPA[idx] < 0.05:
                ageMatA[n2,n1] = ageMatA[n1,n2]
            idx += 1
            
            
    idx = 0
    for n1 in range(nNodes4):
        for n2 in range(n1+1,nNodes4):
            if corrP[idx] < 0.05:
                ageMat[n2,n1] = ageMat[n1,n2]
            idx += 1
    
    
    
    gp.plot_fc_mat(ageMat, outFile=outFileLg4M, atlas='localGlobal', figDpi=72)
    plt.close('all')
    
    gp.plot_fc_mat(ageMatA, outFile=outFileAalM, atlas='aal', figDpi=72)
    plt.close('all')
    
    #gp.plot_fc_mat(meanSAal, outFile=outFileAalMS, atlas='aal', figDpi=72, vmin=0, vmax=0.1, cmap='hot')
    #plt.close('all')
    
    gp.plot_fc_mat(meanLg4, outFile=outFileLg4, atlas='localGlobal', figDpi=72)
    plt.close('all')
    
    
    #ms = np.reshape(meanSAal, (1,nNodesA*nNodesA))
    #ma = np.reshape(meanAal, (1,nNodesA*nNodesA))
       
    #mDist = np.corrcoef( ma, ms )
    #mDist = mutual_info_score(np.reshape(ms, 8836), np.reshape(ma,8836))
    #gp.plot_fc_mat(meanAal, outFile=outFileAal, atlas='aal', figDpi=72, title='R_structural = {:.03f}'.format(mDist))
    #plt.close('all')
    
    
    
    gp.plot_fc_mat(sdLg4, outFile=outFileLg4S, atlas='localGlobal', vmin=0.1, vmax=0.25, cmap='hot', figDpi=72)
    plt.close('all')

    gp.plot_fc_mat(sdAal, outFile=outFileAalS, atlas='aal', vmin=0.1, vmax=0.25, cmap='hot', figDpi=72)
    plt.close('all')


    #%%
    fig = plt.figure(figsize=(20,20), dpi=72, facecolor='w', edgecolor='k')
    if project == 'CRUNCH':
        normAge = incAges/np.max(incAges)
    idx=0
    colors=[(1,0,0,0.3),(0,1,0,0.3),(0,0,1,0.3)]
    for h in histLg4:
        #plt.scatter( incGroups[idx], h, color=colors[incGroups[idx]-1] )
        if project == 'RepImpact':
            plt.plot( h[1][1:], h[0], color=colors[incGroups[idx]-1] )
        if project == 'CRUNCH':
            plt.plot(h[1][1:], h[0], color=(normAge[idx],normAge[idx],normAge[idx]))
        idx+=1
    
    plt.savefig(outFileLg4H)
    plt.close('all')
    #%%
    fig = plt.figure(figsize=(20,20), dpi=72, facecolor='w', edgecolor='k')
    
    idx=0
    for h in histAal:
        if project == 'RepImpact':
            plt.plot( h[1][1:], h[0], color=colors[incGroups[idx]-1] )
        if project == 'CRUNCH':
            plt.plot(h[1][1:], h[0], color=(normAge[idx],normAge[idx],normAge[idx]))
        idx += 1
    
    plt.savefig(outFileAalH)
    plt.close('all')
    print("Done with " + mdl)
    

    '''