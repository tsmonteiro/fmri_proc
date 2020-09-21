#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 13:20:23 2019

@author: u0101486
"""
import numpy as np
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt

from matplotlib import cm
from matplotlib import colors
import matplotlib.colors as mclr


import nilearn.plotting as nlp

import matplotlib
matplotlib.use('Agg')

#from adjustText import adjust_text


def nudge_text(textList):
    '''
    Helps with overlapping labels. Mostly useful in the scatter plot in order
    to identify outliers in teh data.
    
    '''
    yLim = plt.gca().get_ylim()
    yRange = yLim[1]-yLim[0]
    
    xLim = plt.gca().get_xlim()
    xRange = xLim[1]-xLim[0]
    
    its      = 100
    dThr     = yRange * 0.03
    dThrX    = xRange * 0.05
    txtShift = yRange * 0.001
    
    oldPositions = []
    newPositions = []
    for text in textList:
        oldPositions.append(text.get_position())
    
    for it in range(its):
        for text in textList:
            pos1 = text.get_position()
            lbl1 = text.get_text()
            for otherText in textList:
                pos2 = otherText.get_position()
                lbl2 = otherText.get_text()    
                
                
                dy = abs( pos1[1]-pos2[1] )
                dx = abs( pos1[0]-pos2[0] )
                
                if dy < dThr and lbl1 != lbl2 and dx < dThrX:
                    #print(lbl1 + ' - ' + lbl2)
                    if pos1[1] > pos2[1]:
                        text.set_position((pos1[0], pos1[1]+txtShift))
                        otherText.set_position((pos2[0], pos2[1]-txtShift))
                    else:
                        text.set_position((pos1[0], pos1[1]-txtShift))
                        otherText.set_position((pos2[0], pos2[1]+txtShift))
                        
                    pos1 = text.get_position()
            
        
        
    for text in textList:
        newPositions.append(text.get_position())
        
        
    for (pos1,pos2) in zip(oldPositions, newPositions):
        if pos1[1] != pos2[1]:
            plt.arrow(pos2[0], pos2[1], -1, pos1[1]-pos2[1], 
                      head_width=None, linestyle=':',aa=True)
            

def scatter_plot(x, y, labels=None, regLine=True, ofile=None, 
                 figDpi=150, figSize=(8,6), xlabel=None, ylabel=None,
                 resOut=False, correctLabels=False):
    
    '''
    Scatter plot with some extra bling.
    Note: Font and marker size are optimized for file output, but I suppose
          it could be made screen-readable by adjusting the figSize and figDpi parameters.

    Parameters
	--------------------------------

	[REQUIRED]
	   x: 1D vector float
	   y: 1D vector float
       
    [OPTIONAL]
    '''
    plt.figure(figsize=figSize, dpi=figDpi, facecolor='w', edgecolor='k')
    plt.rcParams.update({'font.size': 20, 'font.weight':'bold'} )
    
    if x == None or x == []:
        x = np.array(range(len(y)))
    else:
        xPred = np.linspace(np.min(x), np.max(x), len(x)).reshape(-1,1)
        
        linearRegressor = LinearRegression()    # create object for the class
        linearRegressor.fit(x, y)               # perform linear regression
        yPred = linearRegressor.predict(xPred)  # make predictions
        
        plt.plot(xPred, yPred, color='k', linewidth=2)
    
    yz = y
    
    if resOut == True:
        yz = y - linearRegressor.predict(x)
    
    c = np.zeros((len(x), 3))
    i=0
    
    meanY = np.mean(yz)
    stdY  = np.std(yz)
    
    l1 = meanY - 1*stdY
    l2 = meanY - 2*stdY
    l3 = meanY - 3*stdY
    
    u1 = meanY + 1*stdY
    u2 = meanY + 2*stdY
    u3 = meanY + 3*stdY
    
    lvls = [l3,l2,l1,meanY,u1,u2,u3]
    
    normLvls = np.linspace(0,1, len(lvls))
    
    cmap = cm.get_cmap('nipy_spectral')
    texts = []
    for z in yz:
        
        lvl = np.argmin(abs(z - lvls))
                
        if (lvl <= 2 or lvl >= 5) and labels != None:
            texts.append( plt.text( x[i]+1, y[i][0], labels[i], va='center', fontdict={'size':12} ) )
        
        
        rgba = cmap(int( normLvls[lvl]*255 ))
        c[i,:] = rgba[0:3]
            
        i +=1
    
    c = np.array(c)
    plt.scatter(x, y, 80, c=c, edgecolor='k')

    ax = plt.gca()
    
    for axis in ['left', 'bottom']:
    	ax.spines[axis].set_linewidth(3)
    
    for axis in ['top', 'right']:
    	ax.spines[axis].set_linewidth(0)
        
    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)


    norm = mclr.BoundaryNorm([-0.5,normLvls[0],normLvls[1],normLvls[2],
                                    normLvls[4],normLvls[5],normLvls[6],1.5], cmap.N)
    cbar = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ticks=normLvls,
                        boundaries=[-0.5,normLvls[0],normLvls[1],normLvls[2],
                                    normLvls[4],normLvls[5],normLvls[6],1.5],
                                    extend='both')
    
    cbar.ax.set_yticklabels( ['<= Mu - 3SD','Mu - 2SD','Mu - 1SD',
                              'Mu','Mu + 1SD','Mu + 2SD',
                              '>= Mu + 3SD'], fontdict={'size':16} )
    
    if correctLabels == True:
        nudge_text(texts)

    
    if ofile != None:
        plt.savefig(ofile)
    
    
    
def create_aal_labels(aalLabels='/home/luna.kuleuven.be/u0101486/workspace/fmri_proc/atlas/aal/aal2/aal_labels.txt'):
    labsAal = []
    with open(aalLabels) as fp:
       line = fp.readline()
       cnt = 1
       while line:
           lblLine = str.split(line)
           if cnt%2 == 1:
               l = str.split(lblLine[0], '_')
               labsAal.append(l[0])
           else:
               labsAal.append('')
           
           line = fp.readline()
           cnt += 1
    
    return labsAal

def create_lg_labels():
    lgNetLim = [12,24,43,59,72,85,100,108,113,120,133,143,148,166,187,194,200,
                212,223,243,258,272,284,303,312,318,324,334,350,357,373,384,390,400]
    
    
    netNames = ['L Visual Central', 'L Visual Lat', 'L Motor A', 'L Motor B', 
                'L DAN A', 'L DAN B', 'L SAL A', 'L SAL B', 'L Limbic B', 'L Limbic A',
                'L Control A', 'L Control B', 'L Control C',
                'L DMN A', 'L DMN B', 'L DMN C', 'L TPN',
                'R Visual Central', 'R Visual Lat', 'R Motor A', 'R Motor B', 
                'R DAN A', 'R DAN B', 'R SAL A', 'R SAL B', 'R Limbic B', 'R Limbic A',
                'R Control A', 'R Control B', 'R Control C',
                'R DMN A', 'R DMN B', 'R DMN C', 'R TPN'       ]
    
    labsLg = ['']*400
    prevLim = 0
    i=0
    for lim in lgNetLim:
        #plt.plot( [0, 400], [prevLim,prevLim], color='k' )
        #plt.plot( [0, 400], [lim,lim], color='k' )
        idx = prevLim + round((lim - prevLim)/2) - 5
        prevLim = lim
        labsLg[idx] = netNames[i]
        
        i += 1
        
    return labsLg,lgNetLim

def plot_fc_mat(corrMat, outFile=None, figDpi=300, atlas='localGlobal', vmin=-0.5, vmax=0.5, cmap='jet',
                title=None):
    plt.rcParams.update({'font.size': 20, 'font.weight':'bold'} )
    if atlas == 'aal':
        labs          = create_aal_labels()
        showGrid      = True
        
    if atlas == 'localGlobal':    
        labs,lgNetLim = create_lg_labels()
        showGrid      = False
    

    
    fig = plt.figure(figsize=(20,20), dpi=figDpi, facecolor='w', edgecolor='k')
    
    
    
    nlp.plot_matrix(corrMat, labels=labs, vmin=vmin, vmax=vmax, grid=showGrid, 
                    auto_fit=True, colorbar=True, figure=fig, cmap=cmap)
    
    if not title == None:
        plt.title(title)
    
    if atlas == 'localGlobal':
        for lim in lgNetLim:
            plt.plot( [0, 400], [lim-.5,lim-.5], color='k' )
            plt.plot( [lim-.5,lim-.5],[0, 400], color='k' )
    
    if outFile != None:
        plt.savefig(outFile)

    
    
def histogram_2d(X,Y):
    subDvars = dvars[sidx][0]
    subFd = fd[sidx]
    h, x, y, p = ax.hist2d(X, Y, bins=30, cmap='inferno', range=[[0,2],[0,2.5]],
                  cmax=len(X)/2)
    
    
    plt.cla()
    ax.imshow(np.transpose(h), origin='lower', interpolation='gaussian', extent=[0,800,0,800], 
               vmin=0, vmax=10, cmap='inferno', norm=colors.PowerNorm(0.5))
    
    
    
          
    ax = plt.gca()
    yticksp  = np.linspace(0,800,6)
    xticksp  = np.linspace(0,800,5)
    
    yticks = np.linspace(0,2.5,6)
    xticks = np.linspace(0,2,5)
    
    
    ax.set_xticks(xticksp)
    ax.set_yticks(yticksp)
    
    ax.set_xticklabels(xticks, fontdict={'size':9, 'weight':'normal'})
    ax.set_yticklabels(yticks, fontdict={'size':9, 'weight':'normal'})
    
    ax.set_title(incSubs[sidx], fontdict={'size':9, 'weight':'bold'})