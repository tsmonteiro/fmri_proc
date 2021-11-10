# -*- coding: utf-8 -*-
"""
v 0.1

Very rough visualizer for relevant info regarding melodic/ICA-FIX
component classification

Most pressing, UI updates are painfully slow at the moment, but it does the job

Created on Tue Oct  8 14:32:45 2019

@author: u0101486
"""
import os 
import matplotlib.pyplot as plt
import matplotlib.widgets as wgt
import numpy as np
#import nilearn.plotting as nlp
import nilearn.signal as sgn
import nibabel as nb
import nilearn.image as nimg
import scipy.signal as ssgn


#=======================================
# Event callbacks
#=======================================

# Called when any of the ortho displays is clicked upon
def onclick(event):
    figSize = fig.get_size_inches()*fig.dpi
    axList = []
    
    for ax in [axX,axY,axZ]:
        axPos = ax.get_position()
        
        x0 = axPos.x0 * figSize[0]
        x1 = axPos.x1 * figSize[0]
        
        y0 = axPos.y0 * figSize[1]
        y1 = axPos.y1 * figSize[1]
        axList.append( (x0,x1,y0,y1) )
    
    
    axIdx = get_ortho_clicked(event.x, event.y, axList)
   
    wasAxes = update_position_clicked(event.x, event.y, axIdx, axList)
    if wasAxes == True:
        update_ortho_display(upd_x=True, upd_y=True, upd_z=True, refresh=True)
    
    
    
# Called when textbox with current IC volume is called
def submit(text):
    global curIc
    global maxIc
    
    global text_box
    
    try:
        newIc = int(text)-1
        
        if newIc < (maxIc-1):
            curIc = newIc
            update_ortho_to_max()
            refresh_view()    
            update_tc_pw()
            text_box.set_val(str(curIc+1))
        
    except ValueError:
        print('[WARNING] Invalid IC number = {}'.format(text))
        
    
# Button to increase current viewed IC
def plus_ic(event):
    global curIc
    global text_box
    
    curIc += 1 
    update_ortho_to_max()
    update_ortho_display(upd_x=True, upd_y=True, upd_z=True, refresh=True)
    update_button_colors()
    update_tc_pw()
    text_box.set_val(str(curIc+1))
    
    
    
    

# Button to decrease increase current viewed IC
def minus_ic(event):
    global curIc
    global text_box
    
    curIc -= 1 
    update_ortho_to_max()
    update_ortho_display(upd_x=True, upd_y=True, upd_z=True, refresh=True)
    update_button_colors()
    update_tc_pw()
    text_box.set_val(str(curIc+1))
       


# Button to switch between smoothed and non-smoothed IC display
def switch_ortho(event):
    global isSmooth
    
    if isSmooth == 1:
        isSmooth = 0
    else:
        isSmooth = 1
        
    refresh_view()
    update_tc_pw()

    
def mark_as_artifact(event):
    global curIc
    global icClass
    
    icClass[curIc] = 1
    update_button_colors()
    
    
def mark_as_not_artifact(event):
    global curIc
    global icClass
    
    icClass[curIc] = 0
    update_button_colors()
    
    
    
def save_classification(event):
    global icClass
    global classFile
    global subDir
    np.savetxt(classFile, icClass)
    
    idx = 1
    classString = '['
    for c in icClass:
        if c == 1:
            classString += str(idx) 
            if idx < len(icClass):
                classString += ','
            else:
                classString += ']'
        
        idx += 1
    
    if classString[-1] == ',':
        classString = classString[:-1]
        classString += ']'
    
    classString = list(classString)
    #classString[-2] = ']'
    classString = "".join(classString).strip()
        
    
    
    text_file = open(subDir + '/FIX/hand_labels_noise.txt', "w+")
    text_file.write(classString)
    # Last line must be an empty string
    text_file.write('\n\n')
    text_file.close()
    
    #text_file = open(subDir + '/FIX/hand_labels_noise.txt', "a")
    #    text_file.close()
    
    print('Classification file saved')
    
        
def brain_scroll(event):

    global axX, axY, axZ
    global curZ
    
    
    figSize = fig.get_size_inches()*fig.dpi
    axList = []
    
    for ax in [axX,axY,axZ]:
        axPos = ax.get_position()
        
        x0 = axPos.x0 * figSize[0]
        x1 = axPos.x1 * figSize[0]
        
        y0 = axPos.y0 * figSize[1]
        y1 = axPos.y1 * figSize[1]
        axList.append( (x0,x1,y0,y1) )
    
    
    #axIdx = get_ortho_clicked(event.x, event.y, axList)
    
    print(event.button)
    if event.button == 'up':
        curZ += 1
    else:
        curZ -= 1
        
    #update_ortho(event.x, event.y, axIdx, axList)
    update_ortho_display(upd_x=True, upd_y=True, upd_z=True)
        
        
    
    
    
    
    
#=======================================
# END OF Event callbacks
#=======================================


def refresh_view(waitfor=True):
    # This is now rather outdated...
    update_ortho_display(upd_x=True, upd_y=True, upd_z=True, refresh=False)
    #plt.gcf().canvas.draw()

def update_tc_pw():
    global axTc
    global axPw
    global tcDat
    global pwDat
    
    global curIc
    global mp
    global tr
    global allTcs
    global secondPass
    
    global icClass
    global icClust
    global spatialCorr
    
    
    tcFile = subDir + '/FIX/filtered_func_data.ica/report/t' + str(curIc+1) + '.txt' 
    
    
    tc =np.loadtxt(tcFile)
    
    
    
    if secondPass == 1:
        
        currTc = allTcs[curIc,:]
        tc = []
        
        
        localIcClass = list(icClass)
        for i in range(len(localIcClass)):
            if i != curIc:
                r = np.corrcoef(currTc, allTcs[i,:])[0,1]
                tc.append( r )
            else:
                tc.append(0)
            
        for mr in range(6):    
            r = np.corrcoef(currTc, np.transpose(mp[:,mr]))[0,1]
            tc.append(r)
            localIcClass.append(2)
            
            
        #localIcClass = np.array(localIcClass)
        
        
        tc1 = allTcs[np.where(np.array(localIcClass)==0)[0],:]
        r = np.corrcoef(currTc, np.mean(tc1,axis=0))[0,1]
        tc.append(r)
        localIcClass.append(3)
        
        
        tc1 = allTcs[np.where(np.array(localIcClass)==1)[0],:]
        r = np.corrcoef(currTc, np.mean(tc1,axis=0))[0,1]
        tc.append(r)
        localIcClass.append(4)
        
        tc = np.array(tc)
     
        
        
        
        axTc.cla()
        #axTc.bar(bc1,h1, linewidth=0.5, color='r', edgecolor='k', width=0.05)
        #print(tc[np.where(np.array(icClass)==1)[0]])
        
        for i in range(len(tc)):
            print('.')
            if localIcClass[i] == 0 and tc[i] != 0:
                axTc.plot(tc[i], 1, 'o', markeredgecolor='k', markerfacecolor='b', markersize=4, lw=0.25, alpha=0.85)  
                
            if localIcClass[i] == 1 and tc[i] != 0:
                axTc.plot(tc[i], 0.8, 'o', markeredgecolor='k', markerfacecolor='r', markersize=4, lw=0.25, alpha=0.85)  
                
            if localIcClass[i] == 2 and tc[i] != 0:
                axTc.plot(tc[i], 0.4, 'o', markeredgecolor='k', markerfacecolor=(0,1,0), markersize=4, lw=0.25, alpha=0.85)  
            
            if localIcClass[i] == 3 and tc[i] != 0:
                axTc.plot(tc[i], 0.1, 'o', markeredgecolor='k', markerfacecolor=(.4,.4,1), markersize=4, lw=0.25, alpha=0.85)  
                
            if localIcClass[i] == 4 and tc[i] != 0:
                axTc.plot(tc[i], 0.1, 'o', markeredgecolor='k', markerfacecolor=(1,0.4,.4), markersize=4, lw=0.25, alpha=0.85)  
        
        axTc.vlines(0, 0, 1.2, color='k', lw=1)
        axTc.vlines(0.25, 0, 1.2, color=(0.5,0.5,0.5), lw=0.3)
        axTc.vlines(-0.25, 0, 1.2, color=(0.5,0.5,0.5), lw=0.3)
        
        axTc.vlines(0.5, 0, 1.2, color=(0.5,0.5,0.5), lw=0.3)
        axTc.vlines(-0.5, 0, 1.2, color=(0.5,0.5,0.5), lw=0.3)
        
        axTc.vlines(0.75, 0, 1.2, color=(0.5,0.5,0.5), lw=0.3)
        axTc.vlines(-0.75, 0, 1.2, color=(0.5,0.5,0.5), lw=0.3)
        
        axTc.vlines(1, 0, 1.2, color=(0.0,0.0,0.0), lw=0.7)
        axTc.vlines(-1, 0, 1.2, color=(0.0,0.0,0.0), lw=0.7)
        
        axTc.set_xlim(-1.5, 1.5)
        axTc.set_ylim(-0.1, 1.2)
    elif secondPass == 2:
        axTc.cla()
        
        if (np.nanmax( spatialCorr[:,curIc]) > 0.05 and icClass[curIc] == 1): 
            axTc.bar([0,1], [np.nanmean( spatialCorr[:,curIc]),np.nanmax( spatialCorr[:,curIc])], color=(.8,0,0))
        elif (np.nanmax( spatialCorr[:,curIc]) < 0.05 and icClass[curIc] == 0): 
            axTc.bar([0,1], [np.nanmean( spatialCorr[:,curIc]),np.nanmax( spatialCorr[:,curIc])], color=(.8,0,0))
        else:
            axTc.bar([0,1], [np.nanmean( spatialCorr[:,curIc]),np.nanmax( spatialCorr[:,curIc])], color=(0.4,1,0.4))
        
        axTc.hlines(0.05, -0.4, 1.5, color=(0.0,0.0,0.0), lw=0.7)
        
        axTc.set_ylim(0, 0.1)
        axTc.set_xlim(-0.5, 1.5)
        
    else:
        tc = allTcs[curIc,:]
        #if tcDat == None:
        axTc.cla()
        if icClust[curIc] == 1:
            tcDat = axTc.plot(tc, linewidth=0.5, color='r')
        elif icClust[curIc] == 2:
            tcDat = axTc.plot(tc, linewidth=0.5, color='y')
        else:
            tcDat = axTc.plot(tc, linewidth=0.5, color='k')
        axTc.set_ylim(-6, 6)
        #else:
        #    tcDat[0].set_ydata(tc)
        


    #rs = []
    #for p in range(6):
    #    rs.append( np.corrcoef(tc,mp[:,p])[0,1] )
    tc =np.loadtxt(tcFile)
    tc = allTcs[curIc,:]
    if pwDat == None:
        f, pw= ssgn.welch(tc, fs=1/tr, scaling='spectrum', nperseg=161)
        pw = pw/np.max(pw)
        axPw.cla()
        pwDat = axPw.plot(f, pw, linewidth=0.5)
        xticks_labs=[]
        xticks = np.linspace(f[0], f[-1], 5)
        for l in xticks:
            xticks_labs.append('{:02}'.format(round(l,2)))
        
        
        axPw.set_xticks(xticks)
        axPw.set_xticklabels(xticks_labs)
        axPw.set_ylim(-.05, 1.2)
        #axPw.set_ylim(0, 75)
        
    else:
        f, pw= ssgn.welch(tc, fs=1/tr, scaling='spectrum', nperseg=161)
        pw = pw/np.max(pw)
        #pwFile = subDir + '/melodic.ic/report/f' + str(curIc+1) + '.txt' 
        #pw =np.loadtxt(pwFile)
        
        pwDat[0].set_ydata(pw)
    
    
    



def update_ortho_to_max():
    global melNii4d
    global melSmNii4d
    
    global isSmooth
    global curIc
    
    global curX, curY, curZ
    global nX, nY, nZ
    
    
    if isSmooth == 0:
        comp = melNii4d[:,:,:,(curIc)]
    else:
        comp = melSmNii4d[:,:,:,(curIc)]
     
        update_ortho_display
    mx = np.argmax(comp)
    curX, curY, curZ = np.unravel_index(mx,(nX,nY,nZ))
    curZ=nZ-curZ

    

def update_ortho_display(upd_x=True, upd_y=True, upd_z=True, refresh=True):
    global curX
    global curY
    global curZ
    
    global axX, axY, axZ
    global imX, imY, imZ
    global ovX, ovY, ovZ
    global novX, novY, novZ
    global lX, lY, lZ
    global nX, nY, nZ
    
    global bgImg3d
    
    global melNii4d
    global melSmNii4d
    
    global isSmooth
    global curIc
    
    global brainBbox
    
    if isSmooth == 0:
        mel = melNii4d
        
                
        vmin=2
        vmax=15
        
        #vmaxN=-1.25
        #vminN=-5
        
    else:
        mel = melSmNii4d
        vmin=0
        vmax=8
        
        #vmaxN=-0.5
        #vminN=-3
    #melIca = mel[:,:,:,curIc]
    #vs = np.quantile(melIca[np.where(melIca>0)], [.9,.95])
    #vmin = vs[0]
    #vmax = vs[1]
    
    vminN = -vmax
    vmaxN = -vmin
    
    
    if upd_x == True:
        maskedIc  = np.ma.masked_where((mel[:,:,nZ-curZ,curIc]) < vmin, mel[:,:,nZ-curZ,curIc])
        maskedIcN = np.ma.masked_where((mel[:,:,nZ-curZ,curIc]*-1) < vmin, (mel[:,:,nZ-curZ,curIc]))
        
        if imX == None:
            axX.clear()
            if isSmooth == 0:
                imX = axX.imshow(bgImg3d[:,:,nZ-curZ], cmap='gray',  aspect=None, interpolation=None, alpha=None, 
                   vmin=0.1, vmax=1.35)
                
                ovX = axX.imshow(maskedIc, cmap='hot',  aspect=None, interpolation=None, alpha=0.8, vmin=vmin, vmax=vmax)
                
                novX = axX.imshow(maskedIcN*-1, cmap='Blues',  aspect=None, interpolation=None, alpha=0.8, vmin=vmin, vmax=vmax)
            else:
                imX = axX.imshow(mel[:,:,nZ-curZ,curIc], cmap='gray',  aspect=None, interpolation=None, alpha=None, 
                   vmin=-6, vmax=6)
        else:
            lines = lX
            for line in lines:
                line.remove()

            if isSmooth == 0:            
                imX.set_data(bgImg3d[:,:,nZ-curZ])
                
                imX.set_clim(vmin=0.1, vmax=1.35)
                
                ovX.set_data(maskedIc)
                ovX.set_clim(vmin=vmin, vmax=vmax)
                
                novX.set_data(maskedIcN)
                novX.set_clim(vmin=vminN, vmax=vmaxN)
                
                ovX.set_alpha(0.8)
                novX.set_alpha(0.8)
                
            else:
                imX.set_data(mel[:,:,nZ-curZ,curIc])
                imX.set_clim(vmin=-5, vmax=5)
                
                ovX.set_alpha(0)
                novX.set_alpha(0)
                
                #ovX.set_clim(vmin=10000, vmax=20000)
                #novX.set_clim(vmin=-20000, vmax=-10000)

        lX = []
        lX.append( axX.vlines(curY, 0,nY-1, color=(1,0,1), linewidth=0.5) )
        lX.append( axX.hlines(curX, 0,nX-1, color=(1,0,1), linewidth=0.5) )
    

    if upd_y == True:
        #ic = np.flip(np.transpose(mel[:,curY,:,curIc]), axis=0)
        ic = np.transpose(mel[:,curY,:,curIc])
        maskedIc =np.ma.masked_where((ic) < vmin, ic)
        maskedIcN =np.ma.masked_where((ic*-1) < vmin, ic)
        if imY == None:
            axY.clear()
            if isSmooth == 0:
                imY = axY.imshow(np.transpose(bgImg3d[:,curY,:]), cmap='gray',  alpha=None,
                           vmin=0.1, vmax=1.35,  aspect='auto' )
    
                
                ovY = axY.imshow(np.flip(maskedIc,axis=1), cmap='hot',  interpolation=None, alpha=0.8,  aspect='auto', vmin=vmin, vmax=vmax)
                novY = axY.imshow(np.flip(maskedIcN,axis=1)*-1, cmap='Blues',  interpolation=None, alpha=0.8,  aspect='auto', vmin=vmin, vmax=vmax)
            else:
                imX = axX.imshow(ic, cmap='gray',  aspect=None, interpolation=None, alpha=None, 
                   vmin=-6, vmax=6)
        else:
            lines = lY
            for line in lines:
                line.remove()
            if isSmooth == 0:            
                imY.set_data(np.flip(np.transpose(bgImg3d[:,curY,:]),axis=1))
                imY.set_clim(vmin=0.1, vmax=1.35)
                ovY.set_data(np.flip(maskedIc,axis=1))
                ovY.set_clim(vmin=vmin, vmax=vmax)
                
                novY.set_data(np.flip(maskedIcN,axis=1))
                novY.set_clim(vmin=vminN, vmax=vmaxN)
                
                ovY.set_alpha(0.8)
                novY.set_alpha(0.8)
            else:
                imY.set_data(np.flip(ic,axis=1))
                imY.set_clim(vmin=-5, vmax=5)
                
                ovY.set_alpha(0)
                novY.set_alpha(0)
        
        
        lY = []
        lY.append( axY.hlines(nZ-curZ, 0,nY-1, color=(1,0,1), linewidth=0.5) )
        lY.append( axY.vlines(nX-curX, 0,nZ-1, color=(1,0,1), linewidth=0.5) )
        
        
    
    if upd_z == True:
        #ic = np.flip(np.transpose(mel[curX,:,:,curIc]), axis=0)
        ic = np.transpose(mel[curX,:,:,curIc])
        maskedIc =np.ma.masked_where((ic) < vmin, ic)
        maskedIcN =np.ma.masked_where((ic*-1) < vmin, ic)
        
        if imZ == None:
            axZ.clear()
                
            if isSmooth == 0:
                imZ = axZ.imshow(np.transpose(bgImg3d[curX,:,:]), cmap='gray',  alpha=None,
                           vmin=0.1, vmax=1.35,  aspect='auto' )
                
                ovZ = axZ.imshow(maskedIc, cmap='hot',  interpolation=None, alpha=0.8,  aspect='auto', vmin=vmin, vmax=vmax)
                novZ = axZ.imshow(maskedIcN, cmap='Blues',  interpolation=None, alpha=0.8,  aspect='auto', vmin=vminN, vmax=vmaxN)
            
        else:
            lines = lZ
            for line in lines:
                line.remove()
            
            if isSmooth == 0:
                imZ.set_data(np.transpose(bgImg3d[curX,:,:]))
                imZ.set_clim(vmin=0.1, vmax=1.35)
                ovZ.set_data(maskedIc)
                ovZ.set_clim(vmin=vmin, vmax=vmax)
                
                novZ.set_data(maskedIcN)
                novZ.set_clim(vmin=vminN, vmax=vmaxN)
                
                ovZ.set_alpha(0.8)
                novZ.set_alpha(0.8)
            else:
                imZ.set_data(ic)
                imZ.set_clim(vmin=-5, vmax=5)
                
                ovZ.set_alpha(0)
                novZ.set_alpha(0)
        
        lZ = []
        lZ.append( axZ.vlines(curY, 0,nZ-1, color=(1,0,1), linewidth=0.5) )
        lZ.append( axZ.hlines(nZ-curZ, 0,nY-1, color=(1,0,1), linewidth=0.5) )
        


    
    if refresh == True:
        plt.gcf().canvas.draw()
        #global brainBbox
        #fig = plt.gcf()
        #background = fig.canvas.copy_from_bbox(brainBbox)
        #fig.canvas.blit(brainBbox)
        
    

def get_ortho_clicked(x, y, rect_list):

    idx = 0
    inIdx = -1
    for rect in rect_list:
        x0 = rect[0]
        x1 = rect[1]
        y0 = rect[2]
        y1 = rect[3]

        inX = x >= x0 and x <= x1
        inY = y >= y0 and y <= y1
        
        if inX and inY:
            inIdx = idx
        
        idx += 1
    return inIdx


def update_position_clicked(x, y, clickedAx, axList):
    global bgNii
    global melNii
    
    global curX
    global curY
    global curZ
    
    global nX, nY, nZ
    
    rect = axList[clickedAx]
    x0 = rect[0]
    x1 = rect[1]
    y0 = rect[2]
    y1 = rect[3]
        
    
    xa = (x - x0)/(x1 - x0)
    ya = (y - y0)/(y1 - y0)


    axClicked = False

    if clickedAx == 0:
        curY = int(np.round(xa*nX))
        curX = int(np.round(nY-ya*nY))
        axClicked = True
        
    if clickedAx == 1:
        curX = int(nX-np.round(xa*nX))
        curZ = int(np.round(ya*nZ))
        axClicked = True
    
    if clickedAx == 2:
        curY = int(np.round(xa*nY))
        curZ = int(np.round(ya*nZ))
        axClicked = True
    
    return axClicked
        
def update_button_colors():
    global wbtn
    global wbta
    global curIc
    global icClass
    global btnBbox

    if icClass[curIc] == 0:
        # Not artifact
        wbtn.color = (.4,1,.4)
        wbtn.hovercolor=(.4,1,.4)
        
        wbta.color = (.7,.7,.7)
        wbta.hovercolor=(.7,1,.7)

        # Color and hovercolor are only called on mouse motion
        # This is a harmless hack to overcome that        
        wbtn.ax.set_facecolor(wbtn.color)
        wbta.ax.set_facecolor(wbta.color)
    else:
        # Artifact
        wbta.color = (.4,1,.4)
        wbta.hovercolor=(.4,1,.4)
        
        wbtn.color = (.7,.7,.7)
        wbtn.hovercolor=(.7,1,.7)
        
        wbtn.ax.set_facecolor(wbtn.color)
        wbta.ax.set_facecolor(wbta.color)
    
    
    background = wbta.canvas.copy_from_bbox(btnBbox)
    wbta.canvas.restore_region(background)
    wbta.canvas.blit(btnBbox)
    

    
    
    
    
    


# ================================================================
# END OF function definitions
# ================================================================
    
    
    
    
    
    
# ================================================================
#
#                          MAIN PROGRAM
#
# ================================================================
    
    
# Globals definition. I am certain there is a better way to handle this, but
# for dev time reasons, I'm going with the quick and dirty approach
global classFile
global subDir

global curIc
global maxIc

global curX
global curY
global curZ


global nX, nY, nZ
global imX, imY, imZ
#global bgNii
#global melNii
#global melSmNii

global bgImg3d
global melNii4d
global melSmNii4d


global isSmooth

global axTc
global axPw
global tcDat
global pwDat

global text_box
global icClass
global icClust

global wbtn
global wbta

global mp

global btnBbox
global brainBbox

global tr
global secondPass
global allTcs

#1 - RS005 OK
#2 - RS006 OK
#3 - RS010 OK
#4 - RS015 OK
#5 - RS016 OK
#6 - RS018 OK
#7 - RS020 OK
#8 - RS023 Ok
#9 - RS026 OK
#10 - RS027 OK
#11 - RS029 OK
#12 - RS032 OK
#13 - RS036 OK
#14 - RS050 OK
#15 - RS061 OK
#16 - RS075 OK
#17 - RS081 OK
#18 - RS086 OK
#19 - RS087 OK
#20 - RS103


#4
#B1_07  OK23
#B1_09  OK23
#B1_18  OK23
#B1_23  EXC [too few components to proper classify]
#B1_24  OK23
#B1_28  OK23
#B1_30  OK23
#B1_55  OK23 [*  Difficult participant]
#B1_63  OK23


#B2_06  OK23
#B2_09  OK23
#B2_11  OK23
#B2_13  OK23
#B2_18  OK23
#B2_24  OK23

#B2_30  OK23 [*]
#B2_53  OK23 [*]
#B2_56  OK23
#B2_61  OK23
#B2_63  EXC [LOTS OF MOV]

#B3_06 OK23
#B3_25 OK23 
#B3_27 OK23 [NOT A GOOD SET]
#B3_30 OK23
#B3_32 OK23
#B3_37 OK23
#B3_48 OK23
#B3_53 OK23
#B3_58 OK23
#B3_60 OK23


#2
# N1_02 OK2
# N1_06 OK2
# N1_08 OK2
# N1_12 OK2
# N1_18 OK2
# N1_21 OK2
# N1_27 OK2
# N1_36 OK2
# N1_52 OK2
# N1_67 OK2

# N2_03 OK2
# N2_05 OK2
# N2_12 OK2
# N2_21 OK2
# N2_25 OK2

# N2_37 OK2
# N2_38 OK2
# N2_40 OK2
# N2_42 OK2
# N2_43 OK2
# N2_62 OK2

# N3_06 OK2
# N3_09 OK2
# N3_12 OK2
# N3_25 OK2
# N3_26 OK2
# N3_32 OK2

# N3_34 OK2
# N3_40 OK2
# N3_44 OK2
# N3_58 OK2
# N3_71 OK2



# G1_01 OK2
# G1_11 OK2
# G1_15 OK2
# G1_16 OK2
# G1_30 OK2


# G2_01 OK2
# G2_12 OK2
# G2_14 OK2
# G2_15 OK2
# G2_17 OK2

# G3_09 OK2
# G3_12 OK2
# G3_14 OK2
# G3_15 OK2
# G3_16 OK2
# G3_17 OK2


#O01 OK
#O04 Ok
#O06 OK
#O10 OK
#O13 OK
#O15 OK
#Y01 OK
#Y02 OK
#Y04 OK
#Y07 OK
#Y10 OK
#Y12 OK
#Y14 OK

# Run the 2 lines below if python complains about non-GUI backend
#import matplotlib
#matplotlib.use('Qt5Agg')
tr = 2.1
#baseDir       = '/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/'
baseDir       = '/media/thiago/EXTRALINUX/fmri_proc/DATA/proc/'
baseDir       = '/media/thiago/Data/Leuven/data/RepImpact/tmp/'
#baseDir       = '/home/luna.kuleuven.be/u0101486/workspace/data/RSPET/tmp/'
#baseDir       = '/home/luna.kuleuven.be/u0101486/workspace/data/ConnectEx/tmp/'
#baseDir       = '/home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/tmp/'
#subDir        = baseDir + '/RS103/'
#subDir        = baseDir + '/Y15/'
#subDir        = baseDir + '/O15/'
subDir        = baseDir + '/N1_72/'
#subDir        = baseDir + '/A001_2/'
maskFile        = subDir + '/nat_mask.nii'
#subDir        = baseDir + '/sub124/'
melodicFile   = subDir + '/FIX/filtered_func_data.ica/melodic_IC.nii.gz'
melodicSmFile = subDir + '/FIX/filtered_func_data.ica/melodic_IC_smooth.nii'
groupIcNativeFile = subDir + '/Group_Comps_Native.nii'
bgImg = subDir + '/anat/' + 'anat_native_fsl_func.nii.gz'




imX = None
imY = None
imZ = None

tcDat = None
pwDat = None

classFile = subDir + '/FIX/filtered_func_data.ica/comp_class.txt' 



motionFile = subDir + '/motion_estimate.par'
try:
    mp = np.loadtxt(motionFile)
except:
    # Motion was estimated using INRIA
    mp = np.loadtxt(motionFile, delimiter=',')
    mp = mp[1:,:]

#
plt.rcParams.update({'font.size': 5, 'font.weight':'bold'} )


isSmooth = 0
curIc = 0
curZ = 20

tcFile = subDir + '/FIX/filtered_func_data.ica/report/t' + str(curIc+1) + '.txt' 
pwFile = subDir + '/FIX/filtered_func_data.ica/report/f' + str(curIc+1) + '.txt' 

pw =np.loadtxt(pwFile)



bgNii = nb.load(bgImg)
nX, nY, nZ = bgNii.shape

melNii = nb.load(melodicFile)
melSmNii = nb.load(melodicSmFile)

global spatialCorr



icSize = melNii.shape
maxIc = icSize[3]

bgNii = nb.load(bgImg)

bgImg3d = bgNii.get_fdata()
bgImg3d = bgImg3d / np.quantile(bgImg3d, .98)
melSmNii4d = melSmNii.get_fdata()
melNii4d = melNii.get_fdata()



if os.path.isfile(classFile):
    icClass = np.loadtxt(classFile)
if os.path.isfile(classFile):# and 1 == 2:
    secondPass = 0
    #if os.path.exists(groupIcNativeFile):
    #    
    #    groupIcNii = nb.load(groupIcNativeFile)
    #    mask  = nb.load(maskFile)
    #    mask  = np.reshape(mask.get_fdata(), (nX*nY*nZ,1))
    #    groupIcData = groupIcNii.get_fdata()
    #    nIcGroup = groupIcNii.shape[3]
    #    
    #    
    #    spatialCorr = np.zeros((nIcGroup, maxIc))
    #    
    #    for ig in range(nIcGroup):
    #        gIc = np.reshape( groupIcData[:,:,:,ig], (nX*nY*nZ,1) )
    #        gIc = gIc / np.nanmax(gIc)
    #        
    #        for isub in range(maxIc):
    #            sIc = np.reshape( melNii4d[:,:,:,isub], (nX*nY*nZ,1) )
    #            sIc = sIc/np.nanmax(sIc)
    #            
    #            spatialCorr[ig, isub] = (np.corrcoef( gIc[np.where(mask==1)], sIc[np.where(mask==1)] )[0,1])**2
    #            #spdist.dice( gIc[np.where(mask==1)].astype(bool), sIc[np.where(mask==1)].astype(bool) )
    #            
    #    #spatialCorr = np.max( spatialCorr, axis=0 )
    #    secondPass = 1
    #else:
    #    spatialCorr = None
    icClass = np.loadtxt(classFile)
    #secondPass += 1
    
    
        
    #allTcs = []
    #for i in range(maxIc):
    #    tc = np.loadtxt(subDir + '/FIX/filtered_func_data.ica/report/t' + str(i+1) + '.txt') 
    #    tc = sgn.clean(np.reshape(tc, [-1,1]), t_r=tr, high_pass=0.01, low_pass=None)
    #    allTcs.append( tc )
    #    
    #   
    #   
    #
    #allTcs = np.array(allTcs)
    #allTcs = allTcs[:,:,0]
else:
    icClass = np.ones((maxIc,1))
    
secondPass = 0
allTcs = []
for i in range(maxIc):
    tc = np.loadtxt(subDir + '/FIX/filtered_func_data.ica/report/t' + str(i+1) + '.txt') 
    tc = sgn.clean(np.reshape(tc, [-1,1]), t_r=tr, high_pass=None, low_pass=None)
    allTcs.append( tc )
    
allTcs = np.array(allTcs)
allTcs = allTcs[:,:,0]


# Try to help a bit classifying
physioFile = subDir + '/physio.txt' 
motionFile = subDir + '/motion_estimate.par' 
csfFile    = subDir + '/anat/gm_tpm_native.nii' 


motion = mp

if os.path.isfile(physioFile):
    motion = np.concatenate((motion, np.loadtxt(physioFile)), 1)



import scipy.stats as sst

csfNii = nb.load(csfFile).get_fdata()

csfNii[np.where( np.isnan( csfNii ))] = 0

icClust = np.zeros((maxIc,1))
for i in range(maxIc):
    r = np.zeros((motion.shape[1], 1))
    for ri in range(motion.shape[1]):
        r_ = np.corrcoef(np.transpose(allTcs[i,:]), motion[:,0] )
        r[ri] = r_[0,1]
    r[np.where((r)>=0.98)] = np.nan
    r = abs(r)
    
    
    v1 = np.nansum( (csfNii > 0.66).astype(int) * (melNii4d[:,:,:,i] > 2.31).astype(int) )
    v2 = np.nansum( (melNii4d[:,:,:,i] > 2.31).astype(int) )
    
    
    #WARNING
    # this is deprecated in favor of median_abs_deviation
    # Need to update
    mad = abs(1-sst.median_absolute_deviation(allTcs[i,:]))
    
    vad = np.var(np.diff(allTcs[i,:]))
    print('{}: {:.3f} / {:.3f} / {:.3f} / {:.3f} / {:.3f}'.format(i,np.nanmax(r),np.nanmean(r**2), mad, np.var(np.diff(allTcs[i,:])), v1/v2))
    
    
    if (np.nanmax(r) > 0.2 or np.nanmean(r**2) > 0.05 or mad > 0.33 ):
        icClust[i] = 1
    elif (v1/v2 < 0.5 or vad > 1):
        icClust[i] = 2
    else:
        icClust[i] = 0

##%%
#rMat = np.corrcoef(allTcs)
#db = KMeans(init='k-means++', n_clusters=2, n_init=10).fit(np.mean(rMat, axis=1).reshape(-1,1))
#db.labels_
##%%

fig = plt.figure(figsize=(4,3), dpi=300, facecolor='w', edgecolor='k')
axs = fig.subplots(ncols=3, nrows=5)

gs = axs[0, 0].get_gridspec()

for ax in axs[0:2, 0]:
    ax.remove()

for ax in axs[0:2, 1]:
    ax.remove()
    
for ax in axs[0:2, 2]:
    ax.remove()
    
for ax in axs[2, 0:]:
    ax.remove()
    
for ax in axs[3, 0:]:
    ax.remove()
axTc = fig.add_subplot(gs[2, 0:])
axMp = fig.add_subplot(gs[3, 0:])
axPw = axs[4,0]

axX = fig.add_subplot(gs[0:2, 0])
axY = fig.add_subplot(gs[0:2, 1]) 
axZ = fig.add_subplot(gs[0:2, 2]) 

curX = 40
curZ = 20
curY = 40

brainBbox = axY.bbox
brainBbox = brainBbox.expanded(4,1)



axX.axes.get_xaxis().set_visible(False) 
axX.axes.get_yaxis().set_visible(False) 

axY.axes.get_xaxis().set_visible(False) 
axY.axes.get_yaxis().set_visible(False) 

axZ.axes.get_xaxis().set_visible(False) 
axZ.axes.get_yaxis().set_visible(False) 


  
fig.canvas.mpl_connect('scroll_event', brain_scroll)


axMp.plot(mp, linewidth=0.5)
axMp.set_xticks([])

cid = fig.canvas.mpl_connect('button_press_event', onclick)


axs[4,2].remove()
axs[4,1].remove()
axTxt = plt.axes([0.5, 0.15, 0.05, 0.1])
axBtp = plt.axes([0.57, 0.15, 0.04, 0.1])
axBtm = plt.axes([0.62, 0.15, 0.04, 0.1])
text_box = wgt.TextBox(axTxt, 'IC [' +  str(maxIc) + ']', initial='1', label_pad=0.2)
wbtp = wgt.Button(axBtp, '+')
wbtm = wgt.Button(axBtm, '-')


axBtn = plt.axes([0.41, 0.02, 0.15, 0.1])
wbtn = wgt.Button(axBtn, 'Not Artifact', color=(.4,1,.4), hovercolor=(.4,1,.4))

axBta = plt.axes([0.6, 0.02, 0.15, 0.1])
wbta = wgt.Button(axBta, 'Artifact', color=(.7,.7,.7), hovercolor=(.7,1,.7))


btnBbox = axBtn.bbox
wdt = btnBbox.width
btnBbox = btnBbox.translated(wdt/2,0).expanded(2.5,1)


axBtv = plt.axes([0.8, 0.02, 0.15, 0.1])
wbtv = wgt.Button(axBtv, 'Save', color=(.7,.7,.7), hovercolor=(.6,.6,.6))
wbtv.on_clicked(save_classification)


axBts = plt.axes([0.7, 0.15, 0.18, 0.1])
wbts = wgt.Button(axBts, 'View Smooth')

wbts.on_clicked(switch_ortho)
wbtp.on_clicked(plus_ic)
wbtm.on_clicked(minus_ic)

wbtn.on_clicked(mark_as_not_artifact)
wbta.on_clicked(mark_as_artifact)

#text_box.on_submit(submit)


#

plt.show()

# ===================================
# Initializing ortho brain display
# ===================================

update_ortho_to_max()
update_button_colors()
refresh_view()
update_ortho_display()


# ===================================
# Loading time series
# ===================================
tc = np.loadtxt(tcFile)

update_tc_pw()

#disp.close()