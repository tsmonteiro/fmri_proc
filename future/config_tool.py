#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 11:01:59 2020

@author: u0101486
"""

import guidata
_app = guidata.qapplication() # not required if a QApplication has already been created

import guidata.dataset.datatypes as dt
import guidata.dataset.dataitems as di

from guidata.dataset.dataitems import (FloatItem, IntItem, BoolItem, ChoiceItem,
                             MultipleChoiceItem, ImageChoiceItem, FilesOpenItem,
                             StringItem, TextItem, ColorItem, FileSaveItem,
                             FileOpenItem, DirectoryItem, FloatArrayItem,
                             DateItem, DateTimeItem)

from guidata.dataset.datatypes import DataSetGroup, BeginGroup, EndGroup, ValueProp

import os

class PathForm(dt.DataSet):
    g1 = BeginGroup("Tools")
    homeDir = os.path.expanduser("~")
    fmriProcDir = DirectoryItem("fMRI Proc", homeDir + '/workspace/fmri_proc/', help='This is a help message')
    antsDir = DirectoryItem("ANTs", homeDir + '/ANTs/bin/')
    c3dDir = DirectoryItem("C3d", homeDir + '/Software/itk/ITK/c3d/bin/')
    fixDir = DirectoryItem("FIX", homeDir + '/Software/fix/')
    dicerDir = DirectoryItem("DiCER", homeDir + '/workspace/fmri_proc/ext/DiCER/')
    _g1 = EndGroup("Tools:")
    
    g2 = BeginGroup("Templates:")
    mniRefDir = FileOpenItem("MNI", homeDir + '/emplates/Template_T1_IXI555_MNI152_GS_brain.nii')
    mniRefLowDir = FileOpenItem("MNI Low Res", homeDir + '/emplates/Template_T1_IXI555_MNI152_GS_brain_3mm.nii')
    _g2 = EndGroup("Templates:")
    
    g3 = BeginGroup("Other:")
    outputDir = StringItem('Out Dir', homeDir + '/workspace/data/CAI_China/tmp/@SUB@')
    finalDir = StringItem('Final Dir', homeDir + '/workspace/data/CAI_China/proc/@SUB@')
    _g3 = EndGroup("Other:")
    #a = di.FloatItem("Parameter #1", default=2.3)
    #b = di.IntItem("Parameter #2", min=0, max=10, default=5)
    #type = di.ChoiceItem("Processing algorithm",
    #                     ("type 1", "type 2", "type 3"))
    
    
doCopy       = ValueProp(True)
doReg       = ValueProp(True)
doFMap      = ValueProp(False)
doIng      = ValueProp(False)


class FuncProcForm(dt.DataSet):

    homeDir = os.path.expanduser("~")
    
    tr   = FloatItem('TR', 3)
    acquisitionType = ChoiceItem("Acquisition Type",
                                 ('Sequential', 'Interleaved', 'Parallel'))
    
    doSliceCorr = BoolItem("Correct For Acquisition Timing?")
    
    despike = ChoiceItem("Voxelwise Despike", 
                          ("Do Not Run", "Before MOCO", "After MOCO", "Last Step"))
    
    g1 = BeginGroup("")
    doCopyBtn = BoolItem("Copy Files",
                      help="If disabled, files will not be copied.",
                      default=True).set_prop("display", store=doCopy)
    
    importScript = FileOpenItem("Import Script", 
                                homeDir + '/workspace/fmri_proc/import_scripts/importdata.sh').set_prop("display", active=doCopy)
    
    
    _g1 = EndGroup("")
    
    g2 = BeginGroup("")
    doRegBtn = BoolItem("Process Func Data",
                      help="If disabled, assumes proc_data_native.nii is already present.",
                      default=True).set_prop("display", store=doReg)
    
    moco = ChoiceItem("Motion Correction Method",
                         ("3dvolreg", "SLOMOCO")).set_prop("display", active=doReg)
    
    
    
    _g2 = EndGroup("")
    
    g3 = BeginGroup("")
    doFMapBtn = BoolItem("Distortion Correction",
                      help="",
                      default=False).set_prop("display", store=doFMap)
    
    
    fMapType = ChoiceItem("Correction Type", 
                          ("RevPhase", "Fieldmap", "Anatomical")).set_prop("display", active=doFMap)
    
    ees  = FloatItem("Effective Echo Spacing [s]", 0.0005, 
                     help="EES is calculated by the following formular TODO").set_prop("display", active=doFMap)
    
        
    
    _g3 = EndGroup("").set_prop("display", active=doReg)
    
    g4 = BeginGroup("")
    doIngBtn = BoolItem("Global Mean Intensity Scaling",
                      help="",
                      default=False).set_prop("display", store=doIng)
    
    ingVal  = FloatItem("Global Intensity", 1000, 
                     help="").set_prop("display", active=doIng)
    

    _g4 = EndGroup("")
    
    
doPrior       = ValueProp(False)


class AnatProcForm(dt.DataSet):

    homeDir = os.path.expanduser("~")
    
    extractSurf = BoolItem("Extract cortical surface?", help="NOTE: This takes quite some time.")
    
    extractBrain = ChoiceItem("Skull Stripping Method",
                                 ('DeepBrain', 'BET'))
    
    extraDen = BoolItem("Perform extra denoise steps?", 
                        help="Removes Gaussian noise and attempts to make intensity/contrast better (see 3dUnifize)")
    
    
    g1 = BeginGroup("")
    doSegBtn = BoolItem("Use segmentation priors?",
                      help="",
                      default=False).set_prop("display", store=doPrior)
        
    gmPrior = FileOpenItem("GM Tissue", default='').set_prop("display", active=doPrior)
    wmPrior = FileOpenItem("WM Tissue", default='').set_prop("display", active=doPrior)
    csfPrior = FileOpenItem("CSF Tissue", default='').set_prop("display", active=doPrior)
    _g1 = EndGroup("")
    
    
    
    
# Add the other tabs... ICA, QA, DiCER.. Possibly add them all in a single nuisance tab
pathForm = PathForm(title='PATH')
funcProcForm = FuncProcForm(title='fMRI Proc')
anatProcForm = AnatProcForm(title='T1 Proc')
g = DataSetGroup( [pathForm, funcProcForm, anatProcForm], title='Configuration File Creation Utility v0.1' )
g.edit()

