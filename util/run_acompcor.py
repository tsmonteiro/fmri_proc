#!/usr/bin/env python

# Import required modules
from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from builtins import str
import nipype.algorithms.confounds as npalg
import os
import argparse

parser = argparse.ArgumentParser(description='Running aCompCor')

# Required options
reqoptions = parser.add_argument_group('Required arguments')
reqoptions.add_argument('-d', '-dir', dest="dir", required=True, help='Directory name' )
reqoptions.add_argument('-t', '-tr', dest="tr", required=True, help='Repetition Time' )
reqoptions.add_argument('-i', '-in', dest="inFile", required=True, help='Repetition Time' )

reqoptions.add_argument('-aout',  dest="acomp_out", default='acompcor', required=False, help='Fname acompcor' )
reqoptions.add_argument('-tout',  dest="tcomp_out", default='tcompcor', required=False, help='Fname acompcor' )

reqoptions.add_argument('-n', '-nuis', dest="nuisFile", required=True, help='Nuisance tissue mask [same space as functional, CSFeroded+WMeroded]' )
reqoptions.add_argument('-b', '-bmsk', dest="mskFile", required=True, help='Brain mask [same space as functional]' )

reqoptions.add_argument('-v', '-var', dest="vart", required=False, help='Variance threshold to select components (default 25%)' )

args = parser.parse_args()

sessDir = args.dir
tr = float(args.tr)

vart = args.vart

if not vart:
	vart = 0.25
else:
    vart = float(vart)


inFile   = args.inFile
nuisFile  =  args.nuisFile
mskFile  =  args.mskFile

ofileN = sessDir + '/' + args.acomp_out + '.txt'
ometafileN = sessDir + '/' + args.acomp_out + '_meta.txt'

ofilet  = sessDir + '/' + args.tcomp_out + '.txt'
ometafilet  = sessDir + '/' + args.tcomp_out + '_meta.txt'


ccinterface = npalg.CompCor()
ccinterface.inputs.realigned_file = inFile
ccinterface.inputs.mask_files = mskFile
ccinterface.inputs.variance_threshold = vart
ccinterface.inputs.pre_filter = 'polynomial'
ccinterface.inputs.save_metadata = ometafileN
ccinterface.inputs.regress_poly_degree = 1
ccinterface.inputs.components_file = ofileN
ccinterface.inputs.repetition_time = tr
ccinterface.inputs.high_pass_cutoff = 100

ccinterface.run()




tccinterface = npalg.TCompCor()
tccinterface.inputs.realigned_file = inFile
tccinterface.inputs.percentile_threshold = 0.05
tccinterface.inputs.mask_files = mskFile
tccinterface.inputs.variance_threshold = vart
tccinterface.inputs.pre_filter = 'polynomial'
tccinterface.inputs.regress_poly_degree = 1
tccinterface.inputs.components_file = ofilet
tccinterface.inputs.save_metadata = ometafilet
tccinterface.inputs.repetition_time = tr
tccinterface.inputs.high_pass_cutoff = 100

tccinterface.run()
