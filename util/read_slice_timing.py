#!/usr/bin/env python

# Import required modules
from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from builtins import str


import argparse

import numpy as np

parser = argparse.ArgumentParser(description='Running aCompCor')

# Required options
reqoptions = parser.add_argument_group('Required arguments')
reqoptions.add_argument('-json', dest="jsonFile", required=True, help='Directory name' )
reqoptions.add_argument('-out', dest="outFname", required=True, help='Repetition Time' )

args = parser.parse_args()

jsonFile = args.jsonFile

readingTiming = False
timingData = []
with open(jsonFile) as fp:
   line = fp.readline()
   cnt = 1
   while line:
       line = str.replace(line, '\"', '')
       entry = str.split( line, ':' )
       entry = entry[0].strip()
       #print(entry)
       
       if readingTiming == True:
           if entry[-2] == ']':
               readingTiming = False
               timingData.append( float(entry[0:-2]) )
           else:
               timingData.append( float(entry[0:-1]) )
    
       
       if entry == 'SliceTiming':
           readingTiming = True

       line = fp.readline()
       cnt += 1

timingData = np.array(timingData)
print(timingData)

np.savetxt( args.outFname, timingData, fmt='%.5f' )
