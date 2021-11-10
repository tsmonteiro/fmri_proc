import os
import argparse
import nibabel as nib
import re
import numpy as np

parser = argparse.ArgumentParser(description='Deepbrain skull stripping')

reqoptions = parser.add_argument_group('Required arguments')
reqoptions.add_argument('-i', '-in', dest="inDir", required=True, help='T1 Image' )
reqoptions.add_argument('-b', '-block', dest="block", required=True, help='Desired block' )
reqoptions.add_argument('-p', '-pattern', dest="pattern", required=True, help='Filename pattern' )
reqoptions.add_argument('-t', '-type', dest="type", required=True, help='Filename ending' )
reqoptions.add_argument('-o', '-out', dest="out_path", required=True, help='Brain Mask Path' )

args = parser.parse_args()

inDir  = args.inDir
outFile = args.out_path

pattern = re.compile(args.pattern)

# Get order of files
blocks = []
fnames = []
for root,  dir, files in os.walk(inDir):
	for fname in files:

		if pattern.search(fname) and fname.endswith(args.type):
			fname = ''.join([str(el) for el in fname])

			fnames.append(fname)
			blockNumber  = str.split(fname, '.')[0]
			blockNumber  = str.split(blockNumber, '_')[-2]

			blocks.append(int(blockNumber))



blocksOrd = np.argsort(blocks)

#print('Will import block {}({}) - [{}]'.format(args.block, blocksOrd[int(args.block)-1], os.path.join(inDir, fnames[blocksOrd[int(args.block)-1]])) )


# Load a nifti as 3d numpy image [H, W, D]
img = nib.load(os.path.join(inDir, fnames[blocksOrd[int(args.block)-1]]))
nib.save(img, outFile)
