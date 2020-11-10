#!/usr/bin/env python



import os
import argparse
import nibabel as nib
from deepbrain import Extractor

parser = argparse.ArgumentParser(description='Deepbrain skull stripping')

reqoptions = parser.add_argument_group('Required arguments')
reqoptions.add_argument('-i', '-in', dest="img_path", required=True, help='T1 Image' )
reqoptions.add_argument('-o', '-out', dest="out_path", required=True, help='Brain Mask Path' )

args = parser.parse_args()

img_path  = args.img_path
outFile = args.out_path



# Load a nifti as 3d numpy image [H, W, D]
img = nib.load(img_path)
imgData = img.get_fdata()
ext = Extractor()

# `prob` will be a 3d numpy image containing probability
# of being brain tissue for each of the voxels in `img`

prob = ext.run(imgData)

# mask can be obtained as:
mask = prob > 0.5

mask_img =  nib.Nifti1Image(mask, img.affine, img.header)

nib.save(mask_img, outFile)
