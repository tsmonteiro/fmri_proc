#!/bin/bash
# Clear terminal screen
#printf "\033c"

# I'd recommend running this only between images which possess similar contrast
# For instance, this aligns among T1 images very well, but it is not robust
# when registering T1 to EPI data
ABIN=/home/luna.kuleuven.be/u0101486/ANTs/bin/ # PATH where ANTs is installed

ABIN=$1
REF=$2
MOVING=$3
OUTPREF=$4



${ABIN}/antsRegistration -d 3 -r [$REF,$MOVING,0] -v 1 \
            -m MI[$REF,$MOVING,1,64] -t translation[0.2] -c [500x50,1.e-8,20] \
            -s 8x4vox -f 3x2 -l 1 -n BSpline \
            -m MI[$REF,$MOVING,1,64,Regular,0.5] -t rigid[0.2] -c [500x50,1.e-8,20] \
            -s 8x4vox -f 3x2 -l 1 -n BSpline \
            -m MI[$REF,$MOVING,1,64,Regular,0.5] -t affine[0.2] -c [500x200x200,1.e-8,10] \
            -s 8x6x4vox -f 4x3x2 -l 1 -n BSpline \
            -m CC[$REF,$MOVING,1,3] -t SyN[0.15,3] -c [50x20,1.e-7,10] \
            -s 2x0vox -f 2x1 -l 1 -n BSpline \
            -o [${OUTPREF},${OUTPREF}.nii]
