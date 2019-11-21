#!/bin/tcsh

set inp = $1

set qq  = ( `3dAttribute DATASET_DIMENSIONS ${inp}+orig` )
echo $#qq
if ( $#qq == 0 ) then
  echo "Dataset ${inp}+orig missing DATASET_DIMENSIONS attribute - exiting"
  exit 1
endif

set nz  = $qq[3]
@ nz1   = $nz - 1

### find orientation of dataset

set qq = ( `3dAttribute ORIENT_SPECIFIC ${inp}+orig` )
if ( $#qq == 0 ) then
  echo "Dataset ${inp}+orig missing ORIENT_SPECIFIC attribute - exiting"
  exit 1
endif

switch ( $qq[1] )
  case "0":
  case "1":
    set xxor = "R"
    breaksw
  case "2":
  case "3":
    set xxor = "A"
    breaksw
  case "4":
  case "5":
    set xxor = "I"
    breaksw
  default:
    echo 'Illegal value in ORIENT_SPECIFIC[1] - exiting'
    exit 1
endsw

switch ( $qq[2] )
  case "0":
  case "1":
    set yyor = "R"
    breaksw
  case "2":
  case "3":
    set yyor = "A"
    breaksw
  case "4":
  case "5":
    set yyor = "I"
    breaksw
  default:
    echo 'Illegal value in ORIENT_SPECIFIC[2] - exiting'
    exit 1
endsw

switch ( $qq[3] )
  case "0":
  case "1":
    set zzor = "R" ; set orient = "sagittal"
    breaksw
  case "2":
  case "3":
    set zzor = "A" ; set orient = "coronal"
    breaksw
  case "4":
  case "5":
    set zzor = "I" ; set orient = "axial"
    breaksw
  default:
    echo 'Illegal value in ORIENT_SPECIFIC[3] - exiting'
    exit 1
endsw

echo "Detected slice orientation: $orient"

switch( $zzor )
  case "R":
    set shift = 1 ; set rota = 4 ; set rotb = 6 ; set scala = 7 ; set shra = 10 ; set shrb = 11
    breaksw

  case "A":
    set shift = 2 ; set rota = 4 ; set rotb = 5 ; set scala = 8 ; set shra = 10 ; set shrb = 12
    breaksw

  case "I":
    set shift = 3 ; set rota = 5 ; set rotb = 6 ; set scala = 9 ; set shra = 11 ; set shrb = 12
    breaksw

  default:
    echo "Illegal value of zzor = ${zzor} - exiting"
    exit 1
endsw

echo "Freezing parameters $shift $rota $rotb $scala $shra $shrb in 3dWarpDrive"

### extract each slice in turn, and register it in 2D only;
### suppressed parameters (-parfix) are
###   #shift = shifting along 3rd dimension
###   #rota  = rotation about 1st inplane axis
###   #rotb  = rotation about 2nd inplane axis
###   #scala = dilation factor along 3rd dimension
###   #shra  = shear factor involving 3rd dimension
###   #shrb  = shear factor involving 3rd dimension

echo "-parfix $shift 0 -parfix $rota   0 -parfix $rotb 0 -parfix $scala  1 -parfix $shra   0 -parfix $shrb 0"
