# edit the hard-coded base directory containing the matlab code and the averaged volumes
export PESTICA_DIR=/home/luna.kuleuven.be/u0101486/workspace/fmri_proc/ext/pestica4    # set the one you use pestica4
export SLOMOCO_DIR=$PESTICA_DIR/slomoco
export PESTICA_VOL_DIR=$PESTICA_DIR/template

export MATLAB_AFNI_DIR=$PESTICA_DIR/afni_matlab
export MATLAB_PESTICA_DIR=$PESTICA_DIR/pestica_matlab
export MATLAB_EEGLAB_DIR=$PESTICA_DIR/eeglab

export PATH=$PESTICA_DIR:$PATH
export PATH=$SLOMOCO_DIR:$PATH

# set this to the number of processors/cores you're willing to give to PESTICA at the same time
export PESTICA_MATLAB_POOLSIZE=4
# default code is to take all but one - reserve one core for non-PESTICA work, otherwise comment these lines out and set it above
#a=`cat /proc/cpuinfo  | grep processor | nl | tail -n 1 | gawk '{print $1}'`
#a=`expr $a \- 1`
#export PESTICA_MATLAB_POOLSIZE=$a

# add anything you need to add to the matlab command line here
# this first line seeds the random number generator by a function of the current time
#export MATLABLINE='-nojvm -nosplash -r "c=clock; c=c(3:6); c=round(10*sum(reshape(c(randperm(4)), 2, 2))); normrnd(0,1,c); clear c;"'
# this is if you need a custom matlab license file (use the snippets you need)
#export MATLABLINE="-nojvm -nosplash -c /etc/matlab_license.dat"
#export MATLABLINE="-nojvm -nosplash"
# for recent versions (e.g. R2015a definitely needs this), switch -nojvm to -nodesktop, as you need java to use matlab's Handle Graphics functionality
export MATLABLINE="-nodesktop -nosplash -softwareopengl "
#export MATLABLINE="--eval --no-gui"
