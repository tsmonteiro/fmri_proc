# rs_proc_fmri
Scripts to process and perform QC on resting-state data


# Required software

GNU Parallel [http://www.gnu.org/software/parallel]

Python (Installed as part of Anaconda)
  plus a number of modules: numpy, nibabel, nilearn

AFNI

Matlab

FSL v5
  It works so far with FSL6, though filtershift will not compile

ICA-FIX


ANTS

DeepBrain Extractor (included from https://pypi.org/project/deepbrain/   https://github.com/iitzco/deepbrain]
  As of Ubuntu 18, some changes are necessary for it to work (compatibility with tensorflow module v2)


ITK (following the instructions from here https://itk.org/Wiki/ITK/Getting_Started/Build/Linux)
  plus the Convert3D tool (https://sourceforge.net/projects/c3d/files/c3d/Nightly/)
