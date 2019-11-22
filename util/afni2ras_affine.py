import numpy as np
import argparse
import nibabel as nib


parser = argparse.ArgumentParser(description='Convert AFNI to RAS')

reqoptions = parser.add_argument_group('Required arguments')
reqoptions.add_argument('-d', '-dir', dest="base_dir", required=True, help='Dir' )

args = parser.parse_args()



baseDir = args.base_dir #'/mnt/hgfs/ssd_tmp/ASL/056/'

afni_vec = np.loadtxt(baseDir + '/casl2anat_afni.1D')

ras_vec = np.zeros((4,4))
ras_vec[0,0] = afni_vec[0]
ras_vec[0,1] = afni_vec[1]
ras_vec[0,2] = -afni_vec[2]
ras_vec[0,3] = -afni_vec[3]
ras_vec[1,0] = afni_vec[4]
ras_vec[1,1] = afni_vec[5]
ras_vec[1,2] = -afni_vec[6]
ras_vec[1,3] = -afni_vec[7]
ras_vec[2,0] = -afni_vec[8]
ras_vec[2,1] = -afni_vec[9]
ras_vec[2,2] = afni_vec[10]
ras_vec[2,3] = afni_vec[11]
ras_vec[3,0] = 0
ras_vec[3,1] = 0
ras_vec[3,2] = 0
ras_vec[3,3] = 1

print(afni_vec)
print(ras_vec)

np.savetxt(baseDir + '/casl2anat.mat', ras_vec, fmt='%0.10f')

