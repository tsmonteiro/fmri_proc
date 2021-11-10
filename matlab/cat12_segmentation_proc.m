function cat12_segmentation_proc(ANATDIR)
matlabbatch = {};

SPM12DIR = which('spm');
SPM12DIR = SPM12DIR(1:end-5);


addpath('/media/thiago/EXTRALINUX/fmri_proc/matlab/toolbox/spm12/toolbox/DEM')
addpath('/media/thiago/EXTRALINUX/fmri_proc/matlab/toolbox/spm12/matlabbatch/cfg_basicio')
addpath('/media/thiago/EXTRALINUX/fmri_proc/matlab/toolbox/spm12/toolbox/Shoot')
addpath('/media/thiago/EXTRALINUX/fmri_proc/matlab/toolbox/spm12/toolbox/cat12')
%addpath(genpath('/media/thiago/EXTRALINUX/tmp/cat/cat12_r1109/cat12/'))
addpath('/media/thiago/EXTRALINUX/fmri_proc/matlab/toolbox/spm12/toolbox/TSSS')
addpath('/media/thiago/EXTRALINUX/fmri_proc/matlab/toolbox/spm12/toolbox/SRender')
addpath('/media/thiago/EXTRALINUX/fmri_proc/matlab/toolbox/spm12/toolbox/OldSeg')
addpath('/media/thiago/EXTRALINUX/fmri_proc/matlab/toolbox/spm12/toolbox/OldNorm')
addpath('/media/thiago/EXTRALINUX/fmri_proc/matlab/toolbox/spm12/toolbox/Longitudinal')
addpath('/media/thiago/EXTRALINUX/fmri_proc/matlab/toolbox/spm12/toolbox/DAiSS')
addpath('/media/thiago/EXTRALINUX/fmri_proc/matlab/toolbox/spm12/toolbox/FieldMap')
addpath('/media/thiago/EXTRALINUX/fmri_proc/matlab/toolbox/spm12/config')
addpath('/media/thiago/EXTRALINUX/fmri_proc/matlab/toolbox/spm12/toolbox/DARTEL')
addpath('/media/thiago/EXTRALINUX/fmri_proc/matlab/toolbox/spm12/matlabbatch')



matlabbatch{1}.spm.tools.cat.estwrite.data = {cat(2, ANATDIR, '/anat.nii,1')};
matlabbatch{1}.spm.tools.cat.estwrite.nproc = 0;
matlabbatch{1}.spm.tools.cat.estwrite.opts.tpm = {cat(2, SPM12DIR, '/tpm/TPM.nii')};
matlabbatch{1}.spm.tools.cat.estwrite.opts.affreg = 'mni';
matlabbatch{1}.spm.tools.cat.estwrite.opts.biasstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.opts.accstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.APP = 1070;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.LASstr = 0.75;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.gcutstr = 2;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.dartel.darteltpm = {cat(2, SPM12DIR, '/toolbox/cat12/templates_1.50mm/Template_1_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.cat.estwrite.extopts.vox = 1.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.restypes.fixed = [1 0.1];
matlabbatch{1}.spm.tools.cat.estwrite.output.surface = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.neuromorphometrics = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.lpba40 = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.cobra = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.hammers = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.labelnative = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.warped = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.jacobianwarped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.warps = [0 0];


spm_jobman('run', matlabbatch);


labels = load_nii(cat(2, ANATDIR, '/mri/p0anat.nii'));

csf = labels;
csf.img( csf.img >= 1.2 ) = 0;
csf.img(csf.img>0)=1;
save_nii(csf, cat(2, ANATDIR, '/csf_tpm_anat.nii') );


wm = labels;
wm.img( wm.img <= 2.8 ) = 0;
wm.img(wm.img>0)=1;
save_nii(wm, cat(2, ANATDIR, '/wm_tpm_anat.nii') );

gm = labels;
gm.img( gm.img <= 1.5 | gm.img >= 2.5 ) = 0;
gm.img(gm.img>0)=1;
save_nii(gm, cat(2, ANATDIR, '/gm_tpm_anat.nii') );

end
