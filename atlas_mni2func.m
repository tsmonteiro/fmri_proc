function atlas_mni2func(anat2funcWarp, atlas, meanFunc)
matlabbatch = {};
matlabbatch{1}.spm.spatial.normalise.write.subj.def = {anat2funcWarp};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {atlas};
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                          78 76 85];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 0;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
matlabbatch{2}.spm.spatial.coreg.write.ref = {meanFunc};
matlabbatch{2}.spm.spatial.coreg.write.source(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{2}.spm.spatial.coreg.write.roptions.interp = 0;
matlabbatch{2}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{2}.spm.spatial.coreg.write.roptions.prefix = 'r';

spm_jobman('run', matlabbatch);


% matlabbatch{1}.spm.spatial.normalise.write.subj.def = {'/home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/tmp/sub001/anat/iy_anat_proc.nii'};
% matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {'/home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/tmp/sub001/anat/aal_mni.nii,1'};
% matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
%                                                           78 76 85];
% matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
% matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 0;
% matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
% matlabbatch{2}.spm.spatial.coreg.write.ref = {'/home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/tmp/sub001/anat/mean_func_data_nds.nii,1'};
% matlabbatch{2}.spm.spatial.coreg.write.source(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
% matlabbatch{2}.spm.spatial.coreg.write.roptions.interp = 0;
% matlabbatch{2}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
% matlabbatch{2}.spm.spatial.coreg.write.roptions.mask = 0;