function meanfunc2anat(anatFile, funcFile)
matlabbatch = {};
matlabbatch{1}.spm.spatial.coreg.write.ref = {anatFile};
matlabbatch{1}.spm.spatial.coreg.write.source = {funcFile};
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';

spm_jobman('run', matlabbatch);