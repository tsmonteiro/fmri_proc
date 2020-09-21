function reslice_images(base, imgs)

if size(imgs,1) < size(imgs,2)
    imgs=imgs';
end
matlabbatch={};
matlabbatch{1}.spm.spatial.coreg.write.ref = {base};
matlabbatch{1}.spm.spatial.coreg.write.source = imgs;
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 5;
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
spm_jobman('run', matlabbatch);