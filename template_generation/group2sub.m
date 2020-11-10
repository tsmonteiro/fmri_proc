function group2sub( uFile, funcRef, fileList, wFile)
matlabbatch = {};
matlabbatch{1}.spm.tools.dartel.crt_iwarped.flowfields = {uFile};
matlabbatch{1}.spm.tools.dartel.crt_iwarped.images = cellstr(fileList);
matlabbatch{1}.spm.tools.dartel.crt_iwarped.K = 6;
matlabbatch{1}.spm.tools.dartel.crt_iwarped.interp = 0;

matlabbatch{2}.spm.spatial.coreg.write.ref = {funcRef};
matlabbatch{2}.spm.spatial.coreg.write.source = {wFile};
matlabbatch{2}.spm.spatial.coreg.write.roptions.interp = 0;
matlabbatch{2}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{2}.spm.spatial.coreg.write.roptions.prefix = 'r';

spm_jobman('run', matlabbatch);
end % function atlas2sub
