function coreg_same_image(fixed, moving, other, nVols)
if ~iscell(other)
coregFiles = cell(nVols,1);
for v = 1:nVols
   coregFiles{v} = cat(2, other, ',', num2str(v));
end


matlabbatch = {};

matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {cat(2, fixed)};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {cat(2, moving)};
if length(coregFiles) > 0
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = cellstr(coregFiles);
else
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
end

matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

else
matlabbatch = {};

matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {cat(2, fixed)};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {cat(2, moving)};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = cellstr(other);

matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
end

spm_jobman('run', matlabbatch);
