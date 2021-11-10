function reslice_images(base, imgs, varargin)

if size(imgs,1) < size(imgs,2)
    imgs=imgs';
end

inter = 5;
isSpm = 1;

if ~isempty(varargin)
  inter = varargin{1};

  if length(varargin) > 1
      isSpm = varargin{2};
  end
end


if isSpm == 1
  matlabbatch={};
  matlabbatch{1}.spm.spatial.coreg.write.ref = {base};
  matlabbatch{1}.spm.spatial.coreg.write.source = imgs;
  matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = inter;
  matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
  matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
  matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
  spm_jobman('run', matlabbatch);
else
  voxSize = [];

  if length(varargin) > 2 && ~isempty(varargin{3})
      voxSize = varargin{3};
  end

  for im = 1:length(imgs)
      reslice_nii(  imgs{im}, strrep( imgs{im}, '.nii', '_r.nii'), voxSize, 0, [], 1);
  end
end
