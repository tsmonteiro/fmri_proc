function p_reorient(i_file, varargin )
% %%
%  i_file='/home/luna.kuleuven.be/u0101486/workspace/data/Tests/t1_M12_s1_STM0_2.nii';
%   function reorient(i_file)
%
%   i_file: [string]
%
%   bore: 17 Septembre 2015
%       - creation of reorient
%
% Function to automatically reorient a T2 structural image on the template
% in order further proceed with the segmentation/normalisation.
% very useful as segmentation/normalisation is rather sensitive on the
% starting orienation of the image

% changed 15/02/2019 kirstin
% including another step before the reorient, which sets the origin of the
% T1 using http://www.nemotos.net/?p=177 approach

if nargin<1
    i_file=spm_select(inf,'image');
end
tmpl=fullfile(spm('dir'),'canonical\avg152T1.nii');
vg=spm_vol(tmpl);
flags.regtype='rigid';

if ~isempty(varargin)
  i_file_ = i_file;
  for i = 1:varargin{2}
      % this part is written by Fumio Yamashita
      i_file = cat(2, i_file_, sprintf(varargin{1}, i));
      for i=1:size(i_file,1)
          file = deblank(i_file(i,:));
          st.vol = spm_vol(file);
          vs = st.vol.mat\eye(4);
          vs(1:3,4) = (st.vol.dim+1)/2;
          spm_get_space(st.vol.fname,inv(vs));
      end



      for i=1:size(i_file,1)
          f=strtrim(i_file(i,:));

          fTemp = strrep(f, '.nii', '_temp.nii');
          spm_smooth(f,fTemp,[12 12 12]);
          vf=spm_vol(fTemp);

          [M,scal] = spm_affreg(vg,vf,flags);
          M3=M(1:3,1:3);
          [u s v]=svd(M3);
          M3=u*v';
          M(1:3,1:3)=M3;
          N=nifti(f);
          N.mat=M*N.mat;
          create(N);

          delete(fTemp);
      end
  end
else
  % this part is written by Fumio Yamashita

  for i=1:size(i_file,1)
      file = deblank(i_file(i,:));
      st.vol = spm_vol(file);
      vs = st.vol.mat\eye(4);
      vs(1:3,4) = (st.vol.dim+1)/2;
      spm_get_space(st.vol.fname,inv(vs));
  end



  for i=1:size(i_file,1)
      f=strtrim(i_file(i,:));

      fTemp = strrep(f, '.nii', '_temp.nii');
      spm_smooth(f,fTemp,[12 12 12]);
      vf=spm_vol(fTemp);

      [M,scal] = spm_affreg(vg,vf,flags);
      M3=M(1:3,1:3);
      [u s v]=svd(M3);
      M3=u*v';
      M(1:3,1:3)=M3;
      N=nifti(f);
      N.mat=M*N.mat;
      create(N);

      delete(fTemp);
  end
end
