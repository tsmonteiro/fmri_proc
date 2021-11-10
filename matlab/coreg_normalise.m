function coreg_normalise( refNative, funcNative, nVols, anat, normalise, templates, smoothFwhm, varargin )
rmpath(genpath('/home/luna.kuleuven.be/u0101486/workspace/matlab/toolbox/BrainNetViewer_20191031/'))

    %%
%      baseDir = '/home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/tmp/sub005/anat/';
%      funcNative='func_data.nii';
%       anat = 'anat_proc_brain.nii';
%       refNative = 'mean_func_data_nd.nii';
%       nVols = 200;
%       funcNative = '/home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/tmp/sub002/proc_data_native_fix.nii'
%%

SPM12DIR = which('spm');
SPM12DIR = SPM12DIR(1:end-5);

if iscell(funcNative)
  coregFiles = funcNative';
else
  if length(nVols) > 1
    coregFiles = {};
    for i = 1:length(nVols)
      for v = 1:nVols{i}
         coregFiles{end+1} = cat(2, funcNative{i}, ',', num2str(v));
      end
    end
  else
    coregFiles = cell(nVols,1);
    for v = 1:nVols
       coregFiles{v} = cat(2, funcNative, ',', num2str(v));
    end
  end

end

%smoothFwhm=[5 5 5];
interpMethod = 1;

if ~isempty(varargin)
  interpMethod = varargin{1};
end

[anatPath, anatName, anatExt] = fileparts(anat);
matlabbatch={};
step = 1;
if normalise(1) == 1

    matlabbatch{step}.spm.spatial.coreg.estimate.ref = {cat(2, anat)};
    matlabbatch{step}.spm.spatial.coreg.estimate.source = {cat(2, refNative)};
    if length(coregFiles) > 0
      matlabbatch{step}.spm.spatial.coreg.estimate.other = cellstr(coregFiles);
    else
      matlabbatch{step}.spm.spatial.coreg.estimate.other = {''};
    end
    % or ecc
    matlabbatch{step}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{step}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{step}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{step}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

    step = step + 1;
end

if normalise(1) == 2

    matlabbatch{step}.spm.spatial.coreg.estwrite.ref = {cat(2, anat)};
    matlabbatch{step}.spm.spatial.coreg.estwrite.source = {cat(2, refNative)};
    if length(coregFiles) > 0
      matlabbatch{step}.spm.spatial.coreg.estwrite.other = cellstr(coregFiles);
    else
      matlabbatch{step}.spm.spatial.coreg.estwrite.other = {''};
    end
    % or ecc
    matlabbatch{step}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{step}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{step}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{step}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];

    matlabbatch{step}.spm.spatial.coreg.estwrite.roptions.interp = 5;
    matlabbatch{step}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{step}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{step}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

    step = step + 1;
end



if normalise(2) == 1
    matlabbatch{step}.spm.spatial.preproc.channel.vols = {cat(2, anat)};
    matlabbatch{step}.spm.spatial.preproc.channel.biasreg = 0.001;
    matlabbatch{step}.spm.spatial.preproc.channel.biasfwhm = Inf;
    matlabbatch{step}.spm.spatial.preproc.channel.write = [0 0];
    matlabbatch{step}.spm.spatial.preproc.tissue(1).tpm = {cat(2, SPM12DIR, '/TPM.nii,1')};
    matlabbatch{step}.spm.spatial.preproc.tissue(1).ngaus = 1;
    matlabbatch{step}.spm.spatial.preproc.tissue(1).native = [1 1];
    matlabbatch{step}.spm.spatial.preproc.tissue(1).warped = [0 0];
    matlabbatch{step}.spm.spatial.preproc.tissue(2).tpm = {cat(2, SPM12DIR, '/TPM.nii,2')};
    matlabbatch{step}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch{step}.spm.spatial.preproc.tissue(2).native = [1 1];
    matlabbatch{step}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{step}.spm.spatial.preproc.tissue(3).tpm = {cat(2, SPM12DIR, '/TPM.nii,3')};
    matlabbatch{step}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{step}.spm.spatial.preproc.tissue(3).native = [1 1];
    matlabbatch{step}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{step}.spm.spatial.preproc.tissue(4).tpm = {cat(2, SPM12DIR, '/TPM.nii,4')};
    matlabbatch{step}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{step}.spm.spatial.preproc.tissue(4).native = [1 0];
    matlabbatch{step}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{step}.spm.spatial.preproc.tissue(5).tpm = {cat(2, SPM12DIR, '/TPM.nii,5')};
    matlabbatch{step}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{step}.spm.spatial.preproc.tissue(5).native = [1 0];
    matlabbatch{step}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{step}.spm.spatial.preproc.tissue(6).tpm = {cat(2, SPM12DIR, '/TPM.nii,6')};
    matlabbatch{step}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{step}.spm.spatial.preproc.tissue(6).native = [1 0];
    matlabbatch{step}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{step}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{step}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{step}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{step}.spm.spatial.preproc.warp.affreg = 'mni';

    matlabbatch{step}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{step}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{step}.spm.spatial.preproc.warp.write = [1 1];
    step = step + 1;

end

if normalise(3) == 1
    if ~exist(cat(2, anatPath, filesep, 'u_rc1', anatName, anatExt), 'file')

    matlabbatch{step}.spm.tools.dartel.warp1.images = {
                                                    {cat(2, anatPath, filesep, 'rc1', anatName, '_dartel', anatExt)}
                                                    {cat(2, anatPath, filesep, 'rc2', anatName, '_dartel', anatExt)}
                                                    }';
    matlabbatch{step}.spm.tools.dartel.warp1.settings.rform = 0;
    matlabbatch{step}.spm.tools.dartel.warp1.settings.param(1).its = 3;
    matlabbatch{step}.spm.tools.dartel.warp1.settings.param(1).rparam = [4 2 1e-06];
    matlabbatch{step}.spm.tools.dartel.warp1.settings.param(1).K = 0;
    matlabbatch{step}.spm.tools.dartel.warp1.settings.param(1).template = cellstr(templates{1});
    matlabbatch{step}.spm.tools.dartel.warp1.settings.param(2).its = 3;
    matlabbatch{step}.spm.tools.dartel.warp1.settings.param(2).rparam = [2 1 1e-06];
    matlabbatch{step}.spm.tools.dartel.warp1.settings.param(2).K = 0;
    matlabbatch{step}.spm.tools.dartel.warp1.settings.param(2).template = cellstr(templates{2});
    matlabbatch{step}.spm.tools.dartel.warp1.settings.param(3).its = 3;
    matlabbatch{step}.spm.tools.dartel.warp1.settings.param(3).rparam = [1 0.5 1e-06];
    matlabbatch{step}.spm.tools.dartel.warp1.settings.param(3).K = 1;
    matlabbatch{step}.spm.tools.dartel.warp1.settings.param(3).template = cellstr(templates{3});
    matlabbatch{step}.spm.tools.dartel.warp1.settings.param(4).its = 3;
    matlabbatch{step}.spm.tools.dartel.warp1.settings.param(4).rparam = [0.5 0.25 1e-06];
    matlabbatch{step}.spm.tools.dartel.warp1.settings.param(4).K = 2;
    matlabbatch{step}.spm.tools.dartel.warp1.settings.param(4).template = cellstr(templates{4});
    matlabbatch{step}.spm.tools.dartel.warp1.settings.param(5).its = 3;
    matlabbatch{step}.spm.tools.dartel.warp1.settings.param(5).rparam = [0.25 0.125 1e-06];
    matlabbatch{step}.spm.tools.dartel.warp1.settings.param(5).K = 4;
    matlabbatch{step}.spm.tools.dartel.warp1.settings.param(5).template = cellstr(templates{5});
    matlabbatch{step}.spm.tools.dartel.warp1.settings.param(6).its = 3;
    matlabbatch{step}.spm.tools.dartel.warp1.settings.param(6).rparam = [0.25 0.125 1e-06];
    matlabbatch{step}.spm.tools.dartel.warp1.settings.param(6).K = 6;
    matlabbatch{step}.spm.tools.dartel.warp1.settings.param(6).template = cellstr(templates{6});
    matlabbatch{step}.spm.tools.dartel.warp1.settings.optim.lmreg = 0.01;
    matlabbatch{step}.spm.tools.dartel.warp1.settings.optim.cyc = 3;
    matlabbatch{step}.spm.tools.dartel.warp1.settings.optim.its = 3;
    step = step + 1;
    end
    matlabbatch{step}.spm.tools.dartel.mni_norm.template =  cellstr(templates{6});
    matlabbatch{step}.spm.tools.dartel.mni_norm.data.subj.flowfield = {cat(2, anatPath, filesep, 'u_rc1', anatName, '_dartel',anatExt)};

    if ~iscell(funcNative)
    matlabbatch{step}.spm.tools.dartel.mni_norm.data.subj.images = {funcNative};
    else
    matlabbatch{step}.spm.tools.dartel.mni_norm.data.subj.images = (funcNative)';
    end
    matlabbatch{step}.spm.tools.dartel.mni_norm.vox = [3 3 3];
    matlabbatch{step}.spm.tools.dartel.mni_norm.bb = [NaN NaN NaN
                                                NaN NaN NaN];
    matlabbatch{step}.spm.tools.dartel.mni_norm.preserve = 0;
    matlabbatch{step}.spm.tools.dartel.mni_norm.fwhm = smoothFwhm;


    %matlabbatch{step}.spm.tools.dartel.crt_warped.flowfields = {cat(2, anatPath, filesep, 'u_rc1', anatName, '_dartel', anatExt)};
    %matlabbatch{step}.spm.tools.dartel.crt_warped.images = (funcNative)';
    %matlabbatch{step}.spm.tools.dartel.crt_warped.jactransf = 0;
    %matlabbatch{step}.spm.tools.dartel.crt_warped.K = 6;
    %matlabbatch{step}.spm.tools.dartel.crt_warped.interp = interpMethod;
end


if normalise(3) == 2
      if ~exist(cat(2, anatPath, filesep, 'u_rc1', anatName, anatExt), 'file')
            matlabbatch{step}.spm.tools.dartel.warp1.images = {
                                                            {cat(2, anatPath, filesep, 'rc1', anatName, '_dartel', anatExt)}
                                                            {cat(2, anatPath, filesep, 'rc2', anatName, '_dartel', anatExt)}
                                                            }';



            matlabbatch{step}.spm.tools.dartel.warp1.settings.rform = 0;
            matlabbatch{step}.spm.tools.dartel.warp1.settings.param(1).its = 3;
            matlabbatch{step}.spm.tools.dartel.warp1.settings.param(1).rparam = [4 2 1e-06];
            matlabbatch{step}.spm.tools.dartel.warp1.settings.param(1).K = 0;
            matlabbatch{step}.spm.tools.dartel.warp1.settings.param(1).template = cellstr(templates{1});
            matlabbatch{step}.spm.tools.dartel.warp1.settings.param(2).its = 3;
            matlabbatch{step}.spm.tools.dartel.warp1.settings.param(2).rparam = [2 1 1e-06];
            matlabbatch{step}.spm.tools.dartel.warp1.settings.param(2).K = 0;
            matlabbatch{step}.spm.tools.dartel.warp1.settings.param(2).template = cellstr(templates{2});
            matlabbatch{step}.spm.tools.dartel.warp1.settings.param(3).its = 3;
            matlabbatch{step}.spm.tools.dartel.warp1.settings.param(3).rparam = [1 0.5 1e-06];
            matlabbatch{step}.spm.tools.dartel.warp1.settings.param(3).K = 1;
            matlabbatch{step}.spm.tools.dartel.warp1.settings.param(3).template = cellstr(templates{3});
            matlabbatch{step}.spm.tools.dartel.warp1.settings.param(4).its = 3;
            matlabbatch{step}.spm.tools.dartel.warp1.settings.param(4).rparam = [0.5 0.25 1e-06];
            matlabbatch{step}.spm.tools.dartel.warp1.settings.param(4).K = 2;
            matlabbatch{step}.spm.tools.dartel.warp1.settings.param(4).template = cellstr(templates{4});
            matlabbatch{step}.spm.tools.dartel.warp1.settings.param(5).its = 3;
            matlabbatch{step}.spm.tools.dartel.warp1.settings.param(5).rparam = [0.25 0.125 1e-06];
            matlabbatch{step}.spm.tools.dartel.warp1.settings.param(5).K = 4;
            matlabbatch{step}.spm.tools.dartel.warp1.settings.param(5).template = cellstr(templates{5});
            matlabbatch{step}.spm.tools.dartel.warp1.settings.param(6).its = 3;
            matlabbatch{step}.spm.tools.dartel.warp1.settings.param(6).rparam = [0.25 0.125 1e-06];
            matlabbatch{step}.spm.tools.dartel.warp1.settings.param(6).K = 6;
            matlabbatch{step}.spm.tools.dartel.warp1.settings.param(6).template = cellstr(templates{6});
            matlabbatch{step}.spm.tools.dartel.warp1.settings.optim.lmreg = 0.01;
            matlabbatch{step}.spm.tools.dartel.warp1.settings.optim.cyc = 3;
            matlabbatch{step}.spm.tools.dartel.warp1.settings.optim.its = 3;
            step = step + 1;
      end


    matlabbatch{step}.spm.tools.dartel.crt_iwarped.flowfields = {cat(2, anatPath, filesep, 'u_rc1', anatName, '_dartel', anatExt)};
    matlabbatch{step}.spm.tools.dartel.crt_iwarped.images =  (funcNative)';
    matlabbatch{step}.spm.tools.dartel.crt_iwarped.K = 6;
    matlabbatch{step}.spm.tools.dartel.crt_iwarped.interp = interpMethod;

end

spm_jobman('run', matlabbatch);



%%

%V = spm_vol(cat(2, baseDir, filesep, refNative));
%dlmwrite(cat(2, baseDir, filesep, 'anat2func.mat'), V.mat);
