function upsample(dataIn, dataOut, lf)
  % Superresolution upsampling
% dataIn='/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/tmp/B1_63/anat/mean_func_data_ndsc.nii';
% dataOut='/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/tmp/B1_63/anat/mean_func_data_ndsu.nii';
% lf = [2 2 1];

addpath(genpath('/home/luna.kuleuven.be/u0101486/workspace/matlab/toolbox/nonlocal_superresolution'));
% read LR data


V1=spm_vol(dataIn);
ima=spm_read_vols(V1);
ima=double(ima);

[lima]=NLMUpsample2(ima,lf);

lima(lima<0)=0;

%%

%%
VA=spm_vol(dataIn);
Vtmp = V1;
l = lima;

%%

% Write results
V1 = Vtmp;
lima=l;




lima = lima(:,end:-1:1, end:-1:1);
V1.dim=size(lima);

rat = [V1.mat(1,1)/lf(1) V1.mat(2,2)/lf(2) V1.mat(3,3)/lf(3)]';
V1.mat(1:3,4) = V1.mat(1:3,4).*rat;

V1.mat(1,1)=V1.mat(1,1)/lf(1);
V1.mat(2,2)=V1.mat(2,2)/lf(2);
V1.mat(3,3)=V1.mat(3,3)/lf(3);


V1.fname= dataOut;
spm_write_vol(V1,lima);
