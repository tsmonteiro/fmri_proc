function apply_phycaa_rsn(TR, OUTDIR, FNAME, MASK, CSFMASK, OUTFILE )
%%
addpath(genpath('/home/luna.kuleuven.be/u0101486/workspace/matlab/toolbox/phycaa_plus_2014_09_11'));
% addpath '/home/fsluser/Documents/phycaa/extra_CAA_scripts';
% addpath '/home/fsluser/Documents/phycaa/NIFTI_tools';
%%
%TR=2.5;
%OUTDIR='/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/RS001/';
%FNAME='tmp_smooth_PHY';
%MASK=cat(2, OUTDIR, filesep, 'nat_mask.nii');
%CSFMASK=cat(2, OUTDIR, filesep, 'anat/csf_tpm_native.nii');



%run_PHYCAA_plus( input_cell, mask_name, prior_name, task_SPM_names, dataInfo, num_steps, out_prefix )
run_PHYCAA_plus(cat(2, OUTDIR, filesep, FNAME,'.nii'), ...
    MASK, ...
    CSFMASK,...
    [], struct('FreqCut', 0.1, 'keepmean', 1, 'TR', TR, 'make_output', 1), ...
    2, [OUTDIR filesep 'phyca_data_']);


movefile(cat(2, OUTDIR, filesep, FNAME,'_PHYCAA_step1+2.nii'),...
       cat(2, OUTDIR, filesep, 'p', FNAME, '.nii'));
%%%


mask     = load_untouch_nii(cat(2, OUTDIR, filesep, 'phyca_data__NN_map.nii' ) );
mask1     = 1 - (mask.img == 0);

func     = load_untouch_nii(cat(2, OUTDIR, filesep, FNAME, '.nii' ) );

[x,y,z,t] = size(func.img);
mask2  = std(func.img,0,4) > 0;
func = reshape( func.img, x*y*z, t );



func   = func( mask1(:)== 1 & mask2(:) == 1, :);

%%%

% clf; clc;
[~, score,~,~,expl] = pca(zscore(func'));


nc = max(find(cumsum(expl) < 80));
nc = min([6 nc]);
dlmwrite( cat(2, OUTDIR, filesep, OUTFILE), score(:,1:nc), 'delimiter', '\t' );
%dlmwrite( cat(2, OUTDIR, filesep, OUTFILE), mean(func), 'delimiter', '\t' );

% size(score)

% plot(cumsum(expl), 'sq-k');
