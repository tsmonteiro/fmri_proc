function apply_phycaa_rsn(TR, OUTDIR, FNAME, MASK, CSFMASK )
%%
addpath '/home/fsluser/Documents/phycaa';
addpath '/home/fsluser/Documents/phycaa/extra_CAA_scripts';
addpath '/home/fsluser/Documents/phycaa/NIFTI_tools';
%%
%TR=2.5;
%OUTDIR='/mnt/hgfs/ssd_tmp/RS017/';
%FNAME='proc_data_mni_fix_am';
%MASK=cat(2, OUTDIR, filesep, 'mni_mask.nii');
%CSFMASK=cat(2, OUTDIR, filesep, 'csf_mask_mni.nii');




run_PHYCAA_plus(cat(2, OUTDIR, filesep, FNAME,'.nii'), ...
    MASK, ...
    CSFMASK,...
    [], struct('FreqCut', 0.1, 'keepmean', 1, 'TR', TR, 'make_output', 1), ...
    2, [OUTDIR filesep 'phyca_data_'],2);


movefile(cat(2, OUTDIR, filesep, FNAME,'_PHYCAA_step1+2.nii'),...
       cat(2, OUTDIR, filesep, FNAME, 'p.nii'));
