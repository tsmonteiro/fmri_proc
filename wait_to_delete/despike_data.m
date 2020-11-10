function despike_data(OUTDIR, FNAME, FOUT, RP_FILE, MASK_FILE)
%%
% addpath(genpath('/home/fsluser/Documents/fMRIData'))
% addpath('/home/fsluser/Documents/ArtRepair_v5b3/ArtRepair')
addpath(genpath('/home/fsluser/Documents/SPIKECOR'))
%OUTDIR='/mnt/hgfs/ssd_tmp/TASK_001_1/';
%RUN=1;
%FNAME = 'mfunc_data.nii';
%RP_FILE='motion_estimate.par';
%MASK_FILE='nat_mask.nii';
%FOUT='dmfunc_data.nii';
DATA=cat(2, OUTDIR, filesep, FNAME);
MPE=cat(2, OUTDIR, filesep, RP_FILE);
DATA_OUT=cat(2, OUTDIR, filesep, FOUT);
MASK=cat(2,OUTDIR, filesep, MASK_FILE); 


spikecor(DATA,MASK, MPE, ...
    cat(2, OUTDIR, filesep, 'spk_', FNAME), ...
    'volume',DATA_OUT);
system(['rm ' OUTDIR filesep 'spk*.mat']);



