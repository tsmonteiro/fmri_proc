% function generate_dartel_template(anatDir)


% CHANGE this line to wherever the imported T1's are stored
anatDir = '/media/thiago/EXTRALINUX/fmri_proc/DATA/Template/';
ABIN='/media/thiago/EXTRALINUX/ANTs/abin/bin/';

matlabbatch = {};
c1 = dir(cat(2,anatDir, filesep, 'rc1*.nii'));
c2 = dir(cat(2,anatDir, filesep, 'rc2*.nii'));
c3 = dir(cat(2,anatDir, filesep, 'rc3*.nii'));

t1 = dir(cat(2,anatDir, filesep, 't1*.nii'));

filesC1 = {};
filesC2 = {};
filesT1 = {};

gms = [];
wms = [];

for f = 1:length(c1)
    disp(f);
    
    [~, tid] = strtok( c1(f).name, '_');
%     [tid,~] = strtok(tid, '_');
%     tid = str2double(tid(2:end));
    
    
    filesC1{end+1} = cat(2, anatDir, filesep, c1(f).name );
    filesC2{end+1} = cat(2, anatDir, filesep, c2(f).name );
    filesT1{end+1} = cat(2, anatDir, filesep, 't1', tid);
    
    
    nii            = load_nii(cat(2, anatDir, filesep, c1(f).name ));
    gms(end+1)     = sum(nii.img(:)>0.1);
    
    if gms(end) <= 2e5
        disp(filesC1{end});
    end
    
    nii            = load_nii(cat(2, anatDir, filesep, c2(f).name ));
    wms(end+1)     = sum(nii.img(:)>0.1);
end

%% If not using all subjects, the code below can help in selecting
% some sort of representative sample (in terms of brain size)
[sg,si] = sort(gms+wms);
N       = length(si);
% idx = cat(2, 2:2:20, 21:6:N-21, N-20:2:N );
idx     = cat(2, 2:4:N, N );
filesC1 = filesC1(idx);
filesC2 = filesC2(idx);
filesT1 = filesT1(idx);

clf;

plot(sort(gms+wms), 'sq'); hold on;
plot(idx, sg(idx), 'sqk', 'MarkerSize', 7, 'MArkerFaceColor', 'k');
%%
filesT1=filesT1';
filesC1=filesC1';
filesC2=filesC2';


matlabbatch{1}.spm.tools.dartel.warp.images = {
                                               filesC1
                                               filesC2
                                               }';

matlabbatch{1}.spm.tools.dartel.warp.settings.template = 'Template';
matlabbatch{1}.spm.tools.dartel.warp.settings.rform = 0;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).rparam = [4 2 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).K = 0;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).slam = 16;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).rparam = [2 1 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).K = 0;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).slam = 8;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).rparam = [1 0.5 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).K = 1;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).slam = 4;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).rparam = [0.5 0.25 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).K = 2;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).slam = 2;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).rparam = [0.25 0.125 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).K = 4;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).slam = 1;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).rparam = [0.25 0.125 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).K = 6;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).slam = 0.5;
matlabbatch{1}.spm.tools.dartel.warp.settings.optim.lmreg = 0.01;
matlabbatch{1}.spm.tools.dartel.warp.settings.optim.cyc = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.optim.its = 3;

spm_jobman( 'run', matlabbatch );

%%

for i = 1:length(filesT1)

    matlabbatch = {};

    matlabbatch{1}.spm.tools.dartel.crt_warped.flowfields = {strrep( strrep( filesC1{i}, 'rc1', 'u_rc1'), '.nii', '_Template.nii') };
    matlabbatch{1}.spm.tools.dartel.crt_warped.images     = {cellstr(filesT1{i})};
    matlabbatch{1}.spm.tools.dartel.crt_warped.jactransf  = 0;
    matlabbatch{1}.spm.tools.dartel.crt_warped.K          = 6;
    matlabbatch{1}.spm.tools.dartel.crt_warped.interp     = 5;

    disp(strrep( strrep( filesC1{i}, 'rc1', 'u_rc1'), '.nii', '_Template.nii') );
    disp(cellstr(filesT1{i}))
    disp('-------');
    
    spm_jobman( 'run', matlabbatch );
end


%%

wt1 = dir(cat(2,anatDir, filesep, 'w*.nii'));

niiT1 = [];


for f = 1:length(wt1)
    
    filename = cat(2, anatDir, filesep, wt1(f).name );
    wt1Nii   = load_nii(filename);
    
    wt1Nii_   = wt1Nii.img;
    wt1Nii_(wt1Nii_ < max(wt1Nii_(:))*0.05) = 0;
    spm_smooth(wt1Nii_, wt1Nii_, [1 1 1]);
    
    niiT1    = cat(4, niiT1, wt1Nii_);
    
end


niiT1 = nanmedian(niiT1,4);
wt1Nii.img = niiT1;
save_nii(wt1Nii, cat(2, anatDir, filesep, 'Group_T1_Avg_Brain.nii'));

%%
REF = cat(2, anatDir, filesep, 'Group_T1_Avg_Brain.nii');
OUT = cat(2, anatDir, filesep, 'mni2group');
MOV = '/media/thiago/EXTRALINUX/fmri_proc/atlas/MNI152_T1_1mm.nii';



CMD='/home/luna.kuleuven.be/u0101486/workspace/fmri_proc/coreg_ants.sh ';
CMD=cat(2, CMD, ABIN, ' ');
CMD=cat(2, CMD, REF, ' ');
CMD=cat(2, CMD, MOV, ' ');
CMD=cat(2, CMD, OUT );
clc;

% RUn this on the terminal...
disp(CMD);
