function realign_inria( tmpDir, prefixFunc, refVol, nVols )
%%
addpath('/media/thiago/EXTRALINUX/fmri_proc/matlab/toolbox/spm12/compat');
%addpath('/home/luna.kuleuven.be/u0101486/workspace/matlab/toolbox/spm12/toolbox/INRIAlign');
%tmpDir = '/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/tmp/B1_07';
%prefixFunc = 'asbfunc_data';
%refVol = 150;
refVol = refVol +1;
%nVols = 300;

nii = load_untouch_nii( cat(2, tmpDir, filesep, prefixFunc, '.nii'));
nii.img = cat(4, nii.img(:,:,:,refVol), nii.img(:,:,:,1:end) );
nii.hdr.dime.dim(5) = nii.hdr.dime.dim(5) + 1;
save_untouch_nii(nii,cat(2, tmpDir, filesep, prefixFunc, '.nii') );


files = {cat(2, tmpDir, filesep, prefixFunc, '.nii')};
%         rho_func - A string indicating the cost function that
%                    will be used to caracterize intensity errors.
%                    Possible choices are :
%
%                      'quadratic' -  fast, but not robust...
%                      'absolute'  -  quite slow, not very robust
%                        'huber'   -  Huber function
%                        'cauchy'  -  Cauchy function
%                        'geman'   -  Geman-McClure function
%                       'leclerc'  -  Leclerc-Welsch function
%                        'tukey'   -  Tukey's biweight function
%
%                    DEFAULT: 'geman', usually a good trade-off
%                    between robustness and speed.

flags = struct;
flags.rho_func = 'geman';
flags.cut_off = 2.5;
flags.quality = 0.9;
flags.fwhm = 3;
flags.sep  = 5;



tic;
inria_realign(cellstr(files),flags);
toc;
rmpath('/media/thiago/EXTRALINUX/fmri_proc/matlab/toolbox/spm12/compat');
% %%
flags = struct;
flags.mask = 0;
flags.mean = 0;
flags.interp = 5;
flags.prefix = 'r';

spm_reslice(cellstr(files), flags);

% %%
motParams = zeros(nVols, 6);

fname = strrep( cat(2, tmpDir, filesep, 'realignment_params_', prefixFunc, '.nii'), '.nii', '.txt');

motParams = dlmread(fname);

motParams = motParams(2:end, [4 5 6 1 2 3]);
motParams(1:end,[1 2 3]) = rad2deg( motParams(:,[1 2 3]) );
motParams(isnan(motParams)) = 0;




dlmwrite(cat(2, tmpDir, filesep, 'motion_estimate.par'), motParams, 'delimiter', ' ', 'precision', '%.5f');

fig = figure();

plot( motParams );

saveas(fig,cat(2, tmpDir, filesep, 'motion_estimate_0.png') );

% Remove reference volume
nii     = load_untouch_nii( cat(2, tmpDir, filesep,'r',  prefixFunc, '.nii'));
nii.img = nii.img(:,:,:,2:end);
nii.img(isnan(nii.img)) = 0;
nii.img(nii.img < 0) = 0;

nii.hdr.dime.dim(5) = nii.hdr.dime.dim(5) - 1;
save_untouch_nii(nii,cat(2, tmpDir, filesep, 'r', prefixFunc, '.nii') );

% %%

%plot(motParams)
