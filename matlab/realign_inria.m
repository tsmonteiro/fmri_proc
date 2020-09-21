function realign_inria( tmpDir, prefixFunc )
%%
addpath('/home/luna.kuleuven.be/u0101486/workspace/matlab/toolbox/spm12/compat');
addpath('/home/luna.kuleuven.be/u0101486/workspace/matlab/toolbox/spm12/toolbox/INRIAlign');
tmpDir = '/home/luna.kuleuven.be/u0101486/workspace/data/ConnectEx/tmp/A003';
prefixFunc = 'func_data';


rFiles = dir(cat(2, tmpDir, filesep, prefixFunc, '*.nii'));

files = cell(length(rFiles),1);

for f = 1:length(rFiles)
    files{f} = cat(2, tmpDir, filesep, rFiles(f).name);
end
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
flags.cut_off = 1.5;
flags.quality = 0.9;
flags.fwhm = 0;
flags.sep  = 3;

clc;

inria_realign(cellstr(files),flags);
rmpath('/home/luna.kuleuven.be/u0101486/workspace/matlab/toolbox/spm12/compat');

flags = struct;
flags.mask = 0;
flags.mean =0;
flags.interp = 3;
flags.prefix = 'r2';

spm_reslice(cellstr(files), flags);

%%
motParams = zeros(length(rFiles)-1, 6);

for f = 2:length(rFiles)
    fname = strrep( cat(2, tmpDir, filesep, 'realignment_params_', rFiles(f).name), '.nii', '.txt');
    
    motParams(f-1,:) = dlmread(fname);
    
end

dlmwrite(cat(2, tmpDir, filesep, 'motion_parameters.par'), motParams);


%%

plot(motParams)


