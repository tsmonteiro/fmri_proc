function denoise_mri(inFile, mask, outFile)
%%
    addpath(genpath('/home/luna.kuleuven.be/u0101486/workspace/matlab/toolbox/BM3D'));
    addpath(genpath('/home/luna.kuleuven.be/u0101486/workspace/matlab/toolbox/BM3D_MRI_toolbox'));
    addpath(genpath('/home/luna.kuleuven.be/u0101486/workspace/matlab/toolbox/NIfTI_20140122/'));
    clc
    
    inFile = '/home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/tmp/sub001/t1_brain.nii';
    maskFile = '/home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/tmp/sub001/AnatMask.nii';
    outFile = '/home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/tmp/sub001/ANAT.nii';

    t1 = load_nii(inFile);
    mask = load_nii(maskFile);
    
    for sl = 1:size(t1.img,3)
        fprintf('Doing slice %d / %d\n', sl, size(t1.img,3));
        t = t1.img(:,:,sl);
        if sum(t(:)) < 5
            continue
        end
        [~,t1.img(:,:,sl)]=BM3D_MRI_function(double(t1.img(:,:,sl)), double(mask.img(:,:,sl)));
    end
    
    save_nii(t1, outFile);

