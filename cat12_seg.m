function cat12_seg(SUBDIR, T1)
matlabbatch = {};
matlabbatch{1}.spm.tools.cat.estwrite.data = {T1};
matlabbatch{1}.spm.tools.cat.estwrite.nproc = 0;
matlabbatch{1}.spm.tools.cat.estwrite.opts.tpm = {'/home/luna.kuleuven.be/u0101486/workspace/matlab/toolbox/spm12/tpm/TPM.nii'};
matlabbatch{1}.spm.tools.cat.estwrite.opts.affreg = 'eastern';
matlabbatch{1}.spm.tools.cat.estwrite.opts.biasstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.opts.accstr = 1;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.APP = 2;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.LASstr = 0.75;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.gcutstr = 2;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.dartel.darteltpm = {'/home/luna.kuleuven.be/u0101486/workspace/matlab/toolbox/spm12/toolbox/cat12/templates_1.50mm/Template_1_IXI555_MNI152.nii'};
matlabbatch{1}.spm.tools.cat.estwrite.extopts.vox = 1;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.restypes.fixed = [1 0.1];
matlabbatch{1}.spm.tools.cat.estwrite.output.surface = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.neuromorphometrics = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.lpba40 = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.cobra = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.hammers = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.labelnative = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.warped = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.jacobianwarped = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.warps = [0 0];
matlabbatch{2}.spm.tools.cat.tools.calcvol.data_xml(1) = cfg_dep('CAT12: Segmentation: CAT Report', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','catreport', '()',{':'}));
matlabbatch{2}.spm.tools.cat.tools.calcvol.calcvol_TIV = 0;
matlabbatch{2}.spm.tools.cat.tools.calcvol.calcvol_name = 'TIV.txt';
matlabbatch{3}.spm.tools.cat.tools.calcroi.roi_xml(1) = cfg_dep('CAT12: Segmentation: CAT Report', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','catreport', '()',{':'}));
matlabbatch{3}.spm.tools.cat.tools.calcroi.point = '.';
matlabbatch{3}.spm.tools.cat.tools.calcroi.outdir = {SUBDIR};
matlabbatch{3}.spm.tools.cat.tools.calcroi.calcroi_name = 'ROI';


spm_jobman('run', matlabbatch);
