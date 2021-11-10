function anat_proc_surf(BASEDIR, T1_FILE)

cat_path = which('cat12');
cat_parts = strsplit(cat_path, filesep);

cat_path = [];
for i = 1:length(cat_parts)-1
    cat_path = cat(2, cat_path, filesep, cat_parts{i});
end


spm_path = which('spm');
spm_parts = strsplit(spm_path, filesep);

spm_path = [];
for i = 1:length(spm_parts)-1
    spm_path = cat(2, spm_path, filesep, spm_parts{i});
end



matlabbatch = {};
matlabbatch{1}.spm.tools.cat.estwrite.data = cellstr(cat(2, BASEDIR,filesep,  T1_FILE, ',1'));
matlabbatch{1}.spm.tools.cat.estwrite.nproc = 12;
matlabbatch{1}.spm.tools.cat.estwrite.opts.tpm = cellstr(cat(2, spm_path, filesep, 'tpm/TPM.nii'));
matlabbatch{1}.spm.tools.cat.estwrite.opts.affreg = 'mni';
matlabbatch{1}.spm.tools.cat.estwrite.opts.biasstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.opts.accstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.APP = 1070;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.LASstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.gcutstr = 2;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shooting.shootingtpm = cellstr(cat(2, cat_path, filesep, 'templates_1.50mm/Template_0_IXI555_MNI152_GS.nii'));
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shooting.regstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.vox = 1.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.restypes.fixed = [1 0.1];
matlabbatch{1}.spm.tools.cat.estwrite.output.surface = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.neuromorphometrics = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.lpba40 = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.cobra = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.hammers = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.labelnative = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.jacobianwarped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.warps = [1 0];


matlabbatch{2}.spm.tools.cat.stools.surfresamp.data_surf = cellstr(cat(2, BASEDIR, filesep, 'surf/lh.thickness.t1_f'));
matlabbatch{2}.spm.tools.cat.stools.surfresamp.merge_hemi = 1;
matlabbatch{2}.spm.tools.cat.stools.surfresamp.mesh32k = 0;
matlabbatch{2}.spm.tools.cat.stools.surfresamp.fwhm_surf = 15;
matlabbatch{2}.spm.tools.cat.stools.surfresamp.nproc = 0;

matlabbatch{3}.spm.tools.cat.stools.surfresamp.data_surf = cellstr(cat(2, BASEDIR, filesep, 'surf/lh.thickness.t1_f'));
matlabbatch{3}.spm.tools.cat.stools.surfresamp.merge_hemi = 1;
matlabbatch{3}.spm.tools.cat.stools.surfresamp.mesh32k = 0;
matlabbatch{3}.spm.tools.cat.stools.surfresamp.fwhm_surf = 25;
matlabbatch{3}.spm.tools.cat.stools.surfresamp.nproc = 0;

matlabbatch{4}.spm.tools.cat.stools.surfextract.data_surf = cellstr(cat(2, BASEDIR, filesep, 'surf/lh.central.t1_f.gii'));
matlabbatch{4}.spm.tools.cat.stools.surfextract.GI = 1;
matlabbatch{4}.spm.tools.cat.stools.surfextract.FD = 0;
matlabbatch{4}.spm.tools.cat.stools.surfextract.SD = 1;
matlabbatch{4}.spm.tools.cat.stools.surfextract.nproc = 0;

matlabbatch{5}.spm.tools.cat.stools.surf2roi.cdata = {cellstr(cat(2, BASEDIR, filesep, 'surf/lh.thickness.t1_f'))};


%     matlabbatch{2}.spm.tools.cat.stools.vol2surf.data_vol = '<UNDEFINED>';
%     matlabbatch{2}.spm.tools.cat.stools.vol2surf.data_mesh_lh = '<UNDEFINED>';
%     matlabbatch{21.spm.tools.cat.stools.vol2surf.sample = {'maxabs'};
%     matlabbatch{2}.spm.tools.cat.stools.vol2surf.datafieldname = 'intensity';
%     matlabbatch{2}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.class = 'GM';
%     matlabbatch{2}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.startpoint = -0.5;
%     matlabbatch{2}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.steps = 7;
%     matlabbatch{2}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.endpoint = 0.5;

matlabbatch{4}.spm.tools.cat.tools.calcvol.data_xml =  cellstr(cat(2, BASEDIR, filesep, 'report/cat_t1_f.xml'));
matlabbatch{4}.spm.tools.cat.tools.calcvol.calcvol_TIV = 1;
matlabbatch{4}.spm.tools.cat.tools.calcvol.calcvol_name =  cellstr(cat(2, BASEDIR, filesep, 'TIV.txt'));

spm_jobman('run', matlabbatch);
