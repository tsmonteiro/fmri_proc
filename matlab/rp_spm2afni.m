function rp_spm2afni(infile, outfile)

motParams = dlmread(infile);

motParams = motParams(2:end, [4 5 6 1 2 3]);
motParams(1:end,[1 2 3]) = rad2deg( motParams(:,[1 2 3]) );
motParams(isnan(motParams)) = 0;



dlmwrite(outfile, motParams, 'delimiter', ' ', 'precision', '%.5f');
