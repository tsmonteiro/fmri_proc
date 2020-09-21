function affine_spm2fsl(baseDir, ref, moving)
% WORKING for the Chinese dataset, but needs checking
% Based on http://eeg.sourceforge.net/MJenkinson_coordtransforms.pdf
% However, some changes are necessary as I force everything to be RAS+
% space
% The exception to this is when some FSL function is called, becaused it
% messes thing up.
Vnii = load_nii_hdr(ref);

V = spm_vol(ref);
K = spm_vol(moving);

volCent = [Vnii.hist.qoffset_x Vnii.hist.qoffset_y Vnii.hist.qoffset_z];

VS = [1 0 0 +1; ...
      0 1 0 +1; ...
      0 0 1 +1; ...
      0 0 0 1];

SCdest = [1 0 0 volCent(1); ...
          0 1 0 volCent(2); ...
          0 0 1 volCent(3); ...
          0 0 0 1];

Kmat = K.private.mat ;


P  = spm_imatrix(Kmat );
PT=P;
PT(1)=-P(1);
PT(2)=-P(2);
PT(3)=-P(3);
matOrder  = 'T*R*Z*S';
Amat = spm_matrix(PT, matOrder);
Aspm = Amat;%M.mat(:,:,end);

Ddest = zeros(4,4);
Dsource = zeros(4,4);

for k = 1:3
    Ddest(k,k) = V.private.mat(k,k);
    Dsource(k,k) = K.private.mat(k,k);
end
Ddest(4,4)=1;
Dsource(4,4)=1;

Afsl = inv(Ddest * VS * inv(SCdest) * inv(Ddest) * Aspm  * inv(VS) * inv(Dsource));


dlmwrite(cat(2, baseDir, filesep, 'anat2func_fsl.mat'), Afsl, 'delimiter', '\t');
