function tranmat_zt = convert_1dmat_into_tranmat(transmat1d_zt)

[zdim tdim veclen] = size(transmat1d_zt);
tranmat_zt = zeros(zdim, tdim, 3, 4);
for n = 1:zdim
  for t = 1:tdim
    tmpv = squeeze(transmat1d_zt(n,t,:)); % [12 x 1]
    tranmat_zt(n,t,:,:) = [tmpv(1:4)' ;tmpv(5:8)' ; tmpv(9:12)' ];
%   tranmat = [mat1d(1) mat1d(2)  mat1d(3)  mat1d(4); ...
%              mat1d(5) mat1d(6)  mat1d(7)  mat1d(8); ...
%              mat1d(9) mat1d(10) mat1d(11) mat1d(12)] 
  end
end
