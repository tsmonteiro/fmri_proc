function tranmat_zt = convert_1dmat_into_tranmatarray(transmat1d_zt)

[zdim tdim veclen] = size(transmat1d_zt);
for n = 1:zdim
  for t = 1:tdim
    tmpv = squeeze(transmat1d_zt(n,t,:)); % [12 x 1]
    tranmat_zt(n,t).R = [tmpv(1) tmpv(2)  tmpv(3); ...
                         tmpv(5) tmpv(6)  tmpv(7); ...
                         tmpv(9) tmpv(10) tmpv(11) ];
    
    tranmat_zt(n,t).T = [tmpv(4); tmpv(8); tmpv(12);];
  end
end
