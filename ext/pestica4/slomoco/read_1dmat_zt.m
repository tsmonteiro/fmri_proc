function [transmat_12dof,fpparams_6dof]=read_1dmat_zt(filestr);
% [rotmat_12dof,fpparams_6dof]=read_rotmat_zt(filestr);
%
% Note that read_rotmat_zt reads time series of 12 dof transformation matrix 

mocoparams1=textread(sprintf('%s.0000.1D',filestr),'','headerlines',1);
tdim=size(mocoparams1,1);
[w,s]=system(sprintf('ls %s.????.1D | grep 1D | nl | tail -n 1 | gawk ''{print $1}''',filestr));
zdim=str2num(s);

% rotmat_12dof = [3x3 transformation matrix ; 3x1 shift vector ]
transmat_12dof = zeros(zdim,tdim,12);
for z=1:zdim
  mocoparams=textread(sprintf('%s.%04d.1D',filestr,z-1),'','headerlines',1);
  transmat_12dof(z,:,:) = mocoparams;
end

fpparams_6dof = zeros(zdim,tdim,6);
xyzr = zeros(zdim,tdim,3);

for z = 1:zdim
  for t = 1:tdim
    rotmat_1d = squeeze(transmat_12dof(z,t,:));
    rotmat_3x3 = [rotmat_1d(1) rotmat_1d(2)  rotmat_1d(3) ; ...
                  rotmat_1d(5) rotmat_1d(6)  rotmat_1d(7) ; ...
                  rotmat_1d(9) rotmat_1d(10) rotmat_1d(11)];
    xyzr = convert_rotmat_into_rots(rotmat_3x3); % radian
    fpparams_6dof(z,t,1:3) = [rotmat_1d(4) rotmat_1d(8)  rotmat_1d(12)];
    [xr,yr,zr]=convert_rotmat_into_rots(rotmat_3x3);
    fpparams_6dof(z,t,4:6) = [xr yr zr].*180/pi; % degree
  end
end
