function [fpparams,fpparams_6dof]=read_motion_newslicealg(filestr,slice_timing);
%function [fpparams,fpparams_6dof]=read_motion_newslicealg(filestr,slice_timing);
% return slcdet moco params in same temporal acquisition order as fpfile_S40vol (i.e. timepoint 1 is acquired first, then timepoint 2, then timepoint 3)
% whereas slices may be acquired in some other fashion (typically interleaved), so this output should line up with injected motion (if used) as used
% in SimPACE, but regression corrections must sub-sample this according to the true slice acquisition order 
%   (e.g. slice 1=time1, slice 2=time17, slice 3=time2, etc)
% cleanup ideas: remove periodic (volumetric) motion signals that drift over time, some due to contrast changes (moving average coregistration)
mocoparams1=textread(sprintf('%s.0000.1D',filestr),'','headerlines',1);
tdim=size(mocoparams1,1)
[w,s]=system(sprintf('ls %s.????.1D | grep 1D | nl | tail -n 1 | gawk ''{print $1}''',filestr));
zdim=str2num(s)

[sorted_slice_timing zdims] = sort(slice_timing);
disp('read slice timing info')

for z=1:zdim
 mocoparams1=textread(sprintf('%s.%04d.1D',filestr,z-1),'','headerlines',1);
 % re-interleave to make full z*t vectors in slice-acquisition timing order, so 1st timepoint is acquired first, 2nd acquired 2nd, and so on
 zstart=find(zdims==z);
 % make 12dof matrix for regression
 fpparams(zstart:zdim:zdim*tdim,:)=mocoparams1(:,[4 8 12 1 2 3 5 6 7 9 10 11]);
end

rotmats=fpparams(:,4:12);
for i=1:tdim*zdim
  [xr(i),yr(i),zr(i)]=convert_rotmat_into_rots(reshape(rotmats(i,:),[3 3]));
end
fpparams_6dof=[fpparams(:,1:3) [xr' yr' zr']*180/pi];

