function [inplane,inplane6dof,filt_inplane]=read_slicemocoxy_AFNI_files_auto(moco_folder,slice_timing)
% script is simply called from current workign dir, assumes files are in subdir tempslmoco
% slice_timing is directoy read

if (exist('moco_folder')==0)
  moco_folder='tempslmoco';
end
if (exist(moco_folder)~=7)
  disp(sprintf('moco file folder %s does not exist',moco_folder));
end

[w,s]=system(sprintf('ls %s/motion.allineate.slicewise_inplane.0000.????.1D | grep motion | nl | tail -n 1 | gawk ''{print $1}''',moco_folder));
tdim=str2num(s)
[w,s]=system(sprintf('ls %s/motion.allineate.slicewise_inplane.????.0000.1D | grep motion | nl | tail -n 1 | gawk ''{print $1}''',moco_folder));
zdim=str2num(s)

[sorted_slice_timing slice_acq_order] = sort(slice_timing);

%disp(sprintf('reading in files starting from: %s/motion.allineate.slicewise_inplane.%04d.%04d.1D',moco_folder,0,0));
for t=1:tdim
 for k=1:zdim
  % convert inplane.txt to slicewise motion params for second order motion correction
  motmat=textread(sprintf('%s/motion.allineate.slicewise_inplane.%04d.%04d.1D',moco_folder,k-1,t-1),'','headerlines',1);
  if (isempty(motmat))
    motmat=[1 0 0 0 0 1 0 0 0 0 1 0];
  end
%  motmat=motmat(1:4,1:4);
  % find the temporal location to insert this motion for this slice
  zd=find(slice_acq_order==k);
  % inplane is xyz trans and then the 3x3 rotmat
  inplane(((t-1)*zdim)+zd,:)=[motmat(4) motmat(8) motmat(12) motmat(1:3) motmat(5:7) motmat(9:11)];
 end
end

na=inplane(:,[1:2]);
na=reshape(na,[zdim tdim 2]);
% first remove any trend for non-motion noise
na(find(abs(na)>0.25))=0;
na=na-repmat(mean(na,2),[1 tdim 1]);
na=reshape(inplane(:,1:2),[zdim tdim 2])-repmat(mean(na,2),[1 tdim 1]);
% now remove really big spikes (better way is to check if there are any neighboring outliers - usually slicemocoxy fails catastrophically on one slice, not two)
na(find(abs(na)>4))=0;
% and remove noise motions below 75 microns
na(find(abs(na)<0.075))=0;
filt_inplane=inplane;
filt_inplane(:,1:2)=reshape(na,[zdim*tdim 2]);


rotmats=inplane(:,4:12);
for i=1:tdim*zdim
  [xr(i),yr(i),zr(i)]=convert_rotmat_into_rots(reshape(rotmats(i,:),[3 3]));
end
inplane6dof=[inplane(:,1:3) [xr' yr' zr']*180/pi];

inplane=inplane6dof(:,[1 2 6]);

