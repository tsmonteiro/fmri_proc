function qa_slomoco(ep2d_filename,filestr_out,filestr_in,filter_width,sub_xy_offsets)
%function qa_slomoco(ep2d_filename,filestr_out,filestr_in,slice_timing,filter_width,sub_xy_offsets)
% script reads in SLOMOCO files and fit data in local directory (currently inside pestica/ subdirectory)
% plot motion parameters, histograms of excessive motion, histograms of motion coupling t-score (sum across model)
% clear all
% ep2d_filename='S42vol.slicemocoxy_afni+orig'
% filestr_out='tempslmoco_volslc_alg_vol_S42vol.slicemocoxy_afni/motion.wholevol_zt'
% filestr_in='tempslmocoxy_afni_S42vol'
if (exist('filter_width')==0)
  filter_width=5;
end
if (exist('sub_xy_offsets')==0)
  sub_xy_offsets=1;
end

% BrikInfo only works on AFNI BRIK format files
[err,ainfo] = BrikInfo(ep2d_filename);
xdim=ainfo.DATASET_DIMENSIONS(1);
ydim=ainfo.DATASET_DIMENSIONS(2);
zdim=ainfo.DATASET_DIMENSIONS(3);
tdim=ainfo.DATASET_RANK(2);
dx=ainfo.DELTA(1);
dy=ainfo.DELTA(2);
dz=ainfo.DELTA(3);
TR=double(ainfo.TAXIS_FLOATS(2));
slice_timing=1000*double(ainfo.TAXIS_OFFSETS);  % ms

% check time unit
[TRsec TRms] = TRtimeunitcheck(TR);
[slice_timing_sec slice_timing_ms] = TRtimeunitcheck(slice_timing);
[MB zmbdim uniq_slice_timing_ms uniq_acq_order] = SMSacqcheck(TRms, zdim, slice_timing_ms);

% get scan orientation from header
[rx,cx]=find(ainfo.Orientation=='R');  % for axial == 1
[ry,cy]=find(ainfo.Orientation=='A');  % for axial == 2
[rz,cz]=find(ainfo.Orientation=='I');  % for axial == 3

% apparently its not uncommon for DELTA to be negative on one or more axes, but don't know why that would be...
voxsize=abs(prod(ainfo.DELTA));

% read volumetric params first (assume AFNI 3dvolreg)
motion=textread('mocoafni.txt'); % n roll(I-S) pitch(R-L) yaw(A-P) dS dL dP

% re-order volumetric params
motion=motion(:,[6 7 5 3 4 2]); % dL dP dS pitch(R-L) yaw(A-P) roll(I-S)

% make TD metric
meanvox=parallelepiped_jiang(motion);

% find most corrupted third of volumetric indices
volinds=find(meanvox>2*median(meanvox));

% read out-of-plane motion parameters
% [outofplane,outofplane6]=read_motion_newslicealg(filestr_out,slice_timing_ms);
[rotmat_zt,outofplane6new]=read_1dmat_zt(filestr_out); % rotmat_zt = [zmbdim tdim 12], outofplane6new = [zmbdim tdim 6]
outofplane=zeros(tdim*zmbdim,12);
outofplane6=zeros(tdim*zmbdim,6);
for t=1:tdim
  for z = 1:zmbdim
    acqtp = find(uniq_acq_order==z) + (t-1)*zmbdim;
    outofplane(acqtp,:) = rotmat_zt(z,t,:);      % [tdim*zmbdim [dL dP dS rotmatvector (1:9)]]
    outofplane6(acqtp,:) = outofplane6new(z,t,:);% [tdim*zmbdim [dL dP dS pitch raw roll]]
  end
end

% read in-plane motion parameters
[inplane,inplane6]=read_slicemocoxy_AFNI_files_auto(filestr_in,uniq_slice_timing_ms); 
% inplane = [tdim*zmbdim [dL dP dS rotmatvector (1:9)]]
% inplane6 = [tdim*zmbdim [dL dP dS pitch raw roll]]

% correct for scan axis orientation: 
% axial [rx ry rz]=[1 2 3], sagittal=[3 1 2], coronal=[1 3 2]
axisz=find([rx ry rz]==3);
if (axisz==2)
  inplane6=inplane6(:,[3 1 2 6 4 5]);
  outofplane6=outofplane6(:,[3 1 2 6 4 5]);
  disp('Swapping axes for coronal acquisition');
elseif (axisz==1)
  inplane6=inplane6(:,[3 2 1 6 5 4]);
  outofplane6=outofplane6(:,[3 2 1 6 5 4]);
  % 1, 4 are inverted
  inplane6(:,[1 4])=-1*inplane6(:,[1 4]);
  outofplane6(:,[1 4])=-1*outofplane6(:,[1 4]);
  disp('Swapping axes for sagittal acquisition');
end

% Step 1: remove mean from out-of-plane params (offset removal)
outofplane6=outofplane6-repmat(mean(outofplane6),[length(outofplane6) 1]);
% detrend slicewise trends after skipping large volume motions
    % NOTE: detrending _SHOULD_ be done after identifying large motion areas and excluding them from the detrend
    % but does not make a large difference in most cases - note, large motion at one end of scan may cause this problem
    % on the other hand, if the motion is not volumetric, this won't help at all, and may make it worse. 
    % Conclusion: none yet, perhaps include as an option. For now it doesn't seem to make any difference, possibly
    % only helps with made-up sequence of motion where some is slicewise and some is volumetric.
    % either way, detrending is critical
for i=1:zmbdim
  a=outofplane6(i:zmbdim:end,:);
  for q=1:6
    a(:,q)=pchip(setxor(volinds,1:tdim),a(setxor(volinds,1:tdim),q),1:tdim);
    p=polyfit(1:tdim,a(:,q)',1);
    outofplane6(i:zmbdim:end,q)=(outofplane6(i:zmbdim:end,q) - polyval(p,1:tdim)');
  end
end


% identify outer-most slices according to slice timing
% [tmp acq_order] = sort(slice_timing);
exclude_slices_one=[];
exclude_slices_two=[];
% identify outer-most slices according to slice timing
% [tmp acq_order] = sort(slice_timing);
if MB == 1
  endslices_two = [find(uniq_acq_order == 1) find(uniq_acq_order == zmbdim)  find(uniq_acq_order == 2)  find(uniq_acq_order == zmbdim-1)]; 
  endslices_one = [find(uniq_acq_order == 1) find(uniq_acq_order == zmbdim)];
  for i=1:length(outofplane6)/zmbdim
    exclude_slices_one=[exclude_slices_one [((zmbdim*(i-1))+endslices_one)]];
    exclude_slices_two=[exclude_slices_two [((zmbdim*(i-1))+endslices_two)]];
  end
else
  endslices_two = []; endslices_one = [];
end

% Step 2: normalize stddev across slices for out-of-plane parameters
for i=1:zmbdim
  stdz(i,:)=std(outofplane6(i:zmbdim:end,:));
end
meanstdz=mean(stdz)/max(std(inplane6));
stdz=stdz./repmat(meanstdz,[zmbdim 1]);
for i=1:zmbdim
  % if excluding this slice, don't make it infinity (not going to matter, but for debugging its nice to be able to see it)
  if (length(find(exclude_slices_two==i))>0)
    outofplane6(i:zmbdim:end,:)=detrend(outofplane6(i:zmbdim:end,:)./repmat(mean(stdz(setxor(1:zmbdim,endslices_two),:)),[tdim 1]));
  else
    outofplane6(i:zmbdim:end,:)=detrend(outofplane6(i:zmbdim:end,:)./repmat(stdz(i,:),[tdim 1]));
  end
end


% Step 3: interpolate over outer two end slices (two on each end) for in-planes
inputmesh_two=setxor(1:length(outofplane6),exclude_slices_two);
inputmesh_one=setxor(1:length(outofplane6),exclude_slices_one);
for i=1:6
  % can't trust the in-plane or out-of-plane motion in outer two slices - this may be dependent on # of voxels in those slices
  % this is worst when its slice #2 (half-way thru a stack of odd # of slices) and slice #30 (last even in stack of odds)
  % but the first and last odd is also modestly bad. This is entirely from out-of-plane motion
  % unfortunately, we cannot be sure whether a given motion is really in-plane or just apparent in-plane
  % so we have no choice unless we can obtain some other information
  inplane6(:,i)=pchip(inputmesh_two,inplane6(inputmesh_two,i),1:length(inplane6));
  outofplane6(:,i)=pchip(inputmesh_two,outofplane6(inputmesh_two,i),1:length(outofplane6));
end
% and re-normalize by slice stddev again after the interpolation
newmot=outofplane6;
for i=1:zmbdim
  stdz(i,:)=std(outofplane6(i:zmbdim:end,:));
end
%meanstdz=mean(stdz);
meanstdz=mean(stdz)/max(std(inplane6));
stdz=stdz./repmat(meanstdz,[zmbdim 1]);
for i=1:zmbdim
  % if excluding this slice, don't make it infinity (not going to matter, but for debugging its nice to be able to see it)
  if (length(find(exclude_slices_two==i))>0)
    newmot(i:zmbdim:end,:)=detrend(outofplane6(i:zmbdim:end,:)./repmat(mean(stdz(setxor(1:zmbdim,endslices_two),:)),[tdim 1]));
  else
    newmot(i:zmbdim:end,:)=detrend(outofplane6(i:zmbdim:end,:)./repmat(stdz(i,:),[tdim 1]));
  end
end



% Step 4: normalize out-of-plane, fit the volume-averaged motion to volumetric paramter motion
% important key: average across slices to get a volumetric measure of motion from the out-of-planes
for i=1:6
  volslo(:,i)=mean(reshape(newmot(:,i),[zmbdim tdim]));
  [p(i),stand_err(i),mse] = lscov(volslo(:,i)-mean(volslo(:,i)),motion(:,i)-mean(motion(:,i)));
end
normfactor=std(motion)./std(volslo);
normfactor=p;
% do not alter in-planes, these are to be trusted, as-is
normfactor([1 2 6])=1;
% polarity can be inverted w.r.t. volumetric coregistration, but this is not necessarily important for metrics
% be careful if using these for a grid resampling later
%normfactor([1 2 3])=-1*normfactor([1 2 3]); % 3dvolreg x,y,z is inverted w.r.t. 3dWarpDrive in-plane parameters
newmot=newmot.*repmat(normfactor,[tdim*zmbdim 1]);


% apply a Savitsky-Golay filter to adjacent two slices 
% (this, and other steps, assumes that slices have been read in and re-ordered by temporal acquisition order)
% simple hard-coded SG filter for 3 points
% this should be turned off for data with really fast motion (like SimPACE data with motion on only one slice)
for q=1:6
  if (filter_width==5)
    for i=3:tdim*zmbdim-2
      % for 5-point quadratic, coeffs are -3, 12, 17, 12, -3 (norm = 35)
      newmot(i,q)=(-3*newmot(i-2,q)+12*newmot(i-1,q)+17*newmot(i,q)+12*newmot(i+1,q)-3*newmot(i+2,q))/35;
    end
    newmot(2,q)=(newmot(1,q)+2*newmot(2,q)+newmot(3,q))/4;
    newmot(end-1,q)=(newmot(end-2,q)+2*newmot(end-1,q)+newmot(end,q))/4;
  elseif (filter_width==3)
    for i=2:tdim*zmbdim-1
      % for 3-wide, coeffs are 1,2,1 (norm=4)
      newmot(i,q)=(newmot(i-1,q)+2*newmot(i,q)+newmot(i+1,q))/4;
    end
  else
   continue
  end
  newmot(1,q)=(2*newmot(1,q)+newmot(2,q))/3;
  newmot(end,q)=(2*newmot(end,q)+newmot(end-1,q))/3;
end


% Step 5: combine in-plane and out-of-plane measures and remove slice-dependent x,y-shifts introduced by out-of-plane rotations
combined=inplane6;
combined(:,[3 4 5])=newmot(:,[3 4 5]);

slice_order = uniq_acq_order;
slice_pos=((slice_order-0.5)*dz-dz*zmbdim/2);
for i=1:length(combined)
  % temporal slice acquisition number
  z=mod(i-1,zmbdim)+1;
  zpos=slice_pos(z);
  rotmat=convert_rots_into_rotmat(combined(i,4),combined(i,5),0);
  temptrans(i,:)=[0 0 zpos]'-rotmat(1:3,1:3)*[0 0 zpos]';
  % threshold, so correction is not done unless x-rot,y-rot are above 0.2 degrees
  if (abs(combined(i,4))<0.2 && abs(combined(i,5))<0.2)
    temptrans(i,:)=0*temptrans(i,:);
  end
end
% subtract offsets
combined_rotoffsets=combined;
%combined_rotoffsets(:,1:3)=combined_rotoffsets(:,1:3)-temptrans;

% Step 6: fit these to the in-plane data by slice to improve the slice stddev normalization
% offsetfits should be similar for both x and y
if (sub_xy_offsets==0)
 for i=1:zmbdim
  warning off all
  offsetfits(i)=lscov([temptrans(i:zmbdim:end,1); temptrans(i:zmbdim:end,2)],[combined(i:zmbdim:end,1); combined(i:zmbdim:end,2)]);
  warning on all
 end
 % interpolate over poorly-fitted slices
 inds=find(abs(offsetfits)<eps);
 % but only if less than half were ill-fitted
 if (length(inds)<round(zmbdim/2))
  offsetfits=pchip(setxor(inds,1:zmbdim),offsetfits(setxor(inds,1:zmbdim)),1:zmbdim);
 else
  % don't do an offset fit
  %disp('Not doing an offset fit, falling back to trigonometric derived');
  offsetfits=ones(zmbdim,1);
 end
 for i=1:zmbdim
  combined_rotoffsets(i:zmbdim:end,1)=(combined_rotoffsets(i:zmbdim:end,1)-temptrans(i:zmbdim:end,1)*offsetfits(i));
  combined_rotoffsets(i:zmbdim:end,2)=(combined_rotoffsets(i:zmbdim:end,2)-temptrans(i:zmbdim:end,2)*offsetfits(i));
 end
else
  % skip step 6: this may only be useful for highly synthetic data like simPACE data
  % offset correction affects x, y-translations. use in real data is not justified at present
  combined_rotoffsets=combined;
end

% finally, interpolate over out-of-plane parameters again, this can be important if spikes 
% were re-introduced in the fitting, on those untrustworthy slices alone
for i=3:5
  combined_rotoffsets(:,i)=pchip(inputmesh_two,combined_rotoffsets(inputmesh_two,i),1:length(combined_rotoffsets));
end

% method for generating TD-0D-TRU metric: for each volume, take maximum motion across slices
% we don't do that here, opting to save the data output for every slice and let the user use them
% as appropriate. For example, for censoring one would want to censor slices with motion above some threshold,
% but also adjacent slices, if the present slice moved into them over a threshold. I would also use
% out-of-plane motion only, because in-plane is almost completely corrected. There will be small effects on RX
% and B0 field inhomogeneity in less than half the brain, while the spin history effects will be everywhere a voxel
% went out-of-plane by a few hundred microns
[td_slomoco,td_slomocoz]=parallelepiped_jiang(combined_rotoffsets);
fp=fopen('slomoco.TDmetric.txt','w'); fprintf(fp,'%g\n',td_slomoco); fclose(fp);
fp=fopen('slomoco.TDzmetric.txt','w'); fprintf(fp,'%g\n',td_slomocoz); fclose(fp);
% for a volumetric metric of motion corruption, use the max across slices within a volume
fp=fopen('slomoco.volumetric.TDzmetric.txt','w'); fprintf(fp,'%g\n',max(reshape(td_slomocoz,[zmbdim tdim]))); fclose(fp);
fp=fopen('slomoco.volumetric.TDmetric.txt','w'); fprintf(fp,'%g\n',max(reshape(td_slomoco,[zmbdim tdim]))); fclose(fp);
% 3dvolreg motion x,y,z trans are inverted w.r.t. 3dWarpDrive
motion(:,1:3)=-1*motion(:,1:3);
[td_volmoco,tdz_volmoco]=parallelepiped_jiang(motion);
fp=fopen('volmotion.TDmetric.txt','w'); fprintf(fp,'%g\n',td_volmoco); fclose(fp);
fp=fopen('volmotion.TDzmetric.txt','w'); fprintf(fp,'%g\n',tdz_volmoco); fclose(fp);
for i=1:6
  volmotion(:,i)=reshape(repmat(motion(:,i)',[zmbdim 1]),[zmbdim*tdim 1]);
end
% save the 3dvolreg volumetric motion, repeated over slices
fp=fopen('volmotion.repslices.txt','w'); fprintf(fp,'%g\t%g\t%g\t%g\t%g\t%g\n',volmotion'); fclose(fp);
% to see the difference between volumetric motion and slice motion, plot(slomoco-volmotion)
% this is the residual motion left after volumetric correction, or what volmoco misses
% finally, if prospective motion is turned on (Thesen et al 2002), the volumetric motion is essentially
% subtracted from the data, delayed by one volume, so there will be sharp disruptions at the volume boundary
% we could obtain the true (free-space) motion by shifting the volumetric motion by one volume, and adding to slomoco
% this could also be used to improve the edge slice interpolations and could be done earlier, but many sites do not
% use PACE (the Siemens name, each vendor has their own, as far as I know), and is outside the scope of this work
% NOTE, it will be important to know the real free-space motion for applying RX field and B0 inhomogeneity motion corrections
pacemotion=volmotion;
pacemotion(zmbdim+1:end,:)=volmotion(1:end-zmbdim,:);


%fp=fopen('slomoco.combined.6dof.txt','w'); fprintf(fp,'%g\t%g\t%g\t%g\t%g\t%g\n',combined'); fclose(fp);
% convert back to 12-dof
%for i=1:length(combined)
%  rotmat=convert_rots_into_rotmat(-1*combined(i,4),combined(i,5),combined(i,6)); combined12dof(i,:)=[combined(i,1:3) rotmat(:)'];
%end
%fp=fopen('slomoco.combined.12dof.txt','w'); fprintf(fp,'%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n',combined12dof'); fclose(fp);
fp=fopen('slomoco.combined.reg6dof.txt','w'); fprintf(fp,'%g\t%g\t%g\t%g\t%g\t%g\n',combined_rotoffsets'); fclose(fp);
% convert back to 12-dof for use potentially in voxel-specific regression correction
for i=1:length(combined)
  rotmat=convert_rots_into_rotmat(-1*combined_rotoffsets(i,4),combined_rotoffsets(i,5),combined_rotoffsets(i,6)); combined12dof(i,:)=[combined_rotoffsets(i,1:3) rotmat(:)'];
end
fp=fopen('slomoco.combined.reg12dof.txt','w'); fprintf(fp,'%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n',combined12dof'); fclose(fp);

 figure
 subplot(3,1,1);
 plot(meanvox)
 xlim([0 tdim]);
 legend('avg voxel displacement on cube');
 title('volumetric TD (Jiang parallelepiped method)');
 subplot(3,1,2);
 plot(motion)
 xlim([0 tdim]);
 legend('x-trans','y-trans','z-trans','x-rot','y-rot','z-rot');
 title('volumetric params');
 subplot(3,1,3);
 plot(td_slomoco)
 xlim([0 tdim*zmbdim]);
 title('SLOMOCO TD-0D motion estimator - slicewise motion, converted to TD')
 ylabel('TD')
 xlabel('slice*vol number');
 saveas(gcf,'qa_slomoco_metrics.jpg');

 figure
 subplot(2,1,1);
 plot(combined_rotoffsets(:,[3 4 5]))
 xlim([0 tdim*zmbdim]);
 %plot(newmot)
 title(sprintf('out-of-plane params for %s',ep2d_filename));
 legend('z-trans','x-rot','y-rot');
 subplot(2,1,2);
 plot(combined_rotoffsets(:,[1 2 6]))
 xlim([0 tdim*zmbdim]);
 legend('x-trans','y-trans','z-rot');
 title('in-plane params');
 saveas(gcf,'qa_slomoco_motionvectors.jpg');

% figure
% plot(volmotion(:,3),'b')
% hold on
% plot(combined_rotoffsets(:,3),'r')
% set(gca,'fontsize',24)
% ylabel('z-trans (mm)')
% xlabel('slice #')
% title('Motion at slice and volume level, z-translation')
% xlim([0 tdim*zdim])
% saveas(gcf,'volmotion_discrepancy.jpg')


% plot stddev across slices for out-of-plane z-motion. Outermost slices will be bad
%if (slice_timing==1)
%  if (mod(zdim,2))
%    zdims=[[1:2:zdim] [2:2:zdim]];
%  else
%    zdims=[[2:2:zdim] [1:2:zdim]];
%  end
%elseif (slice_timing==2)
%  zdims=[1:zdim];
%elseif (slice_timing==3)
%  zdims=[[1:2:zdim] [2:2:zdim]];
%else
%  disp('unsupport slice timing')
% return
%end
disp('QA script currently assumes we acquire interleaved asc, odds, then evens');
%for i=1:zdim
%  slnum(i)=find(zdims==i);
%  si(i)=std(combined_rotoffsets(i:zdim:end,3));
%end
%figure
%plot(si(slnum))
%title(sprintf('%s z-motion stddev vs spatial slice number',ep2d_filename));
%saveas(gcf,'qa_slomoco_zmot_by_slice.jpg');



