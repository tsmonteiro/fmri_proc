function matlab_retroicor(ep2d_filename,card,resp,M,mask_filename)
% function matlab_retroicor(ep2d_filename,card,resp,M,slice_timing,mask_filename)
% This function performs the 2nd-order RETROICOR algorithm on an epi dataset, using
% input vectors for physiologic data, and returns coupling coefficients and the corrected image
% created for PESTICA distribution
% variance normalization prior to calculation results in statistic values
% use stat=sqrt(sum(im_ca.^2+im_cb.^2,4)), which follows a
% 2*M-th order Chi-square distribution (where M is order of correction, here M=2)

%for simulated data: cardph=unifrnd(0,3.14159*2,zdim,tdim); card=sin(cardph);
if (exist('card','var')==0 | exist('resp','var')==0 | exist('ep2d_filename','var')==0)
  disp('must input three parameters: filename of EPI data, cardiac, and respiratory traces');
  disp('     fourth optional paramter = order of correction (default=2, which is sin(ph)+sin(2*ph)+cos...)');
  disp('     fifth optional paramter = slice acquisition timing string (''alt-asc'' is typical)');
  disp('     sixth optional paramter  = 3D mask file for EPI data');
  return;
end

Opt.Format = 'matrix';
[err, ima, ainfo, ErrMessage]=BrikLoad(ep2d_filename, Opt);
xdim=ainfo.DATASET_DIMENSIONS(1);
ydim=ainfo.DATASET_DIMENSIONS(2);
zdim=ainfo.DATASET_DIMENSIONS(3);
tdim=ainfo.DATASET_RANK(2);
TR=1000*double(ainfo.TAXIS_FLOATS(2));
slice_timing=1000*double(ainfo.TAXIS_OFFSETS);  % ms

[TRsec TRms] = TRtimeunitcheck(TR);
[slice_timing_sec slice_timing_ms] = TRtimeunitcheck(slice_timing);

% check SMS acquisition
[MBacc zmbdim uniq_slice_timing_ms uniq_acq_order] = SMSacqcheck(TRms, zdim, slice_timing_ms);

if (exist('mask_filename','var')~=0)
  [err,mask,minfo,ErrMessage]=BrikLoad(mask_filename, Opt);
  mask = mask(:,:,:,1);  mask(find(mask~=0))=1;
else
  mask=ones(xdim,ydim,zdim);
end
ima=double(reshape(ima,[xdim ydim zdim tdim]));
mask=double(reshape(mask,[xdim ydim zdim]));

% variance normalize the timeseries
variance=squeeze(std(ima,0,4));
variance=variance.*mask;
ima=ima./repmat(variance,[1 1 1 tdim]);

if (exist('M','var')==0)
  M=2;
  disp(sprintf('setting default RETROICOR model order to %d',M));
end

retima=zeros(xdim,ydim,zdim,tdim);
tim_c=zeros(xdim,ydim,zdim);
tim_r=zeros(xdim,ydim,zdim);

card_vols=convert_timeseries_to_slicextime(card,uniq_acq_order);
resp_vols=convert_timeseries_to_slicextime(resp,uniq_acq_order);

for z=1:zdim
  A = [card_vols(z,:)' resp_vols(z,:)' ones(tdim,1)];
  for y=1:ydim
    for x=1:xdim
      if mask(x,y,z)
        vox=squeeze(ima(x,y,z,:));
      
        [B,BINT,R] = regress(vox,A);
        [u,std_u] = lscov(A,vox);
      
        tim_c(x,y,z) = B(1)/std_u(1);
        tim_r(x,y,z) = B(2)/std_u(2);
      
        retima(x,y,z,:)=R+B(3);
      end
    end
  end
end
      
% correct for induced negative values
retima=retima-min(retima(:));

% re-introduce the image variance
retima=retima.*repmat(variance,[1 1 1 tdim]);

% sum squares of t-statistics
tim_card=abs(tim_c);
tim_resp=abs(tim_r);

[tdist_card tdist_card_bin] = sort(tim_card(find(mask)));
[tdist_resp tdist_resp_bin] = sort(tim_resp(find(mask)));

tthr_card = min(tdist_card(tdist_card_bin(find(tdist_card_bin>0.95*sum(mask(:))))));
tthr_resp = min(tdist_resp(tdist_resp_bin(find(tdist_resp_bin>0.95*sum(mask(:))))));

cmask=zeros(xdim,ydim,zdim);
cmask(find(tim_card>tthr_card))=1;
rmask=zeros(xdim,ydim,zdim);
rmask(find(tim_resp>tthr_resp))=1;

% keep same format as input data
ainfo.BRICK_TYPES=ainfo.BRICK_TYPES(1)*ones(1,tdim); % 1=short, 3=float
ainfo.BRICK_STATS = []; %automatically set
ainfo.BRICK_FLOAT_FACS = [];%automatically set
ainfo.BRICK_LABS = [];
ainfo.BRICK_KEYWORDS = [];
OptOut.Scale = 0;
OptOut.OverWrite= 'y';
OptOut.Prefix = strcat('../',ep2d_filename(1:strfind(ep2d_filename,'+')-1),'.retroicor');
OptOut.verbose = 0;
[err,ErrMessage,InfoOut]=WriteBrik(retima,ainfo,OptOut);

ainfo.BRICK_TYPES=3; %1=short, 3=float
ainfo.BRICK_STATS = []; %automatically set
ainfo.BRICK_FLOAT_FACS = [];%automatically set
ainfo.DATASET_RANK(2)=1;
ainfo.TAXIS_NUMS(1) = 1;
% write out cardiac coupling maps
OptOut.Prefix = 'coupling_ret_card';
[err,ErrMessage,InfoOut]=WriteBrik(tim_card,ainfo,OptOut);
% write out respiratory coupling maps (2 IRFs)
OptOut.Prefix = 'coupling_ret_resp';
[err,ErrMessage,InfoOut]=WriteBrik(tim_resp,ainfo,OptOut);

save('masks_PESTICA4.mat','cmask','rmask')
