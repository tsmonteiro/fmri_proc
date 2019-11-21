function [CARD RESP] = retroicor_get_irf(ep2d_filename,CARD,RESP,M,mask_filename)
%function [irf_resp,irf_card,coeffs_card,coeffs_resp]=retroicor_get_irf(ep2d_filename,cardph,respph,M,slice_timing,mask_filename)
% new matlab version of retroicor which tries to determine impulse response functions to 
% (1) heart beat and (2) respiration.  It does this by fitting sines and cosines of phase
% to raw data and finding the optimal combination as the limit of abs(tscore) increases
%
%  Output is M*2 timeseries of IRF shapes for card and resp each and the coefficients needed to
%   reproduce these using the parallel measured, RETROICOR physiologic phases
%
% The intent is to run this iteratively after regressing out the detected noise to get smaller 
% and smaller components (which will likely have different i.r.f.'s) of the physiologic noise

%for simulated data: cardph=unifrnd(0,3.14159*2,zdim*tdim,1); card=sin(cardph);
if (exist('CARD','var')==0 | exist('RESP','var')==0 | exist('ep2d_filename','var')==0)
  disp('must input three parameters: filename of EPI data, phases of cardiac and respiratory traces');
  disp('     fourth optional paramter = order of correction (default=5, or sin([1:M]*phase)+cos([1:M]*phase)');
  disp('     fifth optional paramter = slice acquisition timing string (''alt-asc'' is typical)');
  disp('     sixth optional paramter  = mask file for EPI data');
  return;
end

% define order of Fourier-decomposition
if (exist('M','var')==0)
  M=2;
  disp(sprintf('setting default order of decomp to %d',M));
end

Opt.Format = 'matrix';
[err, ima, ainfo, ErrMessage]=BrikLoad(ep2d_filename, Opt);
xdim=ainfo.DATASET_DIMENSIONS(1);
ydim=ainfo.DATASET_DIMENSIONS(2);
zdim=ainfo.DATASET_DIMENSIONS(3);
tdim=ainfo.DATASET_RANK(2);

if (exist('mask_filename','var')~=0)
  [err,mask,minfo,ErrMessage]=BrikLoad(mask_filename, Opt);
  mask = mask(:,:,:,1); mask(find(mask~=0))=1;
else
  mask=ones(xdim,ydim,zdim);
end
ima=double(reshape(ima,[xdim ydim zdim tdim]));
mask=double(reshape(mask,[xdim ydim zdim]));

% variance normalizing is unnecessary, although distributions will be easier to understand
variance=squeeze(std(ima,0,4));
ima=ima./repmat(variance,[1 1 1 tdim]);

im_ca=zeros(xdim,ydim,zdim,M);
im_cb=zeros(xdim,ydim,zdim,M);
im_ra=zeros(xdim,ydim,zdim,M);
im_rb=zeros(xdim,ydim,zdim,M);
tim_ca=zeros(xdim,ydim,zdim,M);
tim_cb=zeros(xdim,ydim,zdim,M);
tim_ra=zeros(xdim,ydim,zdim,M);
tim_rb=zeros(xdim,ydim,zdim,M);

im_card=zeros(xdim,ydim,zdim);
im_resp=zeros(xdim,ydim,zdim);
tim_card=zeros(xdim,ydim,zdim);
tim_resp=zeros(xdim,ydim,zdim);
dof=double(tdim-(M*2*2+1));

disp(['Start voxelwise fitting with Fourier Series physio Model (M=' num2str(M) ')'])
for z=1:zdim
  A=[sin((1:M)'*CARD.phz_slc(:,z)')' cos((1:M)'*CARD.phz_slc(:,z)')' sin((1:M)'*RESP.phz_slc(:,z)')' cos((1:M)'*RESP.phz_slc(:,z)')' ones(tdim,1)];
  % QR-decomposition
  [Q,R] = qr(A,0);
  for x=1:xdim
    for y=1:ydim
      if mask(x,y,z)
        vox=squeeze(ima(x,y,z,:));
        % detrending voxel timecourse prior to regression increases tscore
        vox=detrend(vox);
        % matrix division to get amplitudes of fits
        p = R\(Q'*vox);
        im_ca(x,y,z,:)=p(1:M);
        im_cb(x,y,z,:)=p(M+1:M*2);
        im_ra(x,y,z,:)=p(M*2+1:M*3);
        im_rb(x,y,z,:)=p(M*3+1:M*4);
        % and t-scores
        % residuals
        res=vox-A*p;
        % mean square error
        ms=res'*res./(dof+1);
        Rinv=pinv(R);
        % error covariance matrix
        covb=(Rinv*Rinv')*ms;
        % get standardized error (variance that is not explained by fit to design matrix)
        stand_err=sqrt(diag(covb));
        tim_ca(x,y,z,:)=p(1:M)./stand_err(1:M);
        tim_cb(x,y,z,:)=p(M+1:M*2)./stand_err(M+1:M*2);
        tim_ra(x,y,z,:)=p(M*2+1:M*3)./stand_err(M*2+1:M*3);
        tim_rb(x,y,z,:)=p(M*3+1:M*4)./stand_err(M*3+1:M*4);
      end
    end
  end
end
disp(['fitting has been done. '])

% identify M*2 shapes of IRFs for each of card,resp using the following procedure:
% 1. tscore threshold, the distribution of tim_card,tim_resp look like M*2th order chi
%    in order to identify the threshold, fit to uniform dist of rand phase data and plot these 
%    distributions to determine the null hypothesis and an appropriate threshold for significance
%   RESULTS: Tscore threshold is order dependent, 
%       for M=3, use 4.6, for M=4, use 5.3, for M=5, use 5.6
%       for M between 2 and 10, good threshold approximately equiv to 0.52*M+3.03
% 2. In the thresholded voxels, create shapes based on the fit im_ca, im_cb, im_ra, im_rb values,
%    weighting is unimportant since voxels were all variance normalized, so fit coefficients
%    are the appropriate weights
% 3. PCA separately for each slice
% 4. PCA across the slice-determined shapes
% PCA is appropriate since these are constructed from orthogonal components of Fourier series

tim_card=sqrt(sum(tim_ca.^2+tim_cb.^2,4));
tim_resp=sqrt(sum(tim_ra.^2+tim_rb.^2,4));

[tdist_card tdist_card_bin] = sort(tim_card(find(mask)));
[tdist_resp tdist_resp_bin] = sort(tim_resp(find(mask)));

tthr_card = min(tdist_card(tdist_card_bin(find(tdist_card_bin>0.95*sum(mask(:))))));
tthr_resp = min(tdist_resp(tdist_resp_bin(find(tdist_resp_bin>0.95*sum(mask(:))))));

cmask=zeros(xdim,ydim,zdim);
cmask(find(tim_card>tthr_card))=1;
rmask=zeros(xdim,ydim,zdim);
rmask(find(tim_resp>tthr_resp))=1;

% find the impulse response
cdata=linspace(0,2*pi,100);
rdata=linspace(-pi,pi,100);
pct=zeros(z,M*2,length(cdata));
prt=zeros(z,M*2,length(cdata));
resp_impulse = zeros(zdim,100);
card_impulse = zeros(zdim,100);

for z=1:zdim
  % average the impulse response functions per slice
  clear ct rt
  ct = zeros(sum(sum(cmask(:,:,z))),100);
  rt = zeros(sum(sum(cmask(:,:,z))),100);

  c=0;r=0;
  for x=1:xdim
    for y=1:ydim
      if (cmask(x,y,z)==1)
        c=c+1;
        ct(c,:)=squeeze(im_ca(x,y,z,:))'*sin((1:M)'*cdata)+squeeze(im_cb(x,y,z,:))'*cos((1:M)'*cdata);
      end
      if (rmask(x,y,z)==1)
        r=r+1;
        rt(r,:)=squeeze(im_ra(x,y,z,:))'*sin((1:M)'*rdata)+squeeze(im_rb(x,y,z,:))'*cos((1:M)'*rdata);
      end
    end
  end
  
  % check for underfill (not enough significant voxels on this slice)
  if (r<M*2)
    rt(r+1:M*2,:)=zeros(M*2-r,length(rdata));
  end
  if (c<M*2)
    ct(c+1:M*2,:)=zeros(M*2-c,length(cdata));
  end
  % respiration is simple, use mean and assume one impulse only
  resp_impulse(z,:)=mean(rt);
  card_impulse(z,:)=mean(ct);
  resp_impulseCN(1,z)=r;
  card_impulseCN(1,z)=c;
  
  % resp is made ambiguous by the respiratory phase issue which will have to be dealt with)
  % PCA the impulses across the slice, saving only M*2 comps (orthogonal fourier)
  [eg,ev,comps]=pcsquash(ct,M*2);
  pct(z,:,:)=comps;
  [eg,ev,comps]=pcsquash(rt,M*2);
  prt(z,:,:)=comps;
end

CARD.irf_avg = card_impulseCN*card_impulse/sum(card_impulseCN);
RESP.irf_avg = resp_impulseCN*resp_impulse/sum(resp_impulseCN);

t=reshape(pct,[zdim*M*2 100]);
[eg,ev,comps]=pcsquash(t,M*2);
CARD.irf_comps=comps;
t=reshape(prt,[zdim*M*2 100]);
[eg,ev,comps]=pcsquash(t,M*2);
RESP.irf_comps=comps;

% find coefficients
clear coeffs_resp coeffs_card
A=[sin((1:M)'*cdata)' cos((1:M)'*cdata)'];
[Q,R] = qr(A,0);
for i=1:M*2
  coeffs_card(i,:)=R\(Q'*CARD.irf_comps(i,:)');
end

A=[sin((1:M)'*rdata)' cos((1:M)'*rdata)'];
[Q,R] = qr(A,0);
for i=1:M*2
  coeffs_resp(i,:)=R\(Q'*RESP.irf_comps(i,:)');
end

CARD.irf_coeffs = coeffs_card;
RESP.irf_coeffs = coeffs_resp;

save impulse_responses.mat mask cmask rmask CARD RESP

% Based on group analysis of above (see IRF-RET paper), cardiac 4 impulse responses present, respiration has 2

