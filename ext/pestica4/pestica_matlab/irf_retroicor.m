function [tim_card,tim_resp,retima]=irf_retroicor(ep2d_filename,CARD, RESP, mask_filename, ref_file)
%function [tim_card,tim_resp,retima]=irf_retroicor(ep2d_filename,cardph,respph,coeffs_card,coeffs_resp,slice_timing,mask_filename,ref_file)
% Fits Impulse Response Function (IRF) derived from data to the data, so as to perform
% the RETROICOR method without spurious removal of degrees of freedom
% Output is coupling maps (t-maps) and corrected data

%for simulated data: cardph=unifrnd(0,3.14159*2,zdim,tdim); card=sin(cardph);
if (exist('CARD','var')==0 | exist('RESP','var')==0 | exist('ep2d_filename','var')==0)
  disp('must two three parameters: filename of EPI data, CARD and RESP with IRFs');
  disp('     fourth optional paramter = order of correction (default=5, or sin([1:M]*phase)+cos([1:M]*phase)');
  disp('     fifth optional paramter = slice acquisition timing string (''alt-asc'' is typical)');
  disp('     sixth optional paramter  = ANALYZE format 3D mask file for EPI data');
  return;
end

refflag=0;
if (exist(ref_file,'file')~=0)
  disp(sprintf('Using %s as reference file',ref_file));
  refflag=1;
  ref=textread(ref_file);
  if (size(ref,2)>size(ref,1))
    ref=ref';
  end
end

if (exist('M','var')==0)
  M=2;
  disp(sprintf('setting default RETROICOR model order to %d',M));
end

Opt.Format = 'matrix';
[err, ima, ainfo, ErrMessage]=BrikLoad(ep2d_filename, Opt);
xdim=ainfo.DATASET_DIMENSIONS(1);
ydim=ainfo.DATASET_DIMENSIONS(2);
zdim=ainfo.DATASET_DIMENSIONS(3);
tdim=ainfo.DATASET_RANK(2);

if (exist('mask_filename','var')~=0)
  [err,mask,minfo,ErrMessage]=BrikLoad(mask_filename, Opt);
  mask = mask(:,:,:,1);  mask(find(mask~=0))=1;
else
  mask=ones(xdim,ydim,zdim);
end

ima=double(reshape(ima,[xdim ydim zdim tdim]));
mask=double(reshape(mask,[xdim ydim zdim]));

% define order of Fourier-decomposition
M = size(CARD.irf_coeffs,1)/2;

im_card  = zeros(xdim,ydim,zdim,4);
im_resp  = zeros(xdim,ydim,zdim,2);
tim_card = zeros(xdim,ydim,zdim,4);
tim_resp = zeros(xdim,ydim,zdim,2);
% dof is now 6 for 4 card and 2 resp IRFs

dof=double(tdim-(6+1));
retima = zeros(xdim,ydim,zdim,tdim);

for z=1:zdim
  if (refflag)
    A=[(CARD.irf_coeffs(1,1:M)*sin((1:M)'*CARD.phz_slc(:,z)')+CARD.irf_coeffs(1,M+1:end)*cos((1:M)'*CARD.phz_slc(:,z)'))' ...
       (CARD.irf_coeffs(2,1:M)*sin((1:M)'*CARD.phz_slc(:,z)')+CARD.irf_coeffs(2,M+1:end)*cos((1:M)'*CARD.phz_slc(:,z)'))' ...
       (CARD.irf_coeffs(3,1:M)*sin((1:M)'*CARD.phz_slc(:,z)')+CARD.irf_coeffs(3,M+1:end)*cos((1:M)'*CARD.phz_slc(:,z)'))' ...
       (CARD.irf_coeffs(4,1:M)*sin((1:M)'*CARD.phz_slc(:,z)')+CARD.irf_coeffs(4,M+1:end)*cos((1:M)'*CARD.phz_slc(:,z)'))' ...
       (RESP.irf_coeffs(1,1:M)*sin((1:M)'*RESP.phz_slc(:,z)')+RESP.irf_coeffs(1,M+1:end)*cos((1:M)'*RESP.phz_slc(:,z)'))' ...
       (RESP.irf_coeffs(2,1:M)*sin((1:M)'*RESP.phz_slc(:,z)')+RESP.irf_coeffs(2,M+1:end)*cos((1:M)'*RESP.phz_slc(:,z)'))' ref ones(tdim,1)];
  else
    A=[(CARD.irf_coeffs(1,1:M)*sin((1:M)'*CARD.phz_slc(:,z)')+CARD.irf_coeffs(1,M+1:end)*cos((1:M)'*CARD.phz_slc(:,z)'))' ...
       (CARD.irf_coeffs(2,1:M)*sin((1:M)'*CARD.phz_slc(:,z)')+CARD.irf_coeffs(2,M+1:end)*cos((1:M)'*CARD.phz_slc(:,z)'))' ...
       (CARD.irf_coeffs(3,1:M)*sin((1:M)'*CARD.phz_slc(:,z)')+CARD.irf_coeffs(3,M+1:end)*cos((1:M)'*CARD.phz_slc(:,z)'))' ...
       (CARD.irf_coeffs(4,1:M)*sin((1:M)'*CARD.phz_slc(:,z)')+CARD.irf_coeffs(4,M+1:end)*cos((1:M)'*CARD.phz_slc(:,z)'))' ...
       (RESP.irf_coeffs(1,1:M)*sin((1:M)'*RESP.phz_slc(:,z)')+RESP.irf_coeffs(1,M+1:end)*cos((1:M)'*RESP.phz_slc(:,z)'))' ...
       (RESP.irf_coeffs(2,1:M)*sin((1:M)'*RESP.phz_slc(:,z)')+RESP.irf_coeffs(2,M+1:end)*cos((1:M)'*RESP.phz_slc(:,z)'))' ones(tdim,1)];
  end
  % QR-decomposition
  [Q,R] = qr(A,0);
  
  for x=1:xdim
    for y=1:ydim
      if (mask(x,y,z)==0)
        continue;
      end
      vox=squeeze(ima(x,y,z,:));
      if (mask(x,y,z)==0)
        retima(x,y,z,:)=vox;
        continue;
      end
      % matrix division to get amplitudes of fits
      p = R\(Q'*vox);
      im_card(x,y,z,:)=p(1:4);
      im_resp(x,y,z,:)=p(5:6);
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
      tim_card(x,y,z,:)=p(1:4)./stand_err(1:4);
      tim_resp(x,y,z,:)=p(5:6)./stand_err(5:6);
      % do not remove mean from data
      p(end)=0;
      if (refflag)
        p(end-1)=0;
      end
      % corrected data
      retima(x,y,z,:)=vox-A*p;
    end
  end
end

% correct for induced negative values
retima=retima-min(retima(:));

% use same ainfo as we read in from the original data file
ainfo.BRICK_TYPES=ainfo.BRICK_TYPES(1)*ones(1,tdim); % 1=short, 3=float
ainfo.BRICK_STATS = [];           % automatically set
ainfo.BRICK_FLOAT_FACS = [];      % automatically set
ainfo.BRICK_LABS = [];
ainfo.BRICK_KEYWORDS = [];
OptOut.Scale = 0;                 % no scaling
OptOut.OverWrite= 'y';            % overwrite if exists
OptOut.Prefix = strcat(ep2d_filename(1:strfind(ep2d_filename,'+')-1),'.retroicor_pmu');
OptOut.verbose = 0;
[err,ErrMessage,InfoOut]=WriteBrik(retima,ainfo,OptOut);

ainfo.BRICK_STATS = []; %automatically set
ainfo.BRICK_FLOAT_FACS = [];%automatically set
% write out cardiac coupling maps (4 IRFs), preceeded by sum of squares of all IRFs
ainfo.BRICK_TYPES = 3*ones(1,size(tim_card,4)+1); %1=short, 3=float
ainfo.DATASET_RANK(2)=size(tim_card,4)+1;
ainfo.TAXIS_NUMS(1) = size(tim_card,4)+1;
OptOut.Prefix = 'coupling_irfret_card_pmu';
temp1=sqrt(sum(tim_card.^2,4));
temp2=tim_card;
tim_card=zeros(size(tim_card,1),size(tim_card,2),size(tim_card,3),size(tim_card,4)+1);
tim_card(:,:,:,2:size(tim_card,4))=temp2;
tim_card(:,:,:,1)=temp1;
[err,ErrMessage,InfoOut]=WriteBrik(tim_card,ainfo,OptOut);
% write out respiratory coupling maps (2 IRFs)
ainfo.BRICK_TYPES = 3*ones(1,size(tim_resp,4)+1); %1=short, 3=float
ainfo.DATASET_RANK(2)=size(tim_resp,4)+1;
ainfo.TAXIS_NUMS(1) = size(tim_resp,4)+1;
OptOut.Prefix = 'coupling_irfret_resp_pmu';
temp1=sqrt(sum(tim_resp.^2,4));
temp2=tim_resp;
tim_resp=zeros(size(tim_resp,1),size(tim_resp,2),size(tim_resp,3),size(tim_resp,4)+1);
tim_resp(:,:,:,2:size(tim_resp,4))=temp2;
tim_resp(:,:,:,1)=temp1;
[err,ErrMessage,InfoOut]=WriteBrik(tim_resp,ainfo,OptOut);

