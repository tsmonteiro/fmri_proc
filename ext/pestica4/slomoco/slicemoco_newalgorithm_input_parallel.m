function slicemoco_newalgorithm_input(ep2d_filename,mask_filename,filestr,CARD,RESP)
%function slicemoco_newalgorithm_input(ep2d_filename,mask_filename,filestr,CARD,RESP)
    
%disp('warning, motionparams must be in form x y z then 9-element rotation matrix on oneline');
Opt.Format = 'matrix';
[err, im, ainfo, ErrMessage]=BrikLoad(ep2d_filename, Opt);
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

% read moco parameters
[transmat1d_zt,fpparams_6dof]=read_1dmat_zt(filestr);   % rotmat1d_zt = [zdim tdim 1x12]
% tranmat_zt = convert_1dmat_into_tranmat(transmat1d_zt); % tranmat_zt = [zdim tdim 3 x 4 ]
tranmat_zt = convert_1dmat_into_tranmatarray(transmat1d_zt);% tranmat_zt = reshape(tranmat_zt,[zmbdim tdim 3 4]); % (v1) debugged

% if second col, then mult that coord by -1
% if row rx~=1, then have to swap coord according to the row numbers
[rx,cx]=find(ainfo.Orientation=='R');  % for axial == 1
[ry,cy]=find(ainfo.Orientation=='A');  % for axial == 2
[rz,cz]=find(ainfo.Orientation=='I');  % for axial == 3

% load mask file
[err, im_mask, minfo, ErrMessage]=BrikLoad(mask_filename, Opt);
ep2d_filename_pr=ep2d_filename(1:strfind(ep2d_filename,'+')-1);
im=reshape(im,[xdim ydim zdim tdim]);
im_mask=reshape(im_mask,[xdim ydim zdim]);

% slice timing should be provided with milli second unit
physioflag=0; pmuflag=0; pesticaflag=0;
if exist('CARD','var') && exist('RESP','var')
  physioflag=1;
  if isfield(CARD,'irf_coeffs') && isfield(RESP,'irf_coeffs')
    pmuflag=1;
    % define order of Fourier-decomposition
    M = size(CARD.irf_coeffs,1)/2;
  else
    pesticaflag=1;  
    card_vols=convert_timeseries_to_slicextime(CARD,uniq_acq_order);
    resp_vols=convert_timeseries_to_slicextime(RESP,uniq_acq_order);
  end
end

tic
% pre define variables
im_moco=zeros(xdim,ydim,zdim,tdim);
if pmuflag
  dof = double(tdim-(6+1)-12);
  studt=zeros(xdim,ydim,zdim,12+6);
elseif  pesticaflag
  dof = double(tdim-(2+1)-12);
  studt=zeros(xdim,ydim,zdim,12+2);
else
  dof = double(tdim-(1)-12);
  studt=zeros(xdim,ydim,zdim,12);
end
trend_funct=(1:tdim)/tdim;
trend_funct=trend_funct-mean(trend_funct);

parfor k= 1:zdim
  % z definition
  z=((k-(zdim/2.))*dz);
  
  % slice above
  if (k<zdim) || (MB > 1)
    adj_sup_z=((k+1-(zdim/2.))*dz);
  end
  % slice below
  if (k>1) || (MB > 1)
    adj_inf_z=((k-1-(zdim/2.))*dz);
  end
  
  % MB consideration
  if mod(k,zmbdim); kmb = mod(k-1,zmbdim)+1;  % (v1) debugged
  else; kmb = zmbdim; 
  end
  
  if pmuflag 
  % physio regressors - these have already been disassembled into slice-acquisition order
    physio=[(CARD.irf_coeffs(1,1:M)*sin((1:M)'*CARD.phz_slc(:,k)')+CARD.irf_coeffs(1,M+1:end)*cos((1:M)'*CARD.phz_slc(:,k)'))' ...
            (CARD.irf_coeffs(2,1:M)*sin((1:M)'*CARD.phz_slc(:,k)')+CARD.irf_coeffs(2,M+1:end)*cos((1:M)'*CARD.phz_slc(:,k)'))' ...
            (CARD.irf_coeffs(3,1:M)*sin((1:M)'*CARD.phz_slc(:,k)')+CARD.irf_coeffs(3,M+1:end)*cos((1:M)'*CARD.phz_slc(:,k)'))' ...
            (CARD.irf_coeffs(4,1:M)*sin((1:M)'*CARD.phz_slc(:,k)')+CARD.irf_coeffs(4,M+1:end)*cos((1:M)'*CARD.phz_slc(:,k)'))' ...
            (RESP.irf_coeffs(1,1:M)*sin((1:M)'*RESP.phz_slc(:,k)')+RESP.irf_coeffs(1,M+1:end)*cos((1:M)'*RESP.phz_slc(:,k)'))' ...
            (RESP.irf_coeffs(2,1:M)*sin((1:M)'*RESP.phz_slc(:,k)')+RESP.irf_coeffs(2,M+1:end)*cos((1:M)'*RESP.phz_slc(:,k)'))'];
  elseif pesticaflag
    physio = [card_vols(kmb,:)' resp_vols(kmb,:)'];
  else
    physio=[];
  end
  
  % rotmat_t_sup
  if kmb == zmbdim
    if MB > 1
      k_sup=1;
    else
      k_sup=0;
    end
  else 
    k_sup = kmb + 1;
  end
  
  if find(uniq_acq_order==(k_sup)) > find(uniq_acq_order==(kmb)) % if upper slice is acquired later within TR
    half1st_sup=1;
  else
    half1st_sup = 0;
  end
  
  % rotmat_t_inf
  if kmb == 1
    if MB > 1
      k_inf = zmbdim;
    else
      k_inf = 0;
    end
  else 
    k_inf = kmb-1;
  end
  
  if find(uniq_acq_order==(k_inf)) > find(uniq_acq_order==(kmb)) % if lower slice is acquired later within TR
    half1st_inf=1;
  else
    half1st_inf = 0;
  end                
  
  for i=1:xdim
    x=((i-xdim/2.)*dx);
    for j=1:ydim
      y=((j-ydim/2.)*dy);

      if im_mask(i,j,k)
        vox=squeeze(im(i,j,k,:));
        
        % change coordinate axes based on ainfo.Orientation. The AFNI coords are all relative to the RAI coordinate system
        % while z's is always slice acquisition axis and so on.
        a=[x y z];
        x=a(rx); y=a(ry); z=a(rz);
        % invert any axes needing inversion
        x=x*((-2*cx)+3);  y=y*((-2*cy)+3);  z=z*((-2*cz)+3);

        %calculate motion of this voxel
        xyzprime       = zeros(tdim,3);
        zdelta_adj_sup = zeros(tdim,1); 
        zdelta_adj_inf = zeros(tdim,1);  
        
        for t=1:tdim
          xyzprime(t,:) = tranmat_zt(kmb,t).R*[x y z]'+tranmat_zt(kmb,t).T;
          
          % calculate the motion of the adjacent superior slice's voxel into this voxel
          % using the motion timed appropriately for the adjacent slice
          if k_sup 
            a = tranmat_zt(k_sup,t).R*[x y z]'+tranmat_zt(k_sup,t).T;
            zdelta_adj_sup(t) = a(3) - adj_sup_z;
          end

          % calculate the motion of the adjacent inferior slice's voxel into this voxel
          if k_inf
            a = tranmat_zt(k_inf,t).R*[x y z]'+tranmat_zt(k_inf,t).T;
            zdelta_adj_inf(t) = a(3) - adj_inf_z;
          end
        end

        % rotation matrix should always be "special orthogonal", so RR^T=I (R times transpose of R equals identity matrix, also R'=inv(R)) and det(R)=1.  
        % Also, R is normalized so squares of elements in any row/column sum to 1, and dot product of any pair of rows or pair of columns is 0. Finally, 
        % rows of R represent coordinates in the original space of unit vectors along the coord axes of the rotated space (and vice versus for columns)
        xdelta=xyzprime(:,1)-x;  
        ydelta=xyzprime(:,2)-y;  
        zdelta=xyzprime(:,3)-z;
        zdelta_lag=[0; zdelta(1:end-1)];

        %create array of motion regressors
        z1st=[0; diff(zdelta)];

        % if edge slice, set to 1st derivative
        if MB == 1;
          if k==1;        zdelta_adj_inf=z1st;  end
          if k==zmbdim;   zdelta_adj_sup=z1st;  end
        end

        if half1st_inf
          zdelta_adj_inf = [0; zdelta_adj_inf(1:end-1)];
        end
        
        if half1st_sup
          zdelta_adj_sup = [0; zdelta_adj_sup(1:end-1)];
        end

        zdelta_adj_inf_lag = [0; zdelta_adj_inf(1:end-1)];
        zdelta_adj_sup_lag = [0; zdelta_adj_sup(1:end-1)];

        mr=[xdelta xdelta.^2 ydelta ydelta.^2 zdelta zdelta.^2 zdelta_lag zdelta_lag.^2 ...
            zdelta_adj_inf zdelta_adj_inf.^2 zdelta_adj_sup  zdelta_adj_sup.^2 zdelta_adj_inf_lag zdelta_adj_sup_lag];   
          
        % de-mean moco params before adding into model
        mvox=mean(vox);
        vox=vox-mvox;

        % include trend in fit
        A=[trend_funct(1:tdim)' mr physio];
        
        warning off all
        % due to propensity for regressors to covary (and thus solution be ill-conditioned), use more robust least-squares methods
        [p,stand_err,mse] = lscov(A,vox);
        studt(i,j,k,1:length(p)-1)=p(2:end)./stand_err(2:end);
        warning on all
        % don't remove trend or mean in corrected data
        p(1)=0;
        im_moco(i,j,k,:)=mvox+(vox-A*p);
      end
    end
  end
  if k == zdim
    fprintf([num2str(k) '\n']);
  elseif k==1
    fprintf(['slice ' num2str(k) '.']);
  else
    fprintf([num2str(k) '.']);
  end
end
toc

% use same ainfo as we read in from the original data file
ainfo.BRICK_TYPES=ainfo.BRICK_TYPES(1)*ones(1,tdim); % 1=short, 3=float
ainfo.BRICK_STATS = [];           % automatically set
ainfo.BRICK_FLOAT_FACS = [];      % automatically set
ainfo.BRICK_LABS = [];
ainfo.BRICK_KEYWORDS = [];
OptOut.Scale = 0;                 % no scaling
OptOut.OverWrite= 'y';            % overwrite if exists
if pmuflag
  OptOut.Prefix = strcat(ep2d_filename(1:strfind(ep2d_filename,'+')-1),'.slomoco_pmu');
elseif pesticaflag
  OptOut.Prefix = strcat(ep2d_filename(1:strfind(ep2d_filename,'+')-1),'.slomoco_pestica');
else
  OptOut.Prefix = strcat(ep2d_filename(1:strfind(ep2d_filename,'+')-1),'.slomoco');
end
  
OptOut.verbose = 0;
[err,ErrMessage,InfoOut]=WriteBrik(im_moco,ainfo,OptOut);

ErrMessage
InfoOut
err

% write out noise coupling maps (4 IRFs)
ainfo.BRICK_TYPES = 3*ones(size(studt,4)+1,1); %1=short, 3=float
%ainfo.DATASET_DIMENSIONS(4)=5;
ainfo.DATASET_RANK(2)=size(studt,4)+1;
ainfo.TAXIS_NUMS(1) = size(studt,4)+1;
if pmuflag
  OptOut.Prefix = 'slomoco_pmu_coupling_maps';
elseif pesticaflag
  OptOut.Prefix = 'slomoco_pestica_coupling_maps';
else
  OptOut.Prefix = 'slomoco_coupling_maps';
end

temp1=sqrt(sum(studt.^2,4));
tim_moco=zeros(size(studt,1),size(studt,2),size(studt,3),size(studt,4)+1);
tim_moco(:,:,:,2:size(studt,4)+1)=studt;
tim_moco(:,:,:,1)=temp1;
% since NAN crashs afni display, set NAN to 0 (may change later). Jian
tim_moco(find(isnan(tim_moco)==1)) = 0.0;

[err,ErrMessage,InfoOut]=WriteBrik(tim_moco,ainfo,OptOut);

