function physio_qa(ep2d_filename,pmuflag)
%function physio_qa(ep2d_filename,pmuflag)
% script reads in physio and fit data in local directory and pestica/ subdirectory
% plot histograms of periodicities, plot IRFs, report number of voxels significantly coupled
%
% make AFNI shots of mask, coregs, coupling maps overlain, skip first-pass results unless I include option to only do first-pass irfret
% conclude with eog or other display tool to open jpegs I've included in the PESTICA distribution for comparison purposes
% finally also release a new tutorial

if (exist('pmuflag','var')==0)
  pmuflag=0;
end

% BrikInfo only works on AFNI BRIK format files
Opt.Format = 'matrix';
[err, ima, ainfo, ErrMessage]=BrikLoad(ep2d_filename, Opt);
zdim=ainfo.DATASET_DIMENSIONS(3);
TR=double(ainfo.TAXIS_FLOATS(2));
slice_timing=double(ainfo.TAXIS_OFFSETS);  % ms

% correct milli second unit of TR and slice timing
[TRsec TRms] = TRtimeunitcheck(TR);
[slice_timing_sec slice_timing_ms] = TRtimeunitcheck(slice_timing);

% check SMS acquisition
[MBacc zmbdim] = SMSacqcheck(TRms, zdim, slice_timing_ms);

% apparently its not uncommon for DELTA to be negative on one or more axes, but don't know why that would be...
voxsize=abs(prod(ainfo.DELTA));
% its not uncommon for users to specify milliseconds for their TR, so check it here

[TRsec TRms] = TRtimeunitcheck(TR);

disp(sprintf('using %d slices, TR=%f seconds',zdim,TRsec));

if (pmuflag==1)
  % Siemens PMU is sampled at 50Hz
  f_s=50.0;
else
  f_s=zmbdim/TRsec;
end

couplingfnamec = ['coupling_irfret_card_pmu'];
couplingfnamer = ['coupling_irfret_resp_pmu'];
if (exist([couplingfnamec '+orig.BRIK'],'file')) && (exist([couplingfnamer '+orig.BRIK'],'file'))
  [err, ima, ainfo, ErrMessage]=BrikLoad([couplingfnamec '+orig'], Opt);
  ima=ima(:,:,:,1); ima=ima(find(ima~=0)); [n,x]=hist(ima,200); n=raylfit(ima); ctthresh=round(15*n);
  [err, ima, ainfo, ErrMessage]=BrikLoad([couplingfnamer '+orig'], Opt);
  ima=ima(:,:,:,1); ima=ima(find(ima~=0)); [n,x]=hist(ima,200); n=raylfit(ima); rtthresh=round(15*n);

  fp=fopen(['coupling_pmu_thresholds.txt'],'w');
  fprintf(fp,'%d\n',[ctthresh rtthresh]);
  fclose(fp);
end
  
couplingfnamec = ['coupling_ret_card_pestica'];
couplingfnamer = ['coupling_ret_resp_pestica'];
if (exist([couplingfnamec '+orig.BRIK'],'file')) && (exist([couplingfnamer '+orig.BRIK'],'file'))
  [err, ima, ainfo, ErrMessage]=BrikLoad([couplingfnamec '+orig'], Opt);
  ima=ima(:,:,:,1); ima=ima(find(ima~=0)); [n,x]=hist(ima,200); n=raylfit(ima); ctthresh=round(15*n);
  [err, ima, ainfo, ErrMessage]=BrikLoad([couplingfnamer '+orig'], Opt);
  ima=ima(:,:,:,1); ima=ima(find(ima~=0)); [n,x]=hist(ima,200); n=raylfit(ima); rtthresh=round(15*n);

  fp=fopen(['coupling_pestica_thresholds.txt'],'w');
  fprintf(fp,'%d\n',[ctthresh rtthresh]);
  fclose(fp);
end
  

if (pmuflag==1)
  load impulse_responses.mat
  disp(sprintf('Mean Cardiac Rate: %f, Mean Respiratory Rate: %f',length(CARD.prd)/CARD.t(end),length(RESP.prd)/RESP.t(end)))  ;
  h = figure('visible','off');
  subplot(2,1,1);
  hist(60./(diff(CARD.tptrace)),15)
  set(gca,'fontsize',16)
  ylabel('Count');
  xlabel('Beats per Minute');
  title(sprintf('%d cardiac voxels\n%d respiratory voxels',length(find(cmask==1)), length(find(rmask==1))));
  subplot(2,1,2);
  hist(60./(diff(RESP.tptrace)),15)
  set(gca,'fontsize',16)
  ylabel('Count');
  xlabel('Breaths per Minute');
  saveas(gcf,'pmu_qa_hists.png');
  
  h = figure('visible','off');
  subplot(2,1,1);
  plot(CARD.irf_comps(1:4,:)','linewidth',2)
  set(gca,'fontsize',16)
  ylabel('A.U.'); % Arbitrary Units
  xlabel('Relative % Phase in Cycle');
  subplot(2,1,2);
  plot(RESP.irf_comps(1:2,:)','linewidth',2)
  set(gca,'fontsize',16)
  ylabel('A.U.');
  xlabel('Relative % Phase in Cycle');
  saveas(gcf,'pmu_qa_irf.png');
end
