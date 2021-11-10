function process_physio_data(OUTDIR, nScans, nDummy, TR, nSlices)
% clc; clear all;
% OUTDIR    = '/home/luna.kuleuven.be/u0101486/workspace/data/ConnectEx/tmp/A005/';
fs        = 496; % Wireless physio unit, 500Hz for the wired one [Philips at least]
% nScans    = 600;
% nDummy    = 4;
% TR        = 1;
% nSlices   = 42;
nPrepGrad = 3;

% Read SCANLOG file
fid  = fopen(cat(2, OUTDIR, filesep, 'physio.log'));
data = textscan(fid, '%n %n %n %n %n %n %n %n %n %n %n', 'CommentStyle', '#');
fclose(fid);

ppu     = data{5};
resp    = data{6};

gradX   = data{7};
gradY   = data{8};
gradZ   = data{9};

markers = hex2dec( string(data{10}));

% Those markers sadly are not a perfect match for scanning times
lastScanPhys  = find(markers == 32, 1, 'last');
firstScanPhys = find(markers == 16, 1, 'last');


tf = lastScanPhys;
t0 = tf - ((nScans+nDummy) * fs)+1 - fs*4;


% TODO check that these thresholds are actually good enough
% Possibly do it iteratively, reducing the threshold until the correct
% number of peaks is identified (also the distance between them)
[pksX, locsX] = findpeaks(gradX(t0:tf)./max(gradX(t0:tf)), 'MinPeakProminence', 0.8, 'MinPeakDistance', 15);
[pksY, locsY] = findpeaks(gradY(t0:tf)./max(gradY(t0:tf)), 'MinPeakProminence', 0.999, 'MinPeakDistance', 15);
[pksZ, locsZ] = findpeaks(gradZ(t0:tf)./max(gradZ(t0:tf)), 'MinPeakProminence', 0.8, 'MinPeakDistance', 15);

% %%

gx = gradX(t0:tf)./max( gradX(t0:tf) );
gy = gradY(t0:tf)./max( gradY(t0:tf) );
gz = gradZ(t0:tf)./max( gradZ(t0:tf) );

scanTriggers = [];
sliceOnset   = 0;

% Slice dependent
% Importantlty, this might not exactly match the number in the slices in
% the volume, that is there might be a preparatory pulse across gradients.
p0    = 1;
pf    = 15;

l1    = locsX(1);
l2    = locsY(1);
l3    = locsZ(1);


% TODO Check this threshold
slThr = 1e-6;

while (gx(l1-sliceOnset)^2 > slThr || ...
      gy(l1-sliceOnset)^2 >  slThr || ... 
      gz(l1-sliceOnset)^2 >  slThr) && (l1-sliceOnset) > 0
    sliceOnset = sliceOnset + 1;
end

readoutf = 0;
slZ      = locsZ(3); % Again... prep pulse
while gz(slZ+readoutf)^2 > slThr || ...
      gx(slZ+readoutf)^2 > slThr || ...
      gy(slZ+readoutf)^2 > slThr
    readoutf = readoutf + 1;
end

scanTriggers = locsX(1)-sliceOnset;

for p = 16:15:length(locsZ)
    scanTriggers(end+1) = locsX(p)-sliceOnset;
end
% TODO Get this info from parameters
sliceTiming = [];

for p = 2:15
   sliceTiming = cat(2, sliceTiming, locsY(p)-scanTriggers(1));
end
for p = 2:15
   sliceTiming = cat(2, sliceTiming, locsX(p)-scanTriggers(1));
end
for p = 2:15
   sliceTiming = cat(2, sliceTiming, locsZ(p)-scanTriggers(1));
end

sliceEvents  = [];
spulsePerVol = {};
for p = 1+(15*0):length(locsX)
    if mod(p, 15) ~= 1
        sliceEvents = cat(2, sliceEvents, locsY(p));
        sliceEvents = cat(2, sliceEvents, locsX(p));
        sliceEvents = cat(2, sliceEvents, locsZ(p));
    end
    
    if mod(p, 15) == 0
       spulsePerVol{end+1} = sliceEvents(end-41:end); 
    end
end

dlmwrite(cat(2, OUTDIR, filesep, 'slice_timing_g.txt'), round((sliceTiming/fs)*1000)');
dlmwrite(cat(2, OUTDIR, filesep, 'slice_timing_g0.txt'), round(((sliceTiming-min(sliceTiming))/fs)*1000)');

% Checkpoint for correct scan length detected
ds        = (diff(scanTriggers));
estTR     = mean(ds(2:end));
estLength = (scanTriggers(605)-readoutf-1)-scanTriggers(5);

fprintf('Estimated inter-volume time %.3f\n', estTR/fs);
fprintf('Estimated total acquisition time %.3f\n', estLength/fs);


appu  = ppu(t0:tf);
aresp = resp(t0:tf);

dlmwrite(cat(2, OUTDIR, filesep, 'physio.card'), appu(scanTriggers(5):(scanTriggers(605)-readoutf-1)))
dlmwrite(cat(2, OUTDIR, filesep, 'physio.resp'), aresp(scanTriggers(5):(scanTriggers(605)-readoutf-1)))

% %%


appu  = ppu(t0:tf)./max(ppu(t0:tf));
aresp = resp(t0:tf)./max(resp(t0:tf));

% If needed, use the code below to check if the marks are in the right
% place
% 
% k0 = 800;
% kf = 2200;
% 
% clf; hold on;
% plot( (k0:kf)./500, gx(k0:kf) );
% plot( (k0:kf)./500, gy(k0:kf) );
% plot( (k0:kf)./500, gz(k0:kf) );
% 
% plot( (k0:kf)./500, appu(k0:kf), 'color', [0.5, 0.5, 0.5], 'LineWidth', 2 );
% 
% scanTriggers = [];
% 
% for p = 1:30
%    scatter((locsX(p)-1)/500, pksX(p), 'k');
%    scatter((locsY(p)-1)/500, pksY(p), 'k');
%    scatter((locsZ(p)-1)/500, pksZ(p), 'k');
% end
% 
% plot([scanTriggers(1)-1 scanTriggers(1)-1]./fs,[-1,1], 'b', 'LineWidth', 2);
% plot([locsZ(15)+readoutf,locsZ(15)+readoutf]./fs,[-1,1], 'k', 'LineWidth', 2);
% 
% plot([scanTriggers(2)-1, scanTriggers(2)-1]./fs,[-1,1], 'b', 'LineWidth', 2);
% plot([locsZ(30)+readoutf,locsZ(30)+readoutf]./fs,[-1,1], 'k', 'LineWidth', 2);
% 
% plot([scanTriggers(3)-1, scanTriggers(3)-1]./fs,[-1,1], 'b', 'LineWidth', 2);
% plot([locsZ(45)+readoutf,locsZ(45)+readoutf]./fs,[-1,1], 'k', 'LineWidth', 2);
% 
% 
% 
% % (locsZ(15)+readoutf)-(l1-sliceOnset)
% 
% axis([1.5 4 -1.1 1.1])




% %%

% Use parameters to define first and last scan trigger
timeVec    = (scanTriggers(5):(scanTriggers(605)-readoutf-1))./496;
timeVec    = timeVec - min(timeVec);

timeVecAll = (1:length(aresp))./496;
timeVecAll = timeVecAll - min(timeVecAll);

% %%

% 'outs' contains the a mask (physio sampling) of points where the
% respiration measurement was unreliable (e.g. belt was not in contact with participant)
outs       = sqrt(smooth(aresp.^2,496));
outs       = outs./max(outs);
% plot(outs); hold on;
% TODO Check this threshold here
o          = find(outs < 0.1);
% plot(o, outs(o), 'k');
outs(o)      = 0;
outs(outs > 0) = 1;

% 'timeMask' is the downsampled version of outs
timeMask = [];

for i = 1:(length(scanTriggers))
    t0 = scanTriggers(i);
    tf = scanTriggers(i) + fs;
    
    timeMask(end+1) = mode(outs(t0:tf));
end


a_opts = PhLEM_ana_opts({'RESP', 'PPG', 'PPGhigh', 'PPGlow'}, fs);


% ADD all times
phlem  = PhLEM_constructor('RESP', aresp, ...
                           'PPG', appu, ...
                           'PPGHIGH', appu, ...
                           'PPGLOW', appu, ...
                           'Time', timeVecAll, 'a_opts', a_opts );
                      
phlem  = PhLEM_process(phlem);
phlem  = PhLEM_analysis(phlem);

% Uncomment this to check if events are being correctly detected
%phlem = PhLEM_resultsplot(phlem,'PPGHIGH');
% %%

respSig     = aresp;

[b,a]       = cheby2(1,10, [0.01 40]./fs, 'bandpass');
fResp       = filtfilt(b,a,respSig);


erespSig    = wextend('1d', 'ppd', respSig, 5000);
erespSigOut = medfilt1(erespSig', 30)';
erespSigOut = abs(zscore(erespSig-erespSigOut));

[~,outlier] = findpeaks(erespSigOut, 'MinPeakHeight', 4);
 
for k = 1:length(outlier)
    k0 = outlier(k);
    kf = outlier(k);
 
    while k0 > 0 && erespSigOut(k0) > 0.5  
        k0 = k0 - 1;
    end
    
    while kf <= length(erespSigOut) && erespSigOut(kf) > .5  
        kf = kf + 1;
    end
 
    erespSig(k0:kf) = interp1([k0-1 k0 kf kf+1], erespSig([k0-1 k0 kf kf+1]), [k0:kf], 'spline');
end


[b,a] = butter(2, 2*[0.01, 45]/496);
fResp = filtfilt(b,a,erespSig);
fResp = fResp(5001:end-5000);

% close all;
% hold on;
% plot(respSig, 'k', 'LineWidth', 2); hold on;
% plot(fResp); hold on;
% plot(fCard); hold on;

dt                  = fs / nSlices;
physio              = tapas_physio_new('RETROICOR');
model_rvt           = physio.model.rvt;
model_rvt.include   = 1;
model_rvt.delays    = [0:round(dt*5*nSlices):round(dt*20*nSlices)];

sqpar               = struct();
sqpar.onse_slice    = sliceEvents;
sqpar.Nslices       = nSlices;
sqpar.Nscans        = nScans;
sqpar.Ndummies      = nDummy;
sqpar.TR            = TR * 1000;
sqpar.onset_slice   = 1;
 
spulsePerVolMs      = spulsePerVol;

for i = 1:length(spulsePerVolMs)
    ms = spulsePerVolMs{i};
    ms = ms ./ fs;
    
    spulsePerVolMs{i} = ms;
end

ons_secs                = struct();
ons_secs.fr             = fResp;
ons_secs.t              = timeVecAll;
ons_secs.cpulse         = (find(phlem.PPG.events == 1))./fs;
ons_secs.spulse_per_vol = spulsePerVolMs;
 
samplePoints            = sliceEvents./fs;
[convRVTOut, rvtOut, ~] = tapas_physio_create_rvt_regressors(...
                                ons_secs, sqpar,model_rvt, samplePoints);

% clf;
% subplot(211); plot(convRVTOut);
% subplot(212); plot(rvtOut);

dlmwrite(cat(2, OUTDIR, filesep, 'physio_reg.rvt'), convRVTOut(5:604,:));

% %%
verbose             = [];
verbose.level       = 0;
verbose.process_log = {};

[cardiac_phase, ~] = tapas_physio_get_cardiac_phase( ...
                        (find(phlem.PPG.events == 1))./fs, ...
                        timeVecAll', verbose, scanTriggers'./fs);

[rphase, fh]       = tapas_physio_get_respiratory_phase( ...
                            fResp, 1./fs, 0);

sqpar                = struct();
sqpar.onse_slice     = sliceEvents'./fs;
sqpar.Nslices        = nSlices;
sqpar.Nscans         = nScans;
sqpar.Ndummies       = nDummy;
sqpar.TR             = TR * 1000;
sqpar.onset_slice    = 1;
sqpar.NslicesPerBeat = nSlices;
 

order    = struct();
order.c  = 3;
order.r  = 4;
order.cr = 1;
 
ons_secs                = struct();
ons_secs.fr             = fResp;
ons_secs.t              = timeVecAll'.*fs;
ons_secs.r              = aresp;
ons_secs.cpulse         = (find(phlem.PPG.events == 1))./1;
ons_secs.spulse         = sliceEvents'./1;
ons_secs.svolpulse      = scanTriggers'./fs;
ons_secs.spulse_per_vol = spulsePerVol;

 
[cardiac_sess, respire_sess, mult_sess, ~, order, ~, ...
    c_sample_phase, r_sample_phase] ...
    = tapas_physio_create_retroicor_regressors(ons_secs, sqpar, order, verbose);

cardiac_sess = tapas_physio_scaleorthmean_regressors(cardiac_sess);


for c = 1:size(respire_sess,2)
   respire_sess(timeMask(5:604)==0,c)  = max(respire_sess(timeMask(5:604)==1,c));
end

dlmwrite(cat(2, OUTDIR, filesep, 'physio_reg.cr'), cat(2, cardiac_sess, respire_sess, mult_sess, c_sample_phase, r_sample_phase));


sqpar             = struct();
sqpar.onse_slice  = sliceEvents;
sqpar.Nslices     = nSlices;
sqpar.Nscans      = nScans;
sqpar.Ndummies    = nDummy;
sqpar.TR          = TR * 1000;
sqpar.onset_slice = 1;
 
ons_secs                = struct();
ons_secs.cpulse         = (find(phlem.PPG.events == 1));
ons_secs.spulse_per_vol = spulsePerVol;
 
hrv = tapas_physio_create_hrv_regressors(ons_secs, sqpar);
% Uncomment the code below in case outlier removal is desirable
% which I am not sure it is
%
%
%ahrv = abs(zscore(detrend(hrv)));
%ahrv = wextend('1d', 'ppd', ahrv, 20);
%hrv  = wextend('1d', 'ppd', hrv, 20);
%  
%[~,outHrv] = findpeaks(ahrv, 'MinPeakHeight', 1.5);
%  
%for k = 1:length(outHrv)
%   k0 = outHrv(k);
%   kf = outHrv(k);
%    
%   while ahrv(k0) > 0.5 && k0 > 0
%       k0 = k0 - 1;
%   end
%   while ahrv(kf) > .5 && kf <= length(ahrv)
%       kf = kf + 1;
%   end
%     
%   hrv(k0:kf) = interp1([k0-1 k0 kf kf+1], hrv([k0-1 k0 kf kf+1]), [k0:kf], 'spline');
%end
%  
%hrv = hrv(21:end-20);
%  
%plot(hrv, 'k');
%  
dlmwrite(cat(2, OUTDIR, filesep, 'physio_reg.hrv'), hrv.*1000);


[respReg, namesR]   = make_physio_regressors('RESP', phlem, 'TR', 1, 'DOSAVE', 0, 'ana_type', 'fourier');
[respRReg, namesRR] = make_physio_regressors('RESP', phlem, 'TR', 1, 'DOSAVE', 0, 'ana_type', 'rvhr');


[ppgReg, namesP]    = make_physio_regressors('PPG', phlem, 'TR', 1, 'DOSAVE', 0);
[ppghReg, namesPH]  = make_physio_regressors('PPGHIGH', phlem, 'TR', 1, 'DOSAVE', 0, 'ana_type', 'rvhr');
[ppglReg, namesPL]  = make_physio_regressors('PPGLOW', phlem, 'TR', 1, 'DOSAVE', 0, 'ana_type', 'rvhr');


respReg   = respReg(5:604,:);
respRReg  = respRReg(5:604,:);
ppgReg    = ppgReg(5:604,:);
ppghReg   = ppghReg(5:604,:);
ppglReg   = ppglReg(5:604,:);
timeMaskc = timeMask(5:604);

for c = 1:size(respReg,2)
   cs = respReg(:,c);
   cs(timeMaskc == 0) = max(cs);
   respReg(:,c) = cs;
end

for c = 1:size(respRReg,2)
   cs = respRReg(:,c);
   cs(timeMaskc == 0) = max(cs);
   respRReg(:,c) = cs;
end

dlmwrite(cat(2, OUTDIR, filesep, 'physio_reg.resp'), respReg);
dlmwrite(cat(2, OUTDIR, filesep, 'physio_reg.respr'), respRReg);

dlmwrite(cat(2, OUTDIR, filesep, 'physio_reg.ppg'), ppgReg);
dlmwrite(cat(2, OUTDIR, filesep, 'physio_reg.ppgh'), ppghReg);
dlmwrite(cat(2, OUTDIR, filesep, 'physio_reg.ppgl'), ppglReg);


otherRegs = cat(2, c_sample_phase, r_sample_phase, ...
            respReg, respRReg, ppgReg, ppghReg, ppglReg, convRVTOut(5:604,:), hrv);

        
        
[R, verbose] = tapas_physio_orthogonalise_physiological_regressors(cardiac_sess, respire_sess, mult_sess, otherRegs, 'all', verbose);

R = R(:,~isnan(mean(R)));

dlmwrite(cat(2, OUTDIR, filesep, 'physio_reg.ortho'), R);




