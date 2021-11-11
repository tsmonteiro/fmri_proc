function process_physio_data(OUTDIR, nScans, nScansExtra, TR, nSlices, scanAlign)

%%
%  clc;
%  clear all;
addpath(genpath('/home/luna.kuleuven.be/u0101486/workspace/matlab/toolbox/spm12/toolbox/PhysIO'));
addpath(genpath('/home/luna.kuleuven.be/u0101486/workspace/matlab/toolbox/PhLEM/'));
%TODO add PhLEM and TAPAS to the path here

fs        = 500; %496; % Wireless physio unit, 500Hz for the wired one [Philips at least]

 %OUTDIR    = '/home/luna.kuleuven.be/u0101486/workspace/data/ConnectEx/tmp/A002';
 %nScans    = 600;
 % nScansExtra    = 607; % nScansExtra here represents the number of scans in the Physio...
 % scanAlign = -1;

if scanAlign == -1
 lastScan  = nScansExtra;
 firstScan = nScansExtra - nScans + 1;
else

 firstScan = scanAlign +1 ;
 lastScan  = firstScan + nScans -1;
end
% nii = load_untouch_nii(cat(2,OUTDIR, filesep, 'func_data.nii'));
%  [x,y,z,t] = size(nii.img);

%  gs = mean(reshape(nii.img, x*y*z, t) );
%  TR        = 1;
%  nSlices   = 42;
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
t0 = tf - ((nScansExtra) * fs)+1 - fs*2 + 500;



% clf;
% plot( gradX(t0:tf)./max(gradX(t0:tf)) ); hold on;
% plot( gradY(t0:tf)./max(gradY(t0:tf)) ); hold on;
% plot( gradZ(t0:tf)./max(gradZ(t0:tf)) ); hold on;
% % TODO check that these thresholds are actually good enough
% Possibly do it iteratively, reducing the threshold until the correct
% number of peaks is identified (also the distance between them)
[~, locsX] = findpeaks(gradX(t0:tf)./max(gradX(t0:tf)), 'MinPeakProminence', 0.8, 'MinPeakDistance', 15);
[~, locsY] = findpeaks(gradY(t0:tf)./max(gradY(t0:tf)), 'MinPeakProminence', 0.999, 'MinPeakDistance', 15);
[~, locsZ] = findpeaks(gradZ(t0:tf)./max(gradZ(t0:tf)), 'MinPeakProminence', 0.8, 'MinPeakDistance', 15);

%%

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
dlmwrite(cat(2, OUTDIR, filesep, 'slice_timing_gs.txt'), (((sliceTiming-min(sliceTiming))/fs))');

% Checkpoint for correct scan length detected
ds        = (diff(scanTriggers));
estTR     = mean(ds(2:end));
estLength = (scanTriggers(605)-readoutf-1)-scanTriggers(5);

fprintf('Estimated inter-volume time %.3f\n', estTR/fs);
fprintf('Estimated total acquisition time %.3f\n', estLength/fs);


appu  = ppu(t0:tf);
aresp = resp(t0:tf);


%TODO  Correct atf at

if scanAlign == -1
    atf   = min( [(scanTriggers(lastScan) + fs*TR -1), length(appu) ]);
    at0   = atf - ( fs*nScans );
else
    at0   = scanTriggers( firstScan );
    atf   = scanTriggers( lastScan ) - readoutf -1;
end

% dlmwrite(cat(2, OUTDIR, filesep, 'physio.card'), appu(scanTriggers(5):(scanTriggers(605)-readoutf-1)))
% dlmwrite(cat(2, OUTDIR, filesep, 'physio.resp'), aresp(scanTriggers(5):(scanTriggers(605)-readoutf-1)))

dlmwrite(cat(2, OUTDIR, filesep, 'physio.card'), appu(at0:atf));
dlmwrite(cat(2, OUTDIR, filesep, 'physio.resp'), aresp(at0:atf));

%
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




%

% Use parameters to define first and last scan trigger
% timeVec    = (scanTriggers(5):(scanTriggers(605)-readoutf-1))./496;
% timeVec    = (scanTriggers(5):(scanTriggers(605)-readoutf-1))./496;
% timeVec    = timeVec - min(timeVec);

timeVecAll = (1:length(aresp))./496;
timeVecAll = timeVecAll - min(timeVecAll);

%

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
%    %%
[b,a] = butter(2, 1*[0.1, 3]/fs);
% [b,a]       = cheby2(1,4, [0.01 3]./fs, 'bandpass');
faresp       = filtfilt(b,a,aresp);

respRegD = [];
for d = 0:5:20
    arespDown = [];
    daresp    = circshift(cat(1, 0, diff(faresp)), fs*d);


    for i = 1:(length(scanTriggers))
        t0 = scanTriggers(i);
        tf = scanTriggers(i) + fs;

%         timeMask(end+1) = mode(outs(t0:tf));
        arespDown(end+1) = sqrt(sum(daresp(t0:tf).^2));
    end



    erespSig    = wextend('1d', 'ppd', arespDown, 20);
    erespSigOut = abs(zscore( erespSig ));
%     erespSigOut = abs(zscore(erespSig-erespSigOut));
% break
    [~,outlier] = findpeaks(erespSigOut, 'MinPeakHeight', 3);

    for k = 1:length(outlier)
        k0 = outlier(k);
        kf = outlier(k);

        while k0 > 1 && erespSigOut(k0) > 1
            k0 = k0 - 1;
        end

        while kf < length(erespSigOut) && erespSigOut(kf) > 1
            kf = kf + 1;
        end

        if k0 == 1
            if kf < length(erespSigOut)
                erespSig(k0:kf) = interp1([k0 kf kf+1], erespSig([k0 kf kf+1]), [k0:kf], 'spline');
            else
                erespSig(k0:kf) = interp1([k0 kf ], erespSig([k0 kf ]), [k0:kf], 'spline');
            end
        else
            if kf < length(erespSigOut)
                erespSig(k0:kf) = interp1([k0-1 k0 kf kf+1], erespSig([k0-1 k0 kf kf+1]), [k0:kf], 'spline');
            else
                erespSig(k0:kf) = interp1([k0-1 k0 kf ], erespSig([k0-1 k0 kf ]), [k0:kf], 'spline');
            end
        end
    end

    arespDown = erespSig(21:end-20);
    respRegD = cat(2, respRegD, arespDown');


end
%  clf;clc;
%  plot(zscore(respRegD(8:127,1))); hold on;
%  plot(zscore(gs(1:120)));
%  axis([0 120 -4 4]);
%  corr(respRegD(firstScan:lastScan,:), gs')
%   %%
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
%

respSig     = aresp;

% [b,a]       = cheby2(1,10, [0.01 40]./fs, 'bandpass');
% fResp       = filtfilt(b,a,respSig);


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

    while kf < length(erespSigOut) && erespSigOut(kf) > .5
        kf = kf + 1;
    end

    if k0 == 1
        if kf < length(erespSigOut)
            erespSig(k0:kf) = interp1([k0 kf kf+1], erespSig([k0 kf kf+1]), [k0:kf], 'spline');
        else
            erespSig(k0:kf) = interp1([k0 kf ], erespSig([k0 kf ]), [k0:kf], 'spline');
        end
    else
        if kf < length(erespSigOut)
            erespSig(k0:kf) = interp1([k0-1 k0 kf kf+1], erespSig([k0-1 k0 kf kf+1]), [k0:kf], 'spline');
        else
            erespSig(k0:kf) = interp1([k0-1 k0 kf ], erespSig([k0-1 k0 kf ]), [k0:kf], 'spline');
        end
    end
end

%
[b,a] = butter(2, 2*[0.1, 3]/496);
fResp = filtfilt(b,a,erespSig);
fResp = fResp(5001:end-5000);



%   close all;
%   hold on;
%   plot(respSig(2000:8000), 'k', 'LineWidth', 2, 'Color', [.6 .6 .6]); hold on;
%   plot(fResp(2000:8000), 'LineWidth', 2, 'Color', [0 0 0]); hold on;
% plot(fCard); hold on;
%
dt                  = fs / nSlices;
physio              = tapas_physio_new('RETROICOR');
model_rvt           = physio.model.rvt;
model_rvt.include   = 1;
model_rvt.delays    = [0:round(dt*5*nSlices):round(dt*15*nSlices)];

sqpar               = struct();
sqpar.onse_slice    = sliceEvents;
sqpar.Nslices       = nSlices;
sqpar.Nscans        = nScansExtra;
sqpar.Ndummies      = 0;
sqpar.TR            = TR * 1000;
sqpar.onset_slice   = 10;%cat(2, 1:14, 1:14, 1:14);%1;

spulsePerVolMs      = spulsePerVol;

for i = 1:length(spulsePerVolMs)
    ms = spulsePerVolMs{i};
    ms = ms ./ fs;

    spulsePerVolMs{i} = ms;
end


% fRespO = fResp;
% fRespO(outs==0) = mean(fResp);
ons_secs                = struct();
ons_secs.fr             = fResp;
ons_secs.t              = timeVecAll;
ons_secs.cpulse         = (find(phlem.PPG.events == 1))./1;
ons_secs.spulse_per_vol = spulsePerVolMs;
%
% clf;
% plot(fResp); hold on;
% plot(find(outs==0), fResp(find(outs==0)), 'r.', 'linewidth',2);
%
samplePoints            = sliceEvents./fs;
[convRVTOut, rvtOut, ~] = tapas_physio_create_rvt_regressors(...
                                ons_secs, sqpar,model_rvt, samplePoints);

%
%
for d = 1:size(convRVTOut,2)

    erespSig    = wextend('1d', 'ppd', convRVTOut(:,d), 20);
    erespSigOut = abs(zscore(erespSig));
    % erespSigOut = medfilt1(erespSig', 5)';
    % erespSigOut = abs(zscore(erespSig-erespSigOut));

    [~,outlier] = findpeaks(erespSigOut, 'MinPeakHeight', 2, 'MinPeakDistance', 10);

    for k = 1:length(outlier)
        k0 = outlier(k);
        kf = outlier(k);

        while k0 > 1 && erespSigOut(k0) > 0.5
            k0 = k0 - 1;
        end

        while kf < length(erespSigOut) && erespSigOut(kf) > .5
            kf = kf + 1;
        end

        if k0 > 1
            if kf == length(erespSigOut)
                erespSig(k0:kf) = interp1([k0-1 k0 kf ], erespSig([k0-1 k0 kf ]), [k0:kf], 'spline');
            else
                erespSig(k0:kf) = interp1([k0-1 k0 kf kf+1], erespSig([k0-1 k0 kf kf+1]), [k0:kf], 'spline');
            end
        else
            if kf == length(erespSigOut)
                erespSig(k0:kf) = interp1([ k0 kf ], erespSig([k0 kf ]), [k0:kf], 'spline');
            else
                erespSig(k0:kf) = interp1([ k0 kf kf+1], erespSig([k0 kf kf+1]), [k0:kf], 'spline');
            end
        end
    end
    convRVTOut(:,d) = erespSig(21:end-20);
end
% clf;
% hold on;
% plot(rvtOut);
% plot(convRVTOut);

%
%
%   clf;
%   subplot(211); plot(convRVTOut);
%   subplot(212); plot(rvtOut);
%
dlmwrite(cat(2, OUTDIR, filesep, 'physio_reg.rvt'), convRVTOut(firstScan:lastScan,:));
%


%
verbose             = [];
verbose.level       = 0;
verbose.process_log = {};

% [cardiac_phase, ~] = tapas_physio_get_cardiac_phase( ...
%                         (find(phlem.PPG.events == 1))./fs, ...
%                         timeVecAll', verbose, scanTriggers'./fs);
%
% [rphase, fh]       = tapas_physio_get_respiratory_phase( ...
%                             fResp, 1./fs, 0);

sqpar                = struct();
sqpar.onse_slice     = sliceEvents'./fs;
sqpar.Nslices        = nSlices;
sqpar.Nscans         = nScansExtra;
sqpar.Ndummies       = 0;
sqpar.TR             = TR * 1000;
sqpar.onset_slice    = 1;
sqpar.NslicesPerBeat = nSlices;


order    = struct();
order.c  = 3;
order.r  = 4;
order.cr = 2;

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

% cardiac_sess = tapas_physio_scaleorthmean_regressors(cardiac_sess);
%

for c = 1:size(respire_sess,2)
   respire_sess(timeMask==0,c)  = max(respire_sess(timeMask==1,c));
end

dlmwrite(cat(2, OUTDIR, filesep, 'physio_reg.cr'), ...
    cat(2, cardiac_sess(firstScan:lastScan,:), ...
    respire_sess(firstScan:lastScan,:), mult_sess(firstScan:lastScan,:), ...
    c_sample_phase(firstScan:lastScan,:), r_sample_phase(firstScan:lastScan,:)));
% %
physio              = tapas_physio_new('RETROICOR');
model_rvt           = physio.model.hrv;
model_rvt.include   = 1;
model_rvt.delays    = [0:round(dt*5*nSlices):round(dt*15*nSlices)];


sqpar             = struct();
sqpar.onse_slice  = sliceEvents;
sqpar.Nslices     = nSlices;
sqpar.Nscans      = nScansExtra;
sqpar.Ndummies    = 0;
sqpar.TR          = TR * 1000;
sqpar.onset_slice = 10;%cat(2, 1:14, 1:14, 1:14);


ons_secs                = struct();
ons_secs.cpulse         = (find(phlem.PPG.events == 1));
ons_secs.spulse_per_vol = spulsePerVol;

[hrv, ~] = tapas_physio_create_hrv_regressors(ons_secs, sqpar, model_rvt);

% clf ; hold on;
% plot(hrv);

for d = 1:size(hrv,2)
    % Uncomment the code below in case outlier removal is desirable
    % which I am not sure it is

    ahrv = abs(zscore(detrend(hrv(:,d))));
    ahrv = wextend('1d', 'ppd', ahrv, 20);
    hrv_  = wextend('1d', 'ppd', hrv(:,d), 20);

    [~,outHrv] = findpeaks(ahrv, 'MinPeakHeight', 2);

    if ~isempty(outHrv)
      for k = 1:length(outHrv)
        k0 = outHrv(k);
        kf = outHrv(k);

        while k0 > 2 && ahrv(k0) > 0.5
            k0 = k0 - 1;
        end
        while kf < (length(ahrv)-1) && ahrv(kf) > .5
            kf = kf + 1;
        end

        hrv_(k0:kf) = interp1([k0-1 k0 kf kf+1], hrv([k0-1 k0 kf kf+1]), [k0:kf], 'spline');
      end

      hrv(:,d) = hrv_(21:end-20);
    end
end
% plot(hrv, 'k');
%
dlmwrite(cat(2, OUTDIR, filesep, 'physio_reg.hrv'), hrv(firstScan:lastScan,:).*1000);


[respReg, namesR]   = make_physio_regressors('RESP', phlem, 'TR', 1, 'DOSAVE', 0, 'ana_type', 'fourier');
[respRReg, namesRR] = make_physio_regressors('RESP', phlem, 'TR', 1, 'DOSAVE', 0, 'ana_type', 'rvhr');



[ppgReg, namesP]    = make_physio_regressors('PPG', phlem, 'TR', 1, 'DOSAVE', 0);
[ppghReg, namesPH]  = make_physio_regressors('PPGHIGH', phlem, 'TR', 1, 'DOSAVE', 0, 'ana_type', 'rvhr');
[ppglReg, namesPL]  = make_physio_regressors('PPGLOW', phlem, 'TR', 1, 'DOSAVE', 0, 'ana_type', 'rvhr');


respReg   = respReg(firstScan:lastScan,:);
respRReg  = respRReg(firstScan:lastScan,:);
ppgReg    = ppgReg(firstScan:lastScan,:);
ppghReg   = ppghReg(firstScan:lastScan,:);
ppglReg   = ppglReg(firstScan:lastScan,:);
timeMaskc = timeMask(:,firstScan:lastScan)';

for c = 1:size(respReg,2)
   cs = respReg(:,c);
   cs(timeMaskc == 0) = 0;%max(cs);
   respReg(:,c) = cs;
end

for c = 1:size(respRReg,2)
   cs = respRReg(:,c);
   cs(timeMaskc == 0) = 0;%max(cs);
   respRReg(:,c) = cs;
end

dlmwrite(cat(2, OUTDIR, filesep, 'physio_reg.resp'), respReg);
dlmwrite(cat(2, OUTDIR, filesep, 'physio_reg.respr'), respRReg);

dlmwrite(cat(2, OUTDIR, filesep, 'physio_reg.ppg'), ppgReg);
dlmwrite(cat(2, OUTDIR, filesep, 'physio_reg.ppgh'), ppghReg);
dlmwrite(cat(2, OUTDIR, filesep, 'physio_reg.ppgl'), ppglReg);



respire_sess(timeMaskc == 0,:) = 0;
r_sample_phase(timeMaskc == 0) = 0;

respire_sess(isnan(r_sample_phase) ,:) = 0;
r_sample_phase(isnan(r_sample_phase)) = 0;
%
% co = convRVTOut(firstScan:lastScan,:);
% plot(co(timeMaskc==1,1));
%
% rp = dlmread(cat(2, OUTDIR, filesep, 'motion_estimate.par'));
% rp = cat(2, rp, rp.^2);
% rp = cat(2, rp, cat(1, zeros(1,12),diff(rp)));
% ppgReg, ppghReg, ppglReg, respReg, rp
otherRegs = cat(2, c_sample_phase(firstScan:lastScan,:), ...
                   r_sample_phase(firstScan:lastScan,:), ...
                   respRReg, ppghReg, respReg, ...
                   respRegD(firstScan:lastScan,:), ...
                   convRVTOut(firstScan:lastScan,:), hrv(firstScan:lastScan,:) );
 clf;
%  subplot(121);
% imagesc(corr(cat(2,cardiac_sess(firstScan:lastScan,:), ...
%     respire_sess(firstScan:lastScan,:), ...
%     mult_sess(firstScan:lastScan,:), ...
%     otherRegs)), [-1 1]); colormap gray;

% subplot(122);

[R, verbose] = tapas_physio_orthogonalise_physiological_regressors(...
    cardiac_sess(firstScan:lastScan,:), ...
    respire_sess(firstScan:lastScan,:), ...
    mult_sess(firstScan:lastScan,:), ...
    otherRegs, 'all', verbose);

R = R(:,~isnan(mean(R)));


imagesc(R, [-1 1]); colormap gray;
saveas(gcf, cat(2, OUTDIR, filesep, 'physio_model.png'));
% close all;
%

dlmwrite(cat(2, OUTDIR, filesep, 'physio_reg.ortho'), R);

% %%
%
% %  %%
%  clf;
%  plot( corr(detrend(gs'),R).^2, 'ksq--' );
%  axis([0 size(R,2)+1 0 0.5]);
% %%
% % clf;
% % dgs = detrend(gs);
% % dgs = zscore(dgs);
% % plot(dgs(1:120), 'k'); hold on;
% % plot(zscore(arespDown(7:127)));
% % axis([0 120 -5 5]);
% % respRegD
% clf;
% subplot(211);
% plot(zscore(detrend(gs)),'k'); hold on; plot(zscore(R(:,2)))
% subplot(212);
%  plot(zscore(detrend(gs)), 'k');
%  hold on; plot(zscore(respRegD(firstScan:lastScan,1)))
% % axis([0 120 -4 4]);
% % arespDown
