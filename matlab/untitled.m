fid  = fopen('physio.log');
data = textscan(fid, '%n %n %n %n %n %n %n %n %n %n %n', 'CommentStyle', '#');
fclose(fid);

%%

ppu  = data{5};
resp = data{6};


gradX = data{7};
gradY = data{8};
gradZ = data{9};

markers = hex2dec( string(data{10}));

%%
fs = 500;
nScans = 600;
nDummy = 4;

lastScanPhys = find(markers == 32, 1, 'last');
firstScanPhys = find(markers == 16, 1, 'last');

%%


tf = lastScanPhys;
t0 = tf - ((nScans+nDummy) * fs)+1;

clf; hold on;
plot(gradX(t0:tf));
plot(gradY(t0:tf));
plot(gradZ(t0:tf));

%%

[pksX, locsX] = findpeaks(gradX(t0:tf)./max(gradX(t0:tf)), 'MinPeakProminence', 0.8, 'MinPeakDistance', 15);
[pksY, locsY] = findpeaks(gradY(t0:tf)./max(gradY(t0:tf)), 'MinPeakProminence', 0.999, 'MinPeakDistance', 15);
[pksZ, locsZ] = findpeaks(gradZ(t0:tf)./max(gradZ(t0:tf)), 'MinPeakProminence', 0.8, 'MinPeakDistance', 15);
%%

(length(locsX)+length(locsY)+length(locsZ))/45

%%
gx = gradX(t0:tf)./max(gradX(t0:tf) );
gy = gradY(t0:tf)./max(gradY(t0:tf) );
gz = gradZ(t0:tf)./max(gradZ(t0:tf) );
scanTriggers = [];

sliceOnset   = 0;



p0 = 1;
pf = 15;

l1 = locsX(1);
l2 = locsY(1);
l3 = locsZ(1);

slThr = 1e-6;

while (gx(l1-sliceOnset)^2 > slThr || ...
      gy(l1-sliceOnset)^2 >  slThr || ... 
      gz(l1-sliceOnset)^2 >  slThr) && (l1-sliceOnset) > 0
    sliceOnset = sliceOnset + 1;
end

readoutf = 0;
slZ = locsZ(3);
while gz(slZ+readoutf)^2 > 1e-6 || gx(slZ+readoutf)^2 > 1e-6 || gy(slZ+readoutf)^2 > 1e-6
    readoutf = readoutf + 1;
end

scanTriggers = locsX(1)-sliceOnset;

for p = 16:15:length(locsZ)
    scanTriggers(end+1) = locsX(p)-sliceOnset;
end

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

sliceEvents = [];
spulsePerVol = {};
for p = 1+(15*0):length(locsX)
    if mod(p, 15) ~= 1
        sliceEvents = cat(2, sliceEvents, locsY(p));
        sliceEvents = cat(2, sliceEvents, locsX(p));
        sliceEvents = cat(2, sliceEvents, locsZ(p));
    end
    
    if mod( p, 15) == 0
       spulsePerVol{end+1} = sliceEvents(end-41:end); 
    end
end

dlmwrite('slice_timing_g.txt', round((sliceTiming/496)*1000)');
dlmwrite('slice_timing_g0.txt', round(((sliceTiming-min(sliceTiming))/496)*1000)');

clc
ds = (diff(scanTriggers));
mean(ds(2:end))
(scanTriggers(605)-readoutf-1)-scanTriggers(5)

% st = cat(2, scanTriggers', scanTriggers', scanTriggers');

%dlmwrite('scan_triggers.tsv', st(5:605,:) );
appu  = ppu(t0:tf);
aresp = resp(t0:tf);

dlmwrite('physio_2.card', appu(scanTriggers(5):(scanTriggers(605)-readoutf-1)))
dlmwrite('physio_2.resp', fResp(scanTriggers(5):(scanTriggers(605)-readoutf-1)))

%%
clc;
nii = load_untouch_nii('func_data.nii');
msk = load_untouch_nii('mask.nii');

fData = nii.img;
[x,y,z,t] = size(fData);
fData = reshape(fData, x*y*z, t);
gSig = mean(fData(msk.img(:)==1, :));
%%

sResp = aresp(scanTriggers(5):(scanTriggers(605)-readoutf-1));
clf; hold on;

plot(0:599, zscore(gSig));
plot((0:(length(sResp)-1))./496, sResp);

%%


appu  = ppu(t0:tf)./max(ppu(t0:tf));
aresp = resp(t0:tf)./max(resp(t0:tf));

k0 = 800;
kf = 2200;

clf; hold on;
plot( (k0:kf)./500, gx(k0:kf) );
plot( (k0:kf)./500, gy(k0:kf) );
plot( (k0:kf)./500, gz(k0:kf) );


plot( (k0:kf)./500, appu(k0:kf), 'color', [0.5, 0.5, 0.5], 'LineWidth', 2 );

% scanTriggers = [];


for p = 1:30
   scatter((locsX(p)-1)/500, pksX(p), 'k');
   scatter((locsY(p)-1)/500, pksY(p), 'k');
   scatter((locsZ(p)-1)/500, pksZ(p), 'k');
end

plot([scanTriggers(1)-1 scanTriggers(1)-1]./fs,[-1,1], 'b', 'LineWidth', 2);
plot([locsZ(15)+readoutf,locsZ(15)+readoutf]./fs,[-1,1], 'k', 'LineWidth', 2);
% plot([scanTriggers(1)-1+496 scanTriggers(1)-1+496]./500,[-1,1], 'r', 'LineWidth', 2);


plot([scanTriggers(2)-1, scanTriggers(2)-1]./fs,[-1,1], 'b', 'LineWidth', 2);
plot([locsZ(30)+readoutf,locsZ(30)+readoutf]./fs,[-1,1], 'k', 'LineWidth', 2);
% plot([scanTriggers(2)-1+496 scanTriggers(2)-1+496]./500,[-1,1], 'r', 'LineWidth', 2);

plot([scanTriggers(3)-1, scanTriggers(3)-1]./fs,[-1,1], 'b', 'LineWidth', 2);
plot([locsZ(45)+readoutf,locsZ(45)+readoutf]./fs,[-1,1], 'k', 'LineWidth', 2);



% (locsZ(15)+readoutf)-(l1-sliceOnset)

axis([1.5 4 -1.1 1.1])

