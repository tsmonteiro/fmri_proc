BASEDIR='/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp';


meanMov = [];
stdMov  = [];
maxMov  = [];
maxDist = [];
sids    = {};

for i = 1:106
   
   SID   = sprintf('RS%03d',i);
   SFILE = cat(2, BASEDIR, filesep, SID, filesep, 'maximum_disp.1d_delt' );
   SFILE2 = cat(2, BASEDIR, filesep, SID, filesep,  'maximum_disp.1d' );
   
   if exist(SFILE, 'file')
       
       fd = dlmread(SFILE, '\t', 2, 0);
       fd2 = dlmread(SFILE2, '\t', 2, 0);
       
       meanMov(end+1) = mean(fd);
       stdMov(end+1)  = std(fd);
       maxMov(end+1)  = max(fd);
       maxDist(end+1) = max(fd2);
       sids{end+1}    = SID;
   end
    
end


ages = dlmread('/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/ages.txt');
ages = ages(:,1);


scores = zeros(106,1);
[~,idx1]  = sort(meanMov);
[~,idx2]  = sort(stdMov);
[~,idx3]  = sort(maxMov);
[~,idx4]  = sort(maxDist);

for i = 1:20
   scores(idx1(i)) = scores(idx1(i)) + 5;
   scores(idx2(i)) = scores(idx2(i)) + 5;
   scores(idx3(i)) = scores(idx3(i)) + 5;
   scores(idx4(i)) = scores(idx4(i)) + 5;
end
for i = 21:45
   scores(idx1(i)) = scores(idx1(i)) + 4;
   scores(idx2(i)) = scores(idx2(i)) + 4;
   scores(idx3(i)) = scores(idx3(i)) + 4;
   scores(idx4(i)) = scores(idx4(i)) + 4;
end

for i = 46:75
   scores(idx1(i)) = scores(idx1(i)) + 3;
   scores(idx2(i)) = scores(idx2(i)) + 3;
   scores(idx3(i)) = scores(idx3(i)) + 3;
   scores(idx4(i)) = scores(idx4(i)) + 3;
end
%%
clf;
plot(scores,'ksq');

idx = find( scores>16 );
for i = 1:length(idx);
    fprintf('%d - %s\n', i, sids{idx(i)});
end

%%


nii = load_untouch_nii('/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/RS001_old/proc_data_native.nii');
%%

sig = squeeze(nii.img(45,76,10,:));

clf;
subplot(211); plot(sig);
subplot(212); 
[pxx,freqs] = pwelch(detrend(sig),160,[],[],1/2.5);
plot(freqs, pxx);




%%

clf;


subplot(2,2,1);
plot(meanMov, 'sqk');



subplot(2,2,2);
plot(stdMov, 'sqk');

subplot(2,2,3);
plot(maxMov, 'sqk');


subplot(2,2,4);
plot(maxDist, 'sqk');


dlmwrite('/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/Quality_Control/group_motion.txt', ...
    cat(1, meanMov, stdMov, maxMov, maxDist)');



%%
clc;
BASEDIR='/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp';

% ages = dlmread('/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/ages.txt');
% ages = ages(:,1);

SREGMODEL='FAD_DiC_RP24';
REGMODEL=cat(2,'QA_',SREGMODEL,'_AGG');


fcMats = zeros(800,800,0);
struMat = [];
allNodes = [];
inc = [];
exc = [];
nExc = 0;
atlas = 'yeo17_cobra';
atlas = 'local_global_cobra';
% meanFc =zeros(106,1);
meanFc = [];
% excs = sort(unique([3 8 12 14 20 38 69 75 87 90 92 4 11 24 31 39 49 63 85 97 105 106]));
for i = 1:92
   
%    SID   = sprintf('RS%03d',i);
   SID   = param.Subjects{i};
   SFILE = cat(2, BASEDIR, filesep, SID, filesep, REGMODEL, filesep, ...
       sprintf('FC_%s_%s.txt',  SREGMODEL, atlas) );

   NFILE =  cat(2, BASEDIR, filesep, SID, filesep, REGMODEL, filesep, ...
       sprintf('FC_%s_%s_node.txt',  SREGMODEL, atlas) );
   
   MFILE = cat(2, BASEDIR, filesep, SID, filesep, 'maximum_disp.1d_delt' );
   DFILE = cat(2, BASEDIR, filesep, SID, filesep, 'maximum_disp.1d' );
   
%    AFILE = sprintf('/media/u0101486/Seagate Backup Plus Drive/Crunch_SC_AAL2_dirty/R%03d_muw__connectome_tractogram10MSIFT2_FOD_AVE99_WM_norm_zerodia_sym_aal2.mat',i);
       
   fd = dlmread(MFILE, '\t', 2, 0);
   md = dlmread(DFILE, '\t', 2, 0);
   
   if mean(fd) <= 0.3 && max(md) < 5 && sum(fd>0.5) < 30
%    if ~any(i == excs)
       inc(end+1) = i;
   else
       nExc = nExc + 1;
       exc(end+1) = i;
       fprintf('%d - Excluding RS%03d (Aged %.1f)\n', nExc, i, ages(i));
   end
   
   
   if exist(SFILE, 'file')
%        fprintf('%d - Adding RS%03d (Aged %.1f)\n', length(inc), i, ages(i));
       fcMat = dlmread(SFILE);
       if any(i == excs)
           sMat = zeros(94,94).* NaN;
       else
%           sMat  = load(AFILE);
%           sMat =  sMat.sym_conn_muw;
       end
       nodes = dlmread(NFILE);
       nodes = nodes(2:end);
       allNodes = unique(cat(1, allNodes, nodes));
       fcMats(:,:,i) = NaN;
       fcMats(nodes,nodes,i)=fcMat;
       meanFc(end+1) = sqrt(nanmean(fcMat(:).^2));
       
%        if isempty(struMat)
%            struMat = sMat;
%        else
%            struMat = cat(3, struMat, sMat);
%        end
       
   else
       fprintf('%s does not exist\n', SFILE);
   end
    
end

% %%

% %%
% nodes=1:17;
fcMats   = fcMats(nodes, nodes, :);
% %   %%
netL1    = [1, 12,24, 43,59, 72,85, 100,108, 120, 133,143,148, 166,187,194, 200, ...
               212,223, 243,258, 272,284, 303,312, 324, 335,350,357, 373,384,390, 400 ...
               403, 415, 426, 429, 441, 452];
           

snNodes  = length(netL1);
sfcMats   = zeros(snNodes,snNodes, size(fcMats,3));


for n1 = 2:snNodes
    sN1 = netL1(n1-1):netL1(n1);
    for n2 = n1:snNodes
        sN2 = netL1(n2-1):netL1(n2);
        fc = fcMats(sN1, sN2,:);
        fc(fc==0) = NaN;
        tmp = squeeze(nanmean(squeeze(nanmean( fc ))));
        
        sfcMats(n1,n2,:) = tmp;
        sfcMats(n2,n1,:) = sfcMats(n1,n2,:);
        
    end
end

fcMats   = fcMats(1:400, 1:400, :);
% 
% % for i = 1:106
% %    clf;
% %    imagesc(squeeze(sfcMat(:,:,i)), [-.5 .5]); colormap jet;
% %    drawnow;
% %    pause(0.5);
% % end
% 
% 
% % %%
% sfcMats = sfcMats([2:33 36 ],[2:33 36],:);
 sfcMats = sfcMats([2:33 ],[2:33],:);
% % clf;
% % subplot(121);
% % imagesc(mean(sfcMat(2:end,2:end,:),3), [-.5 .5])
% % subplot(122);
% % imagesc(mean(fcMats,3), [-.5 .5]); 
% % colormap jet;
% % 
% % %%
% 
%  %%

fcm = fcMats;
nNodes = size(fcm,1);
ageMats  = zeros(nNodes); 
pageMats = zeros(nNodes); 

meanFc(exc) = [];

% inc(abs(zscore(meanFc))>2)=[];

pvals = [];
for i = 1:nNodes
    
    
    
    for j = i+1:nNodes
       
        fc = squeeze(atanh(fcm(i,j, : )));
        [r,p] = corr( fc, ages(:), 'rows', 'complete');
        
        ageMats(i,j) = r;
        ageMats(j,i) = r;
        
        pvals(end+1) = p;
        pageMats(j,i) = p < 0.05;
        
    end
end

[~,crit_p,~,pcorr] = fdr_bh(pvals);
idx = 1;
for i = 1:nNodes
    for j = i+1:nNodes
        if pcorr(idx) < 0.05
           pageMats(i,j) = 1; 
           pageMats(j,i) = 1; 
        else
            pageMats(i,j) = 0;
            pageMats(j,i) = 0;
        end
        idx = idx + 1;
    end
end
%%
ifcm = fcm(:,:,inc);
clf;
subplot(221);
imagesc(nanmean(ifcm(:,:,ages(:)<=30),3),[-.3,.3]); colormap jet;
subplot(222);
imagesc(nanmean(ifcm(:,:,ages(:)>=65),3),[-.3,.3]); colormap jet;

subplot(223);
imagesc(ageMats,[-.6,.6]); colormap jet;


subplot(224);
%%
m = atanh(ageMats .* pageMats);
degs = zeros(400,2);
for i = 1:400
    ki = length(netL1)-1;
    for k = 1:length(netL1)
        if i > netL1(k)
            ki = k -2;
        end
    end
    for j = i+1:400
        
         
        kj = length(netL1)-1;
        for k = 1:length(netL1)
            if j > netL1(k)
                kj = k -2;
            end
        end
        
        if ki == kj
            degs(i,1) = degs(i,1) + m(i,j);
        else
            degs(i,2) = degs(i,2) + m(i,j);
        end
    end
    
end
%%
nii = load_nii('/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/Results/lg400_cobra_group.nii');
zdegs = zscore(degs(:,1));


hubList = [27    32    42    45    53    73    75    78    86   149   150   155   175   225   228   232   253   276   282   285   286   287   302   304   312   326   380   384];

for i = 1:400
    if any(i == hubList)
      nii.img(nii.img == i) = degs(i);
   else
      nii.img(nii.img == i) = 0;
   end
%    if zdegs(i) >= -0.5
%       nii.img(nii.img == i) = 0;
%    else
%       nii.img(nii.img == i) = abs(degs(i));
%    end
end

nii.img(nii.img > 400) = 0;

save_nii(nii, '/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/ICAPs/Results/hubs_015.nii');

%%

nii = load_nii('/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/iCAPs_results/dFC_Full3_Alpha_5_95_Fraction_0DOT05/K_15_Dist_cosine_Folds_17/iCAPs_z.nii');

nii.img((nii.img)<3) = 0;
nii.img = nii.img(:,:,:,[7 8 10 12 13]);
nii.img = sum(nii.img,4);

nii.hdr.dime.dim(5) = 1;
nii.hdr.dime.dim(1) = 3;

save_nii(nii, '/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/Results/icaps_ya.nii');


%%
clf
netNames = {'L Visual A', 'L Visual B', 'L Somatomotor A', 'L Somatomotor B', ...
        'L Dorsal Attention A', 'L Dorsal Attention B', ...
        'L Ventral Attention A', 'L Ventral Attention B', ...
        'L Limbic', ...
        'L Control A', 'L Control B', 'L Control C',...
        'L DMN A', 'L DMN B', 'L DMN C',...
        'L Frontoparietal',  ...
        'R Visual A', 'R Visual B', 'R Somatomotor A', 'R Somatomotor B', ...
        'R Dorsal Attention A', 'R Dorsal Attention B', ...
        'R Ventral Attention A', 'R Ventral Attention B', ...
        'R Limbic', ...
        'R Control A', 'R Control B', 'R Control C',...
        'R DMN A', 'R DMN B', 'R DMN C',...
        'R Frontoparietal'};

cmap = cat(2, cat(1, linspace(0,.6,127), linspace(.1,1,127), linspace(1,1,127) ), ...
               [1,1,1]', ...
              cat(1, linspace(1,.3,127), linspace(1,0.1,127), linspace(.1,0,127) ) )';

imagesc(atanh(ageMats .* pageMats),[-.5 .5]); colormap(cmap); hold on;
set(gca, 'ytick', 1:length(netNames), 'yticklabels', netNames);
set(gca, 'xtick', 1:length(netNames), 'xticklabels', netNames, 'xticklabelrotation', 90);
set(gcf, 'Color', 'w');

for i = 1:32
%     for j = 1:32
        plot([i-.5 i-.5], [0 32.5], 'color', [.3 .3 .3 .5]);
        plot([0 32.5], [i-.5 i-.5], 'color', [.3 .3 .3 .5]);
%     end
end

plot([16.5 16.5], [0 33], 'k', 'LineWidth', 2);
plot([0, 33], [16.5 16.5], 'k', 'LineWidth', 2);
colorbar();
%%
pvals = [];
aalLabels = '/home/luna.kuleuven.be/u0101486/workspace/fmri_proc/atlas/aal/aal2/aal_labels_nocrb.txt';
fid = fopen(aalLabels);
C = textscan(fid, '%s %s %s');
fclose(fid);
aalLabels = C{1};

sageMat = zeros(94,94);
for i = 1:94
    for j = i+1:94
        cnn = squeeze(struMat(i,j,inc));
        [r,p] = corr(cnn, ages(inc));
        sageMat(i,j) = r;
        sageMat(j,i) = r;
        pvals(end+1) = p;
    end
end

[~,crit_p,~,pcorr] = fdr_bh(pvals);
idx=1;
for i = 1:94
    for j = i+1:94
        if pcorr(idx) > 0.05
            sageMat(i,j) = 0;
        end
        idx = idx + 1;
    end
end

clf;
imagesc( sageMat, [-.5 .5] ); colorbar
set(gca, 'xtick',1:94, 'xticklabels', aalLabels, 'xticklabelrotation', 90);
set(gca, 'ytick',1:94, 'yticklabels', aalLabels, 'yticklabelrotation', 0);

%%

% size(fcMats)
clf;

eigCentr = zeros(400, 92);
eigCentrP = zeros(400, 1);
% 
% %%
% k = 200;
% clc;
% fcGroup = nanmean( fcMats(:,:, ages > 65), 3 ); 
% 
% 
% fcGroup = threshold_proportional( fcGroup, 0.1 );
% rc = rich_club_wu(fcGroup, k );
% 
% 
% degs = degrees_und(fcGroup);
% nrep = 50;
% rrc = zeros(nrep, k);
% 
% parfor rep = 1:nrep
%     rfcGroup = randmio_und(fcGroup, 50);
%     rrc(rep,:) = rich_club_wu(rfcGroup, k);
% end
% 
% rrc = nanmean(rrc);
% rc(isnan(rc))= 0;
% rrc(isnan(rrc))= 1;
% clf;
% plot(rc); hold on;
% plot(rrc); hold on;
% %%
% %    fcThr = threshold_proportional( fcSub, 0.05 );

parfor i = 1:92
   i
   fcSub = fcMats(:,:,i); 
   fcThr = fcSub;
   fcThr = threshold_proportional( fcSub, 0.25 );
%     fcThr(fcThr<0)=0;

   w1 = eigenvector_centrality_und(fcThr);
   w2 = strengths_und(fcThr);
   w3 = clustering_coef_wu(fcThr);
   
   Ci = community_louvain(fcThr, 2);
   w4 = module_degree_zscore(fcThr, Ci);
%    %%
   w5 = pagerank_centrality(fcThr, 0.8, w1);
%    plot(w5);
%    %%
  
%    L = 1./fcThr;
%    L(isinf(L)) = 1000;
   
%    w4 = nanmean(distance_wei(L));
   
%    w5 = betweenness_wei(L);
   N = 20;
   w = zeros(400,1);
   [~, wi] = sort(w1, 'descend');
   w(wi(1:N)) = w(wi(1:N)) + 1;
   
   [~, wi] = sort(w2, 'descend');
   w(wi(1:N)) = w(wi(1:N)) + 1;
    
   [~, wi] = sort(w3);
   w(wi(1:N)) = w(wi(1:N)) + 1;
   
   [~, wi] = sort(w4, 'descend');
   w(wi(1:N)) = w(wi(1:N)) + 1;
   
   [~, wi] = sort(w5, 'descend');
   w(wi(1:N)) = w(wi(1:N)) + 1;
   
   eigCentr(:,i) = w;
end

for n = 1:400
   [r,p] = corr(ages, eigCentr(n,:)' );
   
   if p < 0.05
       n
       if r > 0
           eigCentrP(n) = r;
       else
           eigCentrP(n) = r;
       end
   end
end

%%

nii = load_nii('/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/Results/lg400_cobra_group.nii');
zdegs = zscore(degs(:,1));


hubList = find(degs>=58 & degs <= 82);
% hubList = find(zscore(rc./rrcs) > 0.5);
 nii.img(nii.img > 400) = 0;
 
 eigCentrP = mean( eigCentr(:,ages>65), 2 );
for i = 1:400
    if eigCentrP(i) ~= 0
      nii.img(nii.img == i) = eigCentrP(i);
   else
      nii.img(nii.img == i) = 0;
   end
%    if any( i == hubList )
%       nii.img(nii.img == i) = degs(i);
%    else
%       nii.img(nii.img == i) = 0;
%    end
end

nii.hdr.dime.glmin = -0.5;
nii.hdr.dime.glmax = -0.5;

nii.hdr.dime.cal_min = -0.5;
nii.hdr.dime.cal_max = -0.5;

save_nii(nii, '/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/ICAPs/Results/hubs_015_ya.nii');
% plot(nanmean(rcCurve,2));