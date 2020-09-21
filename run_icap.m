clc;clear all;close all;
addpath(genpath('/home/luna.kuleuven.be/u0101486/workspace/matlab/toolbox/iCAPs_latest'));

rmpath(genpath('/home/luna.kuleuven.be/u0101486/workspace/matlab/toolbox/CanlabCore-master'));
rmpath(genpath('/home/luna.kuleuven.be/u0101486/workspace/matlab/toolbox/MediationToolbox-master'));

TMPDIR='/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/';
FNAME='proc_data_mni';
TR=2.5;
SID='RS001';



dwtmode('sym') % Default mode

% EO: Set whether to use GPU implementation (= 1) or original Matlab.
param.use_cuda = 0;

% setting up all parameters to run TA

% Path where we have our data stored 
param.PathData = TMPDIR;
% TR of the data
param.TR = TR;


% Links towards the data of all subjects to analyze
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% List of subjects on which to run total activation (must be a cell array
% with all group/subject names). This is where the TA folder will be
% created (or looked for) for each subject
param.Subjects = cell(10,1);
inc  = 1:106;
%3 8 12 14 20 38 69 75 87 90 92 4 11 24 31 39 49 63 85 97 105 106
% 49, 93 --> MoCA <= 24 (cutoff score)
% 106 --> Lesion
% 105  --> Withdrew?
% 75 --> No BTT
% %%
excs = sort(unique([4 11 24 31 39 49 63 85 97 93 105 106 38 75]));

 inc(excs) = [];
idx = 1;
for i = inc
    param.Subjects{idx} = sprintf('RS%03d',i);
    
    motFile = cat(2, TMPDIR, filesep, sprintf('RS%03d/motion_estimate.par',i));
    rp = dlmread(motFile);
    
    rpSpm = cat(2, rp(:,4:6), deg2rad(rp(:,1:3)));

    
    mkdir(cat(2, TMPDIR, filesep, sprintf('RS%03d/motion',i)));
    dlmwrite( cat(2, TMPDIR, filesep, sprintf('RS%03d/motion/rp_motion_estimate.txt',i)), rpSpm );
    
    
    idx = idx + 1;
    
    
end

% %%
% param.Subjects = {SID};
param.n_subjects = length(param.Subjects);
param.title = 'dFC';

% %%
% clc;
% for i = 1:length(inc)
%    s = inc(i);
%    
%    %/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/RS002/TA_results/dFC/inputData/STD_MAP.nii
%    ifile = sprintf('/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/RS%03d/TA_results/dFC/Thresholding/Alpha_5_95_Fraction_0DOT05/SignInnov.nii',s); 
%    
%    if ~exist(ifile, 'file')
%        disp(ifile);
%    end
% end
% 
% %%
% Information about the folders where to retrieve functional and structural
% data of relevance for the pipeline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% In all cases, if only one name is given for any variable, the name is
% assumed to be the same across subjects. If there is as many names as
% there is subjects, then the name is assumed to differ across subjects. If
% a [] is provided (e.g. for a folder where to access data), then it is 
% assumed that the data lies within the main path itself without any
% additional subpath


% Name of the folder containing the functional data (if void, will directly
% look into the subject folder itself)
param.Folder_functional = '';
param.TA_func_prefix = FNAME;

% Folder where we can find the probabilistic Gray Matter maps for each
% subject
param.Folder_GM = 'anat';
param.TA_gm_prefix = 'gm_tpm_gro';


% Threshold at which we want to threshold the probabilistic gray matter map
% (values larger than this only will be included in the mask); has to lie
% between 0 and 1
param.T_gm = 0.15;

% select if morphological operations (opening and closure) should be
% run on the GM mask to remove wholes, and if yes, specify the size (in
% voxels) for opening and closing operators
param.is_morpho=0;
param.n_morpho_voxels=2;



% Number of scans to discard due to T1 equilibration effects
param.skipped_scans = 0;

% select if detrending should be run on the data, if set to 0 the fields
% 'DCT_TS' and 'Covariates' do not need to be set
param.doDetrend=1;
param.doNormalize =0;

% Detrending information: cut-off period for the DCT basis (for example,
% 128 means a cutoff of 1/128 = 0.0078 [Hz], and covariates to add (should
% be provided each as a column of 'Covariates')
param.DCT_TS = 110; 
param.Covariates = [];



% select if scrubbing should be run on the data, if 0 the fields
% 'Folder_motion', 'TA_mot_prefix', 'skipped_scans_motionfile',
% 'FD_method', 'FD_threshold' and 'interType' do not need to be set
param.doScrubbing=0;

% Folder where motion data from SPM realignment is stored, if motion data
% is taken from another programm than SPM, a text file with the 6 motion
% parameters (3 translational in mm + 3 rotational in rad) should be set as
% input here
param.Folder_motion = ['motion'];
param.TA_mot_prefix = 'rp_mot';

% Number of lines to ignore at the beginning of the motion file, if empty
% or not set, this will be equal to param.skipped_scans
param.skipped_scans_motionfile = []; 

% Motion information: type of method to use to quantify motion (choose
% between 'Power' and ), and threshold of displacement to use for each
% frame (in [mm])
param.FD_method = 'Power';
param.FD_threshold = 0.5;

% Interpolation method (i.e. 'spline' or 'linear', see interp1 for all 
% possibilities) - default (if []) is 'spline'
param.interType='spline';






% Total activation-related information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% select if you want to force TA to run, and overwriting existing files if
% it has already been done
param.force_TA_on_real=1;
param.force_TA_on_surrogate=1;


% Type of assumed hemodynamic response function type (can be 'bold' or
% 'spmhrf')
param.HRF ='bold';

% Number of outer iterations for the whole forward-backward scheme:
% this means that at most, the solutions with both regularizer 
% terms will be found and averaged five times. This is the 'k_max'
% in Farouj et al. 2016 (ISBI abstract from Younes)
param.Nit=5;

% Number of iterations for which the temporal regularization scheme
% and the spatial one are run
param.NitTemp = 500;
param.NitSpat=400;

% Tolerance threshold for convergence of the TA methods
param.tol = 1e-5;

% Weighting parameters for the temporal and the spatial TA schemes: because
% the temporal part is more elaborate and more extensively tested, Dr
% Karahanoglu and Dr Farouj opted for a 3/4 to 1/4 trade-off =P
param.weights = [0.75 0.25];

% Coefficient somehow used to multiply the regularization weights of
% temporal regularization
param.LambdaTempCoef = 1/0.8095;


% The noise estimation procedure uses a single scale
% wavelet decomposition. Here Daubechies wavelets with 4 vanishing 
% moments are used. The corresponding high pass filter is given by:

g=[0    -0.12941    -0.22414     0.83652    -0.48296];
g = g';
param.daub = g;


% Weight of spatial regularization
param.LambdaSpat=6;




param.COST_SAVE = 0; % save the costs..






% run TA
% Run_TA(param);


% %%




% if set to 1, Thresholding will be forced to run, even if already has been done
param.force_Thresholding=1;


% Alpha-level at which to look for significance of innovation signal frames
% (first element is the percentile of the lower threshold - negative
% innovations - and second element the upper threshold one - positive
% innovations)
param.alpha = [5 95];

% Fraction of voxels from the ones entering total activation for a given 
% subject that should show an innovation at the same time point, so that
% the corresponding frame is retained for iCAPs clustering
param.f_voxels = 5/100;

% Title used to create the folder where thresholding data will be saved
param.thresh_title = ['Alpha_',strrep(num2str(param.alpha(1)),'.','DOT'),'_',...
    strrep(num2str(param.alpha(2)),'.','DOT'),'_Fraction_',...
    strrep(num2str(param.f_voxels),'.','DOT')];

% Number of neighbours that must also show an innovation for a voxel to be
% retained
param.threshold_minclussize = 6;

% Number of neighbours to consider in the process
param.threshold_interconnectivity = 26;



% run Thresholding
% Run_Thresholding(param);

%  %%
% rmpath(genpath('/home/luna.kuleuven.be/u0101486/workspace/matlab/toolbox/CanlabCore-master/CanlabCore/External/spider/'));

% %%
% name of the iCAPs output for this data
% if only a subset of subjects should be included in the clustering, this
% can be useful to save those different runs in different folders
param.data_title = [param.title '_Full3'];


% iCAPs-related information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% specify if clustering should be done (one may want to only run consensus 
% clustering, then set this to 0), default = 1
param.doClustering=1;

% if set to 1, Clustering will be forced to run, even if already has been done
param.force_Aggregating=0;
param.force_Clustering=0;

% file with external mask to use to define the input voxels,
% if nothing is specified, the mask will be the intersection of the GM
% masks of all subjects that are included in the clustering
param.common_mask_file=[];

% file with additional mask to apply additionally to the intersection of
% all GM masks, it has to be in the same space as all subject's significant
% innovations (I am using this here to exclude the cerebellum)
param.extra_mask_file=[];%'GM_mask_MNI333_AAL.nii';




% Number of iCAPs into which to separate the data
param.K = [8 14 20];
param.K = [15 17];

% Type of distance to use for the k-means clustering process (choose
% between 'sqeuclidean' and 'cosine')
param.DistType = 'cosine';

% Number of times the clustering process is run in a row to extract iCAPs
param.n_folds = 17;

% specify if the result of each replicate should be saved during clustering
param.saveClusterReplicateData=0;

% Maximum number of allowed iterations of the kmeans clustering, the Matlab
% default of 100 is sometimes not enough if many frames are included,
% default = 100
param.MaxIter=300;


% save subject-specific iCAPs maps
param.saveSubjectMaps=1;

% save iCAPs regions tables
param.saveRegionTables=0;
param.regTab_thres=1.5; % z-score at which to threshold map for regions table
param.regTab_codeBook='AALcodeBook.mat'; % file with atlas region names
param.regTab_atlasFile='AAL90_correctLR.nii';%'lg400_cobra_group.nii'; % file with atlas data in MNI (has to be in same space as iCAPs results)


% Title used to create the folder where iCAPs data will be saved

if length(param.K)==1
    nK = 1;
    param.iCAPs_title=['K_',num2str(param.K(nK)),'_Dist_',...
            param.DistType,'_Folds_',num2str(param.n_folds)];
else

    for nK=1:length(param.K)
        param.iCAPs_title{nK} = ['K_',num2str(param.K(nK)),'_Dist_',...
                param.DistType,'_Folds_',num2str(param.n_folds)];
    end
end



% consensus clustering - related information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% select if consensus clustering should be run or not
param.doConsensusClustering=1;

% if set to 1, Consensus Clustering will be forced to run, even if already has been done
param.force_ConsensusClustering=0;

% Subsample Type:
% 'subjects' to subsample all frames from a subject; 
% 'items' to subsample frames without taking into account within- or
%   between-subject information
param.Subsample_type='items';
param.Subsample_fraction=0.8;
param.cons_n_folds=15;
param.cons_title=[num2str(param.K(1)) 'to' num2str(param.K(end)) ...
    '_SubsampleType_' param.Subsample_type ...
    '_Fraction_' strrep(num2str(param.Subsample_fraction),'.','DOT') ...
    '_nFolds_' num2str(param.cons_n_folds) ...
    '_Dist_' param.DistType];


% flag to indicate that cluster consensus should be computed, clustering
% and consensus clustering have to be done to compute this measure
param.computeClusterStability=1;


% run Clustering
Run_Clustering(param);

 %%


% setting up all parameters to run Clustering
%Inputs_TimeCourses_Data_OpenfMRI


% if set to 1, regression will be forced to run, even if already has been done
param.force_Regression=0;

% which type of regression should be done, 
%   'unconstrained' - as in the original paper [Karahanoglu et al., NatComm 2015]
%   'transient-informed' - recommended (default), see [Zoeller et al, IEEE TMI 2018]
param.regType='transient-informed';

% parameter for soft cluster assignment in transient-informed regression,
% can be only one value or a vector of multiple values
% for details see [Zoeller et al., IEEE TMI 2018]
param.softClusterThres=[1:0.25:2];

% choose if for the evaluation of soft cluster assignment factor the
% correlation between measured and estimated amplitudes should be computed
% and evaluated. otherwise, only the BIC and AIC results will be evaluated
param.evalAmplitudeCorrs=0;

% threshold above which a z-scored iCAPs time course will be considered
% "active" - according to Karahanoglu et al, NatComm 2015 and Zoller et
% al., IEEE TMI 2018 we select the default z-score of |1|
param.activityThres=1;

% run Clustering
Run_TimeCourses(param);

%%


addpath(genpath('/home/luna.kuleuven.be/u0101486/workspace/matlab/toolbox/CanlabCore-master'));
addpath(genpath('/home/luna.kuleuven.be/u0101486/workspace/matlab/toolbox/MediationToolbox-master'));

allPerf = load('/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/beahv_results2.mat');
t0 = 1;
tf = 48;
zzPerf = squeeze(median(allPerf.coverageSubs(:,t0:tf,4),2));
anPerf = squeeze(median(allPerf.coverageSubs(:,t0:tf,3),2));
l2Perf = squeeze(median(allPerf.coverageSubs(:,t0:tf,2),2));
l1Perf = squeeze(median(allPerf.coverageSubs(:,t0:tf,1),2));

% zzPerf = -sum( squeeze(mean(  diff(allPerf.coverageSubs(:,1:end,1:4),1,3)  ,2)),2);

zzPerf2 = squeeze(median(allPerf.bhiSubs(:,t0:tf,4),2));
anPerf2 = squeeze(median(allPerf.bhiSubs(:,t0:tf,3),2));
l2Perf2 = squeeze(median(allPerf.bhiSubs(:,t0:tf,2),2));
l1Perf2 = squeeze(median(allPerf.bhiSubs(:,t0:tf,1),2));

%length(param.Subjects)
% %%


%TC=load('/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/iCAPs_results/dFC_Full_Alpha_5_95_Fraction_0DOT05/K_14_Dist_cosine_Folds_10/TCs_1_0DOT2_2/TC.mat');
tempChar=load('/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/iCAPs_results/dFC_Full3_Alpha_5_95_Fraction_0DOT05/K_15_Dist_cosine_Folds_17/TCs_1_0DOT25_2/tempChar.mat');
excs = sort(unique([4 11 24 31 39 49 63 85 97 93 105 106 38 75]));
inc = 1:106;
inc(excs) = [];
% 2 ou 3
tempChar = tempChar.tempChar{2};
clc;
thr = 5;
meanCo  = [];
meanCo2 = [];

ages = dlmread('/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/ages.txt');
ages = ages(inc,1);
% ages(38)=[];


% %%

ya = find(ages <= 35);
ma1 = find(ages > 35 & ages <= 50);
ma2 = find(ages > 50 & ages <= 65);
oa = find(ages > 65);

pvals = [];
% 
% clf; 
% 
% visual    =  [10 14 15];
% cognitive =  [4 6 11 13 8];
% fpn       =  [5 9];
% motor     =  [ 7 12 ];
% DMN       =  [3 3];
% crb       =  [1 2];
% 
% netGroups= {visual, cognitive, fpn, DMN, motor, crb};
% netX = [1.5, 4.5, 9.5, 11.5, 13.5];
% netNums = cat(2, visual, cognitive, fpn, motor, crb);
% % for i = 1:length(netGroups)
% for i = 1:length(netNums)
% %     if i == 1
% %         subplot(1,5,1);
% %     else
% %         subplot(1,5,2:5);
% %     end
%     
%     hold on;
% %      yav = nansum(squeeze(tempChar.innov_counts(netGroups{i},ya)));
% %      ma1v = nansum(squeeze(tempChar.innov_counts(netGroups{i},ma1)));
% %      ma2v = nansum(squeeze(tempChar.innov_counts(netGroups{i},ma2)));
% %      oav = nansum(squeeze(tempChar.innov_counts(netGroups{i},oa)));
% 
%      yav = (squeeze(tempChar.innov_counts(netNums(i),ya)));
%      ma1v = (squeeze(tempChar.innov_counts(netNums(i),ma1)));
%      ma2v = (squeeze(tempChar.innov_counts(netNums(i),ma2)));
%      oav = (squeeze(tempChar.innov_counts(netNums(i),oa)));
% %     
%     h1 = bar(i-.4, mean(yav), 'barwidth', 0.2, 'faceColor', [.4 .9 1]);
%     plot([i-.4 i-.4], [mean(yav) mean(yav)+std(yav)], 'k')
%     scat = (i-.4) + (rand(length(yav),1)-.5)./10;
%     plot(scat, yav, '.', 'Color', [.1 .2 1], 'MarkerSize', 5 );
%     
%     h2 = bar(i-.2, mean(ma1v), 'barwidth', 0.2, 'faceColor', [.4 .9 1]./2);
%     plot([i-.2 i-.2], [mean(ma1v) mean(ma1v)+std(ma1v)], 'k')
%     scat = (i-.2) + (rand(length(ma1v),1)-.5)./10;
%     plot(scat, ma1v, '.', 'Color', [0 .1 .2], 'MarkerSize', 5 );
%     
%     
%     h3 = bar(i+.0, mean(ma2v), 'barwidth', 0.2, 'faceColor', [1 .8 .3]./2);
%     plot([i-.0 i-.0], [mean(ma2v) mean(ma2v)+std(ma2v)], 'k')
%     scat = (i-.0) + (rand(length(ma2v),1)-.5)./10;
%     plot(scat, ma2v, '.', 'Color', [.2 .1 0], 'MarkerSize', 5 );
%     
%     
%     h4 = bar(i+.2, mean(oav), 'barwidth', 0.2, 'faceColor', [1 .8 .3]);
%     plot([i+.2 i+.2], [mean(oav) mean(oav)+std(oav)], 'k')
%     scat = (i+.2) + (rand(length(oav),1)-.5)./10;
%     plot(scat, oav, '.', 'Color', [.2 .1 0], 'MarkerSize', 5 );
% %     %%
%     grps = cat(2, zeros(1, length(yav)), zeros(1, length(ma1v))+1, ...
%         zeros(1, length(ma2v))+2,zeros(1, length(oav))+3);
%     tbl = table(cat(2,yav, ma1v, ma2v, oav)', ages([ya; ma1; ma2; oa]), 'VariableNames', {'Innovation','AgeGroup'});
%     
%     mdl = fitlme(tbl, 'Innovation ~ 1 + AgeGroup', 'FitMethod', 'REML'  );
%     man = mdl.anova();
%     
%    
%     pvals(end+1) = mdl.Coefficients.pValue(2);%man.pValue(2);
%     
%     
%     %[p,tbl,stats] = anova1(tbl);
%     
% end
% 
% [~,~,~,pcorr] = fdr_bh(pvals, 0.05);
% % pcorr = bonf_holm(pvals);
% 
% % for i = 1:length(netGroups)
% for i = 1:length(netNums)
%    if pcorr(i) < 0.05
% %        if i == 1
% %            subplot(1,5,1);
%            plot(i-.2, 150, 'k*', 'MarkerSize', 10, 'LineWidth', 1.5); 
% %            set(gca, 'xtick', 1-.2, 'xticklabels', 1);
% %            grid on;
%            
% %        else
% %            subplot(1,5,2:5);
%            plot(i-.2, 150, 'k*', 'MarkerSize', 10, 'LineWidth', 1.5); 
%            
% %        end
%        
%     end 
% end
% ylabel('Innovation Counts [Volumes, Avg within ICAP group]');
% % set(gca, 'xtick', [1:length(netGroups)], 'LineWidth', 1, ...
% %     'xticklabels',{'Visual', 'Cognitive', 'FPN', 'DMN', 'Motor', 'Other'});
% 
% set(gca, 'xtick',netX, 'LineWidth', 1, ...
%     'xticklabels',{'Visual', 'Cognitive', 'FPN', 'Motor', 'Other'});
% 
% grid on;
% ylim([0 160]);
% xlim([0 15]);
% xlabel('ICAP [Grouped]');
% 
% ya = find(ages <= 35);
% ma1 = find(ages > 35 & ages <= 50);
% ma2 = find(ages > 50 & ages <= 65);
% oa = find(ages > 65);
% % 
% % legend([h1,h2,h3,h4],{'YA [Age <= 35yo]', 'MA1 [35 < Age <= 50]', 'MA2 [50 < Age <= 65]', 'OA [Age > 65]'});
% set(gcf, 'Color', 'w');
% 
% % netNums([1 7 10 11])
% %%
% 
% % lgAtlas = load_nii('/home/luna.kuleuven.be/u0101486/workspace/fmri_proc/atlas/CAT12/lg400_cobra.nii');
% icapZ   = load_nii('/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/iCAPs_results/dFC_Full3_Alpha_5_95_Fraction_0DOT05/K_15_Dist_cosine_Folds_17/iCAPs_z.nii');
% 
% icapZ.img = icapZ.img(:,:,:,[10 8 13 7 12]);
% 
% icapZ.img(abs(icapZ.img) < 3) = 0;
% icapZ.img = sum(icapZ.img, 4);
% 
% icapZ.hdr.dime.dim([1 5]) = [3 1];
% 
% save_nii(icapZ, '/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/ICAPs/Results/icap_thresh.nii');
% 
% 
%  %%


innovList = zeros(15, 92);
innovList2 = zeros(15, 92);
coupMat = zeros(15,15, 92);

mediators = zeros(92, 7);

for s = 1:92

    visual    =  [10 14 15];
    cognitive =  [4 6 11 13 8];
    fpn       =  [5 9];
    DMN       =  [3 3];
    motor     =  [ 7 12 ];
    crb       =  [1 2];
    all = cat(2, crb, visual, cognitive, DMN, motor, fpn);
    all2 = cat(2, visual, cognitive,  motor, fpn);
    all3 = cat(2, crb);

%     %%
    
    
%     %%
     nodes = cat(2, DMN, visual);
     
%      cc = sum(tempChar.duration_avg_pos_counts(cat(2, 3),s));
%         cc = sum(tempChar.duration_avg_counts(cat(2, DMN),s));
%         cc2 = sum(tempChar.duration_avg_counts(cat(2, cognitive),s));
%       cc = sum(tempChar.duration_total_counts(cat(2, cognitive),s));
%         cc = sum(tempChar.duration_total_counts(cat(2, cognitive),s));
%          cc = sum(tempChar.occurrences_neg(cat(2, DMN),s));
%         cc2 = sum(tempChar.occurrences_neg(cat(2, visual(3)),s));
%     cc = nanmean(nansum(tempChar.coupling_counts(8, cognitive, s)));
%     cc = nanmean(nansum(tempChar.coupling_counts(8, [1:7  9:15 ], s)));
%       cc = nanmean(nansum(tempChar.coupling_counts(cat(2, fpn, cognitive, motor), cat(2, fpn, cognitive, motor), s)));
      cc = nanmean(nansum(tempChar.coupling_diffSign_counts( cat(2, visual, DMN), cat(2, visual, DMN), s)));
     
%         cc2 = sum( tempChar.coactiveiCAPs_total{s} >= 3);       
%         cc2 = sum( tempChar.coactiveiCAPs_total{s} <= 1);

        tcAct = tempChar.TC_active{s};
%         coact = sum(tcAct(3:end,:) > 0)';
        coact = sum(tcAct(cat(2, cognitive, fpn, motor),:)>0 )' ;
         cc  = sum(coact >= 2 );
%         for t = 1:162
%             if coact(t) == N
%                 
%                 tcActTP = tcAct(3:15,t);
%                 comb = find(tcActTP > 0)+2;
% 
% 
%                 allCombs(end+1,:) = sprintf( '%02d|', sort(comb) );
%             end
% 
%         end
        
%       % Total innovation
%          cc = sum(tempChar.innov_counts_total(s));
%         cc2 = sum(tempChar.innov_counts_total(s));

%         % Average duration of icap 
%          cc = mean(tempChar.duration_total_counts(cat(2,cognitive),s));


%       % Percentage duration of active block per ICAP
%          cc = mean(tempChar.duration_total_perc(cat(2, 3),s));
      
%       % Total duration of activation per ICAP
%        cc = mean(tempChar.duration_total_pos_perc(cat(2, fpn),s));
%        cc = mean(tempChar.duration_total_pos_perc(cat(2, [7 8 10 12]),s));


%        % Number of frames when Innovation occurs
%           cc = sum(tempChar.innov_counts(cat(2,  fpn, cognitive, motor  ),s));
%          nn = 1:15;
%          nn(cat(2, 1:2,cognitive))=[];
%          cc = sum(tempChar.innov_counts(cat(2, [1 2]),s));
        
       
% %        % Number of frames when Innovation occurs as % of n innovation
           
%             cc = sum(tempChar.innov_counts(cat(2, DMN),s));
%            cc2 = sum(tempChar.innov_counts(cat(2, visual),s));
           
%           cc = sum(tempChar.innov_counts(cat(2,  cognitive),s));
%           cc = sum(tempChar.innov_counts(cat(2, [1 2]),s)) - ...
%                   sum(tempChar.innov_counts(cat(2, [7 8 10 12 ]),s));
%        % Number of activation blocks

% 
%       cc = sum(tempChar.occurrences(cat(2, fpn),s));
       
%     nodes = cat(2,  crb);
%     cc = nansum(nansum(tempChar.coupling_jacc(nodes, nodes, s)));
    
%     cc = (tempChar.TC_norm_thes{s});
%     cc = std(cc,0,2);
%     cc = mean(cc(cat(2, motor, visual, cognitive)));
    
    mediators(s,1) = sum(tempChar.innov_counts(cat(2,  visual),s));
    mediators(s,2) = sum(tempChar.innov_counts(cat(2,  cognitive),s));
    mediators(s,3) = sum(tempChar.innov_counts(cat(2,  fpn),s));
    mediators(s,4) = sum(tempChar.innov_counts(cat(2,  DMN),s))/2;
    mediators(s,5) = sum(tempChar.innov_counts(cat(2,  motor),s));
    mediators(s,6) = sum(tempChar.innov_counts(cat(2,  crb),s));
    mediators(s,7) = tempChar.innov_counts(motor(1), s);
    mediators(s,8) = tempChar.innov_counts(motor(2), s);
    
     meanCo(end+1) = cc;
     meanCo2(end+1) = cc2;
     
       
end
% 
% clf; clc;
% 
% ai = find(ages>60);
% meanCo = meanCo(ai)';
% meanCo2 = meanCo2(ai)';
% % [~,~,resid] = regress(meanCo, ages);
% z1 = abs(zscore(meanCo));
% z2 = abs(zscore(meanCo2));
% z = z1 < 2 & z2 < 2;
% 
% plot(meanCo(z), meanCo2(z), '.k');
% [r,p] = corr( meanCo(z), meanCo2(z) )
% 
% 
% %%
perfNon = (l1Perf + l2Perf);
perfSwi = (anPerf + zzPerf);

perfNon2 = (l1Perf2 + l2Perf2);
perfSwi2 = (anPerf2 + zzPerf2);

perf=(l1Perf + l2Perf) - (anPerf + zzPerf);
perf = perf(inc);

perfNon = perfNon(inc);
perfSwi = perfSwi(inc);
perfNon2 = perfNon2(inc);
perfSwi2 = perfSwi2(inc);


 [~,~,resid] = regress(meanCo', ages);

 ai = find((zscore(detrend(resid)))<3 );
% ai = find(ages<650);
dlmwrite('/home/luna.kuleuven.be/u0101486/workspace/rstudio/CRUNCH/mediators_f.txt', mediators(ai, [1,2,3,4,5,6,7,8]));
dlmwrite('/home/luna.kuleuven.be/u0101486/workspace/rstudio/CRUNCH/ages_f.txt', ages(ai));

dlmwrite('/home/luna.kuleuven.be/u0101486/workspace/rstudio/CRUNCH/perf_coverage_f.txt', cat(2, perfNon(ai), perfSwi(ai), perfNon(ai)-perfSwi(ai)));
dlmwrite('/home/luna.kuleuven.be/u0101486/workspace/rstudio/CRUNCH/perf_bhi_f.txt', cat(2, perfNon2(ai), perfSwi2(ai), perfSwi2(ai)-perfNon2(ai)));



% %%
meanCo = meanCo';
% meanCo2 = meanCo2';
% plot(meanCo)
% [r,p] = corr(ages,meanCo);
% clc;
%  innovMat = zeros(15,15);
% for k1 = 1:15
%     for k2 = k1+1:15
% %            s1 = coupMat(k1,k2,ages<35);
% %            s2 = coupMat(k1,k2,ages>65);
%            
% %            [h,p,ci, stats] = ttest2(s1,s2);
% %            [h,p,ci, stats] = ttest(squeeze(coupMat(k1,k2,:)));
%         [r,p] = partialcorr(innovList(k1,:)',innovList2(k2,:)',ages, 'rows', 'complete', 'type', 'spearman');
%          innovMat(k1,k2) = r;
%          innovMat(k2,k1) = r * (p < 0.05);
%     end
% end
% 
% clf;
% imagesc(innovMat, [-.5 .5]); colormap jet;
% set(gca, 'xtick', 1:15);
% set(gca, 'ytick', 1:15);
% grid on;

% ai = find(ages<30);
% [r,p] = partialcorr(meanCo(ai),meanCo2(ai),ages(ai), 'rows', 'complete')
% [r,p] = corr(meanCo(ai),meanCo2(ai), 'rows', 'complete')
% %%

% ylim([300 650]);
% fprintf( 'Corr = %.3f (p = %.5f)\n', r, p );

% text(20, 310, sprintf( 'r = %.3f (p = %.5f)\n', r, p ), 'fontweight', 'bold', 'fontsize', 12);


%  %%

% addpath(genpath('/home/luna.kuleuven.be/u0101486/workspace/matlab/toolbox/CanlabCore-master'));
% addpath(genpath('/home/luna.kuleuven.be/u0101486/workspace/matlab/toolbox/MediationToolbox-master'));

%
% clf; clc;

perf=(l1Perf + l2Perf) - (anPerf + zzPerf);
% perf=(l1Perf2 + l2Perf2) - (anPerf2 + zzPerf2);

perf = anPerf+zzPerf;
% perf = l1Perf+l2Perf;
% perf=(anPerf2 + zzPerf2);
% clf; plot(perf-perf2);
% %%
% perf=zzPerf;
% perf=l2Perf;
perf2=l1Perf;
% perf2 =  purdueLR';
perf = perf(inc);
perf2 = perf2(inc);
% plot(ages, perf, 'k.', 'MarkerSize', 10); hold on;
% plot(meanCo, meanCo2, 'r.', 'MarkerSize', 20);
% [r,p] = corr( ages, meanCo' );
% fprintf('Age vs dFC: %.3f [%.3f]\n', r,p);
% %%
% xlabel('Age [Years]');
% ylabel('Average coupling counts [Volumes]');
ai = find(ages>60);
 [~,~,resid] = regress(meanCo, ages);

%  ai = find((zscore(detrend(resid)))<3 );

clf; clc;
fo = fit(ages, meanCo, 'poly1');
fo2 = fit(ages(ai), meanCo(ai), 'poly1');
x = (min(perf):0.001:max(perf))';
x = (55:80)';

% 
% plot( x, fo.p1 .* x + fo.p2, 'r', 'LineWidth', 2); hold on;
plot( x, fo2.p1 .* x + fo2.p2, 'b', 'LineWidth', 2); hold on;

% if p < 0.05
%     plot( x, fo.p1 .* x + fo.p2, 'r', 'LineWidth', 2); hold on;
% else
%     plot( x, fo.p1 .* x + fo.p2, 'LineWidth', 2, 'Color', [.6 .6 .6]); hold on;
% end

% plot( ages, meanCo, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', [1 .4 .4] );
plot( ages(ai), meanCo(ai), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', [.4 .8 .8] );

ylabel('Average PCC/Precun ICAP Duration');
xlabel('Age');
set(gca, 'box', 'off', 'linewidth',2);
set(gcf, 'Color', 'w');

% xlim([min(perf)*0.8 max(perf)*1.2]);
xlim([55 80]);

%  %%
% % clc;
% close all;
%   1 a   X -> M relationship
%  2 b   M -> Y relationship
%  3 cp  unmediated X -> Y relationship (residual)
%  4 c   X -> Y relationship
%  5 ab  mediated X -> Y by M (a * b)




mediation(ages(ai), perf(ai), meanCo(ai), 'plots0', 'boottop', 'verbose', 'robust0', 'bootsamples', 5000, 'names', {'Age', 'dPerf', 'Volumes w/o ICAP'});


[r,p] = corr( ages(ai), meanCo(ai) );
fprintf('Age vs dFC: %.3f [%.3f]\n', r,p);
[r,p] = corr( perf(ai), meanCo(ai) );
fprintf('Perf vs dFC: %.3f [%.3f]\n', r,p);

[r,p] = partialcorr( perf(ai), meanCo(ai), ages(ai) );
fprintf('Perf vs dFC [partial Age]: %.3f [%.3f]\n', r,p);
disp(size(ai))

rmpath(genpath('/home/luna.kuleuven.be/u0101486/workspace/matlab/toolbox/CanlabCore-master'));
rmpath(genpath('/home/luna.kuleuven.be/u0101486/workspace/matlab/toolbox/MediationToolbox-master'));

%%

coactPct = zeros( 16, 92 );
for s = 1:92
    tcAct = tempChar.TC_active{s};
    % clf;
    coact = sum(tcAct(3:end,:)~=0);
    for k = 0:15
        coactPct(k+1,s) = 100*sum(coact==k)/162;
    end
%     plot(sum(tcAct~=0))
end
%%

clf;
hold on;
g = find(ages<35);
x = linspace(0.5, 18, 16);
y = nanmean(coactPct(:,g), 2);
ys = nanstd(coactPct(:,g), 1, 2) ./ sqrt(length(g));
h1=bar( x, y, 'barwidth', 0.25, 'FaceColor', [.6 .6 1] ); hold on
for i = 1:16
    plot([x(i) x(i)], [y(i) y(i)+ys(i)], 'k');
    plot([x(i) x(i)], [y(i) y(i)-ys(i)], 'k');
end

g = find(ages>65);
x = linspace(0.8, 18.3, 16);
y = nanmean(coactPct(:,g), 2);
ys = nanstd(coactPct(:,g), 1, 2) ./ sqrt(length(g));
h2=bar( x, y, 'barwidth', 0.25, 'FaceColor', [1 .8 .2] ); hold on
for i = 1:16
    plot([x(i) x(i)], [y(i) y(i)+ys(i)], 'k');
    plot([x(i) x(i)], [y(i) y(i)-ys(i)], 'k');
end

set(gca, 'xtick', linspace(0.65, 18.15, 16), 'xticklabels', 0:15);
set(gcf, 'color', 'w');
xlabel('Number of Co-Active ICAPs');
ylabel('Percentage of Scan Length');

legend([h1,h2], {'Age < 35', 'Age > 65'});
%%

% t = tiledlayout(2,1);
gTitle = {'YA', 'MA1', 'MA2', 'OA'};
% ax1 = nexttile;
% pie(ax1,y2010)
% legend(labels)
% title('2010')
% 
% ax2 = nexttile;
% pie(ax2,y2011)
% legend(labels)
% title('2011')
N = 3;
clc; clf;
for g = 5%1:4
%     ax = subplot(2,2,g);
    allCombs = '';
    for s = 1:92
        if g == 1 && ages(s) > 35, continue; end
        if g == 2 && (ages(s) < 35 || ages(s) > 50), continue; end
        if g == 3 && (ages(s) < 50 || ages(s) > 65), continue; end
        if g == 4 && ages(s) < 65, continue; end
        tcAct = tempChar.TC_active{s};
        coact = sum(tcAct(3:end,:) > 0);
        for t = 1:162
            if coact(t) == N
                
                tcActTP = tcAct(3:15,t);
                comb = find(tcActTP > 0)+2;


                allCombs(end+1,:) = sprintf( '%02d|', sort(comb) );
            end

        end
    end

    uniqueCombs = unique(allCombs, 'rows');
    topCombs = zeros(length(uniqueCombs),1);
    for k = 1:length(uniqueCombs)
        for k2 = 1:length(allCombs)
            topCombs(k,1) = topCombs(k,1) + (count(allCombs(k2,:), uniqueCombs(k,:)) ./ size(allCombs,1) );
        end
    end

    [~, topIdx] = sort(topCombs, 'descend');
    topSum = 0;

    pieX = [];
    labels = {};
    for i = 1:5
        pieX(end+1) = topCombs(topIdx(i))*100;
        topSum = topSum + topCombs(topIdx(i))*100;
        fprintf( '%s - %.1f\n', uniqueCombs(topIdx(i),:), topCombs(topIdx(i))*100 );

        labels{end+1} = cat(2, strrep(uniqueCombs(topIdx(i),:),'|',' '), ...
            sprintf(' [%0.2f%%]',pieX(end)));
    end
    
    ns = [];
    lab = split(labels{1}, ' ');
    
    for w = 1:N
        ns(end+1) = str2double(lab{w});
    end
    
    
    subPatt = zeros(92,1);
    for s = 1:92
        tcAct = tempChar.TC_active{s};
        coact = sum(tcAct(3:end,:) > 0);
        for t = 1:162
            subCombs = 0;
            if coact(t) == N
                tcActTP = sum(tcAct(ns,t));
                if tcActTP == N
                    subPatt(s) = subPatt(s) + 1;
                end
            end

        end
    end

    
    
%     pieX(end+1) = 100-topSum;
%     labels{end+1} = sprintf('Others [%.2f%%]', pieX(end));
    pie( pieX,  labels );

%     legend(labels);
    colormap('parula');
%     title(gTitle{g})
end

set(gcf, 'Color', 'w');

clf;
plot(ages, subPatt, '.');
[r,p] = corr(ages, subPatt )
%%
g = find(ages<35);
cc8 = squeeze(nanmean((tempChar.coupling_diffSign_counts(8, :, g)),3));
cc3 = squeeze(nanmean((tempChar.coupling_diffSign_counts(3, :, g)),3));
cc3s = squeeze(nanstd((tempChar.coupling_diffSign_counts(3, :, g)),1,3)) ./ sqrt(length(g));
cc8s = squeeze(nanstd((tempChar.coupling_diffSign_counts(8, :, g)),1,3)) ./ sqrt(length(g));

cc3(3) = 0;
cc3s(3) = 0;

cc8(8) = 0;
cc8s(8) = 0;

clf;
x = linspace(0.5, 17, 15);
bar( x, cc3, 'barwidth', 0.25, 'FaceColor', [.6 .6 1] ); hold on
for i = 1:15
    plot([x(i) x(i)], [cc3(i) cc3(i)+cc3s(i)], 'k');
    plot([x(i) x(i)], [cc3(i) cc3(i)-cc3s(i)], 'k');
end

x = linspace(0.8, 17.3, 15);
bar( x, cc8, 'barwidth', 0.25, 'FaceColor', [1 .8 .2] ); hold on
for i = 1:15
    plot([x(i) x(i)], [cc8(i) cc8(i)+cc8s(i)], 'k');
    plot([x(i) x(i)], [cc8(i) cc8(i)-cc8s(i)], 'k');
end

legend({'DMN', 'Precun'});
set(gca, 'xtick', linspace(0.65, 17.15, 15), 'xticklabels', 1:15);
% bar( 0.5:14.5, cc8 );

%%
clf;
ymat = mean(squeeze(coupMat(1:15,1:15,ages<35)),3);
omat = mean(squeeze(coupMat(1:15,1:15,ages>65)),3);
% subplot(121);

imagesc(omat-ymat, [-2 2]); colormap jet; colorbar;
%%

genders = [1 1 1 1 1  1 1 -1 1 1 1 -1 -1 -1 -1 -1 -1 -1 1 -1 1 -1 1 1 1 -1 -1 -1 -1 1 1 1 -1 -1 -1 -1 -1 ...
    1 -1 1 -1 -1 -1 -1 1 -1 1 -1 -1 1 -1 -1 -1 1 -1 -1 -1 -1 1 1 1 -1 -1 1 ...
    1 -1 -1 1 1 -1 1 -1 1 -1 1 1 1 1 1 -1 1 -1 1 1 1 1 -1 1 1 -1 1 -1 -1 -1 -1 -1 -1 -1 1 -1 1 -1 -1 1 -1 -1];

dlmwrite('/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/genders.txt', genders);


allNiis = [];
TIVs = [];

delete('/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/VBM/*.gii');
delete('/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/VBM/*.nii');
delete('/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/VBM/SPM.mat');
files = cell(92,1);

measureName = 'gyrification';
smoo = 20;

qualExc = [];
for i_ = 1:92
   i = inc(i_);
   
   if any(i == [48 87])
%        qualExc(end+1) = i_;
   end
   
   i_
%    fIn   = sprintf('/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/RS%03d/anat_cat12/mri/mwp1t1.nii', i);
   fMeas = sprintf('/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/RS%03d/anat_cat12/report/cat_t1.mat', i);
   fIn = sprintf('/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/RS%03d/REHO.nii', i);
   mIn = sprintf('/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/RS%03d/anat/gm_mask_group.nii', i);
%    fIn   = sprintf('/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/RS%03d/anat_cat12/surf/s%d.mesh.%s.resampled_32k.t1.gii', i, smoo, measureName);
   fOut = sprintf('/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/VBM/%03d.nii', i);
%     fOut = sprintf('/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/VBM/%03d_s%d.mesh.%s.resampled_32k.t1.gii', i, smoo, measureName);
    
%    catRep = load(fMeas);
   
%    TIVs(end+1) = catRep.S.subjectmeasures.vol_TIV;
   
    nii = load_nii(fIn);
    mii = load_nii(mIn);
    nii.img = nii.img(:,:,:,1);
    nii.img(isnan(nii.img)) = 0;
    nii.img(isinf(nii.img)) = 0;
%    
    tmp = nii.img;
% % %    
    spm_smooth(tmp, tmp, [4, 4, 4]);
%    
    if i_==1
        allNiis = mii.img;
    else
        allNiis = cat(4, allNiis, mii.img);
    end
%    
    nii.img = tmp;
    nii.hdr.dime.dim([1 5]) = [3 1];
    save_nii(nii, fOut);
%     copyfile(fIn, fOut);
    
    
    files{i_} = fOut;
end

nii.img = nanmean(allNiis > 0,4) > 0.8;
nii.hdr.dime.dim([1 5]) = [3 1];

save_nii(nii, '/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/VBM/Reho_Mask.nii');
%%
qualInc = 1:92;
qualInc(qualExc) = [];
% nii.img = allNiis;
% nii.hdr.dime.dim([1 5]) = [4 92];
% %%
% save_nii(nii, '/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/VBM/AllRehs.nii');
delete('/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/VBM/SPM.mat');
matlabbatch = {};
matlabbatch{1}.spm.stats.factorial_design.dir = {'/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/VBM'};

matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = cellstr(files(qualInc));
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.c = ages(qualInc);

matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.cname = 'Age';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

% matlabbatch{2}.spm.tools.cat.tools.check_SPM.spmmat = {'/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/VBM/SPM.mat'};
% matlabbatch{2}.spm.tools.cat.tools.check_SPM.check_SPM_cov.do_check_cov.use_unsmoothed_data = 1;
% matlabbatch{2}.spm.tools.cat.tools.check_SPM.check_SPM_cov.do_check_cov.adjust_data = 1;
% matlabbatch{2}.spm.tools.cat.tools.check_SPM.check_SPM_ortho = 1;

spm_jobman('run', matlabbatch);
%  %%
SPM = load('/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/VBM/SPM.mat');
SPM.SPM.swd ='/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/VBM';
cat_stat_spm(SPM.SPM);

matlabbatch = {};
matlabbatch{1}.spm.stats.con.spmmat = {'/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/VBM/SPM.mat'};
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Age Decrease';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [0 -1];
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;
spm_jobman('run', matlabbatch);

% %%
matlabbatch = {};
matlabbatch{1}.spm.tools.cat.tools.T2x_surf.data_T2x = {'/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/VBM/spmT_0001.gii'};
matlabbatch{1}.spm.tools.cat.tools.T2x_surf.conversion.sel = 2;
matlabbatch{1}.spm.tools.cat.tools.T2x_surf.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{1}.spm.tools.cat.tools.T2x_surf.conversion.inverse = 1;
matlabbatch{1}.spm.tools.cat.tools.T2x_surf.conversion.cluster.none = 1;

matlabbatch{2}.spm.tools.cat.tools.T2x_surf.data_T2x = {'/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/VBM/spmT_0001.gii'};
matlabbatch{2}.spm.tools.cat.tools.T2x_surf.conversion.sel = 2;
matlabbatch{2}.spm.tools.cat.tools.T2x_surf.conversion.threshdesc.fwe.thresh05 = 0.05;
matlabbatch{2}.spm.tools.cat.tools.T2x_surf.conversion.inverse = 1;
matlabbatch{2}.spm.tools.cat.tools.T2x_surf.conversion.cluster.none = 1;
spm_jobman('run', matlabbatch);
%%


for i_ = 1:92
   i = inc(i_);
   
   if any(i == [48 87])
%        qualExc(end+1) = i_;
   end
   
   i_
   fIn   = sprintf('/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/RS%03d/anat_cat12/mri/mwp1t1.nii', i);
   fMeas = sprintf('/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/RS%03d/anat_cat12/report/cat_t1.mat', i);
   
   break;
end

%%

nii = load_nii('/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/iCAPs_results/dFC_Full3_Alpha_5_95_Fraction_0DOT05/K_15_Dist_cosine_Folds_17/iCAPs_z.nii');

for i = 1:15
   inii = nii;
   inii.hdr.dime.dim([1 5]) = [3 1];
   
   inii.img = nii.img(:,:,:,i);
   
   iniip = inii;
   iniin = inii;
   iniin.img = iniin.img  * -1;
   
   
   iniip.img(iniip.img < 1.5) = 0;
   iniin.img(iniin.img < 1.5) = 0;
   
   if sum(iniip.img(:) > 0.5)
      fout = sprintf('/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/iCAPs_results/dFC_Full3_Alpha_5_95_Fraction_0DOT05/K_15_Dist_cosine_Folds_17/ICAPS_3d/ICAP_p_%d.nii', i);
      save_nii(iniip, fout);
   end
  
   if sum(iniin.img(:) > 0.5)
      fout = sprintf('/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/iCAPs_results/dFC_Full3_Alpha_5_95_Fraction_0DOT05/K_15_Dist_cosine_Folds_17/ICAPS_3d/ICAP_n_%d.nii', i);
      save_nii(iniin, fout);
   end
    
end

%%
cd '/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/GABA';
clc
baseDir = '/media/u0101486/MONTEIRO/Day_A/';
baseDirT1 = '/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/';
outDir = '/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/GABA/';
subs = dir(cat(2, baseDir, '*'));

gabaFilesAct = {};
gabaFilesRef = {};
t1s = {};
fsizes = [];
skipG = [];

load('workspace.mat');
for i = 75:length(subs)
    
    baseSub = cat(2, baseDir, subs(i).name);
    
    gabaFiles = dir(cat(2, baseSub, filesep, '*NSA_PRESS_68*.SDAT'));
    
    t1 = sprintf( 'RS%03d', i-2 );
    
    if strcmp(t1, 'RS004' ) || strcmp(t1, 'RS011' ) || strcmp(t1, 'RS024' ) || ...
       strcmp(t1, 'RS031' ) || strcmp(t1, 'RS038' ) || strcmp(t1, 'RS039' ) || ...
       strcmp(t1, 'RS049' ) || strcmp(t1, 'RS063' ) || strcmp(t1, 'RS075' ) || ...
       strcmp(t1, 'RS085' ) || strcmp(t1, 'RS093' ) || strcmp(t1, 'RS097' ) || ...
       strcmp(t1, 'RS105' ) || strcmp(t1, 'RS106' ) 
        continue;
    end
    
    t1 = cat(2, baseDirT1, t1, filesep, 'anat/anat_proc.nii');
    
    
    m1Act = cat(2, baseSub, filesep, gabaFiles(1).name);
    m1ActP = cat(2, baseSub, filesep, strrep(gabaFiles(1).name, 'SDAT', 'SPAR'));
    
    m1Ref = cat(2, baseSub, filesep, gabaFiles(2).name);
    m1RefP = cat(2, baseSub, filesep, strrep(gabaFiles(2).name, 'SDAT', 'SPAR'));
    
    occAct = cat(2, baseSub, filesep, gabaFiles(3).name);
    occActP = cat(2, baseSub, filesep, strrep(gabaFiles(3).name, 'SDAT', 'SPAR'));
    
    occRef = cat(2, baseSub, filesep, gabaFiles(4).name);
    occRefP = cat(2, baseSub, filesep, strrep(gabaFiles(4).name, 'SDAT', 'SPAR'));
    
   
%     copyfile(m1Act, sprintf('%s/%03d_m1_act.SDAT', outDir, i-2));
%     copyfile(m1ActP, sprintf('%s/%03d_m1_act.SPAR', outDir, i-2));
%     
%     copyfile(m1Ref, sprintf('%s/%03d_m1_ref.SDAT', outDir, i-2));
%     copyfile(m1RefP, sprintf('%s/%03d_m1_ref.SPAR', outDir, i-2));
%     
%      copyfile(occAct, sprintf('%s/%03d_occ_act.SDAT', outDir, i-2));
%      copyfile(occActP, sprintf('%s/%03d_occ_act.SPAR', outDir, i-2));
%      
%      copyfile(occRef, sprintf('%s/%03d_occ_ref.SDAT', outDir, i-2));
%      copyfile(occRefP, sprintf('%s/%03d_occ_ref.SPAR', outDir, i-2));
%     
%     copyfile(t1, sprintf('%s/%03d_t1.nii', outDir, i-2));

    s=dir(sprintf('%s/%03d_occ_ref.SDAT', outDir, i-2));
    

%     if s.bytes > 500000
        gabaFilesAct = sprintf('%s/%03d_occ_act.SDAT', outDir, i-2);
        gabaFilesRef = sprintf('%s/%03d_occ_ref.SDAT', outDir, i-2);
    
%         t1s{end+1} = sprintf('%s/%03d_t1.nii', outDir, i-2);
     try
     MRS_struct2 = GannetLoad({gabaFilesAct}, {gabaFilesRef});
     MRS_struct2 = GannetFit(MRS_struct2);
     save(sprintf('Occ_Gaba_%03d.mat', i-2), 'MRS_struct2');
     catch err
     end
        
%     else
%         skipG(end+1) = i-2;
%     end
end
%  %%

% GABA
% MRS_struct2 = GannetLoad(gabaFilesAct, gabaFilesRef);
% MRS_struct2 = GannetFit(MRS_struct2);

% save('Occ_Gaba.mat', 'MRS_struct2');
%MRS_struct = GannetCoRegister(MRS_struct, t1s);
%MRS_struct = GannetSegment(MRS_struct);
%%



%%
gabaOcc = zeros( 1, 106) .* NaN;
glxOcc = zeros( 1, 106) .* NaN;

parfor i = 1:106
    fin = sprintf('Occ_Gaba_%03d.mat', i);
    
    if exist(fin, 'file')
        MRS = load(fin);
        MRS_struct = MRS.MRS_struct2;
        if isreal(MRS_struct.out.vox1.GABA.ConcIU)
            gabaOcc(i) = MRS_struct.out.vox1.GABA.ConcIU;
        end
        
        if isreal(MRS_struct.out.vox1.Glx.ConcIU)
            glxOcc(i) = MRS_struct.out.vox1.Glx.ConcIU;
        end
    end
end

glxOcc = glxOcc(inc);
gabaOcc = gabaOcc(inc);


%%
MRS = load('M1_Gaba.mat');
MRS_struct = MRS.MRS_struct;



%%

addpath(genpath('/home/luna.kuleuven.be/u0101486/workspace/matlab/toolbox/CanlabCore-master'));
addpath(genpath('/home/luna.kuleuven.be/u0101486/workspace/matlab/toolbox/MediationToolbox-master'));
%     mediators(s,1) = sum(tempChar.innov_counts(cat(2,  visual),s));
%     mediators(s,2) = sum(tempChar.innov_counts(cat(2,  cognitive),s));
%     mediators(s,3) = sum(tempChar.innov_counts(cat(2,  fpn),s));
%     mediators(s,4) = sum(tempChar.innov_counts(cat(2,  DMN),s))/2;
%     mediators(s,5) = sum(tempChar.innov_counts(cat(2,  motor),s));
%     mediators(s,6) = sum(tempChar.innov_counts(cat(2,  crb),s));
M = dlmread('/home/luna.kuleuven.be/u0101486/workspace/rstudio/CRUNCH/mediators_f.txt');

glxM1  = MRS_struct.out.vox1.Glx.ConcIU(inc);
gabaM1 = MRS_struct.out.vox1.GABA.ConcIU(inc);

gaba = glxM1;
gaba(isnan(gaba)) = 0;
zgaba = MRS_struct.out.vox1.Glx.FitError(inc);% abs(zscore(gaba));
% zgaba = (abs(zscore(zgaba)));
% zgaba = abs(gaba);
thr = 4;
thr = 10;
ai = find(abs(zgaba)<thr & zgaba > 0);
n = 5;

M = sum(tempChar.innov_counts(motor, :))';

clf;
x = (15:80)';
subplot(121);
scatter( ages(ai), gaba(ai) );

fo = fit(ages(ai), gaba(ai)', 'poly1'); hold on;
plot(x, x.*fo.p1 + fo.p2, 'k');
% 

xlabel('Age');
ylabel('Glx Concentration');

subplot(122);
x = (0:80)';
scatter( M((ai)), gaba(ai) );

fo = fit(M(ai), gaba(ai)', 'poly1'); hold on;
plot(x, x.*fo.p1 + fo.p2, 'k');

xlabel('Innovation Counts (Motor)');
ylabel('Glx Concentration');
% [r,p] = corr( ages(ai), gaba(ai)' )
% plot(zgaba)
% close all

mediation(ages(ai), M(ai), gaba(ai)', 'plots0', 'boottop', 'verbose', 'robust0', 'bootsamples', 5000, 'names', {'Age', 'Innov', 'GABA'});

% [r,p] = corr(ages(ai), gaba(ai)')
rmpath(genpath('/home/luna.kuleuven.be/u0101486/workspace/matlab/toolbox/CanlabCore-master'));
rmpath(genpath('/home/luna.kuleuven.be/u0101486/workspace/matlab/toolbox/MediationToolbox-master'));