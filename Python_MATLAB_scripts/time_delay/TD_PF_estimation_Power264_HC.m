%% Main lag script used to compute time delay matrix and lag projection
% Version 2:
%   In this version, I modified the computation of probability flow
%   projection map. The positive and negative value in the map will be
%   normalized to 1, respectively.
%
% This script was orignally obtained from Github (https://github.com/ryraut/lag-code)
% and customized by Liang Qunjun.
%
% I added the computation of probability flow based on TD and Peak
% correlation. 
% 
% Before running, please make sure the listed scripts are placed in the
% same folder:
%   1. lagged_cov.m
%   2. parabolic_interp.m
%   3. create_blocks.m 
%   4. probability_flow_esti.m 
%   5. probability_flow_norm.m
%
%  Qunjun Liang 2023/01/16

clear; 
clc ;
%% Setup
% Set parameters
num_nodes      = 214;    % number of time series, also, number of ROIs`
wkDir          = '/home/lqj/MDD_patient/NIfTI_convert/derivatives/timeseries_for_TD_Power264/';
outdir         = 'time_lag_estimation';    % set directory for saving out images here
mkdir([wkDir outdir])

lag_lim        = 4;    % lag limit (in seconds)
lags           = -3:3;    % range of TR shifts; max(lags) = round(lag_lim/tr + 1)

% Specify data parameters
sbjNamePath    = '/home/lqj/MDD_patient/R_project/tde_gradient_coupling_in_-mdd/inputs/';
subject_list   = importdata([sbjNamePath 'HC_match_data.tsv']); %!!!!!!!!!!!!!!!!!!!!
subjects       = subject_list.textdata(2:end,1);

% The Yeo 7 network identification for each ROIs
Yeo_7_path     = 'Power264_net_identity.csv';

% how much subject you want to use for estimation, assign 1:numel(subj1ects) for alll subjects
useSubject     = 1:numel(subjects); % 1:91 for MDD and 92:182 for HC
outFile_prefix = 'HC'; % add string to the output csv to identify the analysis
tr             = 2;
motion_thresh  = .2;    % important: must match motion criteria used during preproc
% parallel_core  = 7; % use parallel 
min_block_durn = (max(lags)+1)*tr;   % min. block duration (in seconds)

%% Loop over subjects
% initialize group matrices for running sum
grp_lags = single(nan(num_nodes));    % peak lags
grp_ZL = grp_lags;      % zero-lag correlation, A.K.A., functional connectivity
grp_peak = grp_lags;    % peak correlation
grp_prob_flow = grp_lags; % probability flow 

% collect the NaNs for normalizing the result and compute the mean value
grp_lags_nans = single(zeros(num_nodes));
grp_ZL_nans = grp_lags_nans;
grp_peak_nans = grp_lags_nans;
grp_prob_flow_nans = grp_lags_nans;


%p = parpool(parallel_core);
for s = useSubject
    tic
    subj = subjects{s};
    disp(['Processing ' subj]);
    
    % initialize subject matrices
    subj_lags = single(nan(num_nodes)); % peak lags
    subj_ZL = subj_lags;   % zero-lag correlation
    subj_peak = subj_lags; % peak correlation (correlation at optimal lag)
    
    
    BOLD = load([wkDir subj '_timeseries_power264.csv']); % read in time series matrix
    %BOLD = BOLD.data;
    %good = importdata(''); % read in spatial mask if desired
    good = true(1,num_nodes);
    
    % read in temporal mask/motion time series (e.g., FD or DVARS)
    format = dlmread([wkDir subj '_FD.csv']) <= motion_thresh;
    
    % ignore pre-steady-state frames, commoent out!
    format(1:2) = false; % ignore first X frames
   
    FORMAT = create_blocks(format,min_block_durn,tr);
    
    %% Do the lagged correlation/covariance computation of TD matrices
    Cov = zeros([sum(good) sum(good) numel(lags)]);
    nblocks = numel(FORMAT);
    nframes = 0;
    
    % De-mean time series
    run_mean = nanmean(BOLD(format,:),1);
    BOLD = bsxfun(@minus,BOLD,run_mean);
    
    % Loop over blocks of contiguous frames
    for j = 1:numel(FORMAT)
        nframes = nframes + numel(FORMAT{j});
        FHCR = false(1,numel(format));
        FHCR(FORMAT{j}) = true;
        Cov = Cov + lagged_cov(BOLD(FHCR,good),BOLD(FHCR,good),max(lags));
    end
    
    % Normalize pairwise cross-covariance functions based on entire run
    for k = 1:numel(lags)
        Cov(:,:,k) = Cov(:,:,k)/(nframes - abs(lags(k))*nblocks);
    end
    
    % Parabolic interpolation to get peak lag/correlation
    [pl,pc] = parabolic_interp(Cov,tr);
    pl(abs(pl) > lag_lim) = nan; % Exclude long lags (generally occur when CCF is flat)
    
    % Get zero-lag correlation
    temp = Cov(:,:,lags==0);  % zero-lag correlation
    d = zeros(size(temp));
    d(logical(eye(length(temp)))) = sqrt(diag(temp));
    temp = d^(-1)*temp/d;
    temp = atanh(temp); % Fisher z transform
    temp(isnan(pl)) = nan;
    
    % obtain subject-level estimations and export to CSV
    subj_lags(good,good) = pl;
    subj_ZL(good,good) = temp;
    subj_peak(good,good) = pc;
    
    csvwrite([wkDir outdir '/' subj '_timeDelay.csv'], subj_lags);
    csvwrite([wkDir outdir '/' subj '_peakCorrelation.csv'], subj_peak);
    csvwrite([wkDir outdir '/' subj '_zeroLag.csv'], subj_ZL);
    
    %% Weight subject-level lag projection map
    % Unweighted lag projection
    subj_lags_proj_unweighted = nanmean(subj_lags);

    % Weighted lag projection (inversely weight lags by correlation magnitude
    % to reduce sampling error)
    lag_weights = tan((pi/2)*(1-abs(subj_ZL))).^(-2);    % weighted by 1/f^2(r); f(r) = tan[(pi/2)(1-|r|)]
    lag_weights(logical(eye(size(lag_weights)))) = 0;   % zero-out diagonal weights
    lag_weights(isnan(subj_lags)) = nan;
    subj_lags_mean_wghtd = subj_lags.*lag_weights;
    subj_lags_proj = nansum(subj_lags_mean_wghtd)./nansum(lag_weights);
    
    csvwrite([wkDir outdir '/' subj '_projection_map_weighted.csv'], subj_lags_proj);
    csvwrite([wkDir outdir '/' subj '_projection_map_unweighted.csv'], subj_lags_proj_unweighted);
    
    %% Estimate probability flow
    subj_prob_flow = probability_flow_esti(subj_lags, subj_peak);
    subj_prob_flow_projMap = nanmean(subj_prob_flow);
    csvwrite([wkDir outdir '/' subj '_probaFlow_projectionMap.csv'], subj_prob_flow_projMap);
    
    % normalize the value in probability flow projection map
    subj_prob_flow_projMap = probability_flow_norm(subj_prob_flow_projMap);
    
    csvwrite([wkDir outdir '/' subj '_probaFlow.csv'], subj_prob_flow);
    csvwrite([wkDir outdir '/' subj '_probaFlow_projectionMap_normalization.csv'], subj_prob_flow_projMap);
    %% Add to group running sum
    grp_lags = cat(3,grp_lags,subj_lags); % 'cat' adds to the the third dimention of the matrix
    grp_lags = nansum(grp_lags,3); % 'nansum' computes the sum ignoring the NaNs
    
    grp_ZL = cat(3,grp_ZL,subj_ZL);
    grp_ZL = nansum(grp_ZL,3);
    
    grp_peak = cat(3,grp_peak,subj_peak);
    grp_peak = nansum(grp_peak,3);
    
    grp_prob_flow = cat(3,grp_prob_flow,subj_prob_flow);
    grp_prob_flow = nansum(grp_prob_flow,3);
    
    % running sum of nans
    grp_lags_nans      = grp_lags_nans + isnan(subj_lags);
    grp_ZL_nans        = grp_ZL_nans + isnan(subj_ZL);
    grp_peak_nans      = grp_peak_nans + isnan(subj_peak);
    grp_prob_flow_nans = grp_prob_flow_nans + isnan(subj_prob_flow);
    
    toc
    
end
%delete(p)
%%  Compute group averages
grp_lags_mean = grp_lags ./ (numel(subjects) - grp_lags_nans);
grp_peak_mean = grp_peak ./ (numel(subjects) - grp_peak_nans);
grp_prob_flow_mean = grp_prob_flow ./ (numel(subjects) - grp_prob_flow_nans);
grp_ZL_mean = grp_ZL ./ (numel(subjects) - grp_ZL_nans);
grp_ZL_mean = tanh(grp_ZL_mean); % un-fisher z transform

%% Sort group matrices

assns = importdata(Yeo_7_path); % import ROI network assignments for sorting TD matrix

% Sort by matrices by lag
[M,sorted_inds1] = sort(nanmean(grp_lags_mean));
assns_sort = assns(sorted_inds1);

grp_lags_temp = grp_lags_mean(sorted_inds1,sorted_inds1);
grp_peak_temp = grp_peak_mean(sorted_inds1,sorted_inds1);
grp_ZL_temp = grp_ZL_mean(sorted_inds1,sorted_inds1);
grp_prob_flow_tmp = grp_prob_flow_mean(sorted_inds1,sorted_inds1);

% Sort by network
[N,sorted_inds2] = sort(assns_sort);
sorted_inds2 = sorted_inds2(find(N,1):end);

grp_lags_mat = grp_lags_temp(sorted_inds2,sorted_inds2);
grp_peak_mat = grp_peak_temp(sorted_inds2,sorted_inds2);
grp_ZL_mat = grp_ZL_temp(sorted_inds2,sorted_inds2);
grp_prob_flow_mat = grp_prob_flow_tmp(sorted_inds2,sorted_inds2);

figure;imagesc(grp_lags_mat,[-.6,.6]);colormap(jet);colorbar;title('group time delay matrix');
figure;imagesc(grp_peak_mat);colormap(jet);colorbar;title('group peak time delay matrix');
figure;imagesc(grp_ZL_mat);colormap(jet);colorbar;title('group zero-lag matrix');
figure;imagesc(grp_prob_flow_mat,[-.6,.6]);colormap(jet);colorbar;title('group probability flow matrix');

% save the plots
saveas(1, [wkDir outdir '/' outFile_prefix '_group_lags_plot.png'])
saveas(2, [wkDir outdir '/' outFile_prefix '_group_peak_plot.png'])
saveas(3, [wkDir outdir '/' outFile_prefix '_group_zeroLag_plot.png'])
saveas(4, [wkDir outdir '/' outFile_prefix '_group_probabilityFLow_plot.png'])

%% Make group-level lag and probability flow projection maps
% Unweighted lag projection
grp_lags_proj_unweighted = nanmean(grp_lags_mean);
grp_prob_flow_proj = nanmean(grp_prob_flow_mean);
csvwrite([wkDir outdir '/' outFile_prefix '_group_probability_flow_projection_map.csv'], grp_prob_flow_proj);

grp_prob_flow_proj = probability_flow_norm(grp_prob_flow_proj); % normalization

% Weighted lag projection (inversely weight lags by correlation magnitude
% to reduce sampling error)
lag_weights = tan((pi/2)*(1-abs(grp_ZL_mean))).^(-2);    % weighted by 1/f^2(r); f(r) = tan[(pi/2)(1-|r|)]
lag_weights(logical(eye(size(lag_weights)))) = 0;   % zero-out diagonal weights
lag_weights(isnan(grp_lags_mean)) = nan;
grp_lags_mean_wghtd = grp_lags_mean.*lag_weights;
grp_lags_proj = nansum(grp_lags_mean_wghtd)./nansum(lag_weights);

%% Export the results 

disp('Export results to CSV.');

% export group-mean results
csvwrite([wkDir outdir '/' outFile_prefix '_group_mean_lags.csv'], grp_lags_mean);
csvwrite([wkDir outdir '/' outFile_prefix '_group_mean_peak.csv'], grp_peak_mean);
csvwrite([wkDir outdir '/' outFile_prefix '_group_mean_zeroLag.csv'], grp_ZL_mean);
csvwrite([wkDir outdir '/' outFile_prefix '_group_probability_flow.csv'], grp_prob_flow_mean);

% export group-level projection maps
csvwrite([wkDir outdir '/' outFile_prefix '_group_projection_map_unweighted.csv'], grp_lags_proj_unweighted);
csvwrite([wkDir outdir '/' outFile_prefix '_group_projection_map_weighted.csv'], grp_lags_proj);
csvwrite([wkDir outdir '/' outFile_prefix '_group_probability_flow_projection_map_normalization.csv'], grp_prob_flow_proj);

disp('Pipeline finished without error!');
