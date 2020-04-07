%%% REMEMBER: I was given this code for the purpose of analyzing the SWMRI
%%% data set on 9/11/2018 by Anastasia Nikoulina who, herself received it
%%% from Rick Betzel.
% 0.7895    0.8947    1.0000
%  0.4211    0.8421
% 
%     14

%     0.3684    0.7895    0.8947    1.0000
  

% -----Selecting gamma = 0.7895 and using Power communities as prior to genlouvaion().

clear all; close all; clc

rootDir = '/Volumes/SWMRI/';
fcDir = [rootDir 'bin/fcAnalysis/'];
inoutDir = [fcDir '/communityDetection/'];
supportDir = [rootDir 'bin/fcAnalysis/supportFiles/'];

addpath(genpath(fcDir));

atlas = 'power13'; %power13, yeo7, yeo17
runType = 'LOC';
contrast = 'Letters-Fixation';
subIDs = {'101', '102', '104', '105', '106', '107', '108', '109', '110', '111', '112', '113', '114', '115', '116', '117', '118', '119', '120', '121'};

gammas = [0.8421];
%[0.7895 0.8947 1.0000];

% Load in the MAT file containing the partial correlation matrices for each run of each participant.
load([inoutDir 'SWMRI_' runType '_' contrast '_' atlas '_mcdgsd_rTC_ppiStdBeta_amax.mat']);% created by loadPowerParts 12.m, matrix of the connectivity strenghts

% Load in the MAT file containint the Power ROI community assignments.
load([supportDir 'powerYeoIndex2.mat']);

% Load in the MAT file that codes which runs belong to each participant.
load([supportDir 'SWMRI_LOC_run_sub_idx.mat']);

% Get number of nodes.
N = size(CorMatROIz, 1);

%% AVERAGE ACROSS ALL LOC RUNS FOR ALL SUBS.
% Use CorMatROI because this is the partial correlation matrix using the 264 Power parcellation.
% Use the Fisher r-to-z transformed matrix.

% Get average correlation matrix for this subject (average across their runs).
rho = mean(CorMatROIz, 3);

% Set diagonal to zero, thereby removing Inf from diagonal.
rho(1:(N + 1):end) = 0;

%% PERFORM MODULARITY MAXIMIZATION NREPS NUMBER OF TIMES FOR THIS GAMMA FOR THIS SUB.

for m = 1:size(gammas, 2)
    
    % SET PARAMETERS.

    % Set gamma based on results from communityDetection_static_SWMRI_selectGamma.m.
    gamma = gammas(m);
    
    % Modularity matrix: the matrix of observed functional connections minus what would be expected by chance (in this case chance = 0).
    b = (rho - gamma).* ~eye(N);
    
    % Symmetrize (just in case there's numerical issues).
    b = (b + b')/2;
    
    % Number of times to repeat the algorithm.
    nreps = 1000;
    
    % Initialize ci.
    ci = zeros(N, nreps);
    
    % Loop through reps.
    for irep = 1:nreps
        
        % Generalized Louvain algorithm.
        [ci(:, irep), q(:, irep)] = genlouvain(b, [], [], [], [], powerIndex13);
%         [ci(:, irep), q(:, irep)] = genlouvain(b);
        %         ci(:, irep) = genlouvain(b);
        
    end
    
    %% CONSENSUS PARTITION.
    
    % Get agreement matrix: number of times pairs of nodes were assigned to the same cluster.
    agree = agreement(ci);
    
    % Get consensus clusters.
    [cicon, qcon] = fcn_consensus_communities(ci, nreps, false);
    
    % Reorder network nodes based on consensus clusters.
    [idx, cisort] = fcn_order_partition(rho, cicon);
    
    % Get boundaries of communities based on the sorting of the consensus clusters.
    [gx, gy] = fcn_plot_blocks(cisort);
    
    % Get number of communities for quick reference.
    nci = max(cicon);
    disp(nci)
    
    %% VISUALIZE COMMUNITIES.
    
    f = figure('units', 'inches', 'position', [2, 2, 8, 4]);
    
    s(1) = subplot(1, 2, 1);
    imagesc(rho(idx, idx));   % plot correlation matrix for consensus partition
    hold(gca, 'on');
    plot(gx, gy, 'k');
    xticks([15, 40, 75, 108, 128, 148, 166, 190, 217, 235, 254, 262])
    xticklabels({'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'})
    yticks([15, 40, 75, 108, 128, 148, 166, 190, 217, 235, 254, 262])
    yticklabels({'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'})
    cb = colorbar;

    s(2) = subplot(1, 2, 2);
    imagesc(agree(idx, idx));        % plot agreement matrix
    hold(gca, 'on');
    plot(gx, gy,'k');
    xticks([15, 40, 75, 108, 128, 148, 166, 190, 217, 235, 254, 262])
    xticklabels({'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'})
    yticks([15, 40, 75, 108, 128, 148, 166, 190, 217, 235, 254, 262])
    yticklabels({'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'})
    cb = colorbar;
    
    title(s(1), {['Gamma = ' num2str(gamma) ', ' num2str(nci) ' Communities, Q = ' num2str(qcon)]; 'Correlation Matrix for Consensus Partition'})
    title(s(2), 'Agreement Matrix Across Iterations')
    
%     print([inoutDir 'agreementMatrix_gamma' num2str(gamma) '_group_noprior.eps'], '-depsc2'); 
%     print([inoutDir 'agreementMatrix_gamma' num2str(gamma) '_group_noprior.png'], '-dpng');
    
    hold(gca, 'off')
    
    %% SAVE COMMUNITY ASSIGNMENTS
    
%     save([supportDir 'SWMRI_' runType '_allTRs_' atlas '_mcdgsd_rTC_pc_amax_groupConsensusPartitionBasedOnAveragedCorrelationMatrix_gamma' num2str(gamma) '_noprior.mat'], 'cicon');

end % for gammas

