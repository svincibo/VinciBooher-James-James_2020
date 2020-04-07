clear all; close all; clc;

% User input... only one session at a time, for ease.
for iSes = 1:3
    
    tic
    
    rootDir = '/N/dc2/projects/lifebid/development/';
    inoutDir = [rootDir 'SWMRI_fc/'];
    supportDir = [rootDir 'SWMRI_fc/supportFiles/'];
    subDir = [rootDir 'SWMRI_subjects/'];
    
    addpath(genpath(supportDir));
    addpath(subDir);
    addpath(inoutDir);
    
    % Matrix of MNI coordinates, in the order of ROIs.
    load([supportDir 'MNI_XYZ_Power.mat']);
    TAL_XYZ_Power = mni2tal(MNI_XYZ_Power);
    
    % Load in the MAT file containing the SWMRI LOC ROI community assignments.
    load([supportDir 'SWMRI_LOC_allTRs_power13_mcdgsd_rTC_pc_amax_groupConsensusPartitionBasedOnAveragedCorrelationMatrix_gamma0.8421_Powerprior.mat']);
    
    % Load in the MAT file containing the PPI Correlation matrices from loadPowerParts...
    load([inoutDir 'SWMRI_EXP_sppi_interaction(DI>DnI)>(WD>WS)_SWMRILOC2_r=2_mcdgsd_rTC_ppiStdBeta_amax.mat']);
    
    % Load in the MAT file that codes which runs belong to each participant.
    load([supportDir 'SWMRI_EXP_run_sub_idx.mat']);
    
    % Set general Ns and parameters.
    % netnamesshort = [{'1'}, {'2'}, {'3'}, {'4'}, {'5'}, {'6'}, {'7'}, {'8'}, {'9'}, {'10'}, {'11'}, {'12'}];
    NofNets = length(unique(cicon));
    NSub = 20;
    %     NRand = 10000; % Set number of randomizations requested. Eventually, make NRand = 1000;
    sz = 2;
    
    %     %% Calculate the Mean Beta of the Interaction Contrast for Each Session.
    %
    %     % Initialize.
    %     temp = NaN(NofNets, NofNets);
    %     diff_real = NaN(NofNets, NofNets, NSub);
    %     mu_real = zeros(NofNets, NofNets);
    %     s_real = zeros(NofNets, NofNets);
    %
    %     % Calculate the mean beta of the interaction contrast for each subject.
    %     for iSub = 1:NSub
    %
    %         % Find runs that belong with this subject.
    %         if iSes == 3
    %             % For Session 3, sub105, sub111, and sub116 have only 2 runs.
    %             run_idx = find(EXP_run_sub_indx_3 == iSub);
    %         else
    %             % For Sessions 1 and 2, all subjects have 3 runs.
    %             run_idx = find(EXP_run_sub_indx_1and2 == iSub);
    %         end
    %
    %         % Average over the runs that belong to this subject.
    %         temp = squeeze(mean(CorMatNETppi(:, :, run_idx, iSes), 3));
    %
    %         % Set diagonal to zero. Note: It's effectively zero to begin with (e.g., 1e-15).
    %         temp(1:(NofNets + 1):end) = 0;
    %
    %         % Collect.
    %         diff_real(:, :, iSub) = temp;
    %
    %         % Reinitialize temp matrices.
    %         clear temp
    %         temp = NaN(NofNets, NofNets);
    %
    %     end % end iSub
    %
    %     % Calculate the mean beta for this session, averaging over subjects.
    %     mu_real = nanmean(diff_real, 3);
    %
    %     % Calculate the standard deviation for this session.
    %     s_real = nanstd(diff_real, [], 3);
    
    %% Create Null Distribution of the Betas for Interaction Contrast for Each Session.
    
    % Initialize.
    CorMatNETppi_null = NaN(NofNets, NofNets);
    diff_null = NaN(NofNets, NofNets, NSub);
    %     mu_null = zeros(NofNets, NofNets, NRand);
    %     s_null = zeros(NofNets, NofNets, NRand);
    
    % Create all possible permutations of the interaction
    all_perm = perms([2, 3, 4, 5]); % for interaction
    
    for p = 1:size(all_perm, 1)
        
        disp(['******************************--- PERMUTATION ' num2str(p) ' ---******************************'])
        
        % Create permuted contrast vector for this permutation.
        contrast = all_perm(p, :);
        
        for iSub = 1:NSub
            
            % Get permuted interaction beta matrix.
            [CorMatNETppi_null, ~] = function_loadPowerParts_ppiStdBeta_EXP_interaction(iSub, iSes, contrast, cicon, TAL_XYZ_Power, sz, subDir);
            
            % Set diagonal to zero.
            CorMatNETppi_null(1:(NofNets + 1):end) = 0;
            
            % Assign to diff_null.
            diff_null(:, :, iSub) = CorMatNETppi_null;
            
            % Reinitialize temp matrices.
            clear CorMatNETppi_null
            CorMatNETppi_null = NaN(NofNets, NofNets);
            
        end % end sub
        
        % Save this permutation for this subject to call later when making the null distribution.
        all_diff_null(:, :, :, p) = diff_null;
        
    end % end p
    
    %     clear diff_null
    %     diff_null = NaN(NofNets, NofNets, NSub);
    %
    %     % SECOND: Randomly select from all possible permutations to create the null distribution for a comparison within a session.
    %     for iRand = 1:NRand
    %
    %         disp(['******************************--- RANDOM SELECTION ' num2str(iRand) ' ---******************************'])
    %
    %         % For each subject (so that subject is treated as a random effect),
    %         for iSub = 1:NSub
    %
    %             % Randomly select a subject.
    %             is = randi(NSub);
    %
    %             % Randomly select one of the possible permutations for this subject.
    %             ir = randi(size(all_perm, 1));
    %
    %             % Collect.
    %             diff_null(:, :, iSub) = all_diff_null(:, :, is, ir);
    %
    %         end % end iSub
    %
    %         % Calculate the mean interaction beta for the interaction under the null hypothesis, averaging over subjects.
    %         mu_null(:, :, iRand) = nanmean(diff_null, 3);
    %
    %         % Calculate the standard deviation for the interaction under the null hypothesis.
    %         s_null(:, :, iRand) = nanstd(diff_null, [], 3);
    %
    %         % Reinitialize diff_null.
    %         clear diff_null
    %         diff_null = NaN(NofNets, NofNets, NSub);
    %
    %     end % end iRand
    
    %% Save.
    
    save([inoutDir 'samplingDistribution_interaction_session' num2str(iSes) '_SWMRILOC2partition.mat'], ...
        'all_perm', 'all_diff_null');
    
    toc
    
    clear all_diff_null diff_real all_perm mu_real s_real mu_null s_null
    
end % end iSes