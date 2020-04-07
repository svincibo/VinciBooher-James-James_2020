clear all; close all; clc;

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
netnamesshort = [{'1'}, {'2'}, {'3'}, {'4'}, {'5'}, {'6'}, {'7'}, {'8'}, {'9'}, {'10'}, {'11'}, {'12'}, {'13'}];
NofNets = length(unique(cicon));
NSub = 20;
NRand = 10000; % Set number of randomizations requested. Eventually, make NRand = 1000;
sz = 2;
q = 0.12; % 0.12 is highest q value that gives no significance at Session 1 and was, therefore, selected as the threshold for significance.

% Only one session at a time, for ease.
for iSes = 1:3
    
    %% FIRST: Calculate the Real Distribution of the Interaction Contrast.
    
    % Initialize.
    temp = NaN(NofNets, NofNets);
    diff_real = NaN(NofNets, NofNets, NSub);
    mu_real = zeros(NofNets, NofNets);
    s_real = zeros(NofNets, NofNets);
    
    % Calculate the mean beta of the interaction contrast for each subject.
    for iSub = 1:NSub
        
        % Find runs that belong with this subject.
        if iSes == 3
            % For Session 3, sub105, sub111, and sub116 have only 2 runs.
            run_idx = find(EXP_run_sub_indx_3 == iSub);
        else
            % For Sessions 1 and 2, all subjects have 3 runs.
            run_idx = find(EXP_run_sub_indx_1and2 == iSub);
        end
        
        % Average over the runs that belong to this subject.
        temp = squeeze(mean(CorMatNETppi(:, :, run_idx, iSes), 3));
        
        % Set diagonal to zero. Note: It's effectively zero to begin with (e.g., 1e-15).
        temp(1:(NofNets + 1):end) = 0;
        
        % Collect.
        diff_real(:, :, iSub) = temp;
        
        % Reinitialize temp matrices.
        clear temp
        temp = NaN(NofNets, NofNets);
        
    end % end iSub
    
    % Calculate the mean beta for this session, averaging over subjects.
    mu_real = nanmean(diff_real, 3);
    
    % Calculate the standard deviation for this session.
    s_real = nanstd(diff_real, [], 3);
    
    %% SECOND: Create Null Distribution of the Betas for Interaction Contrast.  
        
    % Load the MAT file containing the matrices for each possible permutation.
    load([inoutDir 'samplingDistribution_interaction_session' num2str(iSes) '_SWMRILOC2partition.mat']);
    
    % Initialize.
    diff_null = NaN(NofNets, NofNets, NSub);
    mu_null = zeros(NofNets, NofNets, NRand);
    s_null = zeros(NofNets, NofNets, NRand);
    
    % SECOND: Randomly select from all possible permutations to create the null distribution for a comparison within a session.
    for iRand = 1:NRand
        
        disp(['******************************--- RANDOM SELECTION ' num2str(iRand) ' ---******************************'])
        
        % For each subject (so that subject is treated as a random effect),
        for iSub = 1:NSub
            
            % Randomly select a subject.
            is = randi(NSub);
            
            % Randomly select one of the possible permutations for this subject.
            ir = randi(size(all_perm, 1));
            
            % Collect.
            diff_null(:, :, iSub) = all_diff_null(:, :, is, ir);
            
        end % end iSub
        
        % Calculate the mean interaction beta for the interaction under the null hypothesis, averaging over subjects.
        mu_null(:, :, iRand) = nanmean(diff_null, 3);
        
        % Calculate the standard deviation for the interaction under the null hypothesis.
        s_null(:, :, iRand) = nanstd(diff_null, [], 3);
        
        % Reinitialize diff_null.
        clear diff_null
        diff_null = NaN(NofNets, NofNets, NSub);
        
    end % end iRand
    
    % Clear loaded variables that were specific to this session.
    clear all_diff_null all_perm
    
    %% THIRD: Test for significance by comparing the real distribution to the null distribution.
    
    % ---Mean under the null.
    m = nanmean(mu_null, 3);
    % ---Standard deviation under the null.
    s = nanstd(mu_null, [], 3);
    % ---Convert to z-score ((actual - mean)/standard deviation).
    z = (mu_real - m)./s;
    % ---Set diagonals to zero.
    z(1:(NofNets + 1):end) = 0;
    % ---Correct for multiple comparisons, FDR.
    p = 1 - normcdf(abs(z), 0, 1);
    padj = fcn_linear_step_up(p(triu(ones(NofNets)) > 0), q);
    mask = p <= padj;
    if size(find(p <= padj), 1) == NofNets*NofNets
        mask = zeros([NofNets NofNets]);
    end
    
    %% FOURTH: Plot.
    
    NRsn = NofNets;
    jdx = (1:NofNets)';
    
    f = figure(iSes);
    temp = z.*mask; 
    temp(find(temp(:) == 0)) = NaN;
    imagesc(temp(jdx,jdx), 'AlphaData', ~isnan(temp(jdx, jdx)), [-5 5]);
    f.InvertHardcopy = 'off';
    set(gca,...
        'ytick',1:NRsn,...
        'yticklabel',netnamesshort(jdx),...
        'xtick',1:NRsn,...
        'xticklabel',netnamesshort(jdx));
    cb = colorbar;
    set(get(cb,'ylabel'),'string','z-score');
    xlabel('communities');
    ylabel('communities');
    title(['Session ' num2str(iSes)])
    set(gca,'color', .75*[1 1 1]);
    
    print([inoutDir 'significanceTesting_interaction_session' num2str(iSes) '_sig_q=' num2str(q) '_SWMRILOC2partition.eps'], '-depsc2');
    print([inoutDir 'significanceTesting_interaction_session' num2str(iSes) '_sig_q=' num2str(q) '_SWMRILOC2partition.png'], '-dpng');
    
    %% FIFTH: Save.
    save([inoutDir 'significanceTesting_interaction_session' num2str(iSes) '_q=' num2str(q) '_SWMRILOC2partition.mat'], ...
    'z', 'mask');   
    
end % end iSes