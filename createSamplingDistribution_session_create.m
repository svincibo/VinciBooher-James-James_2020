clear all
close all
clc

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
NofNets = length(unique(cicon));
NSub = 20;
NSes = 3;
NRand = 10000; % Set number of randomizations requested. Eventually, make NRand = 10000;
sz = 2;

tic

% FIRST: create all possible permutations of session 
all_perm_2v1 = perms(1:2);
all_perm_3v2 = perms(2:3);
all_perm_3v1 = perms([1 3]);

for p = 1:size(all_perm_2v1, 1)
    
    disp(['******************************--- PERMUTATION ' num2str(p) ' ---******************************'])
    
    % Contrast of interest.
    r = [2 3 4 5]; % for DI, DnI, WD, WS

    for iSub = 1:NSub
        
        % Get permuted interaction beta matrix: 1v2.
        [CorMatNETppi_null1, ~] = function_loadPowerParts_ppiStdBeta_EXP_interaction(iSub, all_perm_2v1(p, 1), r, cicon, TAL_XYZ_Power, sz, subDir);
        [CorMatNETppi_null2, ~] = function_loadPowerParts_ppiStdBeta_EXP_interaction(iSub, all_perm_2v1(p, 2), r, cicon, TAL_XYZ_Power, sz, subDir);
        
        % Within iSub loop storage of permuted interaction matrix for this subject and this session.
        temp = CorMatNETppi_null1 - CorMatNETppi_null2 ;
        
        % Set diagonal to zero.
        temp(1:(NofNets + 1):end) = 0;
        
        % Assign to diff.
        diff_null_2v1(:, :, iSub) = temp;
        
        % Reinitialize temp.
        clear temp CorMatNETppi_null1 CorMatNETppi_null2
        temp = NaN(NofNets, NofNets);
        
        % Get permuted interaction beta matrix: 2v3.
        [CorMatNETppi_null1, ~] = function_loadPowerParts_ppiStdBeta_EXP_interaction(iSub, all_perm_3v2(p, 1), r, cicon, TAL_XYZ_Power, sz, subDir);
        [CorMatNETppi_null2, ~] = function_loadPowerParts_ppiStdBeta_EXP_interaction(iSub, all_perm_3v2(p, 2), r, cicon, TAL_XYZ_Power, sz, subDir);
        
        % Within iSub loop storage of permuted interaction matrix for this subject and this session.
        temp = CorMatNETppi_null1 - CorMatNETppi_null2 ;
        
        % Set diagonal to zero.
        temp(1:(NofNets + 1):end) = 0;
        
        % Assign to diff.
        diff_null_3v2(:, :, iSub) = temp;
        
        % Reinitialize temp.
        clear temp CorMatNETppi_null1 CorMatNETppi_null2
        temp = NaN(NofNets, NofNets);
        
        % Get permuted interaction beta matrix: 1v3.
        [CorMatNETppi_null1, ~] = function_loadPowerParts_ppiStdBeta_EXP_interaction(iSub, all_perm_3v1(p, 1), r, cicon, TAL_XYZ_Power, sz, subDir);
        [CorMatNETppi_null2, ~] = function_loadPowerParts_ppiStdBeta_EXP_interaction(iSub, all_perm_3v1(p, 2), r, cicon, TAL_XYZ_Power, sz, subDir);
        
        % Within iSub loop storage of permuted interaction matrix for this subject and this session.
        temp = CorMatNETppi_null1 - CorMatNETppi_null2 ;
        
        % Set diagonal to zero.
        temp(1:(NofNets + 1):end) = 0;
        
        % Assign to diff.
        diff_null_3v1(:, :, iSub) = temp;
        
        % Reinitialize temp.
        clear temp CorMatNETppi_null1 CorMatNETppi_null2
        temp = NaN(NofNets, NofNets);
        
    end % end sub
    
    % Save this permutation for this subject to call later when making the null distribution.
    all_diff_null_2v1(:, :, :, p) = diff_null_2v1;
    all_diff_null_3v2(:, :, :, p) = diff_null_3v2;
    all_diff_null_3v1(:, :, :, p) = diff_null_3v1;
    
end % end p

toc

%% Save.

save([inoutDir 'samplingDistribution_interaction_session_SWMRILOC2partition.mat'], ...
    'all_perm_2v1', 'all_perm_3v2', 'all_perm_3v1', ...
    'all_diff_null_2v1', 'all_diff_null_3v2', 'all_diff_null_3v1');