clear all
close all
clc

rootDir = '/N/dc2/projects/lifebid/development/';
inoutDir = [rootDir 'SWMRI_fc/'];
supportDir = [rootDir 'SWMRI_fc/supportFiles/'];
subDir = [rootDir 'SWMRI_subjects/'];

addpath(genpath(supportDir));
addpath(subDir);
addpath(genpath(inoutDir));

% Load in the MAT file containing the SWMRI LOC ROI community assignments.
load([supportDir 'SWMRI_LOC_allTRs_power13_mcdgsd_rTC_pc_amax_groupConsensusPartitionBasedOnAveragedCorrelationMatrix_gamma0.8421_Powerprior.mat']);

% Load in the MAT file containing the PPI Correlation matrices from loadPowerParts...
load([inoutDir 'SWMRI_EXP_sppi_interaction(DI>DnI)>(WD>WS)_SWMRILOC2_r=2_mcdgsd_rTC_ppiStdBeta_amax.mat']);

% Load in the MAT file that codes which runs belong to each participant.
load([supportDir 'SWMRI_EXP_run_sub_idx.mat']);

% Load in the sampling distribution
load([inoutDir 'samplingDistribution_interaction_session_SWMRILOC2partition.mat']);

% Set general Ns and parameters.
netnamesshort = [{'1'}, {'2'}, {'3'}, {'4'}, {'5'}, {'6'}, {'7'}, {'8'}, {'9'}, {'10'}, {'11'}, {'12'}, {'13'}];
NofNets = length(unique(cicon));
NSub = 20;
NSes = 3;
NRand = 10000; % Set number of randomizations requested.
sz = 2;
q=.10;

%% FIRST: Get Real Distribution for Between-session Differences in Interaction.

% Initialize: Session 2 - Session 1
diff_2v1_real = NaN(NofNets, NofNets, NSub);
mu_2v1_real = zeros(NofNets, NofNets, NSes);
s_2v1_real = zeros(NofNets, NofNets, NSes);

% Initialize: Session 3 - Session 2
diff_3v2_real = NaN(NofNets, NofNets, NSub);
mu_3v2_real = zeros(NofNets, NofNets, NSes);
s_3v2_real = zeros(NofNets, NofNets, NSes);

% Initialize: Session 3 - Session 1
diff_3v1_real = NaN(NofNets, NofNets, NSub);
mu_3v1_real = zeros(NofNets, NofNets, NSes);
s_3v1_real = zeros(NofNets, NofNets, NSes);

temp2v1 = NaN(NofNets, NofNets);
temp3v2 = NaN(NofNets, NofNets);
temp3v1 = NaN(NofNets, NofNets);

% Calculate the mean beta of the interaction contrast for each subject.
for iSub = 1:NSub
    
    % Find runs that belong with this subject. For Session 3, sub105, sub111, and sub116, have only 2 runs.
    run_idx1 = find(EXP_run_sub_indx_1and2 == iSub);
    run_idx2 = find(EXP_run_sub_indx_1and2 == iSub);
    run_idx3 = find(EXP_run_sub_indx_3 == iSub);
    
    % ---Session 2 - Session 1
    
    % Average over the runs that belong to this subject.
    temp2v1 = squeeze(mean(CorMatNETppi(:, :, run_idx2, 2), 3) - mean(CorMatNETppi(:, :, run_idx1, 1), 3));
    
    % Set diagonal to zero. Note: It's effectively zero to begin with (e.g., 1e-15).
    temp2v1(1:(NofNets + 1):end) = 0;
    
    % Assign to diff.
    diff_2v1_real(:, :, iSub) = temp2v1;
    
    % ---Session 3 - Session 2
    
    % Average over the runs that belong to this subject.
    temp3v2 = squeeze(mean(CorMatNETppi(:, :, run_idx3, 3), 3) - mean(CorMatNETppi(:, :, run_idx2, 2), 3));
    
    % Set diagonal to zero. Note: It's effectively zero to begin with (e.g., 1e-15).
    temp3v2(1:(NofNets + 1):end) = 0;
    
    % Assign to diff.
    diff_3v2_real(:, :, iSub) = temp3v2;
    
    % ---Session 3 - Session 1
    
    % Average over the runs that belong to this subject.
    temp3v1 = squeeze(mean(CorMatNETppi(:, :, run_idx3, 3), 3) - mean(CorMatNETppi(:, :, run_idx2, 1), 3));
    
    % Set diagonal to zero. Note: It's effectively zero to begin with (e.g., 1e-15).
    temp3v1(1:(NofNets + 1):end) = 0;
    
    % Assign to diff.
    diff_3v1_real(:, :, iSub) = temp3v1;
    
    % Reinitialize temp2v1 and temp3v2,
    clear temp2v1 temp3v2 temp3v1
    temp2v1 = NaN(NofNets, NofNets);
    temp3v2 = NaN(NofNets, NofNets);
    temp3v1 = NaN(NofNets, NofNets);

end % end iSub

% Calculate the mean beta for each session comparison, averaging over subjects.
mu_2v1_real = nanmean(diff_2v1_real, 3);
mu_3v2_real = nanmean(diff_3v2_real, 3);
mu_3v1_real = nanmean(diff_3v1_real, 3);

% Calculate the standard deviation for this session.
s_2v1_real = nanstd(diff_2v1_real, [], 3);
s_3v2_real = nanstd(diff_3v2_real, [], 3);
s_3v1_real = nanstd(diff_3v1_real, [], 3);

%% Create Null Distribution of the Betas for Interaction Contrast for Each Session.

% Initialize.
temp = NaN(NofNets, NofNets);
diff_null = NaN(NofNets, NofNets, NSub);
mu_null = zeros(NofNets, NofNets, NRand);
s_null = zeros(NofNets, NofNets, NRand);

%% SECOND: Get Null Distribution for Between-session Differences in Interaction.
diff_null = NaN(NofNets, NofNets, NSub);

% Randomly select from all possible permutations to create the null distribution for a comparison between two sessions.
for iRand = 1:NRand
    
    % 2v1: For each subject (so that subject is treated as a random effect),
    for iSub = 1:NSub
        
        % Randomly select a subject.
        is = randi(NSub);
        
        % Randomly select one of the possible permutations.
        ir = randi(size(all_diff_null_2v1, 4));
        
        % Collect.
        diff_null(:, :, iSub) = all_diff_null_2v1(:, :, is, ir);
        
    end % end iSub
    
    % Calculate the mean interaction beta for this session under the null hypothesis, averaging over subjects.
    mu_null_2v1(:, :, iRand) = nanmean(diff_null, 3);
    
    % Calculate the standard deviation for this session under the null hypothesis.
    s_null_2v1(:, :, iRand) = nanstd(diff_null, [], 3);
    
    % Reinitialize diff_null.
    clear diff_null is ir iSub
    diff_null = NaN(NofNets, NofNets, NSub);
    
    % 3v2: For each subject (so that subject is treated as a random effect),
    for iSub = 1:NSub
        
        % Randomly select a subject.
        is = randi(NSub);
        
        % Randomly select one of the possible permutations.
        ir = randi(size(all_diff_null_3v2, 4));
        
        % Collect.
        diff_null(:, :, iSub) = all_diff_null_3v2(:, :, is, ir);
        
    end % end iSub
    
    % Calculate the mean interaction beta for this session under the null hypothesis, averaging over subjects.
    mu_null_3v2(:, :, iRand) = nanmean(diff_null, 3);
    
    % Calculate the standard deviation for this session under the null hypothesis.
    s_null_3v2(:, :, iRand) = nanstd(diff_null, [], 3);
    
    % Reinitialize diff_null.
    clear diff_null is ir iSub
    diff_null = NaN(NofNets, NofNets, NSub);
    
    % For each subject (so that subject is treated as a random effect),
    for iSub = 1:NSub
        
        % Randomly select a subject.
        is = randi(NSub);
        
        % Randomly select one of the possible permutations.
        ir = randi(size(all_diff_null_3v1, 4));
        
        % Collect.
        diff_null(:, :, iSub) = all_diff_null_3v1(:, :, is, ir);
        
    end % end iSub
    
    % Calculate the mean interaction beta for this session under the null hypothesis, averaging over subjects.
    mu_null_3v1(:, :, iRand) = nanmean(diff_null, 3);
    
    % Calculate the standard deviation for this session under the null hypothesis.
    s_null_3v1(:, :, iRand) = nanstd(diff_null, [], 3);
    
    % Reinitialize diff_null.
    clear diff_null is ir iSub
    diff_null = NaN(NofNets, NofNets, NSub);
    
end %end iRand

%% THIRD: Test for Significance by Comparing Real Distribution to Null Distribution.

% Session 2 - Session 1
mu_2v1_null = mu_null_2v1;
% ---Mean under the null.
m2v1 = nanmean(mu_2v1_null, 3);
% ---Standard deviation under the null.
s2v1 = nanstd(mu_2v1_null, [], 3);
% ---Convert to z-score ((actual - mean)/standard deviation).
z2v1 = (mu_2v1_real - m2v1)./s2v1;
% ---Set diagonals to zero.
z2v1(1:(NofNets + 1):end) = 0;
% ---Correct for multiple comparisons, FDR.
p2v1 = 1 - normcdf(abs(z2v1), 0, 1);
padj2v1 = fcn_linear_step_up(p2v1(triu(ones(NofNets)) > 0), q);
mask2v1 = p2v1 <= padj2v1;
if size(find(p2v1 <= padj2v1), 1) == NofNets*NofNets
    mask2v1 = zeros([NofNets NofNets]);
end

% Session 3 - Session 2
mu_3v2_null = mu_null_3v2;
% ---Mean under the null.
m3v2 = nanmean(mu_3v2_null, 3);
% ---Standard deviation under the null.
s3v2 = nanstd(mu_3v2_null, [], 3);
% ---Convert to z-score ((actual - mean)/standard deviation).
z3v2 = (mu_3v2_real - m3v2)./s3v2;
% ---Set diagonals to zero.
z3v2(1:(NofNets + 1):end) = 0;
% ---Correct for multiple comparisons, FDR.
p3v2 = 1 - normcdf(abs(z3v2), 0, 1);
padj3v2 = fcn_linear_step_up(p3v2(triu(ones(NofNets)) > 0), q);
mask3v2 = p3v2 <= padj3v2;
if size(find(p3v2 <= padj3v2), 1) == NofNets*NofNets
    mask3v2 = zeros([NofNets NofNets]);
end

% Session 3 - Session 1
mu_3v1_null = mu_null_3v1;   % keep the same null?
% ---Mean under the null.
m3v1 = nanmean(mu_3v1_null, 3);
% ---Standard deviation under the null.
s3v1 = nanstd(mu_3v1_null, [], 3);
% ---Convert to z-score ((actual - mean)/standard deviation).
z3v1 = (mu_3v1_real - m3v2)./s3v1;
% ---Set diagonals to zero.
z3v1(1:(NofNets + 1):end) = 0;
% ---Correct for multiple comparisons, FDR.
p3v1 = 1 - normcdf(abs(z3v1), 0, 1);
padj3v1 = fcn_linear_step_up(p3v1(triu(ones(NofNets)) > 0), q);
mask3v1 = p3v1 <= padj3v1;
if size(find(p3v1 <= padj3v1), 1) == NofNets*NofNets
    mask3v1 = zeros([NofNets NofNets]);
end

%% FIFTH: Plot.

NRsn = NofNets;
jdx = (1:NofNets)';

f1 = figure(1);
temp = -z2v1.*mask2v1; % negative because original contrasts were done 2>1 but for display we want 1>2
temp(find(temp(:) == 0)) = NaN;
imagesc(temp(jdx,jdx), 'AlphaData', ~isnan(temp(jdx, jdx)), [-5 5]);
f1.InvertHardcopy = 'off';
set(gca,...
    'ytick',1:NRsn,...
    'yticklabel',netnamesshort(jdx),...
    'xtick',1:NRsn,...
    'xticklabel',netnamesshort(jdx));
cb = colorbar;
set(get(cb,'ylabel'),'string','z-score');
xlabel('communities');
ylabel('communities');
title('Session 1 - Session 2')
set(gca,'color', .75*[1 1 1]);

print([inoutDir 'significanceTesting_interaction_session1v2_sig_q=' num2str(q) '_SWMRILOC2partition.eps'], '-depsc2');
print([inoutDir 'significanceTesting_interaction_session1v2_sig_q=' num2str(q) '_SWMRILOC2partition.png'], '-dpng');

f2 = figure(2);
temp = -z3v2.*mask3v2;
temp(find(temp(:) == 0)) = NaN;
imagesc(temp(jdx,jdx), 'AlphaData', ~isnan(temp(jdx, jdx)), [-5 5]);
f2.InvertHardcopy = 'off';
set(gca,...
    'ytick',1:NRsn,...
    'yticklabel',netnamesshort(jdx),...
    'xtick',1:NRsn,...
    'xticklabel',netnamesshort(jdx));
cb = colorbar;
set(get(cb,'ylabel'),'string','z-score');
xlabel('communities');
ylabel('communities');
title('Session 2 - Session 3')
set(gca,'color', .75*[1 1 1]);

print([inoutDir 'significanceTesting_interaction_session2v3_sig_q=' num2str(q) '_SWMRILOC2partition.eps'], '-depsc2');
print([inoutDir 'significanceTesting_interaction_session2v3_sig_q=' num2str(q) '_SWMRILOC2partition.png'], '-dpng');

f3 = figure(3);
temp = -z3v1.*mask3v1;
temp(find(temp(:) == 0)) = NaN;
imagesc(temp(jdx,jdx), 'AlphaData', ~isnan(temp(jdx, jdx)), [-5 5]);
f3.InvertHardcopy = 'off';
set(gca,...
    'ytick',1:NRsn,...
    'yticklabel',netnamesshort(jdx),...
    'xtick',1:NRsn,...
    'xticklabel',netnamesshort(jdx));
cb = colorbar;
set(get(cb,'ylabel'),'string','z-score');
xlabel('communities');
ylabel('communities');
title('Session 1 - Session 3')
set(gca,'color', .75*[1 1 1]);

print([inoutDir 'significanceTesting_interaction_session1v3_sig_q=' num2str(q) '_SWMRILOC2partition.eps'], '-depsc2');
print([inoutDir 'significanceTesting_interaction_session1v3_sig_q=' num2str(q) '_SWMRILOC2partition.png'], '-dpng');

%% Save.
save([inoutDir 'significanceTesting_interaction_session_q=' num2str(q) '_SWMRILOC2partition.mat'], ...
    'mask2v1', 'mask3v2', 'mask3v1', 'z2v1', 'z3v2', 'z3v1');
