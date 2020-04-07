%%% REMEMBER: I was given this code for the purpose of analyzing the SWMRI
%%% data set on 9/11/2018 by Anastasia Nikoulina who, herself received it
%%% from Rick Betzel.

% References for gamma selection algorithm:
% Bassett DS, Porter MA, Wymbs NF, Grafton ST, Carlson JM, Mucha PJ. Robust detection of dynamic  community structure in networks. Chaos. 2013; 23(1):013142. [PubMed: 23556979] [An important  technical paper highlighting many useful guidelines for modularity maximization.]
%
% Traud AL, Kelsic ED, Mucha PJ, Porter MA. Comparing community structure to characteristics in  online collegiate social networks. SIAM Rev. 2011; 53(3):526?43.

clear all; close all; clc

rootDir = '/Volumes/swmri/';
fcDir = [rootDir 'bin/fcAnalysis/'];
inoutDir = [fcDir 'communityDetection/'];
supportDir = [rootDir 'bin/fcAnalysis/supportFiles/'];

addpath(genpath(fcDir));

atlas = 'Power13'; %power13, yeo7, yeo17
runType = 'LOC';
contrast = 'Letters-Fixation'; %  NOTE: CorMatROIz is not created based on this contrast. It is only the partial correlation over the entire run after controlling for nuisance variables. 

% Load in the MAT file containing the partial correlation matrices for each run of each participant.
load([inoutDir 'SWMRI_' runType '_' contrast '_' atlas '_mcdgsd_rTC_ppiStdBeta_amax.mat']);% created by loadPowerParts 12.m, matrix of the connectivity strenghts

% Load in the MAT file containint the Power ROI community assignments.
load([supportDir 'powerYeoIndex2.mat']);

% Load in the MAT file that codes which runs belong to each participant.
load([supportDir 'SWMRI_LOC_run_sub_idx.mat']);

% Use CorMatROI because this is the partial correlation matrix using the 264 Power parcellation. Use the Fisher r-to-z transformed matrix.
rho = mean(CorMatROIz, 3);

% Get number of nodes.
N = size(CorMatROIz, 1);

% Set diagonal to zero, thereby removing Inf from diagonal.
rho(1:(N + 1):end) = 0;

% SET/INITIALIZE PARAMETERS FOR SELECTION OF GAMMA LOOP.

% Number of gamma values to try.
ngam = 20;

% Set the range of gamma values based on the number of gamma values to try.
gammavalues = linspace(0, 1, ngam);

% Set how many times the modularity maximization should be performed for each gamma.
nreps = 1000;

% Initialize the array to hold the community clusters generated.
ci = zeros(N, nreps, ngam);

% Initialize the array to hold the z-scores of the rand coefficients.
z = zeros(nreps*(nreps - 1)/2, ngam);

% LOOP THROUGH RESOLUTION PARAMETERS (GAMMA).
for igam = 1:ngam
    
    % SET PARAMETERS.
    
    % Select current gamma.
    gamma = gammavalues(igam);
    
    % Modularity matrix: the matrix of observed functional connections minus what would be expected by chance (in this case chance = 0).
    b = (rho - gamma).* ~eye(N);
    
    % Symmetrize (just in case there's numerical issues).
    b = (b + b')/2;
    
    % PERFORM MODULARITY MAXIMIZATION NREPS NUMBER OF TIMES FOR THIS GAMMA.
    for irep = 1:nreps
        
%         if irep == 1
        
        % Generalized Louvain algorithm.
%         ci(:, irep, igam) = genlouvain(b, [], [], [], [], powerIndex13);
        
%         else 
            
        ci(:, irep, igam) = genlouvain(b);
        
%         end
        
    end % irep
    
    % COMPARE Z-SCORE OF RAND COEFFICIENT FOR THIS GAMMA TO CHANCE. -- why does this not use cicon instead of ci?
    count = 0;
    for irep = 1:(nreps - 1)
        
        for jrep = (irep + 1):nreps
            
            % Update counter.
            count = count + 1;
            
            % CALCULATE SIMILARITY: z-score of the rand coefficient for this gamma
            z(count, igam) = fcn_zrand(ci(:, irep, igam), ci(:, jrep, igam)); % entering the irep partition and the partition immediately following it (i.e., the jrep partition)
            
        end %jrep
        
    end % irep
    
    disp(igam);
    
end % igam

% We want to focus on local maxima, so let's find the peaks of the median similarity.
pks = findpeaks(median(z));
pks = pks.loc;
f = figure('units', 'inches', 'position', [2,2,4,4]);
ph = plot(gammavalues, median(z), gammavalues(pks), median(z(:, pks)), 'ro');
xlabel('gamma values');
ylabel('z-rand');

saveas(f, [inoutDir 'peak_gammas_fromGroupAveragedConnectivityMatrix_noprior.eps'])
saveas(f, [inoutDir 'peak_gammas_fromGroupAveragedConnectivityMatrix_noprior.png'])

% Display the gamma values at each peak identified.
disp(gammavalues(pks));

% Display the average number of clusters in nrep partitions created with each gamma.
disp(mean(mean(squeeze((max(ci(:, :, pks), [], 1))), 1)));

% Save.
out.peakGammas = gammavalues(pks);
out.avgNumClusters = mean(mean(squeeze((max(ci(:, :, pks), [], 1))), 1));

save([inoutDir 'peak_gammas_fromGroupAveragedConnectivityMatrix_noprior.mat'], 'out');