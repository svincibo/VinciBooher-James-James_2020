% passive music listening experiment

% traveling wave

% Creating THE MATRIX

clear all; clc;

tic

rootDir = '/Volumes/swmri/';
fcDir = [rootDir 'bin/fcAnalysis/'];
outDir = [fcDir 'communityDetection/'];
supportDir = [rootDir 'bin/fcAnalysis/supportFiles/'];
bvDir = rootDir;

addpath(genpath(fcDir));

% Tells us to which of the 13 power networks does each ROI belong.
% Assumption: each ROI belongs to only one nework.
load([supportDir 'powerYeoIndex2.mat']);
[py, pyi] = sort(powerIndex13);

% Matrix of MNI coordinates, in the order of ROIs.
load([supportDir 'MNI_XYZ_Power.mat']);
TAL_XYZ_Power = mni2tal(MNI_XYZ_Power);

%subIDs = {'101', '102'}; % Don't forget 106!
subIDs = {'101', '102', '104', '105', '106', '107', '108', '109', '110', '111', '112', '113', '114', '115', '116', '117', '118', '119', '120', '121'};
NofSubs = length(subIDs);

%category = {'B', 'A'};
category = {'B', 'A', 'B', 'A', 'A', 'A', 'A', 'D', 'D', 'D', 'C', 'C', 'B', 'D', 'C', 'D', 'C', 'D', 'C', 'B'};

runType = 'LOC';
NofRuns = 1;

NofVols = 368; % number of volumes per run, 368 for LOC and 452 for EXP
NofRois = size(TAL_XYZ_Power, 1); % number of ROIs
NofNets = length(unique(powerIndex13)); % unique reproduces a vector, but takes out repeats

realIDs = cell(NofSubs, 1);
realSession = NaN(NofSubs*3, 1); % because 3 sessions per subject

CorMatROI = NaN(NofRois, NofRois, NofSubs, NofRuns); % correlation matrix initialization
CorMatNET = NaN(NofNets, NofNets, NofSubs, NofRuns); % correlation matrix initialization
CorMatROIppi = NaN(NofRois, NofRois, NofSubs, NofRuns); % correlation matrix initialization
CorMatNETppi = NaN(NofNets, NofNets, NofSubs, NofRuns); % correlation matrix initialization
CorMatROInonppi = NaN(NofRois, NofRois, NofSubs, NofRuns);
CorMatNETnonppi = NaN(NofNets, NofNets, NofSubs, NofRuns);
RTC = NaN(NofRois, NofRois ,NofSubs, NofRuns);
NTC = NaN(NofNets ,NofNets, NofSubs, NofRuns);
FD = NaN(NofVols, NofSubs, NofRuns);

realIndex = 0;

for subIndex = 1:NofSubs
    
    subID = subIDs{subIndex};
    disp(subID)
    
    for sesNum = 1:3
        
        subDir = [bvDir 'SWMRI_subjects' filesep subID '_' category{subIndex} '_session' num2str(sesNum) '/']
        
        fileExist = NaN(2, NofRuns); % initializing matrix to see if VTCs and SDMs exist
        
        for runNum = 1:NofRuns
            
            vtcFileName = sprintf('%s_session%d_run%d_%s_SCCAI2_3DMCTS_SD3DSS6.00mm_LTR_THP2c_TAL.vtc', subID, sesNum, runNum, runType); % because run 1 is LOC
            fileExist(1, runNum) = exist([subDir vtcFileName], 'file');
            
            sdmFileName = sprintf('%s_session%d_run%d_%s_SCCAI2_3DMCTS_SD3DSS6.00mm_LTR_THP2c_TAL_3DMC_SR.sdm', subID, sesNum, runNum, runType); % because run 1 is LOC
            fileExist(2, runNum) = exist([subDir sdmFileName], 'file');
            
        end %for runNum
        
        clear runNum
        
        %if all are non-zero
        %         if all(fileExist)
        
        realIndex = realIndex + 1; % only update realIndex if a VTC and SDM exists for this subject and session
        realSession(realIndex, :) = sesNum;
        realIDs{realIndex} = [subID '_' num2str(sesNum)];
        
        for runNum = 1:NofRuns
            
            vtcFileName = sprintf('%s_session%d_run%d_%s_SCCAI2_3DMCTS_SD3DSS6.00mm_LTR_THP2c_TAL.vtc', subID, sesNum, runNum, runType); %because run1 is LOC
            [header, V] = readVTC([subDir vtcFileName]); % readVTC takes vtc file and reads it into header info and data variables
            % size(V) = time, 58, 40, 46; TAL dims (V) : [time, Y, Z, X] [time, CR, AX, SG]; Y=coronoal, Z=axial, X=saggital
            
            rTC = zeros(NofVols, NofRois); % mean time courses for ROIs
            nTC = zeros(NofVols, NofNets); % mean time courses for networks
            
            for rIndex = 1:NofRois
                
                vtc_coords = tal2vtc(TAL_XYZ_Power(rIndex,:));
                % makes a "sphere" or cube in this case
                cr = (-1:+1) + vtc_coords(1); ax = (-1:+1) + vtc_coords(2); sg = (-1:+1) + vtc_coords(3);
                % make sure we don't go outside of where the data are supposed to exist
                cr = cr(cr>=1); cr = cr(cr<=58); ax = ax(ax>=1); ax = ax(ax<=40); sg = sg(sg>=1); sg = sg(sg<=46);
                
                roiTC = reshape(V(:, cr, ax, sg), NofVols, length(cr)*length(ax)*length(sg)); % Collapse 4D matrix for this ROI at size(V(:, cr, ax, sg))=452x3x3x3 cube, into 452 x 27 sheet. No operations necesary, just moving cubes of a block into a 1-block thick sheet.
                rTC(:, rIndex) = mean(roiTC, 2); % adding mean BOLD activation for each ROI to 452 by 264 matrix
                
            end %for rIndex
            
            clear V;
            
            for nIndex = 1:NofNets
                
                nTC(:, nIndex) = mean(rTC(:, powerIndex13 == nIndex), 2); % For some network, take mean across columns/rois for that network
                
            end %for nIndex
            
            sdmFileName = sprintf('%s_session%d_run%d_%s_SCCAI2_3DMCTS_SD3DSS6.00mm_LTR_THP2c_TAL_3DMC_SR.sdm', subID, sesNum, runNum, runType);
            [h, S] = readSDM([subDir filesep sdmFileName]);
            S(isnan(S)) = 0;
                      
            %GS(:, realIndex, runNum) = S(:, 15); % I don't think I should use GS because this is task-based, not steady-state or resting-state..
            % fd = sqrt(sum(S(:,3:5).^2, 2));% this doesn't seem to be used in the correlation
            
            %-----------------
            % ROI-level PPI
            %-----------------
            
            % Select task predictors from the full design matrix.
            M2 = S(:, h.firstConfoundPred:end-1);
            
            % Create a contrast vector. This is different for each session because of the way the PRT files were made.
            if sesNum == 1
                
%                                 contrast = 'Letters-Shapes';
%                                 X = S(:, 2) + S(:, 3) - (S(:, 4) + S(:, 5));
                
                contrast = 'Letters-Fixation';
                X = S(:, 2) + S(:, 3) - (2*S(:, 6));
                
            elseif sesNum == 2
                
%                                 contrast = 'Letters-Shapes';
%                                 X = S(:, 7) + S(:, 8) - (S(:, 9) + S(:, 10));
                
                contrast = 'Letters-Fixation';
                X = S(:, 7) + S(:, 8) - (2*S(:, 11));
                
            elseif sesNum == 3
%                 
%                                 contrast = 'Letters-Shapes';
%                                 X = S(:, 12) + S(:, 13) - (S(:, 14) + S(:, 15));
                
                contrast = 'Letters-Fixation';
                X = S(:, 12) + S(:, 13) - 2*(S(:, 16));
                
            end
            
            % Correlate the seed timecourses for all ROIs with the design matrix.
            CorMatROI(:, :, realIndex, runNum) = partialcorr(rTC, M2);
            
            % Create interaction predictor for all ROIs for a particular contrast.
            Rx = rTC .* X;
            
            %Correlate all interaction predictor with the other task predictors.
            CorMatROInonppi(:, :, realIndex, runNum) = partialcorr(Rx, M2);
            
            % Specify the physiological predictor.
            Y = rTC;
            
            [BI, BA] = ppiStdBeta(Y, X, M2,'amax');
            
            CorMatROIppi(:, :, realIndex, runNum) = BI;
            
            %                 for rIndex = 1:NofRois
            %
            %                     seed = rTC(:,rIndex); % a in ppiStdBeta
            %                     M = S(:,3:16); % nuisance variables
            %                     [prsq,rsq] = ppiStdBeta(seed,rTC,X,M);
            %                     CorMatROIppi(:,rIndex,realIndex,runNum) = prsq;
            %
            %                 end
            
            %-----------------
            % Network-level PPI
            %-----------------
            
            Nx = nTC.* X;
            CorMatNETnonppi(:, :, realIndex, runNum) = partialcorr(Nx, M2);
            
            Yn = nTC;
            [BI,BA] = ppiStdBeta(Yn, X, M2, 'amax');
            CorMatNETppi(:, :, realIndex, runNum) = BI;
            CorMatNET(:, :, realIndex, runNum) = partialcorr(nTC, M2);
            %                 for nIndex = 1:NofNets
            %
            %                     seedN = nTC(:,nIndex);
            %                     M = S(:,3:16);
            %                     [prsq,rsq] = ppiStdBeta(seedN,nTC,X,M);
            %                     CorMatNETppi(:,nIndex,realIndex,runNum) = prsq;
            %
            %                 end
            
            %                 FD(:, realIndex, runNum) = fd;
            %                 RTC(:, :, realIndex, runNum) = rTC;
            %                 NTC(:, :, realIndex, runNum) = nTC;
            X(:, realIndex, runNum) = X;
            
        end %for runNum
                
    end %for sesNum
    
    %     end %if fileExist
    
end %for subID

realIDs = realIDs(1:realIndex);
realSession = realSession(1:realIndex, :);
CorMatNET = CorMatNET(:, :, 1:realIndex, :);
CorMatROI = CorMatROI(:, :, 1:realIndex, :);
CorMatROIppi = CorMatROIppi(:, :, 1:realIndex, :);
CorMatROInonppi = CorMatROInonppi(:, :, 1:realIndex, :);
CorMatNETppi = CorMatNETppi(:, :, 1:realIndex, :);
CorMatNETnonppi = CorMatNETnonppi(:, :, 1:realIndex, :);
% GS = GS(:, 1:realIndex, :);
% FD = FD(:, 1:realIndex, :);
% RTC = RTC(:, :, 1:realIndex, :);
% NTC = NTC(:, :, 1:realIndex, :);
% X = X(:, 1:realIndex, :);

CorMatROIz =.5.*log((1+CorMatROI)./(1-CorMatROI));

save([outDir 'SWMRI_LOC_' contrast '_Power13_mcdgsd_rTC_ppiStdBeta_amax.mat'],'realIDs', 'realSession', 'CorMatNET','CorMatROI', 'CorMatROIz', 'CorMatROIppi','CorMatNETppi','CorMatROInonppi','CorMatNETnonppi','FD');
%save([fcDir 'SWMRI_Power13_mcdgsd_rTC_ppiStdBeta_amax.mat'],'realIDs','realGrps','CorMatNET','CorMatROI', 'CorMatROIz', 'CorMatROIppi','CorMatNETppi','CorMatROInonppi','CorMatNETnonppi','GS','FD','RTC','NTC','X');

toc
