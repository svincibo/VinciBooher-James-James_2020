function [CorMatNETppi_out, CorMatNET_out, CorMatROIppi_out, CorMatROIz_out] = function_loadPowerParts_ppiStdBeta_EXP_interaction(iSub, iSes, r_perm, cicon, TAL_XYZ_Power, sz, subDir)

% Other variables specific to SWMRI.
subIDs = {'101', '102', '104', '105', '106', '107', '108', '109', '110', '111', '112', '113', '114', '115', '116', '117', '118', '119', '120', '121'};
NofSubs = length(subIDs); NSes = 3;
category = {'B', 'A', 'B', 'A', 'A', 'A', 'A', 'D', 'D', 'D', 'C', 'C', 'B', 'D', 'C', 'D', 'C', 'D', 'C', 'B'};
runType = 'EXP';
NofVols = 452; % number of volumes per run, 368 for LOC and 452 for EXP
atlas = 'SWMRILOC2';

% Get network indices according to atlas.
[py, pyi] = sort(cicon);
NofNets = length(unique(cicon));

NofRois = size(TAL_XYZ_Power, 1);

% Get subject ID.
subID = subIDs{iSub};

disp(['SUBJECT:' subID])

% Subject 106 only has 2 EXP runs for Session 3.
if strcmp(subID, '106') && iSes == 3
    
    NofRuns = 2;
    
else
    
    NofRuns = 3;
    
end % end if subID 106

CorMatROI = NaN(NofRois, NofRois, NofRuns); % correlation matrix initialization
CorMatNET = NaN(NofNets, NofNets, NofRuns); % correlation matrix initialization
CorMatROIppi = NaN(NofRois, NofRois, NofRuns); % correlation matrix initialization
CorMatNETppi = NaN(NofNets, NofNets, NofRuns); % correlation matrix initialization

realIndex = 0;
for runNum = 1:NofRuns
        
    % these two runs have problems with their signal -- refer to JPGs of rTC for these two runs
    if strcmp(subID, '111') && iSes == 3 && runNum == 2
        
    elseif strcmp(subID, '116') && iSes == 3 && runNum == 1
        
    else 
    
    disp(['RUN: ' num2str(runNum)])
    
    % Update run counter. Note that this is a within-session run counter.
    realIndex = realIndex + 1;
    
    vtcFileName = sprintf('%s_session%d_run%d_%s_SCCAI2_3DMCTS_SD3DSS6.00mm_LTR_THP2c_TAL.vtc', subID, iSes, runNum+1, runType); %because run1 is LOC
    [~, V] = readVTC([subDir vtcFileName]); % readVTC takes vtc file and reads it into header info and data variables
    % size(V) = time, 58, 40, 46; TAL dims (V) : [time, Y, Z, X] [time, CR, AX, SG]; Y=coronoal, Z=axial, X=saggital
    
    % Mean time courses for ROIs.
    rTC = zeros(NofVols, NofRois);
    for rIndex = 1:NofRois
        
        % Convert Power Talairach ROI coordinates to VTC coordinates.
        vtc_coords = tal2vtc(TAL_XYZ_Power(rIndex, :));
        
        % Make a "sphere", or cube in this case.
        cr = (-sz:+sz) + vtc_coords(1); ax = (-sz:+sz) + vtc_coords(2); sg = (-sz:+sz) + vtc_coords(3);
        
        % Make sure we don't go outside of where the data are supposed to exist.
        cr = cr(cr>=1); cr = cr(cr<=58); ax = ax(ax>=1); ax = ax(ax<=40); sg = sg(sg>=1); sg = sg(sg<=46);
        
        % Collapse 4D matrix for this ROI at size(V(:, cr, ax, sg))=452x3x3x3 cube, into 452 x 27 sheet. No operations necesary, just moving cubes of a block into a 1-block thick sheet.
        roiTC = reshape(V(:, cr, ax, sg), NofVols, length(cr)*length(ax)*length(sg));
        
        % Add mean BOLD activation for each ROI to 452 by 264 matrix
        rTC(:, rIndex) = mean(roiTC, 2);
        
    end %for rIndex
    
    clear V;
    
    % Mean time courses for networks.
    nTC = zeros(NofVols, NofNets);
    for nIndex = 1:NofNets
        
        % For some network, take mean across columns/rois for that network
        nTC(:, nIndex) = nanmean(rTC(:, cicon == nIndex), 2);
        
    end %for nIndex
    
    % Get the SDM file name for this subject, session, and run.
    sdmFileName = sprintf('%s_session%d_run%d_%s_SCCAI2_3DMCTS_SD3DSS6.00mm_LTR_THP2c_TAL_3DMC_SR.sdm', subID, iSes, runNum+1, runType);
    
    % Read in the SDM.
    [h, S] = readSDM([subDir filesep sdmFileName]);
    
    % Replace NaN with 0.
    S(isnan(S)) = 0;
    
    % Select confound predictors from the full design matrix.
    M = S(:, h.firstConfoundPred:end-1);
    
    %-----------------
    % Partial Correlation
    %-----------------
    
    % ROI
    CorMatROI(:, :, realIndex) = partialcorr(rTC, M);
    
    % NET
    CorMatNET(:, :, realIndex) = partialcorr(nTC, M);
    
    %-----------------
    % PPI
    %-----------------
    
    % Create Contrast Vector for interaction. Note: PRT includes all sessions.
    if iSes == 1
        
        %         X = S(:, 2) - S(:, 3) - (S(:, 4) - S(:, 5));
                X = S(:, r_perm(1)) - S(:, r_perm(2)) - (S(:, r_perm(3)) - S(:, r_perm(4)));
%         X = 1.5*S(:, r_perm(1)) + 0.5*S(:, r_perm(2)) - 0.5*(S(:, r_perm(3)) - 1.5*S(:, r_perm(4)));
%         X = S(:, r_perm(1)) - S(:, r_perm(2));
        
    elseif iSes == 2
        
        %         X = S(:, 8) - S(:, 9) - (S(:, 10) - S(:, 11));
                X = S(:, r_perm(1)+6) - S(:, r_perm(2)+6) - (S(:, r_perm(3)+6) - S(:, r_perm(4)+6));
%         X = 1.5*S(:, r_perm(1)+6) + 0.5*S(:, r_perm(2)+6) - 0.5*(S(:, r_perm(3)+6) - 1.5*S(:, r_perm(4)+6));
%          X = S(:, r_perm(1)+6) - S(:, r_perm(2)+6);
         
    elseif iSes == 3
        
        %         X = S(:, 14) - S(:, 15) - (S(:, 16) - S(:, 17));
                X = S(:, r_perm(1)+12) - S(:, r_perm(2)+12) - (S(:, r_perm(3)+12) - S(:, r_perm(4)+12));
%         X = 1.5*S(:, r_perm(1)+12) + 0.5*S(:, r_perm(2)+12) - 0.5*(S(:, r_perm(3)+12) - 1.5*S(:, r_perm(4)+12));
%          X = S(:, r_perm(1)+12) - S(:, r_perm(2)+12);
         
    end
    
    % ROI
    [BI, ~] = ppiStdBetaS(rTC, X, M, 'amax');
    
    CorMatROIppi(:, :, realIndex) = BI;
    
    % NET
    [BIn, ~] = ppiStdBetaS(nTC, X, M, 'amax');
    
    CorMatNETppi(:, :, realIndex) = BIn;
    
    end % if 111 or 116 exception
    
end %for runNum

% Fisher r-to-z transform.
CorMatROIz =.5.*log((1+CorMatROI)./(1-CorMatROI));

%% Save.
% Average over all runs for this subject and this session.
CorMatNETppi_out = squeeze(mean(CorMatNETppi, 3));
CorMatNET_out = squeeze(mean(CorMatNET, 3));
CorMatROIppi_out = squeeze(mean(CorMatROIppi, 3));
CorMatROIz_out = squeeze(mean(CorMatROIz, 3));
