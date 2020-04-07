function [BI, BA] = ppiStdBetaS(Y, x, M, method)
% BI = ppiStdBeta(Y,x,M,method)
%
% Y:    t x r matrix, multi-roi timeseries (r rois)
% x:    t x 1 vector, psychological contrast
% M:    t x v matrix, nuisance variables (v variables)
%
% BI:   r x r matrix, ppi interaction term beta
% BA:   r x r matrix, physiological term beta
%
% method: method for dealing with matrix asymmetry
%   'asym': raw, asymmetric
%   'max','min': use maximum or minimum value of upper and lower triangles
%   'amax','amin': use maximum or minimum of the absolute value, but
%                   maintain sign
%
% Betas are standardized: a value of 1 represents a 1-sd change in Y
%

% If you transform all IVs and DVs are transformed to z-scores, then you
% get standardized coefs. So, let's do that here:
YS = ( Y - (ones(size(Y, 1), 1) * mean(Y, 1)) ) ./ (ones(size(Y, 1), 1) * std(Y, 0, 1));

% preassign
BI = NaN(size(YS,2),size(YS,2));
BA = NaN(size(YS,2),size(YS,2));

% For each ROI
for r = 1:size(YS, 2)
    
    % Seect the standardized coefficients for just one ROI.
    as = YS(:, r);
    
    % Only do for ROIs that actually have data.
    if all(isnan(as))
        
    else
        
        % ppi model design matrix: [psychological, physiological, interaction, nuisance variables, constant]
        W = [x, as, x.*as, M, ones(length(as), 1) ];
        
        % fit model
        B = pinv(W'* W) * W' * YS;
        
        BI(r, :) = B(3, :);
        BA(r, :) = B(2, :);
        
    end
    
end

if strcmpi(method,'asym')
    
    return;
    
else
    
    MI = cat(3,triu(BI),tril(BI)');
    MA = cat(3,triu(BA),tril(BA)');
    
    if strcmpi(method,'max')
        
        Mt = max(MI,[],3);
        BI = Mt + Mt';
        Mt = max(MA,[],3);
        BA = Mt + Mt';
        
    elseif strcmpi(method,'min')
        
        Mt = min(MI,[],3);
        BI = Mt + Mt';
        Mt = min(MA,[],3);
        BA = Mt + Mt';
        
    elseif strcmpi(method,'amax')
        
        [Mt,Mi] = max(abs(MI),[],3);
        Mu = sign(MI(:,:,1));
        Mu(Mi==2) = 0;
        Ml = sign(MI(:,:,2));
        Ml(Mi==1) = 0;
        Ms = Mu + Ml;
        
        BI = (Mt + Mt') .* (Ms + Ms');
        
        [Mt,Mi] = max(abs(MA),[],3);
        Mu = sign(MA(:,:,1));
        Mu(Mi==2) = 0;
        Ml = sign(MA(:,:,2));
        Ml(Mi==1) = 0;
        Ms = Mu + Ml;
        
        BA = (Mt + Mt') .* (Ms + Ms');
        
    elseif strcmpi(method,'amin')
        
        [Mt,Mi] = min(abs(MI),[],3);
        Mu = sign(MI(:,:,1));
        Mu(Mi==2) = 0;
        Ml = sign(MI(:,:,2));
        Ml(Mi==1) = 0;
        Ms = Mu + Ml;
        
        BI = (Mt + Mt') .* (Ms + Ms');
        
        [Mt,Mi] = min(abs(MA),[],3);
        Mu = sign(MA(:,:,1));
        Mu(Mi==2) = 0;
        Ml = sign(MA(:,:,2));
        Ml(Mi==1) = 0;
        Ms = Mu + Ml;
        
        BA = (Mt + Mt') .* (Ms + Ms');
        
    else
        
        error('method %s not recognized\n', method);
        
    end
    
    return;
    
end


