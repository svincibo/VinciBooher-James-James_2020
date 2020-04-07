function [BIG, BA] = ppiStdBetaG(Y, X, M, method)
% BI = ppiStdBetaG(Y,X,M,method)
%
% Y:    t x r matrix, multi-roi timeseries (r rois)
% X:    t x p vector, gppi psychological contrasts (p contrasts)
% M:    t x v matrix, nuisance variables (v variables)
%
% BIG:   r x r x p matrix, ppi interaction term betas
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
 
% If all IVs and DVs are transformed to z-scores, then you get standardized coefs. So, let's do that here:
YS = ( Y - (ones(size(Y, 1), 1) * mean(Y, 1)) ) ./ (ones(size(Y, 1), 1) * std(Y, 0, 1));
 
% preassign
BIG = NaN(size(YS, 2), size(YS, 2), size(X, 2));
BA = NaN(size(YS, 2), size(YS, 2));
 
% For each ROI
for r = 1:size(YS, 2)

    % Select the standardized coefficients for just one ROI.
    as = YS(:, r);
    
    if all(isnan(as))
        
    else
    
    % ppi model design matrix: [psychological, physiological, interaction, nuisance variables, constant]
    W = [X, as, X .* (as * ones(1, size(X, 2))), M, ones(length(as), 1) ];
    
    % fit model
    B = pinv(W'* W) * W' * YS;
    
    % resize so that for this r there is a 
    BIG(r, :, :) = reshape(B(5 + (1:size(X, 2)), :)', 1, size(B, 2), size(X, 2));
    
    BA(r, :) = B(2, :);
    
    end
    
end
 
if strcmpi(method,'asym')
    
    return;
    
else
    
    MA = cat(3,triu(BA),tril(BA)');
    
    if strcmpi(method,'max')
        Mt = max(MA,[],3);
        BA = Mt + Mt';
 
    elseif strcmpi(method,'min')
        Mt = min(MA,[],3);
        BA = Mt + Mt';
 
    elseif strcmpi(method,'amax')
        [Mt,Mi] = max(abs(MA),[],3);
        Mu = sign(MA(:,:,1));
        Mu(Mi==2) = 0;
        Ml = sign(MA(:,:,2));
        Ml(Mi==1) = 0;
        Ms = Mu + Ml;
 
        BA = (Mt + Mt') .* (Ms + Ms');
        
    elseif strcmpi(method,'amin')
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

    for p = 1:size(X,2)
        
        MI = cat(3,triu(BIG(:,:,p)),tril(BIG(:,:,p))');
 
        if strcmpi(method,'max')
            Mt = max(MI,[],3);
            BIG(:,:,p) = Mt + Mt';

        elseif strcmpi(method,'min')
            Mt = min(MI,[],3);
            BIG(:,:,p) = Mt + Mt';

        elseif strcmpi(method,'amax')
            [Mt,Mi] = max(abs(MI),[],3);
            Mu = sign(MI(:,:,1));
            Mu(Mi==2) = 0;
            Ml = sign(MI(:,:,2));
            Ml(Mi==1) = 0;
            Ms = Mu + Ml;

            BIG(:,:,p) = (Mt + Mt') .* (Ms + Ms');

        elseif strcmpi(method,'amin')
            [Mt,Mi] = min(abs(MI),[],3);
            Mu = sign(MI(:,:,1));
            Mu(Mi==2) = 0;
            Ml = sign(MI(:,:,2));
            Ml(Mi==1) = 0;
            Ms = Mu + Ml;

            BIG(:,:,p) = (Mt + Mt') .* (Ms + Ms');

        end
        
    end
    
    return;
 
end


