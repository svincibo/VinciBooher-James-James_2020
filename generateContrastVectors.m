function [X] = generateContrastVectors(sesNum, contrast, S)
%GENERATECONTRASTVECTORS generates several contrast vectors at to be
%evaluated -- only here to keep it clear in the main code

% This is different for each session because of the way the PRT files were made.

if sesNum == 1
    
    if strcmp(contrast, 'interaction(DI>DnI)>(WD>WS)')
        
        X = S(:, 2) - S(:, 3) - (S(:, 4) - S(:, 5));
        
    elseif strcmp(contrast, '(DI+DnI)>(WD+WS)')
        
        X = S(:, 2) + S(:, 3) - (S(:, 4) + S(:, 5));
        
    elseif strcmp(contrast, '(DI+WD)>(DnI+WS)')
        
        X = S(:, 2) + S(:, 4) - (S(:, 3) + S(:, 5));
        
    elseif strcmp(contrast, 'DIDnI')
        
        X = S(:, 2) - S(:, 3);
        
    elseif strcmp(contrast, 'DI>WD')
        
        X = S(:, 2) - S(:, 4);
        
    elseif strcmp(contrast, 'WD>WS')
        
        X = S(:, 4) - S(:, 5);
        
    elseif strcmp(contrast,'DI>UL')
        
        X = S(:, 2) - S(:, 6);
        
    elseif strcmp(contrast, 'DnI>UL')
        
        X = S(:, 3) - S(:, 6);
        
    elseif strcmp(contrast, 'WD>UL')
        
        X = S(:, 4) - S(:, 6);
        
    elseif strcmp(contrast, 'WS>UL')
        
        X = S(:, 5) - S(:, 6);
        
    elseif strcmp(contrast, 'DI')
        
        X = S(:, 2);
        
    elseif strcmp(contrast, 'DnI')
        
        X = S(:, 3);
        
    elseif strcmp(contrast, 'WD')
        
        X = S(:, 4);
        
    elseif strcmp(contrast, 'WS')
        
        X = S(:, 5);
        
    elseif strcmp(contrast, 'UL')
        
        X = S(:, 6);
        
    elseif strcmp(contrast, 'DIDnIWDUL')
        
        X = 1.5*S(:, 2) + 0.5*S(:, 3) - 0.5*S(:, 4) - 1.5*S(:, 6);
        
    elseif strcmp(contrast, 'DIDnIWDWS')
        
        X = 1.5*S(:, 2) + 0.5*S(:, 3) - 0.5*S(:, 4) - 1.5*S(:, 5);
        
    elseif strcmp(contrast, 'DIDnIWDFix')
        
        X = 1.5*S(:, 2) + 0.5*S(:, 3) - 0.5*S(:, 4) - 1.5*S(:, 1);
        
    end
    
    
elseif sesNum == 2
    
    if strcmp(contrast, 'interaction(DI>DnI)>(WD>WS)')
        
        X = S(:, 8) - S(:, 9) - (S(:, 10) - S(:, 11));
        
    elseif strcmp(contrast, '(DI+DnI)>(WD+WS)')
        
        X = S(:, 8) + S(:, 9) - (S(:, 10) + S(:, 11));
        
    elseif strcmp(contrast, '(DI+WD)>(DnI+WS)')
        
        X = S(:, 8) + S(:, 10) - (S(:, 9) + S(:, 11));
        
    elseif strcmp(contrast, 'DIDnI')
        
        X = S(:, 8) - S(:, 9);
        
    elseif strcmp(contrast, 'DI>WD')
        
        X = S(:, 8) - S(:, 10);
        
    elseif strcmp(contrast, 'WD>WS')
        
        X = S(:, 10) - S(:, 11);
        
    elseif strcmp(contrast,'DI>UL')
        
        X = S(:, 8) - S(:, 12);
        
    elseif strcmp(contrast, 'DnI>UL')
        
        X = S(:, 9) - S(:, 12);
        
    elseif strcmp(contrast, 'WD>UL')
        
        X = S(:, 10) - S(:, 12);
        
    elseif strcmp(contrast, 'WS>UL')
        
        X = S(:, 11) - S(:, 12);
        
    elseif strcmp(contrast, 'DI')
        
        X = S(:, 8);
        
    elseif strcmp(contrast, 'DnI')
        
        X = S(:, 9);
        
    elseif strcmp(contrast, 'WD')
        
        X = S(:, 10);
        
    elseif strcmp(contrast, 'WS')
        
        X = S(:, 11);
        
    elseif strcmp(contrast, 'UL')
        
        X = S(:, 12);
        
    elseif strcmp(contrast, 'DIDnIWDUL')
        
        X = 1.5*S(:, 8) + 0.5*S(:, 9) - 0.5*S(:, 10) - 1.5*S(:, 12);
        
    elseif strcmp(contrast, 'DIDnIWDWS')
        
        X = 1.5*S(:, 8) + 0.5*S(:, 9) - 0.5*S(:, 10) - 1.5*S(:, 11);
        
    elseif strcmp(contrast, 'DIDnIWDFix')
        
        X = 1.5*S(:, 8) + 0.5*S(:, 9) - 0.5*S(:, 10) - 1.5*S(:, 7);
        
    end
    
elseif sesNum == 3
    
    if strcmp(contrast, 'interaction(DI>DnI)>(WD>WS)')
        
        X = S(:, 14) - S(:, 15) - (S(:, 16) - S(:, 17));
        
    elseif strcmp(contrast, '(DI+DnI)>(WD+WS)')
        
        X = S(:, 14) + S(:, 15) - (S(:, 16) + S(:, 17));
        
    elseif strcmp(contrast, '(DI+WD)>(DnI+WS)')
        
        X = S(:, 14) + S(:, 16) - (S(:, 15) + S(:, 17));
        
    elseif strcmp(contrast, 'DIDnI')
        
        X = S(:, 14) - S(:, 15);
        
    elseif strcmp(contrast, 'DI>WD')
        
        X = S(:, 14) - S(:, 16);
        
    elseif strcmp(contrast, 'WD>WS')
        
        X = S(:, 16) - S(:, 17);
        
    elseif strcmp(contrast,'DI>UL')
        
        X = S(:, 14) - S(:, 18);
        
    elseif strcmp(contrast, 'DnI>UL')
        
        X = S(:, 15) - S(:, 18);
        
    elseif strcmp(contrast, 'WD>UL')
        
        X = S(:, 16) - S(:, 18);
        
    elseif strcmp(contrast, 'WS>UL')
        
        X = S(:, 17) - S(:, 18);
        
    elseif strcmp(contrast, 'DI')
        
        X = S(:, 14);
        
    elseif strcmp(contrast, 'DnI')
        
        X = S(:, 15);
        
    elseif strcmp(contrast, 'WD')
        
        X = S(:, 16);
        
    elseif strcmp(contrast, 'WS')
        
        X = S(:, 17);
        
    elseif strcmp(contrast, 'UL')
        
        X = S(:, 18);
        
    elseif strcmp(contrast, 'DIDnIWDUL')
        
        X = 1.5*S(:, 14) + 0.5*S(:, 15) - 0.5*S(:, 16) - 1.5*S(:, 18);
        
    elseif strcmp(contrast, 'DIDnIWDWS')
        
        X = 1.5*S(:, 14) + 0.5*S(:, 15) - 0.5*S(:, 16) - 1.5*S(:, 17);
        
    elseif strcmp(contrast, 'DIDnIWDFix')
        
        X = 1.5*S(:, 14) + 0.5*S(:, 15) - 0.5*S(:, 16) - 1.5*S(:, 13);
        
    end
    
end

