function zDataNorm = calcNormTPM(zdata,zeroOut,denomFactor,doLog)
    
    if nargin < 4
        doLog = 1;
    end

    zsum = nansum(zdata);
    if ~exist('zeroOut','var')       
        zeroOut = 0;
    end        
    if ~exist('denomFactor','var') || isempty(denomFactor);
        denomFactor = full(median(zsum));
        fprintf(1,'Median factor %d\n',round(denomFactor));
    end
        
    zDataNorm = bsxfun(@rdivide,zdata,zsum)*denomFactor;
    if (zeroOut)
        fprintf('Zeroing output below 1\n');
        zDataNorm(zDataNorm < 1) = 0;
    end        
    
    if doLog
        zDataNorm = log1p(zDataNorm)/log(2);  
    end
end