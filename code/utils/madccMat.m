function rmad = madccMat(x,y,doMean)
% Median Absolute Deviation Correlation Coefficient
%
% Impliments Median Absolute Deviation Correlation Coefficient, (as
% described in Shevlyakov & Smirnov (2011).  The estimate is
%     r = ( mad(Sp)^2 - mad(Sm)^2 ) / ( mad(Sp)^2 + mad(Sm)^2 )
% where
%     Sp = (x-m(x))/mad(x) + (y-m(y))/mad(y);
%     Sm = (x-m(x))/mad(x) - (y-m(y))/mad(y);
% and m() is median and mad() is the median absolute deviation,
%     mad(x) = m(abs(x-m(x)))
%
% For autocorrelation estimation, where we correlate lagged versions
% of a vector, we can assume mad(x)=mad(y), and divisors cancel and we
% can compute Sp & Sm as:
%     Sp = (x-m(x)) + (y-m(y));
%     Sm = (x-m(x)) - (y-m(y));
%
% If we further (reasonably) assume that m(Sp)=m(Sm)=0, then a slight
% further simplification is possible, where mad()'s can be replaced with
% m(abs())'s (Kharin et al., 2011): 
%     r = ( m(abs(Sp))^2 - m(abs(Sm))^2 ) / ( m(abs((Sp))^2 + m(abs((Sm)))^2 )
% This expression is what is computed below.
%
% REFERENCES
%
% Shevlyakov, G., & Smirnov, P. (2011). Robust estimation of a
% correlation coefficient: an attempt of survey. Australian & New
% Zealand Journal of Statistics, 40(1), 147?156.
% 
% Kharin, Y. S., & Voloshko, V. A. (2011). Robust estimation of AR 
% coefficients under simultaneously influencing outliers and missing 
% values. Journal of Statistical Planning and Inference, 141(9),
% 3276?3288.
%
% 2014-07-08
% Thomas Nichols http://warwick.ac.uk/tenichols

    if nargin < 3
        doMean = 1;
    end

    [xN,xD] = size(x);
    [yN,yD] = size(y);
        
    rmad = nan(xD,yD);
    
    for i = 1:xD
        cx = full(x(:,i));
        for j = 1:yD
            cy = full(y(:,j));
            if any(isnan(cx)) || any(isnan(cy))
                rmad(i,j) = nanmadcc(cx,cy,doMean);
            else
                rmad(i,j) = madcc(cx,cy,doMean);
            end
        end
    end
    
    

end 


function rmad = madcc(x,y,doMean)

        if doMean == 1
            mx    = mean(x);
            my    = mean(y);
        else        
            mx    = median(x);
            my    = median(y);
        end

        madx  = median(abs(x-mx));
        mady  = median(abs(y-my));
        Sp    = (x-mx)/madx + (y-my)/mady;
        Sm    = (x-mx)/madx - (y-my)/mady;
        madSp = median(abs(Sp-median(Sp)));
        madSm = median(abs(Sm-median(Sm)));
        rmad = (madSp^2 - madSm^2)/(madSp^2 + madSm^2);
    
end

function rmad = nanmadcc(x,y,doMean)

    rmad=NaN;
    I=find(all(~isnan([x(:) y(:)]),2));

    if ~isempty(I)

        if doMean == 1
            mx    = mean(x(I));
            my    = mean(y(I));
        else        
            mx    = median(x(I));
            my    = median(y(I));
        end

        madx  = median(abs(x(I)-mx));
        mady  = median(abs(y(I)-my));
        Sp    = (x(I)-mx)/madx + (y(I)-my)/mady;
        Sm    = (x(I)-mx)/madx - (y(I)-my)/mady;
        madSp = median(abs(Sp-median(Sp)));
        madSm = median(abs(Sm-median(Sm)));
        rmad = (madSp^2 - madSm^2)/(madSp^2 + madSm^2);
    end
end