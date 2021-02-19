function m = nnzcenter(x,zMean)
%NANMEAN Mean value, ignoring NaNs and zeros.
%   M = NANMEAN(X) returns the sample mean of X, treating NaNs as missing
%   values.  For vector input, M is the mean value of the non-NaN elements
%   in X.  For matrix input, M is a row vector containing the mean value of
%   non-NaN elements in each column.  For N-D arrays, NANMEAN operates
%   along the first non-singleton dimension.
%
%   NANMEAN(X,DIM) takes the mean along dimension DIM of X.
%
%   See also MEAN, NANMEDIAN, NANSTD, NANVAR, NANMIN, NANMAX, NANSUM.

%   Copyright 1993-2004 The MathWorks, Inc.


% Find NaNs and set them to zero

isZero = x == 0;
isNaNvect = isnan(zMean);

if size(zMean,1) == 1    
    m = bsxfun(@minus,x,zMean);
    m(:,isNaNvect) = x(:,isNaNvect); 
elseif  size(zMean,2) == 1
    m = bsxfun(@minus,x,zMean);
    m(isNaNvect,:) = x(isNaNvect,:); 
end
    
m(isZero) = 0;
