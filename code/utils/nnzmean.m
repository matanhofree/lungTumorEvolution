function m = nnzmean(x,dim)
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


% Find zeros and set them to nan;

isZeros = x == 0;
x(isZeros) = nan;

if nargin == 1 % let sum deal with figuring out which dimension to use
    % Count up non-NaNs.
    m = nanmean(x);
else

    m = nanmean(x,dim);
end



