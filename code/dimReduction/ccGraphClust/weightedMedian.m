function [wMed] = weightedMedian(D,W)

% ----------------------------------------------------------------------
% Function for calculating the weighted median
% Input:    D ... matrix of observed values
%           W ... matrix of weights, W = ( w_ij )
% Output:   wMed ... weighted median
% ----------------------------------------------------------------------


if nargin ~= 2
    error('weightedMedian:wrongNumberOfArguments', ...
      'Wrong number of arguments.');
end

if size(D) ~= size(W)
    error('weightedMedian:wrongMatrixDimension', ...
      'The dimensions of the input-matrices must match.');
end

% normalize the weights, such that: sum ( w_ij ) = 1
% (sum of all weights equal to one)

WSum = sum(W(:));
W = W / WSum;

[dSort,zidx] = sort(D);
wSort = W(zidx);


% vector for cumulative sums of the weights
sumVect = cumsum(wSort);


j = find(sumVect >= 0.5,1);
wMed = dSort(j);
% wValMed = wOrg(j);

% jW = (0.5-sumVect(j-1))/((0.5-sumVect(j-1))+(sumVect(j)-0.5));
% jMinus1W = ((sumVect(j)-0.5))/((0.5-sumVect(j-1))+(sumVect(j)-0.5));
% wValMed = dSort(j)*jW + dSort(j-1)*jMinus1W;

end