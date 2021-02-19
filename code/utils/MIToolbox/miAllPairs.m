function output = miAllPairs(X)
%function output = mi(X,Y)
%X & Y can be matrices which are converted into a joint variable
%before computation
%
%expects variables to be column-wise
%
%returns the mutual information between X and Y, I(X;Y)

% Xt = X.';


[output] = miAllPairsMex(1,X);
output(tril(true(size(X,2)))) = nan;