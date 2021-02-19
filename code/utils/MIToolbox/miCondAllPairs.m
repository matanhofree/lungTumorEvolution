function output = miCondAllPairs(X,Z)
%function output = mi(X,Y)
%X & Y can be matrices which are converted into a joint variable
%before computation
%
%expects variables to be column-wise
%
%returns the mutual information between X and Y, I(X;Y|Z)

% Xt = X.';


[output] = miAllPairsMexSub(3,X,Z);
output(tril(true(size(X,2)))) = nan;