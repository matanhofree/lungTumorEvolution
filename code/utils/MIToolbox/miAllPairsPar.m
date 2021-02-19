function output = miAllPairsPar(X,cores)
%function output = mi(X,Y)
%X & Y can be matrices which are converted into a joint variable
%before computation
%
%expects variables to be column-wise
%
%returns the mutual information between X and Y, I(X;Y)


numCols = size(X,2);
stepSize = floor(numCols/cores);
runVector = 1:stepSize:numCols;
%%
[xi,xj] = meshgrid(runVector);
xi = xi(tril(true(length(runVector))));
xj = xj(tril(true(length(runVector))));

%% Collect data
output_t = cell(length(xi));
parfor i = 1:length(xi)
    rowSelect = xi(i):min(xi(i)+stepSize-1,numCols);
    colSelect = xj(i):min(xj(i)+stepSize-1,numCols);
    output_t{i} = miAllPairsMexSub(2,X,int32(rowSelect),int32(colSelect));
end

% Assign data

output = nan(numCols);
for i = 1:length(xi)
    rowSelect = xi(i):min(xi(i)+stepSize-1,numCols);
    colSelect = xj(i):min(xj(i)+stepSize-1,numCols);
    output(rowSelect,colSelect) = output_t{i};
end

