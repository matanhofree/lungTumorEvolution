function output = miCondAllPairsPar(Xmat,Zmat,cores)
%function output = mi(X,Y)
%X & Y can be matrices which are converted into a joint variable
%before computation
%
%expects variables to be column-wise
%
%returns the mutual information between X and Y, I(X;Y)


numCols = size(Xmat,2);
stepSize = floor(numCols/cores);
runVector = 1:stepSize:numCols;
%%
[xi,xj] = meshgrid(runVector);
xi = xi(triu(true(length(runVector))));
xj = xj(triu(true(length(runVector))));

% if (any(size(Xmat) - size(Zmat)))
%     error('Dimensional missmatch');
% end

%% Collect data

output_t = cell(length(xi));

if (all(size(Xmat) == size(Zmat)))
    parfor i = 1:length(xi)
        rowSelect = xi(i):min(xi(i)+stepSize-1,numCols);
        colSelect = xj(i):min(xj(i)+stepSize-1,numCols);
        output_t{i} = miAllPairsMexSub(4,Xmat,Zmat,int32(rowSelect),int32(colSelect));
    end
elseif (size(Zmat,1) == size(Xmat,1) && size(Zmat,2) == 1)
    ZmatIdx = grp2idx(Zmat);
    parfor i = 1:length(xi)
        rowSelect = xi(i):min(xi(i)+stepSize-1,numCols);
        colSelect = xj(i):min(xj(i)+stepSize-1,numCols);
        output_t{i} = miAllPairsMexSub(5,Xmat,ZmatIdx,int32(rowSelect),int32(colSelect));
    end
else
    error('Dimension missmatch')
end
        
    
% Assign data
output = nan(numCols);
for i = 1:length(xi)
    rowSelect = xi(i):min(xi(i)+stepSize-1,numCols);
    colSelect = xj(i):min(xj(i)+stepSize-1,numCols);
    output(rowSelect,colSelect) = output_t{i};
end

