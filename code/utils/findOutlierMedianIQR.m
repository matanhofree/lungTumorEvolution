function [isOutlier,statMedian,statIQR,gNames,isTested] = findOutlierMedianIQR(stat,groupV,rFactor,lowerB,upperB,minGroup)
    
    stat = full(stat(:));
    groupV = groupV(:);
    
    if nargin < 4
        lowerB = 1;
        upperB = 1;
    elseif nargin < 5 
        upperB = 1;        
    end
    
    if nargin < 6 
        minGroup = 0;
    end

    if ~exist('groupV','var') || isempty(groupV)
        groupV = ones(size(stat));        
    end

    if ~exist('rFactor','var') || isempty(rFactor)
        rFactor = 1.5;        
    end

    nD = length(stat);
    
    %[gV,gNames] = grp2idx(groupV);
    %%
    [gNames,~,gV,cntGroup,groupPos] = fastUnique(groupV);
    %%

    % gVlist = unique(gV);
        
    statMedian = cellfun(@(x)quantile(stat(x),[0.25 0.5 0.75]),groupPos,'uniformoutput',0);
    statMedian = cell2mat(statMedian);
    
    statMedian(cntGroup<minGroup,:) = nan;
    isTested = ismember(gV,find(~(cntGroup<minGroup)));
    
    
    % statIQR = arrayfun(@(x)iqr(stat(gV == x)),gVlist);
    % ratioFactor = (stat - statMedian(gV,2))./statIQR(gV);    
    
    statIQR = statMedian(:,3) - statMedian(:,1);    
    
    
    %%
    
    gV = double(gV);   
    if lowerB 
        lowerThr = statMedian(:,1) - statIQR*rFactor;
        
        isOutlier = (stat < lowerThr(gV));        
    else
        isOutlier = false(nD,1);
    end
    
    if upperB
        upperThr = statMedian(:,3) + statIQR*rFactor;
        isOutlier = isOutlier | (stat > upperThr(gV));     
    end
end