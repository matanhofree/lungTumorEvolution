function [recluster,reclustMap] = clusterReorderByVal(inCluster,inRef,refSub,inOpts)

    defaultOpts.equalFreqWeight = 1;
    defaultOpts.weightedMedian = 1;


    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end        
    clear defaultOpts;
    disp(opts);
    
    if ~isempty(refSub)      
        zClust = inCluster(refSub);
        refNum = inRef(refSub);
    else
        zClust = inCluster;
        refNum = inRef;        
    end

       
    if opts.equalFreqWeight == 1     
                        
        % [uniqType,~,~,uniqTypeCnt] = fastUnique(refNum);
        
        zType_weight = ones(size(refNum));
        for i = unique(refNum)'
            zSel = refNum == i;
            zType_weight(zSel) = 1/sum(zSel);
        end

        clustTypeMean = arrayfun(@(x)sum(refNum(zClust == x).*zType_weight(zClust == x))/sum(zType_weight(zClust == x)),unique(zClust));
        clustTypeMed = arrayfun(@(x)weightedMedian(refNum(zClust == x),zType_weight(zClust == x)),unique(zClust));
    else
        clustTypeMed = arrayfun(@(x)median(refNum(zClust == x)),unique(zClust))        
        clustTypeMean = arrayfun(@(x)mean(refNum(zClust == x)),unique(zClust))        
    end
    
    if opts.weightedMedian
        clustTypeMix = clustTypeMean;
    else
        clustTypeMix = clustTypeMed + 0.1*clustTypeMean;
    end        
    
    cV = unique(zClust);
    [~,zOrd] = sort(clustTypeMix);
    
    reclustMap = containers.Map(cV(zOrd),1:length(zOrd));
    recluster = nanvalues(reclustMap,inCluster);
end