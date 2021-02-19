function newCCnmf = updateOutConsensusNMF(ccNMFobj,includeGSEA,includeCC)

    if nargin < 2 || isempty(includeGSEA)
        includeGSEA = 0;
    end    

    if nargin < 3 || isempty(includeCC)
        includeCC = 0;
    end    

    
    if isfield(ccNMFobj,'ccNMFobjVer') && ccNMFobj.ccNMFobjVer == 2
        newCCnmf = ccNMFobj;
        return;
    end
    
    newCCnmf.Wlist = ccNMFobj.expandedW;
    newCCnmf.Hlist = ccNMFobj.extrapH;
    newCCnmf.Kdim = ccNMFobj.testK;
    
    
    newCCnmf.sampleID = ccNMFobj.sampleID;
    newCCnmf.geneID = ccNMFobj.geneID;
    
    newCCnmf.listVar = ccNMFobj.listVar;
    
    if isfield(ccNMFobj,'geneSymbol')
        newCCnmf.geneSymbol = ccNMFobj.geneSymbol;
    end
    
    tK = find(ccNMFobj.testK == size(ccNMFobj.consW,2))
    if length(tK) == 1
        newCCnmf.Kopt = tK;
    else
        [cdiff,zi] = min(cellfun(@(X)norm(ccNMFobj.consW-X,'fro'),ccNMFobj.mergedW(tK)));       
        newCCnmf.Kopt = tK(zi);
        
        if cdiff > 1e-3
            warning('Difference from consW and ref is %f',cdiff)
        end
    end
    
    newCCnmf.KoptMetrics = structSelectField(ccNMFobj,{'reconstError','reconstErrorPostOpt','cMeanSil'});
    
    if includeGSEA
        newCCnmf.geneSetE = structSelectField(ccNMFobj,{'ccEnrichment','ccTF'});    
    end
    
    if includeCC 
        newCCnmf.ccCompStruct = structSelectField(ccNMFobj,{'baseMergeClean','subStability'});    
    end
    
    newCCnmf.ccNMFobjVer = 2;    

end