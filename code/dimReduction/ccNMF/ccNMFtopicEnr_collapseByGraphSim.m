function [enSubCollapse,zidxJ] = ccNMFtopicEnr_collapseByGraphSim(enSub,cSetGeneOverlap,opts)

        nGS = length(enSub.geneSetName);

        cDsim = 1-pdist(cSetGeneOverlap,'jaccard');
        cDsim(isnan(cDsim)) = 0;
        
        
        if opts.jaccardSetThreshold == 0 
           %%
            [cDsimOrd] = sort(cDsim,'descend');

            cjN = floor(nGS*opts.jaccardSetThresholdFactor);            
            if cjN > length(cDsimOrd) 
                fprintf('Choosing threshold failed');                        
                
                jThr = 0;
            else
                jThr = round(cDsimOrd(cjN),2);            
            end
            
            fprintf('Choosing threshold by rank %d heuristic=%f\n',cjN,jThr);                        


        end
        
        %%
        cDsim(cDsim < jThr) = 0;                    
        cSetGraph = graph(squareform(cDsim));
        
        %%
        cSetGraphIdx = conncomp(cSetGraph);
        
        [cIdxList,~,enSubBin,enSubBinCnt,enSubBinPos] = fastUnique(cSetGraphIdx);
        

        enSubCollapse = enSub;                
        enSubCollapse.enSubBinPos = enSubBin;
        enSubCollapse.sortBin = cell(nGS,1);
        
        cidxString = sprintf('%%0%dd_%%0%dd',floor(log10(length(cIdxList))+1),floor(log10(max(double(enSubBinCnt)))+1));
        
        % cidxString = sprintf('%%0%dd_%%0%dd',floor(log10(length(cIdxList))+1),floor(log10(max(double(enSubBinCnt)))+1));
        zidxJ = nan(nGS,1);
        isPrime = false(nGS,1);
        %%
        zi = 1;
        zp = 0;
        enSubBin = double(enSubBin);
        while any(~isnan(enSubBin))
            cList = unique(enSubBin)';            
            cList(isnan(cList)) = [];
            cM = enSubBin == cList;
            
            nL = length(cList);
            cFirst = arrayfun(@(x)find(cM(:,x),1,'first'),1:nL);
            
            zidxJ(cFirst(:)) = (zp + (1:nL))';
            enSubCollapse.sortBin(cFirst) = mergeStringPair(cidxString,cList,zi);                      
            
            if zi == 1
                isPrime(cFirst(:)) = 1;
            end
            
            enSubBin(cFirst) = nan;
            zi = zi + 1;
            zp = zp + nL;
        end      

%         enSubCollapse.sortAll = cell2mat(zidxJ');
        enSubCollapse.sortAll = zidxJ;
        enSubCollapse.isPrime = isPrime;
        
end